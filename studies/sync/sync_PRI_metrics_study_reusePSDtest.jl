
using NCDatasets
using Plots
using Statistics
using JLD2 # note: may have to Pkg.add("JLD2")
using Distributed, SharedArrays

maxprocs = 65 # maximum number of cores to use
curr_procs = nprocs()
if curr_procs < maxprocs
    addprocs(maxprocs - curr_procs)
    
    # not adding the else case where we would remove processes. Possible but not useful
end#if
curr_procs = nprocs()
println("Current procs: " * "$curr_procs")
## includes
@everywhere begin
    include("../../modules/generate_raw_data.jl")
    include("../../modules/process_raw_data.jl")
    include("../../modules/geometry.jl")
    include("../../modules/scene.jl")
    include("../../modules/sync.jl")
    include("../../modules/range_spread_function.jl") # as RSF
    include("../../modules/orbits.jl")
    include("../../modules/error_sources.jl")
    include("../../modules/performance_metrics.jl")
end#begin

## Determining Parameters
c=299792458 # speed of light (m/s)
# planetary shape constants
earth_radius=6378.137e3 # semi-major axis at equator
earth_eccentricity=sqrt(0.00669437999015)
# MIMO parameters
mode=2 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)
# radar parameters
fc=1.25e9 # center frequency (Hz)
fp=10 # pulse repetition frequency (Hz)
SNR=50 # SNR for single platform and single pulse before fast-time processing dB (for additive random noise only) TODO calculate based on sigma-zero (which depends on target type, wavelength, look angle, polarization) and NESZ (which depends on radar specs and processing)
# platform locations in xyz taken from orbits (including slow-time)
orbit_filename="inputs/orbitOutput_082020.nc" # position in km, time in sec

SAR_duration=5 # synthetic aperture duration (s)
SAR_start_time=0 # SAR imaging start time (s)
# target locations and reflectvities

target_pos_mode="CR" #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
ts_coord_sys="SCH" # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene
if ts_coord_sys=="SCH" # if SCH, target and scene locations are defined relative to the point where look angle vector intersects the surface
    look_angle=30 # in cross-track direction, required only if SCH coordinates, using same look angle for targets and scene (deg)
    p_avg_heading=0.1 # average heading of platforms, due North is 0, due East is 90 (deg), required only if SCH coordinates TODO we should get this from orbits!
end
if target_pos_mode=="grid" # target positions are defined as a volumetric grid (useful for distributed target)
    t_loc_1=-10:10:10 # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=-20:20:20 # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=0:30:30 # m  heights if LLH or SCH, Z if XYZ
    t_ref=rand(Float64,length(t_loc_1),length(t_loc_2),length(t_loc_3)) # uniform random reflectivities between 0 and 1, a 3D input array (e.g. 3D image) can be used instead
elseif target_pos_mode=="CR" # ("CR" for corner reflector) target positions are defined as 3xN array (useful for a few discrete targets)
    # length(t_loc_1)==length(t_loc_2)==length(t_loc_3) should hold
    t_loc_1=[0] # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=[0] # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=[40] # m  heights if LLH or SCH, Z if XYZ
    t_ref=  [1] # reflectivities
end
# image/scene pixel coordinates
s_loc_1=-40:.5:40 # deg latitude if LLH, along-track if SCH, X if XYZ
s_loc_2=-60:1:60 # deg longitude if LLH, cross-track if SCH, Y if XYZ
s_loc_3=  0:1:80 # m  heights if LLH or SCH, Z if XYZ
# range spread function (RSF) parameters
pulse_length=10e-6 # s pulse length
Δt=1e-8 # s fast-time resolution (ADC sampling rate effect is excluded for now)
bandwidth=40e6 # bandwidth (Hz)
# performance metrics
res_dB=3 # dB two-sided resolution relative power level (set to 0 for peak-to-null Rayleigh resolution), positive value needed
PSF_peak_target=1 # 1: peak, 2: target
# simulation options
enable_thermal_noise=false # whether to enable or disable random additive noise (e.g. thermal noise)
enable_fast_time=true # whether to enable or disable fast-time axis, 0:disable, 1: enable
display_geometry=false # whether to display geometry plots
display_RSF_rawdata=false # whether to display RSF and rawdata plots
display_tomograms=0 # how to display tomograms, 0: do not display, 1: display only 3 slices at the scene center, 2: display all slices in each dimension, 3: display as 3D scatter plot
disable_freq_offset = false # true = no linear phase ramp (ideal osc frequency), false = linear phase ramp error

sync_processing_time = 0.001 # processing time between stage 1 and stage 2 sync
sync_signal_len = 1024 # waveform length
sync_fc = 1e9 # waveform center frequency
sync_fs = 25e6; # sync receiver sampling rate
sync_fbw = sync_fs # LFM bandwidth

osc_type = "USO" # putting a oscillator type variable here to auto-name save files
# osc_type = "USRP"
#defines oscillator quality. Either leave as single row to use across all platforms, or define values for each platform as a new row
if osc_type == "USO"
    a_coeff_dB = [-95 -90 -200 -130 -155] # [USO: Krieger]
elseif osc_type == "USRP"
    a_coeff_dB = [-28 -40 -200 -130 -155] # [USRP E312]
end
osc_coeffs = 0 # temporary holding value, will get overwritten

if disable_freq_offset == true # option to remove linear phase drift due to osc frequency offset
    sigma_freq_offsets = zeros(1)
else
    sigma_freq_offsets = 1.5e-3 # Hz - std. dev. of the frequency offset of the oscillator. This is the linear phase ramp value
    sigma_freq_offsets = sigma_freq_offsets .* ones(1) # convert to matrix form, one value for each oscillator
end

sync_fmin   = 0.01 # minimum frequency > 0 in Hz to window PSD
f_osc       = 10e6 # local oscillator frequency
sync_clk_fs = 1e3; # sample rate of clock error process
master      = 1 # selection of master transmitter for sync (assumes a simplified communication achitecture- all talking with one master platform)
sync_pri    = 1 # to be overwritten in loop

no_sync_flag = false; # if flag == true, no sync is used. flag == false results in normal sync process estimation
## make a struct of important input parameters
#list key parameters in here, they will get passed to most(?) modules
@everywhere mutable struct keyParameters
    # radar parameters
    mode #1: SAR (ping-pong), 2:SIMO, 3:MIMO
    tx_el # which element transmits for SIMO (max value N)
    fc # center frequency (Hz)
    fp # pulse repetition frequency (Hz)
    # TODO: add others... set defaults

    # clock parameters
    sync_pri
    sync_processing_time
    sync_signal_len
    sync_fc
    sync_fs
    sync_fbw
    sync_fmin
    f_osc
    sync_clk_fs
    master
    osc_coeffs
    sigma_freq_offsets

    no_sync_flag
end

## define instance of parameter structure
parameters = keyParameters(mode,
tx_el,
fc,
fp,
sync_pri,
sync_processing_time,
sync_signal_len,
sync_fc,
sync_fs,
sync_fbw,
sync_fmin,
f_osc,
sync_clk_fs,
master,
osc_coeffs,
sigma_freq_offsets,
no_sync_flag)



Ntrials = 64 # number of trials per SRI in Monte Carlo simulations
sync_PRIs = [.1 1 2 3 4 5]
numSRI = length(sync_PRIs)

## PLATFORM LOCATIONS
orbit_dataset=Dataset(orbit_filename) # Read orbits data in NetCDF format
orbit_dataset=Dataset(orbit_filename) # Read orbits data in NetCDF format

t12_orbits=orbit_dataset["time"][1:2] # first two time samples
dt_orbits=t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
orbit_time_index=(Int(round(SAR_start_time/dt_orbits))+1:1:Int(round((SAR_start_time+SAR_duration)/dt_orbits))+1) # index range for orbit times for time interval of interest
orbit_time=orbit_dataset["time"][orbit_time_index] # read in time data
orbit_pos=orbit_dataset["position"][:,:,orbit_time_index] # read in position data #TODO convert ECI to ECEF?
slow_time=(SAR_start_time:1/fp:SAR_start_time+SAR_duration) # create slow time axis
orbit_pos_interp=Orbits.interp_orbit(orbit_time,orbit_pos,slow_time) # interpolate orbit to slow time
p_xyz=1e3*orbit_pos_interp # convert km to m
Np=size(orbit_pos)[2] # number of platforms
Nst=size(slow_time)[1] # number of slow-time samples (pulses processed)

parameters.osc_coeffs = repeat(a_coeff_dB,Np) # put together all oscillator coeffs
if disable_freq_offset == 1 # option to remove linear phase drift due to osc frequency offset
    sigma_freq_offsets = zeros(Np)
else
    sigma_freq_offsets = 1.5e-3 # Hz - std. dev. of the frequency offset of the oscillator. This is the linear phase ramp value
    sigma_freq_offsets = sigma_freq_offsets .* ones(Np) # convert to matrix form, one value for each oscillator
end
parameters.sigma_freq_offsets=sigma_freq_offsets # overwrite here

orbit_pos_all=reshape(p_xyz,3,Np*Nst) # platform positions in xyz; for each platform, its position at each PRF treated as a different platform; output loops over platforms first, then slow-time

## TARGET/SCENE LOCATIONS
targets,Nt=Scene.construct_targets_str(target_pos_mode,t_loc_1,t_loc_2,t_loc_3,t_ref) # Nt: number of targets, targets: structure array containing target locations and reflectivities
targets_loc=zeros(3,Nt);for i=1:Nt;targets_loc[:,i]=targets[i].loc;end # 3xN
s_loc_3xN=Scene.form3Dgrid_for(s_loc_1,s_loc_2,s_loc_3) # using 3 nested for loops
t_xyz_3xN,s_xyz_3xN=Scene.convert_target_scene_coord_to_XYZ(ts_coord_sys,s_loc_3xN,targets_loc,p_xyz,look_angle,p_avg_heading,earth_radius,earth_eccentricity)
## TARGET REFLECTIVITIES
targets_ref=zeros(1,Nt);for i=1:Nt;targets_ref[i]=targets[i].ref;end

## RANGE SPREAD FUNCTION (matched filter output)
min_range,max_range=Geometry.find_min_max_range(t_xyz_3xN,p_xyz)
Trx=2*(max_range-min_range)/c+2*pulse_length # s duration of RX window
if enable_fast_time # matched filter gain is included in Srx
    Srx,MF,ft,t_rx=RSF.ideal_RSF(pulse_length,Δt,bandwidth,Trx) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
    # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing
end

## GENERATE RAW DATA
ref_range=Geometry.distance(mean(t_xyz_3xN,dims=2),mean(mean(p_xyz,dims=2),dims=3)) # reference range (equal to slant_range in sch?)
if enable_fast_time # with fastime and slowtime; matched filter gain is included in Srx
    rawdata=Generate_Raw_Data.main_RSF_slowtime(t_xyz_3xN,p_xyz,mode,tx_el,fc,Srx,t_rx,ref_range,targets_ref) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
else # without fastime, with slowtime; matched filter gain is included inside the function
    rawdata=Generate_Raw_Data.main_noRSF_slowtime(t_xyz_3xN,p_xyz,mode,tx_el,fc,targets_ref) # rawdata is a: 2D array of size Nst x Np (SAR/SIMO), 3D array of size Nst x Np(RX) x Np(TX) (MIMO)
end
if !enable_fast_time;SNR=SNR*pulse_length*bandwidth;end # SNR increases after matched filter
if enable_thermal_noise;rawdata=Error_Sources.random_noise(rawdata,SNR,enable_fast_time,mode);end # adding random noise based on SNR after range (fast-time) processing

## find Ideal case results first
# PROCESS RAW DATA TO GENERATE IMAGE
if enable_fast_time # with fastime, with slowtime
    image_1xN=Process_Raw_Data.main_RSF_slowtime(rawdata,s_xyz_3xN,p_xyz,mode,tx_el,fc,t_rx,ref_range)
else # without fastime, with slowtime
    image_1xN=Process_Raw_Data.main_noRSF_slowtime(rawdata,s_xyz_3xN,p_xyz,mode,tx_el,fc)
end
Ns_1=length(s_loc_1);Ns_2=length(s_loc_2);Ns_3=length(s_loc_3)
image_3D=Scene.convert_image_1xN_to_3D(image_1xN,Ns_1,Ns_2,Ns_3)

# PERFORMANCE METRICS
(ideal_peak, ideal_idx) = findmax(image_3D) # finds maximum and index of max

if size(t_xyz_3xN,2)==1 # PSF related performance metrics are calculated when there is only one point target
    target_index1=findall(t_loc_1 .==s_loc_1)
    target_index2=findall(t_loc_2 .==s_loc_2)
    target_index3=findall(t_loc_3 .==s_loc_3)
    if isempty(target_index1) || isempty(target_index2) || isempty(target_index3)
        show("PSF related performance metrics cannot be calculated since target is not inside the scene!")
        PSF_metrics=false
    else
        PSF_metrics=true
        target_location=[t_loc_1 t_loc_2 t_loc_3] # point target location
        ideal_res,ideal_PSLR,ideal_ISLR,loc_error=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_loc_1,s_loc_2,s_loc_3,PSF_peak_target) # resolutions in each of the 3 axes
    end
else
    PSF_metrics=false
    show("PSF related performance metrics cannot be calculated since there are more than 1 targets!")
end

# Precalculate all the PSDs, do this in parallel
sync_PSDs = SharedArray{Float64}(numSRI,undef)




#---- Parallel portion -----
## Initialize result vectors
peaks       = SharedArray{Float64}(numSRI,Ntrials)
resolutions = SharedArray{Float64}(3,numSRI,Ntrials)
PSLRs       = SharedArray{Float64}(3,numSRI,Ntrials)
ISLRs       = SharedArray{Float64}(3,numSRI,Ntrials)
loc_errors  = SharedArray{Float64}(3,numSRI,Ntrials)
## run trials
@sync @distributed for ntrial = 1 : Ntrials
# Threads.@threads for ntrial = 1 : Ntrials
    # println("Trial Number: ", ntrial)
    
     for k = 1 : numSRI
        sync_pri = sync_PRIs[k]
        parameters.sync_pri = sync_pri
        println("Starting SRI value: ", sync_pri)
    
        ## add in error sources
        rawdata_sync = Error_Sources.synchronization_errors(rawdata,slow_time,orbit_pos_interp,enable_fast_time,parameters)
    
        ## PROCESS RAW DATA TO GENERATE IMAGE
        if enable_fast_time # with fastime, with slowtime
            image_1xN=Process_Raw_Data.main_RSF_slowtime(rawdata_sync,s_xyz_3xN,p_xyz,mode,tx_el,fc,t_rx,ref_range)
        else # without fastime, with slowtime
            image_1xN=Process_Raw_Data.main_noRSF_slowtime(rawdata_sync,s_xyz_3xN,p_xyz,mode,tx_el,fc)
        end
        Ns_1=length(s_loc_1);Ns_2=length(s_loc_2);Ns_3=length(s_loc_3)
        image_3D=Scene.convert_image_1xN_to_3D(image_1xN,Ns_1,Ns_2,Ns_3)
    
        ## PERFORMANCE METRICS
        # PSF metrics
        if size(t_xyz_3xN,2)==1 # PSF related performance metrics are calculated when there is only one point target
            target_index1=findall(t_loc_1 .==s_loc_1)
            target_index2=findall(t_loc_2 .==s_loc_2)
            target_index3=findall(t_loc_3 .==s_loc_3)
            if isempty(target_index1) || isempty(target_index2) || isempty(target_index3)
                show("PSF related performance metrics cannot be calculated since target is not inside the scene!")
                resolution=[NaN,NaN,NaN]
                PSLR=[NaN,NaN,NaN]
                ISLR=[NaN,NaN,NaN]
                PSF_metrics=false
                loc_error=[NaN,NaN,NaN]
            else
                PSF_metrics=true
                target_location=[t_loc_1 t_loc_2 t_loc_3] # point target location
                try
                    resolution,PSLR,ISLR,loc_error=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_loc_1,s_loc_2,s_loc_3,PSF_peak_target) # resolutions in each of the 3 axes
                catch
                    resolution=[NaN,NaN,NaN]
                    PSLR=[NaN,NaN,NaN]
                    ISLR=[NaN,NaN,NaN]
                    PSF_metrics=false
                    loc_error=[NaN,NaN,NaN]
                    show("PSF related performance metrics cannot be calculated -- error in metric calculation.")
                end#try
            end#if
        else
            resolution=[NaN,NaN,NaN]
            PSLR=[NaN,NaN,NaN]
            ISLR=[NaN,NaN,NaN]
            PSF_metrics=false
            loc_error=[NaN,NaN,NaN]
            show("PSF related performance metrics cannot be calculated since there are more than 1 targets!")
        end
        
        (peak, idx)             = findmax(image_3D) # finds maximum and index of max
        peaks[k,ntrial]         = peak
        loc_errors[:,k,ntrial]  = loc_error
        resolutions[:,k,ntrial] = resolution
        PSLRs[:,k,ntrial]       = PSLR
        ISLRs[:,k,ntrial]       = ISLR
    end#N SRIs
    
end#Ntrials

if disable_freq_offset # test to denote frequency error or not in save file name
    freq_text = "noFreq"
else
    freq_text = "wFreq"
end

outputfilename = "syncModule_MonteCarlo_mode_$mode"*"_$osc_type"*"_sync_pri_sweep_"*freq_text* ".jld2" # this is the output filename that the data is saved to using JLD2
# this saves the data into a JLD2 file. Data includes the estimates
@save outputfilename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
#println(std(resolutions[1,:]))
# Note: JLD2 can be read using "@load filename var1 var2...
println("Run Complete, and file saved to " *outputfilename)


