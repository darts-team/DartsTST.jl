
using NCDatasets
using Plots
using Statistics
using JLD2 # note: may have to Pkg.add("JLD2")
using Distributed, SharedArrays

maxprocs = 81 # maximum number of cores to use
curr_procs = nprocs()
if curr_procs < maxprocs
    addprocs(maxprocs - curr_procs)
    
    # not adding the else case where we would remove processes. Possible but not useful
end#if
curr_procs = nprocs()
println("Current procs: " * "$curr_procs")
## includes
@everywhere begin
    include("../modules/generate_raw_data.jl")
    include("../modules/process_raw_data.jl")
    include("../modules/geometry.jl")
    include("../modules/scene.jl")
    include("../modules/sync.jl")
    include("../modules/range_spread_function.jl") # as RSF
    include("../modules/orbits.jl")
    include("../modules/error_sources.jl")
    include("../modules/performance_metrics.jl")
end#begin
# @everywhere using Orbits, Performance_Metrics, Generate_Raw_Data, Geometry, Sync, Scene, RSF, Process_Raw_Data, Error_Sources

# @everywhere include("../darts-simtool/modules/generate_raw_data.jl")
# @everywhere include("../darts-simtool/modules/process_raw_data.jl")
# @everywhere include("../darts-simtool/modules/geometry.jl")
# @everywhere include("../darts-simtool/modules/scene.jl")
# @everywhere include("../darts-simtool/modules/sync.jl")
# @everywhere include("../darts-simtool/modules/range_spread_function.jl") # as RSF
# @everywhere include("../darts-simtool/modules/orbits.jl")
# @everywhere include("../darts-simtool/modules/error_sources.jl")
# @everywhere include("../darts-simtool/modules/performance_metrics.jl")
## Determining Parameters
c=299792458 # speed of light (m/s)
# planetary shape constants
earth_radius=6378.137e3 # semi-major axis at equator
earth_eccentricity=sqrt(0.00669437999015)
# MIMO parameters
mode=3 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)
# radar parameters
fc=1e9 # center frequency (Hz)
fp=100 # pulse repetition frequency (Hz)
SNR=50 # SNR for single platform and single pulse before fast-time processing dB (for additive random noise only) TODO calculate based on sigma-zero (which depends on target type, wavelength, look angle, polarization) and NESZ (which depends on radar specs and processing)
# platform locations in xyz taken from orbits (including slow-time)
orbit_filename="inputs/orbitOutput_082020.nc" # position in km, time in sec
SAR_duration=3 # synthetic aperture duration (s)
SAR_start_time=0 # SAR imaging start time (s)
# target locations (volumetric grid) defined in geo (θϕh)
# t_θ=0 # deg latitude
# t_ϕ=0 # deg longitude
# t_h=0 # m  heights
# # image/scene pixel coordinates
# s_θ=-0.0001:0.000001:0.0001 # deg latitude
# s_ϕ=-0.0005:0.00001:0.0005 # deg longitude
# s_h=-40:1:40 # m  heights

# target locations (volumetric grid) defined in geo (θϕh) - for slant looking test
t_θ=7 # deg latitude
t_ϕ=0 # deg longitude
t_h=0 # m  heights
# image/scene pixel coordinates
# s_θ=7-0.0003:0.0000025:7+0.0003 # deg latitude
# s_ϕ=-0.001:0.00001:0.001 # deg longitude
# s_h=-35:0.5:35 # m  heights
s_θ=7-0.0005:0.00001:7+0.0005 # deg latitude
s_ϕ=-0.002:0.0001:0.002 # deg longitude
s_h=-25:0.5:25 # m  heights

# range spread function (RSF) parameters
enable_fast_time = true # whether to enable or disable fast-time axis, 0:disable, 1: enable
enable_thermal_noise=false # whether to enable or disable random additive noise (e.g. thermal noise)
disable_freq_offset = false # true = no linear phase ramp (ideal osc frequency), false = linear phase ramp error
Trx=300e-6 # s duration of RX window (may need to be increased if aperture or scene is large) TODO (adjust based on max/min range)
pulse_length=10e-6 # s pulse length
Δt=1e-8 # s fast-time resolution (ADC sampling rate effect is excluded for now)
bandwidth=40e6 # bandwidth (Hz)
# performance metrics
res_dB=3 # dB two-sided resolution relative power level (set to 0 for peak-to-null Rayleigh resolution), positive value needed
PSF_peak_target=2


sync_processing_time = 0.001 # processing time between stage 1 and stage 2 sync
sync_signal_len = 1024 # waveform length
sync_fc = 1e9 # waveform center frequency
sync_fs = 25e6; # sync receiver sampling rate
sync_fbw = sync_fs # LFM bandwidth

# osc_type = "USO" # putting a oscillator type variable here to auto-name save files
osc_type = "USRP"
#defines oscillator quality. Either leave as single row to use across all platforms, or define values for each platform as a new row
if osc_type == "USO"
    a_coeff_dB = [-95 -90 -200 -130 -155] # [USO: Krieger]
elseif osc_type == "USRP"
    a_coeff_dB = [-28 -40 -200 -130 -155] # [USRP E312]
end
osc_coeffs = 0 # temporary holding value, will get overwritten

if disable_freq_offset == 1 # option to remove linear phase drift due to osc frequency offset
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



Ntrials = 80 # number of trials per SRI in Monte Carlo simulations
sync_PRIs = [.1 1 2 5 10]
numSRI = length(sync_PRIs)

## RANGE SPREAD FUNCTION (matched filter output)
if enable_fast_time # matched filter gain is included in Srx
    Srx,MF,ft,t_rx=RSF.ideal_RSF(pulse_length,Δt,bandwidth,Trx) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
end

## PLATFORM LOCATIONS
orbit_dataset=Dataset(orbit_filename) # Read orbits data in NetCDF format
t12_orbits=orbit_dataset["time"][1:2] # first two time samples
dt_orbits=t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
orbit_time_index=(Int(round(SAR_start_time/dt_orbits))+1:1:Int(round((SAR_start_time+SAR_duration)/dt_orbits))+1) # index range for orbit times for time interval of interest
orbit_time=orbit_dataset["time"][orbit_time_index] # read in time data
orbit_pos=orbit_dataset["position"][:,:,orbit_time_index] # read in position data
#TODO convert ECI to ECEF
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

## TARGET LOCATIONS
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,earth_radius,earth_eccentricity)

## GENERATE RAW DATA
ref_range=Generate_Raw_Data.distance(mean(t_xyz_grid,dims=2),mean(mean(p_xyz,dims=2),dims=3)) # reference range

if enable_fast_time # with fastime and slowtime; matched filter gain is included in Srx
    rawdata=Generate_Raw_Data.main_RSF_slowtime(t_xyz_grid,p_xyz,mode,tx_el,fc,Srx,t_rx,ref_range)
else # without fastime, with slowtime; matched filter gain is included inside the function
    rawdata=Generate_Raw_Data.main_noRSF_slowtime(t_xyz_grid,p_xyz,mode,tx_el,fc)
end

if !enable_fast_time
    SNR=SNR*pulse_length*bandwidth # SNR increases after matched filter
end

## IMAGE SCENE
Ns_θ=length(s_θ)
Ns_ϕ=length(s_ϕ)
Ns_h=length(s_h)
s_geo_grid=Scene.form3Dgrid_for(s_θ,s_ϕ,s_h) # using 3 nested for loops
s_xyz_grid=Geometry.geo_to_xyz(s_geo_grid,earth_radius,earth_eccentricity)


## find Ideal case results first
#PROCESS RAW DATA TO GENERATE IMAGE
if enable_fast_time # with fastime, with slowtime
    image_1xN=Process_Raw_Data.main_RSF_slowtime(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc,t_rx,ref_range)
else # without fastime, with slowtime
    image_1xN=Process_Raw_Data.main_noRSF_slowtime(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc)
end
image_3D=Scene.convert_image_1xN_to_3D(image_1xN,Ns_θ,Ns_ϕ,Ns_h)

# PERFORMANCE METRICS
(ideal_peak, ideal_idx) = findmax(image_3D) # finds maximum and index of max

if size(t_xyz_grid)[2]==1 # PSF related performance metrics are calculated when there is only one point target
    target_index1=findall(t_θ .==s_θ)
    target_index2=findall(t_ϕ .==s_ϕ)
    target_index3=findall(t_h .==s_h)
    if isempty(target_index1) || isempty(target_index2) || isempty(target_index3)
        show("PSF related performance metrics cannot be calculated since target is not inside the scene!")
        PSF_metrics=false
    else
        PSF_metrics=true
        target_location=[t_θ t_ϕ t_h] # point target location
        ideal_res,ideal_PSLR,ideal_ISLR,loc_error=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_θ,s_ϕ,s_h,PSF_peak_target) # resolutions in each of the 3 axes
    end
else
    PSF_metrics=false
    show("PSF related performance metrics cannot be calculated since there are more than 1 targets!")
end

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
            image_1xN=Process_Raw_Data.main_RSF_slowtime(rawdata_sync,s_xyz_grid,p_xyz,mode,tx_el,fc,t_rx,ref_range)
        else # without fastime, with slowtime
            image_1xN=Process_Raw_Data.main_noRSF_slowtime(rawdata_sync,s_xyz_grid,p_xyz,mode,tx_el,fc)
        end
        image_3D=Scene.convert_image_1xN_to_3D(image_1xN,Ns_θ,Ns_ϕ,Ns_h)
    
        ## PERFORMANCE METRICS
        # PSF metrics
        if size(t_xyz_grid)[2]==1 # PSF related performance metrics are calculated when there is only one point target
            target_location=[t_θ t_ϕ t_h] # point target location
            resolution,PSLR,ISLR,loc_error=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_θ,s_ϕ,s_h,PSF_peak_target) # resolutions in each of the 3 axes
        else
            resolution=[NaN,NaN,NaN]
            PSLR=[NaN,NaN,NaN]
            ISLR=[NaN,NaN,NaN]
            # println("PSF related performance metrics cannot be calculated since there are more than 1 targets!")
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
@save outputfilename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_θ s_ϕ s_h
#println(std(resolutions[1,:]))
# Note: JLD2 can be read using "@load filename var1 var2...
println("Run Complete, and file saved to " *outputfilename)


