include("modules/generate_raw_data.jl")
include("modules/process_raw_data.jl")
include("modules/geometry.jl")
include("modules/scene.jl")
include("inputs/input_parameters_sync_monte_carlo.jl")
include("modules/range_spread_function.jl") # as RSF
include("modules/orbits.jl")
include("modules/sync.jl")
include("modules/error_sources.jl")
include("modules/Performance_Metrics.jl")

using NCDatasets
using Plots
using Statistics
using JLD2 # note: may have to Pkg.add("JLD2")

#pyplot()
plotly()
#gr()
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
    image_3xN=Process_Raw_Data.main_RSF_slowtime(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc,t_rx,ref_range)
else # without fastime, with slowtime
    image_3xN=Process_Raw_Data.main_noRSF_slowtime(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc)
end
image_3D=Scene.convert_image_3xN_to_3D(image_3xN,Ns_θ,Ns_ϕ,Ns_h)

# PERFORMANCE METRICS
# PSF metrics
if size(t_xyz_grid)[2]==1 # PSF related performance metrics are calculated when there is only one point target
    target_location=[t_θ t_ϕ t_h] # point target location
    ideal_res,ideal_PSLR,ideal_ISLR,loc_errors=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_θ,s_ϕ,s_h,PSF_peak_target) # resolutions in each of the 3 axes
else
    resolution=[NaN,NaN,NaN]
    PSLR=[NaN,NaN,NaN]
    ISLR=[NaN,NaN,NaN]
    println("PSF related performance metrics cannot be calculated since there are more than 1 targets!")
end
(ideal_peak, ideal_idx) = findmax(image_3D) # finds maximum and index of max
ideal_peak_idx1         = Int64(ideal_idx[1])
ideal_peak_idx2         = Int64(ideal_idx[2])
ideal_peak_idx3         = Int64(ideal_idx[3])

## Initialize up result vectors
peaks       = zeros(Ntrials)
peak_idx1   = zeros(Int64,Ntrials) # save the index, then convert to position based on axes
peak_idx2   = zeros(Int64,Ntrials)
peak_idx3   = zeros(Int64,Ntrials)
resolutions = zeros(3,Ntrials)
PSLRs       = zeros(3,Ntrials)
ISLRs       = zeros(3,Ntrials)
loc_errors  = zeros(3,Ntrials)
## run trials
for ntrial = 1 : Ntrials
    println("Trial Number: ", ntrial)
    ## add in error sources, loop over Ntrials and number of SRIs
    if add_phase_errors == true
        include("inputs/input_parameters_sync.jl")
        rawdata_sync = Error_Sources.synchronization_errors(rawdata,slow_time,orbit_pos_interp,enable_fast_time,parameters)
    end

    ## PROCESS RAW DATA TO GENERATE IMAGE
    if enable_fast_time # with fastime, with slowtime
        image_3xN=Process_Raw_Data.main_RSF_slowtime(rawdata_sync,s_xyz_grid,p_xyz,mode,tx_el,fc,t_rx,ref_range)
    else # without fastime, with slowtime
        image_3xN=Process_Raw_Data.main_noRSF_slowtime(rawdata_sync,s_xyz_grid,p_xyz,mode,tx_el,fc)
    end
    image_3D=Scene.convert_image_3xN_to_3D(image_3xN,Ns_θ,Ns_ϕ,Ns_h)

    ## PERFORMANCE METRICS
    # PSF metrics
    if size(t_xyz_grid)[2]==1 # PSF related performance metrics are calculated when there is only one point target
        target_location=[t_θ t_ϕ t_h] # point target location
        resolution,PSLR,ISLR,loc_errors=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_θ,s_ϕ,s_h,PSF_peak_target) # resolutions in each of the 3 axes
    else
        resolution=[NaN,NaN,NaN]
        PSLR=[NaN,NaN,NaN]
        ISLR=[NaN,NaN,NaN]
        println("PSF related performance metrics cannot be calculated since there are more than 1 targets!")
    end
    (peak, idx)             = findmax(image_3D) # finds maximum and index of max
    peaks[ntrial]           = peak  
    peak_idx1[ntrial]       = Int64(idx[1])
    peak_idx2[ntrial]       = Int64(idx[2])
    peak_idx3[ntrial]       = Int64(idx[3])
    resolutions[:,ntrial]   = resolution
    loc_errors[:,ntrial]    = loc_error
    PSLRs[:,ntrial]         = PSLR
    ISLRs[:,ntrial]         = ISLR
end#Ntrials

# convert peak indices to locations
peak_θ = s_θ[peak_idx1]
peak_ϕ = s_θ[peak_idx2]
peak_h = s_h[peak_idx3] 
ideal_peak_θ = s_θ[ideal_peak_idx1]
ideal_peak_ϕ = s_ϕ[ideal_peak_idx2]
ideal_peak_h = s_h[ideal_peak_idx3]

if disable_freq_offset # test to denote frequency error or not in save file name
    freq_text = "noFreq"
else
    freq_text = "wFreq"
end

outputfilename = "syncModule_MonteCarlo_$osc_type"*"_$sync_pri"*"s_"*freq_text* ".jld2" # this is the output filename that the data is saved to using JLD2
# this saves the data into a JLD2 file. Data includes the estimates
@save outputfilename peaks peak_θ peak_ϕ peak_h resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors
#println(std(resolutions[1,:]))
# Note: JLD2 can be read using "@load filename var1 var2...
println("Run Complete")
