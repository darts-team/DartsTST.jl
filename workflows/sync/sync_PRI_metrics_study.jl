
using NCDatasets
using Plots
using Statistics
using Dates
using JLD2 # note: may have to Pkg.add("JLD2")
using Distributed, SharedArrays
using Parameters
using StaticArrays
using .UserParameters

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
    include("../../modules/antenna.jl")
    include("../../modules/simsetup.jl")
    include("../../modules/user_parameters.jl")
end#begin

# Define user parameters
params = UserParameters.inputParameters()

@unpack pulse_length, ts_coord_sys, display_geometry, display_RSF_rawdata, processing_steps,
    display_input_scene, display_tomograms, display_geometry_coord = params

#----- Overwriting default parameter values-----
params.SAR_duration=5 #  set synthetic aperture duration (s)
params.s_loc_1 = -60:1:60 #(m)
params.PSF_image_point = 1 #peak location
params.sync_osc_type = "MicroSemi"
params.display_tomograms = 0 # do not display
# sync_osc_type = "USO" # putting a oscillator type variable here to auto-name save files


# We'll leave this if-else structure here because it's convenient for switching the sync_osc_type. However we will pass the variables into the params struct
#defines oscillator quality. Either leave as single row to use across all platforms, or define values for each platform as a new row
if params.sync_osc_type == "USO"
    params.sync_a_coeff_dB = [-95 -90 -200 -130 -155] # [USO: Krieger]
elseif params.sync_osc_type == "USRP"
    params.sync_a_coeff_dB = [-66 -62 -80 -110 -153] # [USRP E312]
elseif params.sync_osc_type == "Wenzel5MHz"
    params.sync_a_coeff_dB = [-1000 -128 -1000 -150 -178] # [Wenzel 5MHz oscillator] - NOTE: fractional dB values were rounded up for Wenzel oscillators (to keep as Int64 values)
    params.sync_f_osc = 5e6 # local oscillator frequency
elseif params.sync_osc_type == "Wenzel100MHz"
    params.sync_a_coeff_dB = [-1000 -73 -1000 -104 -181] # [Wenzel 100MHz oscillator]
    params.sync_f_osc = 100e6 # local oscillator frequency
elseif params.sync_osc_type == "MicroSemi"
    params.sync_a_coeff_dB = [-120 -114 -999 -134 -166 ] # [Microsemi GPS-3500 oscillator]
    #TODO(correct center freq = 10MHz?)
end

Ntrials = 64 # number of trials per SRI in Monte Carlo simulations
sync_PRIs = [.1 1 2 3 4 5]
numSRI = length(sync_PRIs)

# Compute orbits time, position, and velocity
orbit_time, orbit_pos, orbit_vel = Orbits.computeTimePosVel(params)

# interpolate orbit to slow time, 3 x Np x Nst, convert km to m
p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)

# Create target/scene location 
targets_loc, targets_ref, Nt = Scene.construct_targets_str(params) # Nt: number of targets, targets: structure array containing target locations and reflectivities
s_loc_3xN  = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops
t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions

# Read number of platforms (todo: move into a struct)
const Np  = size(orbit_pos)[2] # number of platforms
sync_osc_coeffs = repeat(params.sync_a_coeff_dB, Np)

# Apply antenna pattern
if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
    Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
end

# Generate range spread function (matched filter output)
min_range, max_range = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
Trx = 2*(max_range-min_range)/c + 5*pulse_length # s duration of RX window
Srx, MF, ft, t_rx = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
# Srx,MF,ft,t_rx=RSF.non_ideal_RSF(pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

# Generate TomoSAR raw data
ref_range = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
rawdata   = Generate_Raw_Data.main_RSF_slowtime(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
if params.enable_thermal_noise # adding random noise based on SNR after range (fast-time) processing
    rawdata = Error_Sources.random_noise(rawdata, params)
end 

## find Ideal case results first
# Process raw data to generate image
if processing_steps === :bp3d # 1-step processing
    image_3D = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
elseif processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
    SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
    image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
end

ideal_image_3D = image_3D # for saving later

# performance metrics for ideal case
(ideal_peak, ideal_idx) = findmax(image_3D) # finds maximum and index of max
if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
    ideal_res, ideal_PSLR, ideal_ISLR, loc_error, scene_axis11, scene_axis22, scene_axis33, PSF_metrics = Performance_Metrics.computePTPerformanceMetrics(image_3D, params)
else
    error("PSF related performance metrics cannot be calculated for more than 1 target.")
end

## Initialize result vectors
peaks       = SharedArray{Float64}(numSRI,Ntrials)
resolutions = SharedArray{Float64}(3,numSRI,Ntrials)
PSLRs       = SharedArray{Float64}(3,numSRI,Ntrials)
ISLRs       = SharedArray{Float64}(3,numSRI,Ntrials)
loc_errors  = SharedArray{Float64}(3,numSRI,Ntrials)
tomo_data   = SharedArray{Float64}(Ns_1,Ns_2,Ns_3,numSRI,Ntrials)

## run trials
@sync @distributed for ntrial = 1 : Ntrials
     for k = 1 : numSRI
        sync_pri = sync_PRIs[k]
        parameters.sync_pri = sync_pri
        println("Starting SRI value: ", sync_pri)

        ## add in error sources
        rawdata_sync = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, sync_osc_coeffs, params)

        # Process raw data to generate image
        if processing_steps === :bp3d # 1-step processing
            image_3D = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
        elseif processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
            SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
            image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
        end

        #store 3D image data into shared array
        tomo_data[:,:,:,k,ntrial] = image_3D

        ## PERFORMANCE METRICS
        try
            resolution, PSLR, ISLR, loc_error, scene_axis11, scene_axis22, scene_axis33, PSF_metrics = Performance_Metrics.computePTPerformanceMetrics(image_3D, params)
        catch
            resolution=[NaN,NaN,NaN]
            PSLR=[NaN,NaN,NaN]
            ISLR=[NaN,NaN,NaN]
            PSF_metrics=false
            loc_error=[NaN,NaN,NaN]
            show("PSF related performance metrics cannot be calculated -- error in metric calculation.")
        end#try
    
        # store metric data into SharedArrays
        (peak, idx)             = findmax(image_3D) # finds maximum and index of max
        peaks[k,ntrial]         = peak
        loc_errors[:,k,ntrial]  = loc_error
        resolutions[:,k,ntrial] = resolution
        PSLRs[:,k,ntrial]       = PSLR
        ISLRs[:,k,ntrial]       = ISLR
    end#N SRIs
end#Ntrials

outputfilename = "syncModule_MonteCarlo_mode_$mode"*"_$sync_osc_type"*"_sync_pri_sweep.jld2" # this is the output filename that the data is saved to using JLD2
# this saves the data into a JLD2 file. Data includes the estimates
@save outputfilename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3

outputfilename_data = "syncModule_MonteCarlo_mode_$mode"*"_$sync_osc_type"*"_sync_pri_sweep_imageData.jld2" # output filename for image data. doesn't save metrics
@save outputfilename_data ideal_image_3D tomo_data sync_PRIs s_loc_1 s_loc_2 s_loc_3

println("Run Complete, and file saved to " *outputfilename)
# Note: JLD2 can be read using "@load filename var1 var2...