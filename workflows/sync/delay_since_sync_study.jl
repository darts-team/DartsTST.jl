
using NCDatasets
using Plots
using Statistics
using Dates
using JLD2
using Distributed, SharedArrays
using Parameters
using StaticArrays
# using UserParameters

maxprocs = 12 # maximum number of cores to use
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
    include("../../inputs/predefined-input-parameters.jl")
    include("../../modules/tomo_workflow_function.jl")
end#begin

#----- Setting parameter values to overwrite defaults-----
sync_osc_type = "Measured"
radar_mode=1
sar_len = 3.0
at_dim = -20:2:20
usr_orbit = 1
use_meas_flag = false # to be overwritten if need be
center_freq = 1e9

# We'll leave this if-else structure here because it's convenient for switching the sync_osc_type. However we will pass the variables into the params struct
#defines oscillator quality. Either leave as single row to use across all platforms, or define values for each platform as a new row
if sync_osc_type == "USO"
    coeffs = [-95 -90 -200 -130 -155] # [USO: Krieger]
elseif sync_osc_type == "USRP"
    coeffs = [-66 -62 -80 -110 -153] # [USRP E312]
elseif sync_osc_type == "Wenzel5MHz"
    coeffs = [-1000 -128 -1000 -150 -178] # [Wenzel 5MHz oscillator] - NOTE: fractional dB values were rounded up for Wenzel oscillators (to keep as Int64 values)
    sync_f_osc = 5e6 # local oscillator frequency
elseif sync_osc_type == "Wenzel100MHz"
    coeffs = [-1000 -73 -1000 -104 -181] # [Wenzel 100MHz oscillator]
    sync_f_osc = 100e6 # local oscillator frequency
elseif sync_osc_type == "MicroSemi"
    coeffs = [-120 -114 -999 -134 -166 ] # [Microsemi GPS-3500 oscillator]
elseif sync_osc_type == "RFSoc"
    coeffs = [-120 -114 -999 -134 -166 ] # [Very rough estimate of measured RFSoC oscillators]
elseif sync_osc_type == "Measured"
    use_meas_flag = true
end

Ntrials = 64 # number of trials per delay in Monte Carlo simulations
delays = [0 .1 .2 .5 1 2 3 4 5 6 7 8 9 10 50 100 ]
numDelays = length(delays)

## find Ideal case results first
params = UserParameters.inputParameters(s_loc_1=at_dim,PSF_image_point=1,PSF_cuts=1,display_tomograms=0,
user_defined_orbit=usr_orbit,SAR_duration=sar_len,mode=radar_mode, fc = center_freq)

# Check consistency of input parameters
paramsIsValid = UserParameters.validateInputParams(params)

# generate tomo using parameters and updated oscillator coefficients
ideal_image_3D = TomoWorkflow.generate_tomo(params)

# Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(ideal_image_3D, params)

(peak, idx) = findmax(ideal_image_3D)
ideal_peak  = peak
# Calculate point target performance metrics
ideal_res, ideal_PSLR, ideal_ISLR, loc_error  = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)

## Initialize result vectors
peaks       = SharedArray{Float64}(numDelays,Ntrials)
resolutions = SharedArray{Float64}(3,numDelays,Ntrials)
PSLRs       = SharedArray{Float64}(3,numDelays,Ntrials)
ISLRs       = SharedArray{Float64}(3,numDelays,Ntrials)
loc_errors  = SharedArray{Float64}(3,numDelays,Ntrials)
tomo_data   = SharedArray{Float64}(params.Ns_1,params.Ns_2,params.Ns_3,numDelays,Ntrials)

## run trials
@sync @distributed for ntrial = 1 : Ntrials
     for k = 1 : numDelays
        
            delay = delays[k]
            println("Starting Delay Since Sync value: ", delay)
            
            if use_meas_flag
                params = UserParameters.inputParameters(PSF_image_point=1, PSF_cuts=1, display_tomograms=0, user_defined_orbit=usr_orbit, s_loc_1 = at_dim, use_measured_psd_flag=use_meas_flag, mode=radar_mode,
                enable_sync_phase_error=true, delay_since_sync=delay, no_sync_flag=true, SAR_duration=sar_len, fc = center_freq)
            else
                params = UserParameters.inputParameters(PSF_image_point=1, PSF_cuts=1, display_tomograms=0, user_defined_orbit=usr_orbit, s_loc_1 = at_dim, sync_a_coeff_dB = coeffs,mode=radar_mode,
                enable_sync_phase_error=true, delay_since_sync=delay, no_sync_flag=true, SAR_duration=sar_len, fc = center_freq)
            end

        image_3D = TomoWorkflow.generate_tomo(params)
        #store 3D image data into shared array
        tomo_data[:,:,:,k,ntrial] = image_3D

        ## PERFORMANCE METRICS
        # try
        resolution=[NaN,NaN,NaN]
        PSLR=[NaN,NaN,NaN]
        ISLR=[NaN,NaN,NaN]
        PSF_metrics=false
        loc_error=[NaN,NaN,NaN]
        scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)
        # Calculate point target performance metrics
        try
            resolution, PSLR, ISLR, loc_error  = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
        catch
        end
        

        # store metric data into SharedArrays
        (peak, idx)             = findmax(image_3D) # finds maximum and index of max
        peaks[k,ntrial]         = peak
        loc_errors[:,k,ntrial]  = loc_error
        resolutions[:,k,ntrial] = resolution
        PSLRs[:,k,ntrial]       = PSLR
        ISLRs[:,k,ntrial]       = ISLR
    end#N delays
end#Ntrials
@unpack mode, s_loc_1, s_loc_2, s_loc_3 = params
outputfilename = "syncModule_MonteCarlo_mode_$mode"*"_$sync_osc_type"*"_sync_DSS_sweep.jld2" # this is the output filename that the data is saved to using JLD2
# this saves the data into a JLD2 file. Data includes the estimates
@save outputfilename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors delays s_loc_1 s_loc_2 s_loc_3 center_freq

outputfilename_data = "syncModule_MonteCarlo_mode_$mode"*"_$sync_osc_type"*"_sync_DSS_sweep_imageData.jld2" # output filename for image data. doesn't save metrics
@save outputfilename_data ideal_image_3D tomo_data delays s_loc_1 s_loc_2 s_loc_3

println("Run Complete, and file saved to " *outputfilename)
