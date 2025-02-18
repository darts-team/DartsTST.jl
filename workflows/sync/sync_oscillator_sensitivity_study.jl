
using NCDatasets
using Plots
using Statistics
using JLD2 # note: may have to Pkg.add("JLD2")
using Distributed, SharedArrays
using Parameters

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
    # include("../../inputs/predefined-input-parameters.jl")
    include("../../modules/tomo_workflow_function.jl")
end#begin


a_coeff_dB = [-1000 -1000 -1000 -1000 -1000] # near zero noise oscillator. Used for sweep values

Ntrials = 64 # number of trials per SRI in Monte Carlo simulations

coeff_number = 1 # which coefficient to sweep over (1 through 5)
# osc_coeff_sweep = [-70 -60 -50 -40 -35 -30 -25 -20]
osc_coeff_sweep = [-100 -90 -80 -70 -60 -50 -40 -30]
numCoeffs = length(osc_coeff_sweep)

# generate set of example parameters to get scene size
params = UserParameters.inputParameters(sync_a_coeff_dB = a_coeff_dB,PSF_image_point=1,display_tomograms=0,user_defined_orbit=0,fc = 1e9,mode=2,s_loc_1=-40:2:40)
# Check consistency of input parameters
paramsIsValid = UserParameters.validateInputParams(params)

image_3D = TomoWorkflow.generate_tomo(params)
scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)
ideal_res, ideal_PSLR, ideal_ISLR, loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
(peak, idx) = findmax(image_3D)
ideal_peak  = peak

@unpack  Ns_1, Ns_2, Ns_3 = params

## Initialize result vectors
peaks       = SharedArray{Float64}(numCoeffs,Ntrials)
resolutions = SharedArray{Float64}(3,numCoeffs,Ntrials)
PSLRs       = SharedArray{Float64}(3,numCoeffs,Ntrials)
ISLRs       = SharedArray{Float64}(3,numCoeffs,Ntrials)
loc_errors  = SharedArray{Float64}(3,numCoeffs,Ntrials)
tomo_data   = SharedArray{Float64}(Ns_1,Ns_2,Ns_3,numCoeffs,Ntrials)
## run trials
@sync @distributed for ntrial = 1 : Ntrials

     for k = 1 : numCoeffs
        test_coeff = osc_coeff_sweep[k]
        a_coeff_dB[coeff_number] = test_coeff

        #TODO: do we want to use sync or not? 
        params = UserParameters.inputParameters(sync_a_coeff_dB = a_coeff_dB,PSF_image_point=1,PSF_cuts=1,display_tomograms=0,mode=3,s_loc_1=-40:2:40,
        no_sync_flag = false,enable_sync_phase_error=true,fc = 1e9)
        #take tomographic cut
        # params = UserParameters.inputParameters(sync_a_coeff_dB = a_coeff_dB,PSF_image_point=1,PSF_cuts=1,display_tomograms=0,user_defined_orbit=0,include_antenna=false)
        # Check consistency of input parameters
        paramsIsValid = UserParameters.validateInputParams(params)

        # generate tomo using parameters and updated oscillator coefficients
        image_3D = TomoWorkflow.generate_tomo(params)

        #store 3D image data into shared array
        tomo_data[:,:,:,k,ntrial] = image_3D

        ## PERFORMANCE METRICS
        # Calculate point target performance metrics
        resolution=[NaN,NaN,NaN]
        PSLR=[NaN,NaN,NaN]
        ISLR=[NaN,NaN,NaN]
        loc_error=[NaN,NaN,NaN]
        try
  	   # Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
           scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)
           resolution, PSLR, ISLR, loc_error = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
        catch            
        end

        (peak, idx)             = findmax(image_3D) # finds maximum and index of max
        peaks[k,ntrial]         = peak
        loc_errors[:,k,ntrial]  = loc_error
        resolutions[:,k,ntrial] = resolution
        PSLRs[:,k,ntrial]       = PSLR
        ISLRs[:,k,ntrial]       = ISLR
    end#N SRIs

end#Ntrials

if params.no_sync_flag # if flag == true, no sync is used. flag == false results in normal sync process estimation
    sync_text = "noSync"
else
    sync_text = "wSync"
end

@unpack  s_loc_1, s_loc_2, s_loc_3, mode = params

outputfilename = "syncModule_MonteCarlo_mode_$mode"*"_coeff_number"*"$coeff_number"*"_sync_osc_sweep_"*sync_text* ".jld2" # this is the output filename that the data is saved to using JLD2
# this saves the data into a JLD2 file. Data includes the estimates
@save outputfilename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors osc_coeff_sweep s_loc_1 s_loc_2 s_loc_3

println("Run Complete, and file saved to " *outputfilename)

# ## Second version looking only at tomo axis
# 
# using NCDatasets
# using Plots
# using Statistics
# using JLD2 # note: may have to Pkg.add("JLD2")
# using Distributed, SharedArrays
# using Parameters
# 
# maxprocs = 12 # maximum number of cores to use
# curr_procs = nprocs()
# if curr_procs < maxprocs
#     addprocs(maxprocs - curr_procs)
#     # not adding the else case where we would remove processes. Possible but not useful
# end#if
# curr_procs = nprocs()
# println("Current procs: " * "$curr_procs")
# ## includes
# @everywhere begin
#     include("../../modules/generate_raw_data.jl")
#     include("../../modules/process_raw_data.jl")
#     include("../../modules/geometry.jl")
#     include("../../modules/scene.jl")
#     include("../../modules/sync.jl")
#     include("../../modules/range_spread_function.jl") # as RSF
#     include("../../modules/orbits.jl")
#     include("../../modules/error_sources.jl")
#     include("../../modules/performance_metrics.jl")
#     include("../../modules/antenna.jl")
#     include("../../modules/simsetup.jl")
#     include("../../modules/user_parameters.jl")
#     # include("../../inputs/predefined-input-parameters.jl")
#     include("../../modules/tomo_workflow_function.jl")
# end#begin
# 
# 
# a_coeff_dB = [-500 -500 -500 -500 -500] # near zero noise oscillator. Used for sweep values
# 
# Ntrials = 64 # number of trials per SRI in Monte Carlo simulations
# 
# coeff_number = 1 # which coefficient to sweep over (1 through 5)
# osc_coeff_sweep = [-100 -90 -80 -70 -60 -50 -40 -30]
# numCoeffs = length(osc_coeff_sweep)
# 
# # generate set of example parameters to get scene size
# params = UserParameters.inputParameters(sync_a_coeff_dB = a_coeff_dB,PSF_image_point=1,display_tomograms=0)
# # Check consistency of input parameters
# paramsIsValid = UserParameters.validateInputParams(params)
# 
# image_3D = TomoWorkflow.generate_tomo(params)
# # Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
# scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)
# ideal_res, ideal_PSLR, ideal_ISLR, loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
# (peak, idx) = findmax(image_3D)
# ideal_peak  = peak
# 
# @unpack  Ns_1, Ns_2, Ns_3 = params
# 
# ## Initialize result vectors
# peaks       = SharedArray{Float64}(numCoeffs,Ntrials)
# resolutions = SharedArray{Float64}(3,numCoeffs,Ntrials)
# PSLRs       = SharedArray{Float64}(numCoeffs,Ntrials)
# ISLRs       = SharedArray{Float64}(numCoeffs,Ntrials)
# loc_errors  = SharedArray{Float64}(3,numCoeffs,Ntrials)
# tomo_data   = SharedArray{Float64}(Ns_1,Ns_2,Ns_3,numCoeffs,Ntrials)
# ## run trials
# @sync @distributed for ntrial = 1 : Ntrials
# 
#      for k = 1 : numCoeffs
#         test_coeff = osc_coeff_sweep[k]
#         a_coeff_dB[coeff_number] = test_coeff
# 
#         #TODO: do we want to use sync or not? 
#         params = UserParameters.inputParameters(sync_a_coeff_dB = a_coeff_dB,PSF_image_point=1,display_tomograms=0, no_sync_flag = true, enable_sync_phase_error=true)
#         #take tomographic cut
#         # params = UserParameters.inputParameters(sync_a_coeff_dB = a_coeff_dB,PSF_image_point=1,PSF_cuts=1,display_tomograms=0,user_defined_orbit=0,include_antenna=false)
#         # Check consistency of input parameters
#         paramsIsValid = UserParameters.validateInputParams(params)
# 
#         # generate tomo using parameters and updated oscillator coefficients
#         image_3D = TomoWorkflow.generate_tomo(params)
# 
#         #store 3D image data into shared array
#         tomo_data[:,:,:,k,ntrial] = image_3D
# 
#         ## PERFORMANCE METRICS
#         # Calculate point target performance metrics
#         resolution=[NaN]
#         PSLR=[NaN]
#         ISLR=[NaN]
#         loc_error=[NaN]
#         try
#   	   # Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
#            scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)
#            resolution, PSLR, ISLR, loc_error = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
#         catch            
#         end
# 
#         (peak, idx)           = findmax(image_3D) # finds maximum and index of max
#         peaks[k,ntrial]       = peak
#         # loc_errors[1,k,ntrial]  = loc_error
#         # resolutions[1,k,ntrial] = resolution
#         PSLRs[k,ntrial]       = PSLR
#         ISLRs[k,ntrial]       = ISLR
#     end#N SRIs
# 
# end#Ntrials
# 
# if params.no_sync_flag # if flag == true, no sync is used. flag == false results in normal sync process estimation
#     sync_text = "noSync"
# else
#     sync_text = "wSync"
# end
# 
# @unpack  s_loc_1, s_loc_2, s_loc_3, mode = params
# 
# outputfilename = "syncModule_MonteCarlo_mode_$mode"*"_coeff_number"*"$coeff_number"*"_sync_osc_sweep_"*sync_text* "_along_n.jld2" # this is the output filename that the data is saved to using JLD2
# # this saves the data into a JLD2 file. Data includes the estimates
# @save outputfilename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors osc_coeff_sweep s_loc_1 s_loc_2 s_loc_3
# 
# println("Run Complete, and file saved to " *outputfilename)
