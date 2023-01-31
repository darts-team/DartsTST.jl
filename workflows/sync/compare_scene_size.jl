# compare_scene_size.jl

# Examine difference in SRI study results for different scene size
# observations: increased scene extend led to significant increase in peak loss, greater PSLR, worse performance generally
using Plots
using Statistics
using Dates
using JLD2
using Distributed, SharedArrays
using Parameters
using StaticArrays
using StatsPlots, PyCall, CurveFit, MAT


# Compare output files from Monte Carlo sim using RFSOC + GPSDO

# filename3 =  "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_3_Measured_sync_pri_sweep_along_n_expandedScene.jld2"
#   @load filename3 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3

filename_reg = "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_pri_sweep_imageData_along_n.jld2"

filename_exp = "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_pri_sweep_imageData_along_n_expandedScene.jld2"

@load filename_exp ideal_image_3D tomo_data sync_PRIs s_loc_1 s_loc_2 s_loc_3

ideal_image_3D_exp = ideal_image_3D; tomo_data_exp = tomo_data;
s_loc_1_exp = s_loc_1; s_loc_2_exp = s_loc_2; s_loc_3_exp = s_loc_3

@load filename_reg ideal_image_3D tomo_data sync_PRIs s_loc_1 s_loc_2 s_loc_3

println(size(ideal_image_3D))
println(size(ideal_image_3D_exp))

## setting up some parameters so we can use plotting.jl to plot tomograms
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
include("../../modules/plotting.jl")

# params for regular scene size
params = UserParameters.inputParameters(s_loc_1=s_loc_1,PSF_image_point=1,PSF_cuts=2,display_tomograms=1,mode=2)
# Check consistency of input parameters
paramsIsValid = UserParameters.validateInputParams(params)

# params for regular scene size
params_exp = UserParameters.inputParameters(s_loc_1=s_loc_1_exp,s_loc_2=s_loc_2_exp,s_loc_3=s_loc_3_exp,PSF_image_point=1,PSF_cuts=2,display_tomograms=1,mode=2)
# Check consistency of input parameters
paramsIsValid = UserParameters.validateInputParams(params_exp)

ts_coord_txt=Plotting.coordinates(params.ts_coord_sys)

#select index of saved tomo data
SRIs = sync_PRIs
sri_num=6; trial_num=2;
image_3D = tomo_data[:,:,:,sri_num,trial_num]
image_3D_exp = tomo_data_exp[:,:,:,sri_num,trial_num]


# plotting ideal tomograms first
scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(ideal_image_3D, params)
Plotting.plot_tomogram(ideal_image_3D, ts_coord_txt, scene_axis11, scene_axis22, scene_axis33, params)

scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(ideal_image_3D_exp, params_exp)
Plotting.plot_tomogram(ideal_image_3D_exp, ts_coord_txt, scene_axis11, scene_axis22, scene_axis33, params_exp)


# plotting tomograms at peak location
scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)
Plotting.plot_tomogram(image_3D, ts_coord_txt, scene_axis11, scene_axis22, scene_axis33, params)

scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D_exp, params_exp)
Plotting.plot_tomogram(image_3D_exp, ts_coord_txt, scene_axis11, scene_axis22, scene_axis33, params_exp)

brightest=maximum(ideal_image_3D)
brightest_exp=maximum(ideal_image_3D_exp)
if brightest != brightest_exp
    error("uh oh")
end

peak = maximum(image_3D)
peak_exp = maximum(image_3D_exp)
println("Ideal peak: $brightest")
SRI = SRIs[sri_num]
println("SRI: $SRI")
println("Peak: $peak")
println("Peak Expanded Scene: $peak_exp")
