include("../modules/generate_raw_data.jl")
include("../modules/process_raw_data.jl")
include("../modules/geometry.jl")
include("../modules/scene.jl")
include("../modules/range_spread_function.jl") # as RSF
include("../modules/orbits.jl")
include("../modules/sync.jl")
include("../modules/error_sources.jl")
include("../modules/performance_metrics.jl")
include("../modules/antenna.jl")
include("../modules/simsetup.jl")
include("../modules/user_parameters.jl")
using NCDatasets
using Statistics
using Parameters
using Dates
using StaticArrays
using .UserParameters
c = 299792458 #TODO does not work without redefining c here

# Define user parameters
#include("../inputs/predefined-input-parameters.jl") TODO gives errors
# params = UserParameters.inputParameters()

# params = UserParameters.inputParameters(PSF_image_point=1,PSF_cuts=1,display_tomograms=1,enable_sync_phase_error=true,use_measured_psd_flag=true,no_sync_flag=false,mode=2,sync_f_osc=12.8e6,
# sync_pri=.1,phase_offset_flag=true,SAR_duration=5,sync_a_coeff_dB=[-120 -114 -999 -134 -166 ],user_defined_orbit=1,osc_psd_meas_filename="inputs/12_8MHz no ref 7M 220321_1531.xlsx")


 # params = UserParameters.inputParameters(PSF_image_point=1,PSF_cuts=2,display_tomograms=1,mode=2,
 # no_sync_flag = true,SAR_duration=3,enable_sync_phase_error=true,sync_pri=1,phase_offset_flag=false)

params = UserParameters.inputParameters(PSF_image_point=1,PSF_cuts=2,display_tomograms=1,enable_sync_phase_error=true,use_measured_psd_flag=true,mode=3,
no_sync_flag=true,osc_psd_meas_filename = "inputs/roseL_osc_specs.jld2",phase_offset_flag=false,sync_a_coeff_dB = [-95 -90 -200 -130 -155],
user_defined_orbit=1)

# params = UserParameters.inputParameters(PSF_image_point=1,PSF_cuts=2,display_tomograms=1,SAR_duration=3,user_defined_orbit=1,mode=2)

# Check consistency of input parameters
paramsIsValid = UserParameters.validateInputParams(params)

# Compute orbits time, position, and velocity
orbit_time, orbit_pos, orbit_vel = Orbits.computeTimePosVel(params)

# interpolate orbit to slow time, 3 x Np x Nst, convert km to m
const p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)

# Create target/scene location
targets_loc, targets_ref, Nt = Scene.construct_targets_str(params) # Nt: number of targets, targets: structure array containing target locations and reflectivities
const s_loc_3xN  = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops
t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions
#t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, orbit_vel, params) ## calculate avg heading from platform positions/velocities

# Read number of platforms (todo: move into a struct)
const Np  = size(orbit_pos)[2] # number of platforms

# Apply antenna pattern
if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
    Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
end

# Generate range spread function (matched filter output)
const min_range, max_range = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
const Trx = 2*(max_range-min_range)/c + 5*params.pulse_length # s duration of RX window
const Srx, MF, ft, t_rx = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
# Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

# Generate TomoSAR raw data
const ref_range = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
const rawdata = Generate_Raw_Data.main_RSF_slowtime(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
if params.enable_thermal_noise # adding random noise based on SNR after range (fast-time) processing
    const rawdata = Error_Sources.random_noise(rawdata, params)
end

# Add phase error
const sync_osc_coeffs = repeat(params.sync_a_coeff_dB, Np)
if params.enable_sync_phase_error
    const rawdata = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, t_xyz_3xN, sync_osc_coeffs, params)
end

# Process raw data to generate image
if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
    const image_3D = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
    const SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
    const image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
end

# Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
if params.left_right_look == "left";PSFcutdir=-1;elseif params.left_right_look == "right";PSFcutdir=1;end
params.PSF_direction[3]=PSFcutdir*params.PSF_direction[3]
scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)

# Calculate point target performance metrics
if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
    const resolutions, PSLRs, ISLRs, loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
    println("Resolutions: ",round.(resolutions,digits=8)," in scene axes units")
    println("Location Errors: ",round.(loc_errors,digits=8)," in scene axes units")
    println("PSLRs: ",round.(PSLRs,digits=2)," dB")
    println("ISLRs: ",round.(ISLRs,digits=2)," dB")
    println("PSF Peak Amplitude: ",round(maximum(20*log10.(image_3D)),digits=2)," dB")
else
    @warn "PSF related performance metrics cannot be calculated for more than 1 target."
end

# Relative Radiometric Accuracy (amplitude difference between input 3D scene and output 3D image, max normalized to 1)
const inputscene_3D = Scene.generate_input_scene_3D(targets_ref, Nt, params)
const diff_image3D, mean_diff_image, std_diff_image = Performance_Metrics.relative_radiometric_accuracy(inputscene_3D, image_3D)
println("Relative Radiometric Accuracy: Mean: ", round(mean_diff_image, digits=2),", Std: ",round(std_diff_image, digits=2)) # mean=0 & std_dev=0 means perfect result


# Plots (1D PSF cuts are displayed by default in the performance.metrics module)
if params.display_geometry || params.display_RSF_rawdata || params.display_input_scene || params.display_tomograms != 0
    include("../modules/plotting.jl")
    const display_geometry_coord_txt=Plotting.coordinates(params.display_geometry_coord)
    const ts_coord_txt=Plotting.coordinates(params.ts_coord_sys)
    if params.display_RSF_rawdata; Plotting.plot_RSF_rawdata(ft, t_rx, MF, Srx, Np, Nst, rawdata, params); end
    if params.display_geometry
        # convert platform and target locations to desired coordinate system
        p_loc,t_loc,s_loc=Geometry.convert_platform_target_scene_coordinates(Np,Nst,Nt,p_xyz,t_xyz_3xN,targets_loc,s_xyz_3xN,s_loc_3xN,avg_peg, params)
        Plotting.plot_geometry(orbit_time,orbit_pos,p_loc,t_loc,s_loc,display_geometry_coord_txt)
    end
    if params.display_input_scene; Plotting.plot_input_scene(inputscene_3D, ts_coord_txt, params);end
    if params.display_tomograms != 0; Plotting.plot_tomogram(image_3D, ts_coord_txt, scene_axis11, scene_axis22, scene_axis33, params);end
    if params.display_input_scene; Plotting.plot_input_scene(diff_image3D, ts_coord_txt, params);end
end
