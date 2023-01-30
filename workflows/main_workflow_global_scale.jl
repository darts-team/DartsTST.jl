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
include("../modules/global_scale_support.jl")
using NCDatasets
using Statistics
using Parameters
using Dates
using StaticArrays
using .UserParameters
using SharedArrays
c = 299792458 #TODO does not work without redefining c here

#Read canopy heights
filepath_GEDIL3 = "/Users/joshil/Documents/GEDI_Data/GEDI_L3_LandSurface_Metrics_V2_1952/data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
grid_res        = 100;
Canopy_heights, Geo_location, size_row, size_col = Global_Scale_Support.read_GEDI_L3_data(filepath_GEDIL3, grid_res)

# Define outout variables
global Output_stat          = SharedArray(zeros(size_row,size_col,4))
global Orbit_index          = SharedArray(zeros(size_row,size_col,1))
global lookang_all          = SharedArray(zeros(size_row,size_col,1))
#Norm baselines
global Norm_baseline_max    = SharedArray(zeros(size_row,size_col,1))
global Norm_baseline_min    = SharedArray(zeros(size_row,size_col,1))
global Norm_baseline_mean   = SharedArray(zeros(size_row,size_col,1))
# perp baselines
global Perp_baseline_max    = SharedArray(zeros(size_row,size_col,1))
global Perp_baseline_min    = SharedArray(zeros(size_row,size_col,1))
global Perp_baseline_mean   = SharedArray(zeros(size_row,size_col,1))
# AT baselines
global Par_baseline_max     = SharedArray(zeros(size_row,size_col,1))
global Par_baseline_min     = SharedArray(zeros(size_row,size_col,1))
global Par_baseline_mean    = SharedArray(zeros(size_row,size_col,1))

#region_xlims = [50,110]
#region_ylims = [15,55]
#region_xlims = 59:62
#region_ylims = 17:19
region_xlims        = 80:81
region_ylims        = 37:37

lat_lon_idx         = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

orbit_dataset       = Dataset("inputs/orbitOutput_082020.nc") # Read orbits data in NetCDF format
global mast_plat            = 1
flag_plat           = 1 #descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, mast_plat, flag_plat)

for i1 = 1:size(lat_lon_idx,1)   
    
    if isnan(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        #Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =0
        global Output_stat[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Output_stat[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = NaN
        global Output_stat[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = NaN
        global Output_stat[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = NaN
        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Norm_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Perp_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Par_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
    
    else

                
        global close_val_lat_lon   = Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = close_val_lat_lon[2]

        slrng_temp2 = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]]; 0]))
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = acosd(orbit_pos_geo_all[:,close_val_lat_lon[2]][3] / slrng_temp2) 

        ref_profile_height, ref_profile_value, NoPeaks  = Global_Scale_Support.constrct_reflectivity_profile_exp(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])

        # Define user parameters
        #include("../inputs/predefined-input-parameters.jl") TODO gives errors
        params = UserParameters.inputParameters(
            target_pos_mode = "layered-grid",
            t_loc_3 = ref_profile_height,
            t_ref = ref_profile_value,
        )

        # Check consistency of input parameters
        paramsIsValid = UserParameters.validateInputParams(params)

        # Compute orbits time, position, and velocity
        orbit_time = orbit_time_all[close_val_lat_lon[2]-5:close_val_lat_lon[2]+5]
        orbit_pos = orbit_pos_all[:,:,close_val_lat_lon[2]-5:close_val_lat_lon[2]+5]
        orbit_vel = orbit_vel_all[:,:,close_val_lat_lon[2]-5:close_val_lat_lon[2]+5]

        SAR_start_time = orbit_time_all[close_val_lat_lon[2]] - (params.SAR_duration / 2)

        ref_plat = 1 #incicate the reference platform
        bperp, b_at, bnorm = Orbits.get_perp_baselines(orbit_pos, orbit_vel, lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], ref_plat)

        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bnorm) ./ 1e3
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bnorm)) ./ 1e3
        global Norm_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,bnorm)) ./ 1e3
        
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bperp) ./ 1e3
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bperp)) ./ 1e3
        global Perp_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,bperp)) ./ 1e3
        
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(b_at) ./ 1e3
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,b_at)) ./ 1e3
        global Par_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,b_at)) ./ 1e3


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
            const rawdata = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, sync_osc_coeffs, params)
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


end