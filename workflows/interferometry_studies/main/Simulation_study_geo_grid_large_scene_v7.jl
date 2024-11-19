
include("../../../modules/generate_raw_data.jl")
include("../../../modules/process_raw_data.jl")
include("../../../modules/geometry.jl")
include("../../../modules/scene.jl")
include("../../../modules/range_spread_function.jl") 
include("../../../modules/orbits.jl")
include("../../../modules/antenna.jl")
include("../../../modules/dem.jl")
include("../../../modules/simsetup.jl")
include("../../../modules/user_parameters.jl")
include("../../../modules/scattering.jl")
include("../../../modules/interferometry.jl")
include("../../../modules/export_output.jl")

using Statistics
using Parameters
using .UserParameters
using TimerOutputs
using JLD2
using SharedArrays
using Interpolations
using DSP

#Global constants
to                          = TimerOutput()
global c                    = 299792458 # Speed of light
global earth_radius         = 6378.137e3 # Earth semi-major axis at equator1

@timeit to "Overall time" begin

##-------------------------------------------------------------------------
#Simulation setup parameters
Sim_idx                 = 1067   # For output file reference
savepath                = "/u/intrepid-z0/joshil/Outputs/TST_sims_geogrid_1/"*string(Sim_idx)*"/"
#savepath2               = "/u/intrepid-z0/joshil/Outputs/TST_sims_geogrid_1/Outputs_sims_all/"

#Target_scene_setup based on scene and look angle
lat_res                 = 0.000033*2
lon_res                 = 0.000033*1.5
lat_lims                = [35.030 35.031] # [35.030 35.042] #[34.431 34.443]
lon_lims                =  [-111.9 -111.68] #[-111.9 -110.6] # #[-115.2450 -115.2300]
#lon_lims                = [-111.685 -111.670] # #[-115.2450 -115.2300]

lat_extent              = abs.(lat_lims[2] - lat_lims[1])
lon_extent              = abs.(lon_lims[2] - lon_lims[1])
NB_lat                  = 1
NB_lon                  = 1
target_mode             = 2     # 1: target fixed in center, 2: Distributed target, 3: Distributed target with 1 dominant scatterer
num_targ_vol            = 10    # number of targets in each voxel
ref_scene_height        = 0.0 #550.0 #0.0 #450.0 #2000.0

platform_heading        = 0.0
platform_look_dir       = "right"
platform_look_angle     = 30.0 
platform_height         = 697.5e3

# Setup scene based on reference target location
trg_ref_point           = [lat_lims[1];lon_lims[1];0.0]
platform_ref_point      = SimSetup.get_platform_ref_from_trg_ref(trg_ref_point, platform_height, platform_heading, platform_look_angle, platform_look_dir, earth_radius)
trg_ref_lat_list, trg_ref_lon_list, B_lat_extent, B_lon_extent = SimSetup.segment_simulation_grid(lat_lims[1] , lon_lims[1], lat_extent, lon_extent, NB_lat, NB_lon)

## Setup scene based on reference platform location
#platform_ref_point      = [lat_lims[1];lon_lims[1];platform_height]
#trg_ref_point           = SimSetup.get_trg_ref_from_platform_ref(platform_ref_point, trg_height, platform_heading, platform_look_angle, platform_look_dir, earth_radius) 
#trg_ref_lat_list, trg_ref_lon_list, B_lat_extent, B_lon_extent = SimSetup.segment_simulation_grid(trg_ref_point[1] , trg_ref_point[1], lat_extent, lon_extent, NB_lat, NB_lon)

# Validation
trg_ref_point_xyz       = Geometry.geo_to_xyz(trg_ref_point)
platform_ref_point_xyz  = Geometry.geo_to_xyz(platform_ref_point)
slrng_temp              = Geometry.distance(trg_ref_point_xyz,platform_ref_point_xyz)
look_angle_test         = Scene.slantrange_to_lookangle(earth_radius,slrng_temp,platform_ref_point[3],trg_ref_point[3])[2]

##-------------------------------------------------------------------------

for B_idx = 1:length(trg_ref_lat_list)  
        
    @timeit to "Block time " begin
        
        @timeit to "Get DEM and slopes " begin
        #DEM_full, Geo_location_lat, Geo_location_lon, ag_geotransform, ag_ref = DEM.read_interp_DEM_from_source("/u/intrepid-z0/joshil/data/nisar-dem-copernicus/EPSG4326.vrt", trg_ref_lat_list[B_idx], trg_ref_lon_list[B_idx], B_lat_extent, B_lon_extent, lat_res, lon_res, 1)
        DEM_full, Geo_location_lat, Geo_location_lon, ag_geotransform, ag_ref = DEM.create_slope_DEM(trg_ref_lat_list[B_idx], trg_ref_lon_list[B_idx], B_lat_extent, B_lon_extent, lat_res, lon_res, 0.0, 60.0) # slope dem for testing

        #DEM_full, Geo_location_lat, Geo_location_lon, ag_geotransform, ag_ref = DEM.custom_DEM3(trg_ref_lat_list[B_idx], trg_ref_lon_list[B_idx], B_lat_extent, B_lon_extent, lat_res, lon_res)
        #DEM_full, Geo_location_lat, Geo_location_lon, ag_geotransform, ag_ref = DEM.custom_DEM4(trg_ref_lat_list[B_idx], trg_ref_lon_list[B_idx], B_lat_extent, B_lon_extent, lat_res, lon_res)

        slope_lat, slope_lon    = DEM.get_slopes_from_DEM(DEM_full, Geo_location_lon, Geo_location_lat)
        DEM_region              = DEM_full[2:end-1,2:end-1]'
        slope_lat_region        = slope_lat[2:end-1,2:end-1]' 
        slope_lon_region        = slope_lon[2:end-1,2:end-1]'
        N_all                   = DEM.get_terrain_norm(DEM_region, slope_lat_region, slope_lon_region)
        end

        t_loc_1_range           = Geo_location_lat[1,2:end-1]   # 1-grid range
        t_loc_2_range           = Geo_location_lon[2:end-1,1]   # 2-grid range
        t_loc_3_range           = 0.0   # 3-grid range

        #t_loc_3xN               = SimSetup.define_target_pixels(t_loc_1_range, t_loc_2_range, t_loc_3_range, DEM_region, num_targ_vol, target_mode, lat_res/2, lon_res/2, 0.5)
        t_loc_3xN               = SimSetup.define_target_pixels(t_loc_1_range, t_loc_2_range, t_loc_3_range, DEM_region, slope_lat_region, slope_lon_region, num_targ_vol, target_mode, lat_res/2, lon_res/2, 0.0)

        t_xyz_3xN               = Geometry.geo_to_xyz(t_loc_3xN,earth_radius,earth_eccentricity)

        scene1_loc_3xN          = SimSetup.define_scene_pixels(t_loc_1_range, t_loc_2_range, 0, DEM_region)
        scene1_xyz_3xN          = Geometry.geo_to_xyz(scene1_loc_3xN,earth_radius,earth_eccentricity)

        scene2_loc_3xN          = SimSetup.define_scene_pixels(t_loc_1_range, t_loc_2_range, 0, ref_scene_height*ones(size(DEM_region)))
        scene2_xyz_3xN          = Geometry.geo_to_xyz(scene2_loc_3xN,earth_radius,earth_eccentricity)

        # Define user parameters
        params = UserParameters.inputParameters(
            mode                = 1, #1:SAR
            processing_mode     = 1, #1:All platforms considered for processing
            pos_n               = [0 2.4]*1e3 , #Platform positions along n
            s_loc_1             = t_loc_1_range, 
            s_loc_2             = t_loc_2_range, #t_loc_2_range[1]-0.02:lon_res:t_loc_2_range[end]+0.02, #t_loc_2_range, 
            s_loc_3             = ref_scene_height, 
            SAR_duration        = 1.3,
            SAR_start_time      = -0.65,
            look_angle          = look_angle_test, #25.5,#2, #platform_look_angle,
            p_heading           = platform_heading, 
            ts_coord_sys        = "LLH", #"LLH",
            #ROSE-L parameters
            fp                  = 200, #1550, # Modified PRF based on scene length in along-track
            p_t0_LLH            = platform_ref_point,
            pulse_length        = 40e-6,
            dt_orbits           = 0.05,
            bandwidth           = 54e6, #80e6,
            user_defined_orbit  = 2,
            fc                  = 1.26e9,
            Δt                  = 0.1e-9,
        )
        # Check consistency of input parameters
        paramsIsValid = UserParameters.validateInputParams(params)

        if params.mode == 1 # SAR
            global p_mode       = 2
        elseif params.mode == 2 # SIMO
            global p_mode       = 1
        elseif params.mode == 3 # MIMO
            global p_mode       = 1.38
        end

        # Compute orbits time, position, and velocity
        orbit_time, orbit_pos, orbit_vel = Orbits.computeTimePosVel(params)

        # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
        p_xyz, Nst, slow_time   = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)
        v_xyz, Nst, slow_time   = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_vel, params)

        # Read number of platforms
        Np                      = size(orbit_pos)[2] # number of platforms
        ref_plat                = 1 #incicate the reference platform

        @timeit to "Compute geometries " begin
            #geometry computations based on scene
            s_slant_range_p1, s_look_angle_p1, s_incidence_angle_p1, s_slant_range_p2, s_look_angle_p2, s_incidence_angle_p2, s_perp_baseline, s_vert_wavnum, s_local_incidence_angle, s_range_slope_angle, s_critical_baseline, s_correlation_theo = Interferometry.get_scene_geometry_values(p_xyz, v_xyz, scene2_xyz_3xN, N_all, 1, 2, p_mode, params, "Flat")
            #geometry computations based on targets
            t_slant_range_p1, t_look_angle_p1, t_incidence_angle_p1, t_slant_range_p2, t_look_angle_p2, t_incidence_angle_p2, t_perp_baseline, t_vert_wavnum, t_local_incidence_angle, t_range_slope_angle, t_critical_baseline, t_correlation_theo = Interferometry.get_scene_geometry_values(p_xyz, v_xyz, scene1_xyz_3xN, N_all, 1, 2, p_mode, params, "")
        end

        # Get target reflectivities
        t_targets_ref_corr      = zeros(size(t_xyz_3xN,2),1)
        pixel_area              = (lat_res*pi/180*earth_radius) * (lon_res*pi/180*(earth_radius+cosd(params.p_heading)))

        # Get local incidence angle in 2D target grid
        t_local_incidence_angle_all = zeros(size(t_xyz_3xN,2),1)
        tk=1
        for tj = 1:size(t_xyz_3xN,2)
            t_local_incidence_angle_all[tj] = t_local_incidence_angle[tk]
            if mod(tj,num_targ_vol)==0
                tk=tk+1
            end
        end

        #Get sigma0 from scattering module
        for ti=1:size(t_xyz_3xN,2)
            t_targets_ref_corr[ti]    =   ( (pixel_area/num_targ_vol) .*  Scattering.TST_surface_brcs(2,params.λ,t_local_incidence_angle_all[ti],0.0,t_local_incidence_angle_all[ti],180.0-0.0,3,1.0) )
        end

        #Matrix dimensions
        dims_s_1                    = length(t_loc_1_range)
        dims_s_2                    = length(t_loc_2_range)
        dims_s_3                    = length(t_loc_3_range)

        # Apply antenna pattern
        #if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
        #    Antenna.applyAntennaPattern!(t_targets_ref_corr, p_xyz, orbit_vel, params)
        #end

        # Generate range spread function (matched filter output)
        min_range, max_range = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
        t_rx_tot            = 2*(max_range-min_range)/c + (1*params.pulse_length) #1
        t_rx                = -(t_rx_tot/2):params.Δt:(t_rx_tot/2)
        
        #truncate the slowtime and generate the RSF function
        t_rsf_tot           = 0.01*params.pulse_length # s duration of RSF #0.01 #0.00064
        S_rsf, t_rsf        = RSF.ideal_RSF_pulsewindow(t_rsf_tot, params) 

        # Generate TomoSAR raw data
        ref_range                       = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)

        @timeit to "rawdata new2" begin # This takes major chunk of simulation time
            rawdata = Generate_Raw_Data.main_gen_rawdata_method1_distributed(t_xyz_3xN, p_xyz, S_rsf, t_rsf, t_rx, ref_range, t_targets_ref_corr, params, 24) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
        end
    
        #=
        @timeit to "processdata data - scene grid" begin
            s_SAR_processed_3D = Process_Raw_Data.SAR_processing(rawdata, scene2_xyz_3xN, p_xyz, t_rx, ref_range, params, dims_s_1, dims_s_2, dims_s_3)
        end
        s_SAR_processed_2D_1p = s_SAR_processed_3D[1,:,:,1]'
        s_SAR_processed_2D_2p = s_SAR_processed_3D[2,:,:,1]'
        =#

        @timeit to "processdata data - target grid" begin 
            t_SAR_processed_3D = Process_Raw_Data.SAR_processing(rawdata, scene1_xyz_3xN, p_xyz, t_rx, ref_range, params, dims_s_1, dims_s_2, dims_s_3)
        end
        t_SAR_processed_2D_1p = t_SAR_processed_3D[1,:,:,1]'
        t_SAR_processed_2D_2p = t_SAR_processed_3D[2,:,:,1]'

        # Write outputs to files
        if ~ispath(savepath)
            mkdir(savepath)
        end

        #=
        # Write Scene geometry parameters to file
        output_data                 = zeros(dims_s_2, dims_s_1,14)
        output_data[:,:,1]          = reshape(s_slant_range_p1, dims_s_2, dims_s_1)
        output_data[:,:,2]          = reshape(s_slant_range_p2, dims_s_2, dims_s_1)
        output_data[:,:,3]          = reshape(s_look_angle_p1, dims_s_2, dims_s_1)
        output_data[:,:,4]          = reshape(s_look_angle_p2, dims_s_2, dims_s_1)
        output_data[:,:,5]          = reshape(s_incidence_angle_p1, dims_s_2, dims_s_1)
        output_data[:,:,6]          = reshape(s_incidence_angle_p2, dims_s_2, dims_s_1)
        output_data[:,:,7]          = reshape(s_perp_baseline, dims_s_2, dims_s_1)
        output_data[:,:,8]          = reshape(s_vert_wavnum, dims_s_2, dims_s_1)
        output_data[:,:,9]          = reshape(s_local_incidence_angle, dims_s_2, dims_s_1)
        output_data[:,:,10]         = reshape(s_range_slope_angle, dims_s_2, dims_s_1)
        output_data[:,:,11]         = reshape(s_critical_baseline, dims_s_2, dims_s_1)
        output_data[:,:,12]         = reshape(s_correlation_theo, dims_s_2, dims_s_1)
        output_data[:,:,14]         = ref_scene_height .* ones(dims_s_2, dims_s_1)
        file_flag = Export_Output.export_tif_file(savepath*string(Sim_idx)*"_scene_geometry_"*string(B_idx)*".tif", output_data, Float64, ag_geotransform, ag_ref)

        # Write SAR processed output data to file
        output_data                 = zeros(ComplexF64, dims_s_2, dims_s_1,2)
        output_data[:,:,1]          = reshape(s_SAR_processed_2D_1p, dims_s_2, dims_s_1)
        output_data[:,:,2]          = reshape(s_SAR_processed_2D_2p, dims_s_2, dims_s_1)
        file_flag                   = Export_Output.export_tif_file(savepath*string(Sim_idx)*"_sim_output_main_scene_"*string(B_idx)*".tif", output_data, ComplexF64, ag_geotransform, ag_ref)
        =#

        # Write Target geometry parameters to file
        output_data                 = zeros(dims_s_2, dims_s_1,14)
        output_data[:,:,1]          = reshape(t_slant_range_p1, dims_s_2, dims_s_1)
        output_data[:,:,2]          = reshape(t_slant_range_p2, dims_s_2, dims_s_1)
        output_data[:,:,3]          = reshape(t_look_angle_p1, dims_s_2, dims_s_1)
        output_data[:,:,4]          = reshape(t_look_angle_p2, dims_s_2, dims_s_1)
        output_data[:,:,5]          = reshape(t_incidence_angle_p1, dims_s_2, dims_s_1)
        output_data[:,:,6]          = reshape(t_incidence_angle_p2, dims_s_2, dims_s_1)
        output_data[:,:,7]          = reshape(t_perp_baseline, dims_s_2, dims_s_1)
        output_data[:,:,8]          = reshape(t_vert_wavnum, dims_s_2, dims_s_1)
        output_data[:,:,9]          = reshape(t_local_incidence_angle, dims_s_2, dims_s_1)
        output_data[:,:,10]         = reshape(t_range_slope_angle, dims_s_2, dims_s_1)
        output_data[:,:,11]         = reshape(t_critical_baseline, dims_s_2, dims_s_1)
        output_data[:,:,12]         = reshape(t_correlation_theo, dims_s_2, dims_s_1)
        output_data[:,:,13]         = reshape(t_targets_ref_corr,num_targ_vol, dims_s_2, dims_s_1)[1,:,:] # reflectivity of 1st target
        output_data[:,:,14]         = DEM_region[:,:]'
        file_flag                   = Export_Output.export_tif_file(savepath*string(Sim_idx)*"_target_geometry_"*string(B_idx)*".tif", output_data, Float64, ag_geotransform, ag_ref)

        output_data                 = zeros(ComplexF64, dims_s_2, dims_s_1, 2)
        output_data[:,:,1]          = reshape(t_SAR_processed_2D_1p, dims_s_2, dims_s_1)
        output_data[:,:,2]          = reshape(t_SAR_processed_2D_2p, dims_s_2, dims_s_1)
        file_flag                   = Export_Output.export_tif_file(savepath*string(Sim_idx)*"_sim_output_main_target_"*string(B_idx)*".tif", output_data, ComplexF64, ag_geotransform, ag_ref)
        
        #if ~ispath(savepath2)
        #    mkdir(savepath2)
        #end
        #@save savepath2*"Simulation_info_"*string(Sim_idx)*"_"*string(B_idx)*".jld" params scene1_xyz_3xN scene2_xyz_3xN t_xyz_3xN p_xyz N_all t_rx ref_range rawdata S_rsf
        end #timeit
    end 

end #timeit
    
println(to)
# END