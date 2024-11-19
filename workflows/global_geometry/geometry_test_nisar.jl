using Distributed
if Sys.islinux()
    addprocs(24) 
    #@everywhere DEPOT_PATH[1]="/u/epstein-z0/wblr/joshil/Julia/.julia" # for Epstein
else
    addprocs(2) 
end

@everywhere include("../../../modules/generate_raw_data.jl")
@everywhere include("../../../modules/process_raw_data.jl")
@everywhere include("../../../modules/geometry.jl")
@everywhere include("../../../modules/scene.jl")
@everywhere include("../../../modules/range_spread_function.jl") 
@everywhere include("../../../modules/orbits.jl")
@everywhere include("../../../modules/antenna.jl")
@everywhere include("../../../modules/dem.jl")
@everywhere include("../../../modules/simsetup.jl")
@everywhere include("../../../modules/user_parameters.jl")
@everywhere include("../../../modules/scattering.jl")
@everywhere include("../../../modules/interferometry.jl")
@everywhere include("../../../modules/export_output.jl")
@everywhere include("../../../modules/data_processing.jl")


@everywhere using Statistics
@everywhere using Parameters
@everywhere using .UserParameters
@everywhere using TimerOutputs
@everywhere using JLD2
@everywhere using SharedArrays
@everywhere using NCDatasets
@everywhere using Dates


@everywhere  to         = TimerOutput()

global c                       = 299792458 # Speed of light
global earth_radius            = 6378.137e3 # Earth semi-major axis at equator
global ref_start_time          = 24188 #84618 #41369  #67437 #101799 #58550 TODO: 24188
global k                       = 0

@timeit to "Overall time" begin

for idx_t = 1:200

    if mod(idx_t-1,5) == 0
        global ref_start_time  = ref_start_time +1
        k = 0
    end

    orbit_dataset               = Dataset("/u/intrepid-z0/joshil/data/orbits/helical/ROSE_L_Helical_4sat.nc")

    Sim_idx                 = 14   # For output file reference
    savepath                = "/u/intrepid-z0/joshil/Outputs/TST_sims_geometry_grid_1/"*string(Sim_idx)*"/"
    savepath2               = "/u/intrepid-z0/joshil/Outputs/TST_sims_geometry_grid_1/Outputs_sims_all/"

    # Define user parameters
    params = UserParameters.inputParameters(
        mode                = 2, #1:SAR
        processing_mode     = 1, #1:All platforms considered for processing
        SAR_duration        = 2,
        look_angle          = 30, #25.5,#2, #platform_look_angle,
        ts_coord_sys        = "LLH", #"LLH",
        #ROSE-L parameters
        fp                  = 1550, #1550, # Modified PRF based on scene length in along-track
        pulse_length        = 40e-6,
        bandwidth           = 54e6, #80e6,
        fc                  = 1.26e9,
        )

        if params.mode == 1 # SAR
            global p_mode       = 2
        elseif params.mode == 2 # SIMO
            global p_mode       = 1
        elseif params.mode == 3 # MIMO
            global p_mode       = 1.38
        end

        t12_orbits 		        = orbit_dataset["time"][1:2] # first two time samples
        orbit_time_index      = Int(1):length(orbit_dataset["time"])
        # index range for orbit times for time interval of interest
        orbit_time1 		      = orbit_dataset["time"][orbit_time_index] # read in time data
        orbit_pos_ECI 	      = 1e3*orbit_dataset["position"][:,:,orbit_time_index] # read in position data, 3 x Np x Nt
        orbit_vel_ECI         = 1e3*orbit_dataset["velocity"][:,:,orbit_time_index] # read in velocity data, 3 x Np x Nt (used optionally in avg peg and heading calculation)
        dv 				            = orbit_dataset.attrib["epoch"];
        if typeof(dv) == String
        temp = dv
        dv= []
        dv = [parse(Int,temp[1:4]);parse(Int,temp[6:7]);parse(Int,temp[9:10]);parse(Int,temp[12:13]);parse(Int,temp[15:16]);parse(Int,temp[18:19])];
        end
        epoch 			          = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
        #global dcm 		        = Orbits.eci_dcm(orbit_time1, epoch);
        global dcm = orbit_dataset["dcm"][:,:,orbit_time_index] 
        orbit_pos1,orbit_vel1 = Orbits.ecef_orbitpos(orbit_pos_ECI,orbit_vel_ECI,dcm)
        orbit_time = orbit_time1
        orbit_pos = orbit_pos1
        orbit_vel = orbit_vel1

        # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
        slow_time = orbit_time[ref_start_time]+k : 1/params.fp : orbit_time[ref_start_time]+k+2 # create slow time axis
        k = k+2

        Nst = size(slow_time)[1] # number of slow-time samples (pulses processed)
    
        if Nst == 1;
            p_xyz = orbit_pos[:,:,ref_start_time]
            v_xyz = orbit_vel[:,:,ref_start_time]
        else
            p_xyz = Orbits.interp_orbit(orbit_time[ref_start_time-5:ref_start_time+5],orbit_pos[:,:,ref_start_time-5:ref_start_time+5],slow_time)
            v_xyz = Orbits.interp_orbit(orbit_time[ref_start_time-5:ref_start_time+5],orbit_vel[:,:,ref_start_time-5:ref_start_time+5],slow_time)
        end # interpolate orbit to slow time, 3 x Np x Nst, convert km to m

        ref_plat = 1
        sec_plat = 2

        platform_ref_point = Geometry.xyz_to_geo(p_xyz[:,ref_plat,1])
        platform_ref_point2 = Geometry.xyz_to_geo(p_xyz[:,ref_plat,2])
        

        platform_heading = Geometry.compute_heading([platform_ref_point[1];platform_ref_point2[1]], [platform_ref_point[2];platform_ref_point2[2]])

        ground_range_ini        = Scene.lookangle_to_range(params.look_angle,platform_ref_point[3],0.0, earth_radius)[2]  
        displacement_N          = ground_range_ini .* cosd(platform_heading[1] + 90) 
        displacement_E          = ground_range_ini .* sind(platform_heading[1] + 90)
        if params.left_right_look == "right"
            trg_ref_lat             = platform_ref_point[1] + (displacement_N/earth_radius*180/pi)
            trg_ref_lon             = platform_ref_point[2] + (displacement_E/(earth_radius+cosd(platform_ref_point[1]))*180/pi) 
        elseif params.left_right_look == "left"
            trg_ref_lat             = platform_ref_point[1] - (displacement_N/earth_radius*180/pi)
            trg_ref_lon             = platform_ref_point[2] - (displacement_E/(earth_radius+cosd(platform_ref_point[1]))*180/pi) 
        end 

        swath_wdith = 240e3
        lon_inc = swath_wdith*180/(pi*earth_radius * cos(pi*trg_ref_lat/180))

        DEM_full, Geo_location_lat_v, Geo_location_lon_v, ag_geotransform, ag_ref = DEM.read_DEM_from_source("/u/intrepid-z0/joshil/data/nisar-dem-copernicus/EPSG4326.vrt", [trg_ref_lat+0.135;trg_ref_lat], [trg_ref_lon;trg_ref_lon+lon_inc],1)

        ag_geotransform[1] = trg_ref_lon
        ag_geotransform[4] = trg_ref_lat

        Geo_location_lat = repeat(Geo_location_lat_v', length(Geo_location_lon_v))
        Geo_location_lon = repeat(Geo_location_lon_v, 1,length(Geo_location_lat_v))

        slope_lat, slope_lon    = DEM.get_slopes_from_DEM(DEM_full, Geo_location_lon, Geo_location_lat)

        DEM_region              = DEM_full[2:end-1,2:end-1]'
        slope_lat_region        = slope_lat[2:end-1,2:end-1]' 
        slope_lon_region        = slope_lon[2:end-1,2:end-1]'

        N_all                   = DEM.get_terrain_norm(DEM_region, slope_lat_region, slope_lon_region)
        
        scene1_loc_3xN          = SimSetup.define_scene_pixels(Geo_location_lat[1,2:end-1] , Geo_location_lon[2:end-1,1] , 0, DEM_region)
        s_xyz_3xN               = Geometry.geo_to_xyz(scene1_loc_3xN,earth_radius,earth_eccentricity)


        # Read number of platforms
        Np                      = size(orbit_pos)[2] # number of platforms

        @timeit to "Compute geometries " begin
            #geometry computations based on scene
            #s_slant_range_p1, s_look_angle_p1, s_incidence_angle_p1, s_slant_range_p2, s_look_angle_p2, s_incidence_angle_p2, s_perp_baseline, s_vert_wavnum, s_local_incidence_angle, s_range_slope_angle, s_critical_baseline, s_correlation_theo = Interferometry.get_scene_geometry_values(p_xyz, v_xyz, scene1_xyz_3xN, N_all, 1, 2, p_mode, params, "Flat")
            @unpack λ, mode, bandwidth, left_right_look = params

            #geometry computations based on scene
            global slant_range_ref                 = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global look_angle_ref                  = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global incidence_angle_ref             = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global slant_range_sec                 = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global look_angle_sec                  = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global incidence_angle_sec             = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global Perp_baseline_ref               = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global AT_baseline_ref                 = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global Vert_wavnum_ref                 = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global local_incidence_angle_ref       = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global range_slope_angle_ref           = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global Critical_baseline_ref           = SharedArray(zeros(size(s_xyz_3xN,2),1))
            global Correlation_theo_ref            = SharedArray(zeros(size(s_xyz_3xN,2),1))
        
            global mean_plats_pos_ref              = mean(p_xyz[:,ref_plat,:], dims=2)
            global mean_plats_pos_sec              = mean(p_xyz[:,sec_plat,:], dims=2)
        
            @sync @distributed for ti = 1:size(s_xyz_3xN,2)
                 slant_range_ref[ti]         = Geometry.distance( mean_plats_pos_ref  , s_xyz_3xN[:,ti] )
                 look_angle_ref[ti]          = Scene.slantrange_to_lookangle(earth_radius,slant_range_ref[ti],Geometry.xyz_to_geo(mean_plats_pos_ref)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3])[2]
                 incidence_angle_ref[ti]     = Scene.lookangle_to_incangle(look_angle_ref[ti],Geometry.xyz_to_geo(mean_plats_pos_ref)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3],earth_radius)
        
                 slant_range_sec[ti]         = Geometry.distance( mean_plats_pos_sec  , s_xyz_3xN[:,ti] )
                 look_angle_sec[ti]          = Scene.slantrange_to_lookangle(earth_radius,slant_range_sec[ti],Geometry.xyz_to_geo(mean_plats_pos_sec)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3])[2]
                 incidence_angle_sec[ti]     = Scene.lookangle_to_incangle(look_angle_sec[ti],Geometry.xyz_to_geo(mean_plats_pos_sec)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3],earth_radius)
        
                bs_perp, bs_at, bs_norm     = Orbits.get_perp_baselines_new(mean(p_xyz[:,:,:],dims=3), mean(v_xyz[:,:,:],dims=3), look_angle_ref[ti], 0.0, left_right_look, 1)
                 Perp_baseline_ref[ti]       = bs_perp[1,sec_plat,1]
                 AT_baseline_ref[ti]         = bs_at[1,sec_plat,1]
        
                Va                          = mean(p_xyz[:, ref_plat, :], dims=2) - s_xyz_3xN[:, ti]
                Vb                          = mean(p_xyz[:, sec_plat, :], dims=2) - s_xyz_3xN[:, ti]
                angle_ip                    = Data_Processing.angle_2vec(Va, Vb) * 1
        
                if mode == 1
                     Vert_wavnum_ref[ti]     = (4 * pi * (angle_ip * pi / 180) ) / (λ * sind(look_angle_ref[ti]))
                elseif mode == 2
                     Vert_wavnum_ref[ti]     = (2 * pi * (angle_ip * pi / 180) ) / (λ * sind(look_angle_ref[ti]))
                end
        
                #plat_pt_xyz                 =  mean(mean(p_xyz,dims=2),dims=3)[:] #??????
                plat_pt_xyz                 =  mean(p_xyz[:,ref_plat,:],dims=2)
                look_vec_xyz                = (plat_pt_xyz - s_xyz_3xN[:,ti])
                look_vec_xyz_norm           = (plat_pt_xyz - s_xyz_3xN[:,ti]) / Geometry.distance( plat_pt_xyz,s_xyz_3xN[:,ti])
        
                Geo_location                = Geometry.xyz_to_geo(s_xyz_3xN[:,ti])
                pegθ                        = Geo_location[1]*π/180
                pegϕ                        = Geo_location[2]*π/180
                #ENU to XYZ transformation matrix
                Menu_xyz                    = [-sin(pegϕ) -sin(pegθ)*cos(pegϕ) cos(pegθ)*cos(pegϕ);
                                            cos(pegϕ) -sin(pegθ)*sin(pegϕ) cos(pegθ)*sin(pegϕ);
                                            0            cos(pegθ)             sin(pegθ)]
                #XYZ to ENU transformation matrix
                Mxyz_enu                    = [-sin(pegϕ)           cos(pegϕ)             0;
                                            -sin(pegθ)*cos(pegϕ) -sin(pegθ)*sin(pegϕ)  cos(pegθ)  ;
                                            cos(pegθ)*cos(pegϕ)   cos(pegθ)*sin(pegϕ)  sin(pegθ)]
                           
                look_vec_enu                = Mxyz_enu * look_vec_xyz
                look_direction_norm         = sqrt( look_vec_enu[1] * look_vec_enu[1]  + look_vec_enu[2] * look_vec_enu[2])
        

                Nxyz                        = Menu_xyz * N_all[ti,:];
                 local_incidence_angle_ref[ti] = Data_Processing.angle_2vec(look_vec_xyz_norm, Nxyz)
                 range_slope_angle_ref[ti]   = atand((N_all[ti,1] * (look_vec_enu[1] / look_direction_norm)) + (N_all[ti,2] * (look_vec_enu[2] / look_direction_norm)))          
   
                #Critical_baseline_ref[ti]   = λ * ((2*bandwidth)/c) * slant_range_ref[ti] * tand(local_incidence_angle_ref[ti] - range_slope_angle_ref[ti]) / p_mode
                 Critical_baseline_ref[ti]   = λ * ((2*bandwidth)/c) * slant_range_ref[ti] * tand(local_incidence_angle_ref[ti]) / p_mode
        
                 Correlation_theo_ref[ti]    = 1 - (Perp_baseline_ref[ti] ./ (Critical_baseline_ref[ti]))
        
            end

        end

        global sigma0_ref      = SharedArray(zeros(size(s_xyz_3xN,2),1))
        pixel_area              = ((Geo_location_lat_v[1]-Geo_location_lat_v[2])*pi/180*earth_radius) * ((Geo_location_lon_v[2]-Geo_location_lon_v[1])*pi/180*(earth_radius+cosd(params.p_heading)))

        for ti=1:size(s_xyz_3xN,2)
            global sigma0_ref[ti]    =  (pixel_area) .*  Scattering.TST_surface_brcs(2,params.λ,local_incidence_angle_ref[ti],0.0,local_incidence_angle_ref[ti],180.0-0.0,3,1.0) 
        end

        dims_s_1                    = length(Geo_location_lat[1,2:end-1] )
        dims_s_2                    = length(Geo_location_lon[2:end-1,1] )
        dims_s_3                    = 1

        # Write outputs to files
        if ~ispath(savepath)
            mkdir(savepath)
        end

        # Write Scene geometry parameters to file
        output_data                 = zeros(dims_s_2, dims_s_1,15)
        output_data[:,:,1]          = reshape(slant_range_ref, dims_s_2, dims_s_1)
        output_data[:,:,2]          = reshape(slant_range_sec, dims_s_2, dims_s_1)
        output_data[:,:,3]          = reshape(look_angle_ref, dims_s_2, dims_s_1)
        output_data[:,:,4]          = reshape(look_angle_sec, dims_s_2, dims_s_1)
        output_data[:,:,5]          = reshape(incidence_angle_ref, dims_s_2, dims_s_1)
        output_data[:,:,6]          = reshape(incidence_angle_sec, dims_s_2, dims_s_1)
        output_data[:,:,7]          = reshape(Perp_baseline_ref, dims_s_2, dims_s_1)
        output_data[:,:,8]          = reshape(AT_baseline_ref, dims_s_2, dims_s_1)
        output_data[:,:,9]          = reshape(local_incidence_angle_ref, dims_s_2, dims_s_1)
        output_data[:,:,10]         = reshape(range_slope_angle_ref, dims_s_2, dims_s_1)
        output_data[:,:,11]         = reshape(Critical_baseline_ref, dims_s_2, dims_s_1)
        output_data[:,:,12]         = reshape(Correlation_theo_ref, dims_s_2, dims_s_1)
        output_data[:,:,13]         = reshape(Vert_wavnum_ref, dims_s_2, dims_s_1)
        output_data[:,:,14]         = reshape(sigma0_ref, dims_s_2, dims_s_1)
        output_data[:,:,15]         = DEM_region[:,:]'

        file_flag = Export_Output.export_tif_file(savepath*string(Sim_idx)*"_scene_geometry_"*string(idx_t)*".tif", output_data, Float64, ag_geotransform, ag_ref)

end

end #timeit
    
[rmprocs(p) for p in workers()]

println(to)
# END
