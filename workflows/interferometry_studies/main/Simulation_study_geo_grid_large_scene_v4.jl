using Distributed
if Sys.islinux()
    addprocs(12) 
    @everywhere DEPOT_PATH[1]="/u/epstein-z0/wblr/joshil/Julia/.julia" # for Epstein
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
@everywhere include("../../../modules/simsetup.jl")
@everywhere include("../../../modules/user_parameters.jl")
@everywhere include("../../../modules/scattering.jl")
@everywhere include("../../../modules/interferometry.jl")
@everywhere include("../../../modules/dem.jl")
@everywhere include("../../../modules/export_output.jl")

@everywhere using Statistics
@everywhere using Parameters
@everywhere using .UserParameters
@everywhere using TimerOutputs
@everywhere using JLD2
@everywhere using SharedArrays

@everywhere  to         = TimerOutput()

c                       = 299792458 # Speed of light
earth_radius            = 6378.137e3 # Earth semi-major axis at equator

@timeit to "Overall time" begin

#platform_ref_point      = [34.2071;-121.83;697.5e3]
platform_ref_point      = [34.8043;-117.8153;697.5e3]
platform_heading        = 0.0
platform_look_dir       = "right" 
platform_look_angle     = 30.0  

lat_res                 = 0.000033
lon_res                 = 0.000033
lat_extent              = 0.000333*2 #0.0108 #0.000333 #0.0108 #0.000333
lon_extent              = 0.02261*1 #0.05
NB_lat                  = 1
NB_lon                  = 1
target_mode             = 2     # 1: target fixed in center, 2: Distributed target, 3: Distributed target with 1 dominant scatterer
num_targ_vol            = 3     # number of targets in each voxel
Sim_idx                 = 2     # For output file reference
savepath                = "/u/epstein-z0/darts/joshil/Outputs/Geo_outputs_set1/"


ground_range_ini        = Scene.lookangle_to_range(platform_look_angle,platform_ref_point[3],0.0, earth_radius)[2]  
displacement_N          = ground_range_ini .* cosd(platform_heading + 90) 
displacement_E          = ground_range_ini .* sind(platform_heading + 90)
trg_ref_lat             = platform_ref_point[1] + (displacement_N/earth_radius*180/pi)
trg_ref_lon             = platform_ref_point[2] + (displacement_E/(earth_radius+cosd(platform_ref_point[1]))*180/pi) 

trg_ref_lat_list, trg_ref_lon_list, B_lat_extent, B_lon_extent = SimSetup.segment_simulation_grid(trg_ref_lat, trg_ref_lon, lat_extent, lon_extent, NB_lat, NB_lon)

for B_idx = 1:length(trg_ref_lat_list)  
        
        @timeit to "Block time " begin
        
        @timeit to "Get DEM and slopes " begin
        #DEM_full, Geo_location_lat, Geo_location_lon, ag_geotransform, ag_ref = DEM.read_interp_DEM_from_source("/home/mlavalle/dat/nisar-dem-copernicus/EPSG4326.vrt", trg_ref_lat_list[B_idx], trg_ref_lon_list[B_idx], B_lat_extent, B_lon_extent, lat_res, lon_res, 1)
        DEM_full, Geo_location_lat, Geo_location_lon, ag_geotransform, ag_ref = DEM.create_slope_DEM(trg_ref_lat_list[B_idx], trg_ref_lon_list[B_idx], B_lat_extent, B_lon_extent, lat_res, lon_res, 0.0, 0.0) # slope dem for testing
 
        slope_lat, slope_lon    = DEM.get_slopes_from_DEM(DEM_full, Geo_location_lon, Geo_location_lat)

        DEM_region              = DEM_full[2:end-1,2:end-1]'
        slope_lat_region        = slope_lat[2:end-1,2:end-1]' 
        slope_lon_region        = slope_lon[2:end-1,2:end-1]'

        N_all                   = DEM.get_terrain_norm(DEM_region, slope_lat_region, slope_lon_region)
        end

        t_loc_1_range           = Geo_location_lat[1,2:end-1]   # 1-grid range
        t_loc_2_range           = Geo_location_lon[2:end-1,1]   # 2-grid range
        t_loc_3_range           = 0.0   # 3-grid range

        t_loc_1, t_loc_2, t_loc_3, t_ref_val = SimSetup.define_target_pixels(t_loc_1_range, t_loc_2_range, t_loc_3_range, num_targ_vol, target_mode, lat_res/2, lon_res/2, 0.5, DEM_region, 1)

        # Define user parameters
        params = UserParameters.inputParameters(
            mode                = 1, #1:SAR
            processing_mode     = 1, #1:All platforms considered for processing
            pos_n               = [0 2.4]*1e3 , #Platform positions along n
            s_loc_1             = t_loc_1_range, 
            s_loc_2             = t_loc_2_range, 
            s_loc_3             = 0.0, #t_loc_3_range, 
            SAR_duration        = 1.3,
            SAR_start_time      = -0.65,
            look_angle          = platform_look_angle,
            p_heading           = platform_heading, 
            target_pos_mode     = "CR",
            t_loc_1             = t_loc_1',
            t_loc_2             = t_loc_2',
            t_loc_3             = t_loc_3',
            t_ref               = t_ref_val',
            ts_coord_sys        = "LLH",
            #ROSE-L parameters
            fp                  = 200, #1550, # Modified PRF based on scene length in along-track
            p_t0_LLH            = platform_ref_point,
            pulse_length        = 40e-6,
            dt_orbits           = 0.05,
            bandwidth           = 54e6,
            user_defined_orbit  = 2,
            fc                  = 1.26e9,
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

        # Create target/scene location
        targets_loc, targets_ref, Nt    = Scene.construct_targets_str(params); # Nt: number of targets, targets: structure array containing target locations and reflectivities
        s_loc_3xN               = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops

        # Target location based only on the reference platform
        t_xyz_3xN, s_xyz_3xN, avg_peg   = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, reshape(orbit_pos[:,ref_plat,:],(size(orbit_pos)[1],1,size(orbit_pos)[3])), params) ## calculate avg heading from platform positions

        @timeit to "Compute geometries " begin
        #geometry computations based on scene
        s_slant_range_p1, s_look_angle_p1, s_incidence_angle_p1, s_slant_range_p2, s_look_angle_p2, s_incidence_angle_p2, s_perp_baseline, s_vert_wavnum, s_local_incidence_angle, s_range_slope_angle, s_critical_baseline, s_correlation_theo = Interferometry.get_scene_geometry_values(p_xyz, v_xyz, s_xyz_3xN, N_all, 1, 2, p_mode, params, "Scene", 3)
        #geometry computations based on targets
        t_slant_range_p1, t_look_angle_p1, t_incidence_angle_p1, t_slant_range_p2, t_look_angle_p2, t_incidence_angle_p2, t_perp_baseline, t_vert_wavnum, t_local_incidence_angle, t_range_slope_angle, t_critical_baseline, t_correlation_theo = Interferometry.get_scene_geometry_values(p_xyz, v_xyz, t_xyz_3xN, N_all, 1, 2, p_mode, params, "Target", 3)
        end

        # Get target reflectivities
        t_targets_ref_corr            = zeros(size(t_xyz_3xN,2),1)
        pixel_area                      = (lat_res*pi/180*earth_radius) * (lon_res*pi/180*(earth_radius+cosd(params.p_heading)))

        for ti=1:size(t_xyz_3xN,2)
            t_targets_ref_corr[ti]    = (pixel_area/num_targ_vol) .*  Scattering.TST_surface_brcs(2,params.λ,t_local_incidence_angle[ti],0.0,t_local_incidence_angle[ti],180.0-0.0,3,1.0)
        end

        # Apply antenna pattern
        if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
            Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
        end

        # Generate range spread function (matched filter output)
        min_range, max_range            = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
        Trx                             = 2*(max_range-min_range)/c + 2*params.pulse_length # s duration of RX window
        Srx, MF, ft, t_rx               = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window

        trunc_fac                       = Int(round(2*(max_range-min_range)/c/params.Δt,digits=-1))+26000
        C,D                             = findmax(Srx)
        Srx                             = Srx[D-Int(trunc_fac/2):D+Int(trunc_fac/2)]
        #MF                              = MF[D-Int(trunc_fac/2):D+Int(trunc_fac/2)]
        #ft                              = ft[D-Int(trunc_fac/2):D+Int(trunc_fac/2)]
        t_rx                            = t_rx[D-Int(trunc_fac/2):D+Int(trunc_fac/2)]

        # Generate TomoSAR raw data
        ref_range                       = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
        #rawdata2                       = Generate_Raw_Data.main_RSF_slowtime_perf_opt(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)

        @timeit to "rawdata " begin # This takes major chunk of simulation time
        # Dist for rawdata
        global Nt=size(t_xyz_3xN,2) # number of targets
        global Np=size(p_xyz,2) # number of platforms
        global Nft=length(t_rx) # number of fast-time samples
        global Nst=size(p_xyz,3) # number of slow-time samples
        global Δt_ft=t_rx[2]-t_rx[1] # fast-time resolution

        global rawdata=SharedArray(zeros(ComplexF64, Nst,Np,Nft) )
        global ref_delay=2*ref_range/c # reference delay

        @sync @distributed for s=1:Nst # slow-time (pulses)
            Srx_shifted = zeros(Nft)
            temp_sum=zeros(ComplexF64,Nft)

            for i=1:Np # RX platform
                temp_sum.=0.0;
            for j=1:Nt # targets
                if t_targets_ref_corr[j]!=0
                    range_tx=Geometry.distance(t_xyz_3xN[:,j],p_xyz[:,i,s])
                    range_rx=Geometry.distance(t_xyz_3xN[:,j],p_xyz[:,i,s])
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=Int(round(rel_delay/Δt_ft))
                        if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                            circshift!(Srx_shifted, Srx, rel_delay_ind)
                            Srx_shifted[1:rel_delay_ind] .= zeros(rel_delay_ind);
                        elseif rel_delay_ind<0
                            circshift!(Srx_shifted, Srx, rel_delay_ind)
                            Srx_shifted[Nft-abs(rel_delay_ind)+1:end] .= zeros(abs(rel_delay_ind))
                        end
                        temp_sum .= temp_sum .+  (t_targets_ref_corr[j].*exp(-im*2*pi/(0.23793052222222222)*(range_tx+range_rx)).*Srx_shifted)

                end
            end
                rawdata[s,i,:].= temp_sum
            end
        end
        end
    
        @timeit to "processdata data - scene grid" begin
            s_SAR_processed_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
        end
        s_SAR_processed_2D_1p = s_SAR_processed_3D[1,:,:,1]'
        s_SAR_processed_2D_2p = s_SAR_processed_3D[2,:,:,1]'

        t_xyz_3xN_2         = reshape(t_xyz_3xN,3,num_targ_vol,length(params.s_loc_2),length(params.s_loc_1))[:,1,:,:]
        t_xyz_3xN_3         = reshape(t_xyz_3xN_2,3,length(params.s_loc_2)*length(params.s_loc_1))

        @timeit to "processdata data - target grid" begin 
            t_SAR_processed_3D = Process_Raw_Data.SAR_processing(rawdata, t_xyz_3xN_3, p_xyz, t_rx, ref_range, params)
        end
        t_SAR_processed_2D_1p = t_SAR_processed_3D[1,:,:,1]'
        t_SAR_processed_2D_2p = t_SAR_processed_3D[2,:,:,1]'

        # Write outputs to files
        if ~ispath(savepath)
            mkdir(savepath)
        end

        # Write Scene geometry parameters to file
        output_data                 = zeros(length(params.s_loc_2),length(params.s_loc_1),14)
        output_data[:,:,1]          = reshape(s_slant_range_p1,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,2]          = reshape(s_slant_range_p2,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,3]          = reshape(s_look_angle_p1,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,4]          = reshape(s_look_angle_p2,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,5]          = reshape(s_incidence_angle_p1,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,6]          = reshape(s_incidence_angle_p2,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,7]          = reshape(s_perp_baseline,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,8]          = reshape(s_vert_wavnum,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,9]          = reshape(s_local_incidence_angle,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,10]         = reshape(s_range_slope_angle,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,11]         = reshape(s_critical_baseline,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,12]         = reshape(s_correlation_theo,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,14]         = DEM_region[:,:]'
        file_flag = Export_Output.export_tif_file(savepath*string(Sim_idx)*"_scene_geometry_"*string(B_idx)*".tif", output_data, Float64, ag_geotransform, ag_ref)

        # Write Target geometry parameters to file based on reference target
        output_data                 = zeros(length(params.s_loc_2),length(params.s_loc_1),14)
        reference_target            = 1
        output_data[:,:,1]          = reshape(t_slant_range_p1,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,2]          = reshape(t_slant_range_p2,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,3]          = reshape(t_look_angle_p1,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,4]          = reshape(t_look_angle_p2,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,5]          = reshape(t_incidence_angle_p1,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,6]          = reshape(t_incidence_angle_p2,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,7]          = reshape(t_perp_baseline,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,8]          = reshape(t_vert_wavnum,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,9]          = reshape(t_local_incidence_angle,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,10]         = reshape(t_range_slope_angle,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,11]         = reshape(t_critical_baseline,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,12]         = reshape(t_correlation_theo,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,13]         = reshape(t_targets_ref_corr,num_targ_vol, length(params.s_loc_2),length(params.s_loc_1))[reference_target,:,:]
        output_data[:,:,14]         = DEM_region[:,:]'
        file_flag = Export_Output.export_tif_file(savepath*string(Sim_idx)*"_target_geometry_"*string(B_idx)*".tif", output_data, Float64, ag_geotransform, ag_ref)

        # Write SAR processed output data to file
        output_data                 = zeros(ComplexF64, length(params.s_loc_2),length(params.s_loc_1),4)
        output_data[:,:,1]          = reshape(s_SAR_processed_2D_1p,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,2]          = reshape(s_SAR_processed_2D_2p,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,3]          = reshape(t_SAR_processed_2D_1p,length(params.s_loc_2),length(params.s_loc_1))
        output_data[:,:,4]          = reshape(t_SAR_processed_2D_2p,length(params.s_loc_2),length(params.s_loc_1))
        file_flag = Export_Output.export_tif_file(savepath*string(Sim_idx)*"_sim_output_main_"*string(B_idx)*".tif", output_data, ComplexF64, ag_geotransform, ag_ref)

        @save savepath*"Simulation_info_"*string(Sim_idx)*"_"*string(B_idx)*".jld" params s_xyz_3xN t_xyz_3xN p_xyz N_all t_rx ref_range rawdata
        
        #@save savepath*"Output12_10km_30la_geo_la_"*string(Sim_idx)*"_"*string(B_idx)*".jld" SAR_images_3D stat_var_all_1p stat_var_all_2p rawdata params slant_range_all look_angle_all incidence_angle_all Critical_baseline_all Correlation_theo_all Perp_baseline_all Vert_wavnum_all local_incidence_angle range_slope_angle trg_slant_range_all trg_look_angle_all trg_incidence_angle_all trg_Critical_baseline_all trg_Correlation_theo_all trg_Perp_baseline_all trg_Vert_wavnum_all trg_local_incidence_angle trg_range_slope_angle trg_targets_ref_corr s_xyz_3xN t_xyz_3xN p_xyz DEM_region N_all t_rx ref_range trg_SAR_images_3D trg_stat_var_all_1p trg_stat_var_all_2p 

        end #timeit

end 

end #timeit
    
[rmprocs(p) for p in workers()]

println(to)
# END

