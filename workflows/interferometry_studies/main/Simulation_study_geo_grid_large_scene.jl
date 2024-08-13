using Distributed
if Sys.islinux()
    addprocs(24) 
    @everywhere DEPOT_PATH[1]="/u/epstein-z0/wblr/joshil/Julia/.julia" # for Epstein
else
    addprocs(2) 
end

@everywhere include("../modules/generate_raw_data.jl")
@everywhere include("../modules/process_raw_data.jl")
@everywhere include("../modules/geometry.jl")
@everywhere include("../modules/scene.jl")
@everywhere include("../modules/range_spread_function.jl") # as RSF
@everywhere include("../modules/orbits.jl")
@everywhere include("../modules/antenna.jl")
@everywhere include("../modules/simsetup.jl")
@everywhere include("../modules/user_parameters.jl")
@everywhere include("../modules/scattering.jl")
@everywhere include("../modules/interferometry.jl")
@everywhere include("../modules/dem.jl")

@everywhere using Statistics
@everywhere using Parameters
@everywhere using .UserParameters
@everywhere using TimerOutputs
@everywhere using JLD2
@everywhere using Distributions
@everywhere using SharedArrays

@everywhere  to = TimerOutput()

c               = 299792458
earth_radius    = 6378.137e3 # Earth semi-major axis at equator

function segment_simulation_grid(trg_ref_lat, trg_ref_lon, lat_extent, lon_extent, NB_lat, NB_lon)

    B_lat_extent = lat_extent / NB_lat
    B_lon_extent = lon_extent / NB_lon

    Add_B_lat_extent = 0 #B_lat_extent*0.1
    Add_B_lon_extent = B_lon_extent*0.1


    trg_ref_lat_list = zeros(NB_lat * NB_lon)
    trg_ref_lon_list = zeros(NB_lat * NB_lon)

    k = 1
    for i=1:NB_lat
        for j=1:NB_lon
            trg_ref_lat_list[k] = (trg_ref_lat - Add_B_lat_extent) + ((i-1) * (B_lat_extent  + Add_B_lat_extent))
            trg_ref_lon_list[k] = (trg_ref_lon - Add_B_lon_extent) + ((j-1) * (B_lon_extent  + Add_B_lon_extent))
            k=k+1
        end
    end

    return trg_ref_lat_list, trg_ref_lon_list, B_lat_extent+(2*Add_B_lat_extent), B_lon_extent+(2*Add_B_lon_extent)

end


@timeit to "Overall time" begin

#platform_ref_point      = [34.2071;-121.83;697.5e3]
platform_ref_point      = [34.8043;-117.8153;697.5e3]
platform_heading_ip     = 0.0
look_angle_ip 		    = 30.0  
lat_res                	= 0.000033
lon_res                 = 0.000033
lat_extent              = 0.0108 #0.0108 #0.000333
lon_extent              = 0.02261*1 #0.05
NB_lat                  = 1
NB_lon                  = 1
target_mode             = 2     # 1: target fixed in center, 2: Distributed target, 3: Distributed target with 1 dominant scatterer
num_targ_vol            = 3     # number of targets in each voxel

ground_range_ini        = Scene.lookangle_to_range(look_angle_ip,platform_ref_point[3],0.0, earth_radius)[2]
displacement_N          = ground_range_ini .* cosd(platform_heading_ip + 90) 
displacement_E          = ground_range_ini .* sind(platform_heading_ip + 90)
trg_ref_lat             = platform_ref_point[1] + (displacement_N/earth_radius*180/pi)
trg_ref_lon             = platform_ref_point[2] + (displacement_E/(earth_radius+cosd(platform_ref_point[1]))*180/pi) 

trg_ref_lat_list, trg_ref_lon_list, B_lat_extent, B_lon_extent = segment_simulation_grid(trg_ref_lat, trg_ref_lon, lat_extent, lon_extent, NB_lat, NB_lon)


for B_idx = 1:length(trg_ref_lat_list)  
        
        @timeit to "Block time " begin
        
        DEM_full, Geo_location_lat, Geo_location_lon, ag_geotransform, ag_ref = DEM.read_interp_DEM_from_source("/home/mlavalle/dat/nisar-dem-copernicus/EPSG4326.vrt", trg_ref_lat_list[B_idx], trg_ref_lon_list[B_idx], B_lat_extent, B_lon_extent, lat_res, lon_res)

        if iseven(size(DEM_full)[2]) || iseven(size(DEM_full)[1])
            B_lat_extent_new = B_lat_extent
            B_lon_extent_new = B_lon_extent
            if iseven(size(DEM_full)[2])
                B_lat_extent_new = B_lat_extent + lat_res
            end
            if iseven(size(DEM_full)[1])
                B_lon_extent_new = B_lon_extent + lon_res
            end
            DEM_full, Geo_location_lat, Geo_location_lon, ag_geotransform, ag_ref = DEM.read_interp_DEM_from_source("/home/mlavalle/dat/nisar-dem-copernicus/EPSG4326.vrt", trg_ref_lat_list[B_idx], trg_ref_lon_list[B_idx], B_lat_extent_new, B_lon_extent_new, lat_res, lon_res)
        end

        #DEM_region_o            = repeat(0:0.1:0.1*(size(DEM_full)[1]-1),1,size(DEM_full)[2]) # This line for testing
        DEM_region              = DEM_full[2:end-1,2:end-1]'

        slope_lat, slope_lon    = DEM.get_slopes_from_DEM(DEM_full, Geo_location_lon, Geo_location_lat)

        slope_lat_region 	    = slope_lat[2:end-1,2:end-1]' 
        slope_lon_region  	    = slope_lon[2:end-1,2:end-1]'

        N_all                   = DEM.get_terrain_norm(DEM_region, slope_lat_region, slope_lon_region)


        t_loc_S_range           = Geo_location_lat[1,2:end-1] #trg_ref_lat:lat_res:(trg_ref_lat+lat_extent)-lat_res	# S-grid range
        t_loc_C_range           = Geo_location_lon[2:end-1,1] #trg_ref_lon:lon_res:(trg_ref_lon+lon_extent)-lon_res	# C-grid range
        t_loc_H_range           = 0.0 	# H-grid range
        # Define targets and targets reflectivity TODO: move this code to a function
        t_loc_S                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_C                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_H                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_ref_val               =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)

        
         m=1
        for i=1:length(t_loc_S_range)
            for j= 1:length(t_loc_C_range)
                #temp_height=temp_height+2
                for k=1:length(t_loc_H_range)
                    for l=1:num_targ_vol
                        if target_mode == 1
                            global t_loc_S[m] = t_loc_S_range[i]
                            global t_loc_C[m] = t_loc_C_range[j]
                            global t_loc_H[m] = t_loc_H_range[k]
                            global t_ref_val[m] = 1
                            m=m+1
                        elseif target_mode == 2
                            global t_loc_S[m] = t_loc_S_range[i] + (rand(Uniform(-1,1)) .* (lat_res/2))
                            global t_loc_C[m] = t_loc_C_range[j] + (rand(Uniform(-1,1)) .* (lon_res/2))
                            global t_loc_H[m] = DEM_region[i,j] + (rand(Uniform(-1,1)) .* 0.5)
			    #global t_loc_H[m] = t_loc_H_range[k] + (rand(Uniform(-1,1)) .* 0.5)
                            global t_ref_val[m] = 1 
                            m=m+1
                        elseif target_mode == 3
                            if l==1
                                global t_loc_S[m] = t_loc_S_range[i]
                                global t_loc_C[m] = t_loc_C_range[j]
                                global t_loc_H[m] = t_loc_H_range[k]
                                global t_ref_val[m] = 1
                                m=m+1
                            else
                                global t_loc_S[m] = t_loc_S_range[i] + (rand(Uniform(-1,1)) .* (lat_res/2))
                                global t_loc_C[m] = t_loc_C_range[j] + (rand(Uniform(-1,1)) .* (lon_res/2))
                                global t_loc_H[m] = DEM_region[i,j]  + (rand(Uniform(-1,1)) .* 0.5)
                                #global t_loc_H[m] = t_loc_H_range[k] + (rand(Uniform(-1,1)) .* 0.5)
                                global t_ref_val[m] = 1
                                m=m+1
                            end
                        end
                    end
                end
            end
        end

        # Define user parameters
        params = UserParameters.inputParameters(
            mode                = 1, #1:SAR
            processing_mode     = 1, #1:All platforms considered for processing
            pos_n               = [0 2.4]*1e3 , #Platform positions along n

            s_loc_1             = t_loc_S_range, #trg_ref_lat:lat_res:(trg_ref_lat+lat_extent)-lat_res, 
            s_loc_2             = t_loc_C_range, #trg_ref_lon:lon_res:(trg_ref_lon+lon_extent)-lon_res,
            s_loc_3             = t_loc_H_range, #0.0,

            SAR_duration        = 1.3,
            SAR_start_time      = -0.65,

            look_angle          = look_angle_ip,
            p_heading           = platform_heading_ip, 

            target_pos_mode     = "CR",
            t_loc_1             = t_loc_S',
            t_loc_2             = t_loc_C',
            t_loc_3             = t_loc_H',
            t_ref               = t_ref_val',
	        ts_coord_sys        = "LLH",

            #ROSE-L parameters
            fp                  = 200, #1550,
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
            global p_mode   = 2
        elseif params.mode == 2 # SIMO
            global p_mode   = 1
        elseif params.mode == 3 # MIMO
            global p_mode   = 1.38
        end

        # Compute orbits time, position, and velocity
        orbit_time, orbit_pos, orbit_vel = Orbits.computeTimePosVel(params)

        # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
        p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)
        v_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_vel, params)

        # Read number of platforms
        Np  = size(orbit_pos)[2] # number of platforms

        if params.processing_mode == 1
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,1:end,:], orbit_vel[:,1:end,:], params.look_angle, 0.0, params.left_right_look,  1)
            avg_sep = maximum(bperp)/(Np - 1)
            ref_plat = 1 #incicate the reference platform
        elseif params.processing_mode == 2
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,2:end,:], orbit_vel[:,2:end,:],  params.look_angle, 0.0, params.left_right_look,  1)
            avg_sep = maximum(bperp)/(Np - 2)
            ref_plat = 2 #incicate the reference platform
        end

        mast_plat = 1

        # theoretical resolution along-n
        range_s, range_g = Scene.lookangle_to_range(params.look_angle, mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
        global res_theory_n = (c/params.fc)*range_s/p_mode/ (maximum(bperp) + avg_sep)

        # theoretical resolution along-track
        mu = 3.986004418e14
        sc_speed = sqrt(mu./(mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:])+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
        Lsa = sc_speed*params.SAR_duration + sc_speed/params.fp
        global res_theory_s = (c/params.fc)*range_s/2/Lsa

        global amb_H = (c/params.fc)*range_s/p_mode/avg_sep*sind(params.look_angle)
        global amb_N = (c/params.fc)*range_s/p_mode/avg_sep

        # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
        #p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, SAR_start_time, params)

        # Create target/scene location
        targets_loc, targets_ref, Nt = Scene.construct_targets_str(params); # Nt: number of targets, targets: structure array containing target locations and reflectivities
        s_loc_3xN  = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops

        # For closed form configuration
        #t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions

        # Target location based only on the reference platform
        t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, reshape(orbit_pos[:,ref_plat,:],(size(orbit_pos)[1],1,size(orbit_pos)[3])), params) ## calculate avg heading from platform positions


        #geometry computations based on scene
        slant_range_all                 = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        look_angle_all                  = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        incidence_angle_all             = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        Critical_baseline_all           = zeros(size(s_xyz_3xN,2),1)
        Correlation_theo_all            = zeros(size(s_xyz_3xN,2),1)
        Perp_baseline_all               = zeros(size(s_xyz_3xN,2),1)
        Vert_wavnum_all                 = zeros(size(s_xyz_3xN,2),1)
        local_incidence_angle           = zeros(size(s_xyz_3xN,2),1)
        range_slope_angle               = zeros(size(s_xyz_3xN,2),1)

        slant_range_all[:,1], look_angle_all[:,1], incidence_angle_all[:,1], slant_range_all[:,2], look_angle_all[:,2], incidence_angle_all[:,2], Perp_baseline_all, Vert_wavnum_all, local_incidence_angle, range_slope_angle, Critical_baseline_all, Correlation_theo_all = Interferometry.get_scene_geometry_values(p_xyz, v_xyz, s_xyz_3xN, N_all, 1, 2, p_mode, params, "Scene", 3)

        #geometry computations based on targets
        trg_slant_range_all             = zeros(size(t_xyz_3xN,2),size(p_xyz,2))
        trg_look_angle_all              = zeros(size(t_xyz_3xN,2),size(p_xyz,2))
        trg_incidence_angle_all         = zeros(size(t_xyz_3xN,2),size(p_xyz,2))
        trg_Critical_baseline_all       = zeros(size(t_xyz_3xN,2),1)
        trg_Correlation_theo_all        = zeros(size(t_xyz_3xN,2),1)
        trg_Perp_baseline_all           = zeros(size(t_xyz_3xN,2),1)
        trg_Vert_wavnum_all             = zeros(size(t_xyz_3xN,2),1)
        trg_local_incidence_angle       = zeros(size(t_xyz_3xN,2),1)
        trg_range_slope_angle           = zeros(size(t_xyz_3xN,2),1)

        trg_slant_range_all[:,1], trg_look_angle_all[:,1], trg_incidence_angle_all[:,1], trg_slant_range_all[:,2], trg_look_angle_all[:,2], trg_incidence_angle_all[:,2], trg_Perp_baseline_all, trg_Vert_wavnum_all, trg_local_incidence_angle, trg_range_slope_angle, trg_Critical_baseline_all, trg_Correlation_theo_all = Interferometry.get_scene_geometry_values(p_xyz, v_xyz, t_xyz_3xN, N_all, 1, 2, p_mode, params, "Target", 3)

        # Get target reflectivities
        trg_targets_ref_corr            = zeros(size(t_xyz_3xN,2),1)
        pixel_area                      = (lat_res*pi/180*earth_radius) * (lon_res*pi/180*(earth_radius+cosd(params.p_heading)))

        for ti=1:size(t_xyz_3xN,2)
            trg_targets_ref_corr[ti]            = (pixel_area/num_targ_vol) .*  Scattering.TST_surface_brcs(2,params.λ,trg_local_incidence_angle[ti],0.0,trg_local_incidence_angle[ti],180.0-0.0,3,1.0)
        end

        # Apply antenna pattern
        if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
            Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
        end

        # Generate range spread function (matched filter output)
        min_range, max_range    = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
        Trx                     = 2*(max_range-min_range)/c + 2*params.pulse_length # s duration of RX window
        Srx, MF, ft, t_rx       = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
        # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

        trunc_fac   = Int(round(2*(max_range-min_range)/c/params.Δt,digits=-1))+26000
        C,D         = findmax(Srx)
        Srx         = Srx[D-Int(trunc_fac/2):D+Int(trunc_fac/2)]
        #MF         = MF[D-Int(trunc_fac/2):D+Int(trunc_fac/2)]
        #ft         = ft[D-Int(trunc_fac/2):D+Int(trunc_fac/2)]
        t_rx        = t_rx[D-Int(trunc_fac/2):D+Int(trunc_fac/2)]

        # Generate TomoSAR raw data
        ref_range               = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
        #rawdata2                 = Generate_Raw_Data.main_RSF_slowtime_perf_opt(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)

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
                    if trg_targets_ref_corr[j]!=0
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
                            temp_sum .= temp_sum .+  (trg_targets_ref_corr[j].*exp(-im*2*pi/(0.23793052222222222)*(range_tx+range_rx)).*Srx_shifted)
    
                    end
                end
                    rawdata[s,i,:].= temp_sum
                end
            end
        end

    
        @timeit to "processdata data" begin # This takes second major chunk of simulation time
            SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
        end
        stat_var_all_1p = SAR_images_3D[1,:,:,1]
        stat_var_all_2p = SAR_images_3D[2,:,:,1]

        t_xyz_3xN_2 = reshape(t_xyz_3xN,3,num_targ_vol,length(params.s_loc_2),length(params.s_loc_1))[:,1,:,:]
        t_xyz_3xN_3 = reshape(t_xyz_3xN_2,3,length(params.s_loc_2)*length(params.s_loc_1))

        @timeit to "processdata data2" begin # This takes second major chunk of simulation time
            trg_SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, t_xyz_3xN_3, p_xyz, t_rx, ref_range, params)
        end
        trg_stat_var_all_1p = trg_SAR_images_3D[1,:,:,1]
        trg_stat_var_all_2p = trg_SAR_images_3D[2,:,:,1]

        savepath                 = "/u/epstein-z0/darts/joshil/Outputs/SC_grid_la/"
        if ~ispath(savepath)
            mkdir(savepath)
        end

        @save savepath*"Output5_10km_30la_geo_la5_"*string(B_idx)*".jld" SAR_images_3D stat_var_all_1p stat_var_all_2p rawdata params slant_range_all look_angle_all incidence_angle_all Critical_baseline_all Correlation_theo_all Perp_baseline_all Vert_wavnum_all local_incidence_angle range_slope_angle trg_slant_range_all trg_look_angle_all trg_incidence_angle_all trg_Critical_baseline_all trg_Correlation_theo_all trg_Perp_baseline_all trg_Vert_wavnum_all trg_local_incidence_angle trg_range_slope_angle trg_targets_ref_corr s_xyz_3xN t_xyz_3xN p_xyz DEM_region N_all t_rx ref_range trg_SAR_images_3D trg_stat_var_all_1p trg_stat_var_all_2p 
         
        end #timeit

end 

end #timeit
    
[rmprocs(p) for p in workers()]

println(to)
            
# END