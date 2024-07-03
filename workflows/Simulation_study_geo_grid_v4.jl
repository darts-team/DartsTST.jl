using Distributed
if Sys.islinux()
    addprocs(4) 
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
@everywhere include("../modules/sync.jl")
@everywhere include("../modules/error_sources.jl")
@everywhere include("../modules/performance_metrics.jl")
@everywhere include("../modules/antenna.jl")
@everywhere include("../modules/simsetup.jl")
@everywhere include("../modules/user_parameters.jl")
@everywhere include("../modules/global_scale_support.jl")
@everywhere include("../modules/data_processing.jl")
@everywhere include("../modules/data_plotting.jl")

@everywhere using NCDatasets
@everywhere using Statistics
@everywhere using Parameters
@everywhere using Dates
@everywhere using StaticArrays
@everywhere using .UserParameters
@everywhere using SharedArrays
@everywhere using Interpolations
@everywhere using Plots
@everywhere using Peaks
@everywhere using TimerOutputs
@everywhere using JLD2
@everywhere using GeoDatasets
@everywhere using Distributions
@everywhere using Measures
@everywhere using ArchGDAL
@everywhere using GeoArrays

@everywhere import ArchGDAL as AG


@everywhere  to = TimerOutput()

c               = 299792458
earth_radius    = 6378.137e3 # Earth semi-major axis at equator

for monte_carlo_idx = 1:1 # Just running one simulation currently

    @timeit to "Overall time" begin

        # Read Copernicus 30m resolution DEM
        dataset         = AG.read("/home/mlavalle/dat/nisar-dem-copernicus/EPSG4326.vrt")

        DEM_orig        = AG.getband(dataset,1)

        gt              = AG.getgeotransform(dataset)
        ulx, uly        = gt[1], gt[4]
        x_pixel_size, y_pixel_size = gt[2], gt[6]
        num_cols        = AG.width(DEM_orig)
        num_rows        = AG.height(DEM_orig)

        Lat_vals_full        = zeros(num_rows,1);
        Lon_vals_full        = zeros(num_cols,1);

        for i=1:num_rows
            Lat_vals_full[i] = uly + i * y_pixel_size;
        end

        for i=1:num_cols
            Lon_vals_full[i] = ulx + i * x_pixel_size;
        end

        target_mode             = 2 #1: target fixed in center, 2: Distributed target, 3: Distributed target with 1 dominant scatterer
        num_targ_vol            = 3 # number of targets in each voxel
        t_loc_S_range           = 0#-0.00033:0.000033:0.00033-0.000033    #-80:8:72 #-516:8:516 #-780:8:772 #-516:8:516 #-516:8:516 # S-grid range
        t_loc_C_range           = 3.59:0.000033:3.64   #-5200:4:5196 #-520:4:516 #-10400:4:10396 #-5200:4:5196 #-520:4:516 # C-grid range
        t_loc_H_range           = 0 # H-grid range

        # Define targets and targets reflectivity TODO: move this code to a function
        t_loc_S                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_C                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_H                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_ref_val               =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)

        #temp_height = 0       
        m=1
        for i=1:length(t_loc_S_range)
            temp_ref = 3.16
	    #temp_height = -262
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
                            global t_loc_S[m] = t_loc_S_range[i] + (rand(Uniform(-1,1)) .* (0.000033/2))
                            global t_loc_C[m] = t_loc_C_range[j] + (rand(Uniform(-1,1)) .* (0.000033/2))
                            #global t_loc_H[m] = temp_height  + (rand(Uniform(-1,1)) .* 0.5)
                            A, B =  findmin(abs.(t_loc_S[m].-Lat_vals_full))
                            C, D =  findmin(abs.(t_loc_C[m].-Lon_vals_full))
                            global t_loc_H[m] = DEM_orig[D[1],B[1]] + (rand(Uniform(-1,1)) .* 0.5)
			    #global t_loc_H[m] = t_loc_H_range[k] + (rand(Uniform(-1,1)) .* 0.5)
                            global t_ref_val[m] = 1 #temp_ref #1
                            m=m+1
                        elseif target_mode == 3
                            if l==1
                                global t_loc_S[m] = t_loc_S_range[i]
                                global t_loc_C[m] = t_loc_C_range[j]
                                global t_loc_H[m] = t_loc_H_range[k]
                                global t_ref_val[m] = 1
                                m=m+1
                            else
                                global t_loc_S[m] = t_loc_S_range[i] + (rand(Uniform(-1,1)) .* 4)
                                global t_loc_C[m] = t_loc_C_range[j] + (rand(Uniform(-1,1)) .* 2)
                                #global t_loc_H[m] = temp_height  + (rand(Uniform(-1,1)) .* 0.5)
                                A, B =  findmin(abs.(t_loc_S[m].-Lat_vals_full))
                                C, D =  findmin(abs.(t_loc_C[m].-Lon_vals_full))
                                global t_loc_H[m] = DEM_orig[D[1],B[1]] + (rand(Uniform(-1,1)) .* 0.5)
                                #global t_loc_H[m] = t_loc_H_range[k] + (rand(Uniform(-1,1)) .* 0.5)
                                global t_ref_val[m] = 1
                                m=m+1
                            end
                        end
                    end
                end
                if mod(j,4)==0
                    temp_ref = temp_ref + 0.011 #0.072
                end
            end
        end

        # Define user parameters
        params = UserParameters.inputParameters(
            mode                = 1, #1:SAR
            processing_mode     = 1, #1:All platforms considered for processing
            pos_n               = [0 2.4]*1e3 , #Platform positions along n

            s_loc_1             = -0.00033:0.000033:0.00033-0.000033, #-80:8:72, #-516:8:516, #-780:8:772, #-516:8:516, #-516:8:516,
            s_loc_2             = 3.59:0.000033:3.64, #-5200:4:5196, #-10000:4:9996, #-5000:4:4996, #-10000:4:9996, #-520:4:516, #-10400:4:10396, #-5200:4:5196, #-520:4:516,
            s_loc_3             = 0,

            SAR_duration        = 1.3,
            SAR_start_time      = -0.65,

            look_angle          = 30,

            target_pos_mode     = "CR",
            t_loc_1             = t_loc_S',
            t_loc_2             = t_loc_C',
            t_loc_3             = t_loc_H',
            t_ref               = t_ref_val',
	    ts_coord_sys        = "LLH",

            #ROSE-L parameters
            fp                  = 1550,
            p_t0_LLH            = [0.0;0.0;697.5e3],
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

        global Norm_baseline_max    = maximum(bnorm) ./ 1e3
        global Norm_baseline_min     = minimum(filter(!iszero,bnorm)) ./ 1e3
        global Norm_baseline_mean    = mean(filter(!iszero,bnorm)) ./ 1e3

        global Perp_baseline_max   = maximum(bperp) ./ 1e3
        global Perp_baseline_min     = minimum(filter(!iszero,bperp)) ./ 1e3
        global Perp_baseline_mean    = mean(filter(!iszero,bperp)) ./ 1e3

        global Par_baseline_max    = maximum(b_at) ./ 1e3
        global Par_baseline_min      = minimum(filter(!iszero,b_at)) ./ 1e3
        global Par_baseline_mean    = mean(filter(!iszero,b_at)) ./ 1e3

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

        # Apply antenna pattern
        if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
            Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
        end

        # Generate range spread function (matched filter output)
        min_range, max_range    = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
        Trx                     = 2*(max_range-min_range)/c + 2*params.pulse_length # s duration of RX window
        Srx, MF, ft, t_rx       = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
        # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

        trunc_fac   = Int(round(2*(max_range-min_range)/c/params.Δt,digits=-1))+16000
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
                    if targets_ref[j]!=0
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
                            temp_sum .= temp_sum .+  (targets_ref[j].*exp(-im*2*pi/(0.23793052222222222)*(range_tx+range_rx)).*Srx_shifted)
    
                    end
                end
                    rawdata[s,i,:].= temp_sum
                end
            end
        end

    
            @timeit to "processdata data" begin # This takes second major chunk of simulation time
                    # Process raw data to generate image
                    if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
                        image_3D        = Process_Raw_Data.main_SAR_tomo_3D_new(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
                    elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
                        if params.processing_mode == 1
                            @timeit to "Step 1" SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
                        elseif params.processing_mode == 2
                            # for co-flyer
                            @timeit to "Step 1" SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata[:,2:size(p_xyz)[2],:], s_xyz_3xN, p_xyz, t_rx, ref_range, params)
                        end
                    @timeit to "Step 2" image_3D           = Process_Raw_Data.tomo_processing_afterSAR_full(SAR_images_3D)
                    end
    
            end

        stat_var_all_1p = SAR_images_3D[1,:,:,1]
        stat_var_all_2p = SAR_images_3D[2,:,:,1]


        #geometry computations
        rng_slope                   = 0
        slant_range_all             = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        look_angle_all              = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        incidence_angle_all         = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        Critical_baseline_all       = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        Correlation_theo_all        = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        Perp_baseline_all           = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        Vert_wavnum_all             = zeros(size(s_xyz_3xN,2),size(p_xyz,2))


        for oi = 1:size(p_xyz,2)
            mean_plats_pos              = mean(p_xyz[:,oi,:], dims=2)
            for ti = 1:size(s_xyz_3xN,2)
                    slant_range_all[ti,oi]       = Geometry.distance( mean_plats_pos  , s_xyz_3xN[:,ti] )
                    look_angle_all[ti,oi]        = Scene.slantrange_to_lookangle(earth_radius,slant_range_all[ti,oi],Geometry.xyz_to_geo(mean_plats_pos)[3],Geometry.xyz_to_geo(t_xyz_3xN[:,ti])[3])[2]
                    incidence_angle_all[ti,oi]   = Scene.lookangle_to_incangle(look_angle_all[ti,oi],Geometry.xyz_to_geo(mean_plats_pos)[3],Geometry.xyz_to_geo(t_xyz_3xN[:,ti])[3],earth_radius)
                    Critical_baseline_all[ti,oi] = params.λ * ((2*params.bandwidth)/c) * slant_range_all[ti,oi] * tand(incidence_angle_all[ti,oi] - rng_slope) / p_mode
            end
        end
      
        for ti=1:size(s_xyz_3xN,2)
            bs_perp2, bs_at2, bs_norm2          = Orbits.get_perp_baselines_new(mean(p_xyz[:,:,:],dims=3), mean(v_xyz[:,:,:],dims=3), look_angle_all[ti,1], 0.0, params.left_right_look, 1)

            Perp_baseline_all[ti]               = bs_perp2[1,2,1]
            Correlation_theo_all[ti]            = 1 - (Perp_baseline_all[ti] ./ (Critical_baseline_all[ti,1]))

            Va                                  = mean(p_xyz[:, 1, :], dims=2) - s_xyz_3xN[:, ti]
            Vb                                  = mean(p_xyz[:, 2, :], dims=2) - s_xyz_3xN[:, ti]
                                    
            angle_ip                            = Data_Processing.angle_2vec(Va, Vb) * 1
     
            if params.mode == 1
                Vert_wavnum_all[ti]            = (4 * pi * (angle_ip * pi / 180) ) / (params.λ * sind(look_angle_all[ti,1]))
            elseif params.mode == 2
                Vert_wavnum_all[ti]            = (2 * pi * (angle_ip * pi / 180) ) / (params.λ * sind(look_angle_all[ti,1]))
            end
        end

        savepath                 = "/u/epstein-z0/darts/joshil/Outputs/SC_grid/"
        if ~ispath(savepath)
            mkdir(savepath)
        end
        
        @save savepath*"Output58_10km_30la_testgeo2.jld" image_3D SAR_images_3D stat_var_all_1p stat_var_all_2p rawdata params slant_range_all look_angle_all incidence_angle_all Critical_baseline_all Correlation_theo_all Perp_baseline_all Vert_wavnum_all rng_slope s_xyz_3xN p_xyz t_rx ref_range t_loc_H t_xyz_3xN

    end
    
    end

    [rmprocs(p) for p in workers()]

println(to)
            

