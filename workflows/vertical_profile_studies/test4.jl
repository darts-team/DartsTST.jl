using Distributed
addprocs(8) 

#@everywhere DEPOT_PATH[1]="/u/epstein-z0/wblr/joshil/Julia/.julia" # for Epstein

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

c = 299792458 
earth_radius = 6378.137e3 # Earth semi-major axis at equator

temp_op=zeros(51,10)
stat_var=zeros(51,10)

for jjj=1:10

        num_targ_vol            = 1
        t_loc_S_range           = 0#-16:2:16#0#-2:0.25:2 #20
        t_loc_C_range           = 0#-30:5:30#0#-2:0.25:2 #40
        t_loc_H_range           = 0:1:22#[1 2 3 4 5 6 7 8 9 10]
        #t_ref_all = [1 1 1]

        t_ref_all = [0.873827  0.949371  0.966266  0.942564  0.896318  0.845578  0.805504  0.779678  0.76879  0.773529  0.794583  0.831239  0.877165  0.924627  0.96589  0.993218  1.0  0.984113  0.944556  0.880329  0.790431  0.675962  0.546409]

       #t_ref_all = [ 0.3 0.8 1 0.8 0.5 0.75 0.9 1 0.8 0.7 0.6 0.3 0.1 0.2 0.6 1 0.8]

        tgrid_midS = 6
        tgrid_midC = 6

        t_loc_S                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_C                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_H                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_ref_val               =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
    #=
        m=1
        for i=1:length(t_loc_S_range)
            for j= 1:length(t_loc_C_range)
                for k=1:length(t_loc_H_range)
                    for l=1:num_targ_vol
                            global t_loc_S[m] = t_loc_S_range[i] + (rand(Uniform(-1,1)) .* 0.5 )#(t_loc_S_range[2]-t_loc_S_range[1]/2) )
                            global t_loc_C[m] = t_loc_C_range[j] + (rand(Uniform(-1,1)) .* 0.5)#(t_loc_C_range[2]-t_loc_C_range[1]/2) )
                            global t_loc_H[m] = t_loc_H_range[k] + (rand(Uniform(-1,1)) .* 0.5)#(0.24/2))
                            global t_ref_val[m] = 1 #itp_ht(t_loc_H[m])
                            m=m+1;
                    end
                end
            end
        end
=#
        
        m=1
        n=1
        for i=1:length(t_loc_S_range)

            for j= 1:length(t_loc_C_range)
                for k=1:length(t_loc_H_range)
                    for l=1:num_targ_vol
                        if l==1
                           global t_loc_S[m] = t_loc_S_range[i] #+ (rand(Uniform(-1,1)) .* 0.5 )#(t_loc_S_range[2]-t_loc_S_range[1]/2) )
                            global t_loc_C[m] = t_loc_C_range[j] #+ (rand(Uniform(-1,1)) .* 0.5)#(t_loc_C_range[2]-t_loc_C_range[1]/2) )
                            global t_loc_H[m] = t_loc_H_range[k] #+ (rand(Uniform(-1,1)) .* 0.5)#(0.24/2))
                            global t_ref_val[m] = t_ref_all[k] #1#itp_ht(t_loc_H[m])
                            #global t_ref_val[m] = t_ref_all[n] #1#itp_ht(t_loc_H[m])
                            
                            m=m+1;
                        else
                            global t_loc_S[m] = t_loc_S_range[i] + (rand(Uniform(-1,1)) .* 1 )#(t_loc_S_range[2]-t_loc_S_range[1]/2) )
                            global t_loc_C[m] = t_loc_C_range[j] + (rand(Uniform(-1,1)) .* 2.5)#(t_loc_C_range[2]-t_loc_C_range[1]/2) )
                            global t_loc_H[m] = t_loc_H_range[k] + (rand(Uniform(-1,1)) .* 0.5)#(0.24/2))
                            global t_ref_val[m] = t_ref_all[k] #1#itp_ht(t_loc_H[m])  
                            #global t_ref_val[m] = t_ref_all[n] #1#itp_ht(t_loc_H[m])
                          
                            m=m+1;
                        end
                    end
                end
            end
            n=n+1
        end



#=        num_targ_vol            = 1
        t_loc_S_range           = -5:1:5
        t_loc_C_range           = -5:1:5
        t_loc_H_range           = 0 #Canopy_profile_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:2:length(0.0:0.5:ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))];

        tgrid_midS = 6
        tgrid_midC = 6

        t_loc_S                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_C                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_H                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_ref_val               =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
    
        m=1
        for i=1:length(t_loc_S_range)
            for j= 1:length(t_loc_C_range)
                for k=1:length(t_loc_H_range)
                    for l=1:num_targ_vol
                        if l==1
                            global t_loc_S[m] = t_loc_S_range[i] #+ (rand(Uniform(-1,1)) .* 0.5 )#(t_loc_S_range[2]-t_loc_S_range[1]/2) )
                            global t_loc_C[m] = t_loc_C_range[j] #+ (rand(Uniform(-1,1)) .* 0.5)#(t_loc_C_range[2]-t_loc_C_range[1]/2) )
                            global t_loc_H[m] = t_loc_H_range[k] #+ (rand(Uniform(-1,1)) .* 0.5)#(0.24/2))
                            global t_ref_val[m] = 1#itp_ht(t_loc_H[m])
                            m=m+1;
                        else
                            global t_loc_S[m] = t_loc_S_range[i] + (rand(Uniform(-1,1)) .* 0.5 )#(t_loc_S_range[2]-t_loc_S_range[1]/2) )
                            global t_loc_C[m] = t_loc_C_range[j] + (rand(Uniform(-1,1)) .* 0.5)#(t_loc_C_range[2]-t_loc_C_range[1]/2) )
                            global t_loc_H[m] = t_loc_H_range[k] + (rand(Uniform(-1,1)) .* 0.5)#(0.24/2))
                            global t_ref_val[m] = itp_ht(t_loc_H[m])  
                            m=m+1;
                        end
                    end
                end
            end
        end

=#

        #=
        profile = reshape(t_ref_val,3,23,11,11)
        profiles_mean = mean(mean(mean(profile,dims=1),dims=3),dims=4)

        heights = reshape(t_loc_H,3,23,11,11)
        heights_mean = mean(mean(mean(heights,dims=1),dims=3),dims=4)

        display(plot(profiles_mean[1,:,1,1,],heights_mean[1,:,1,1]))

        p1=plot()
        for i=1:3
            for j=1:11
                for k=1:11
                    p1 = (plot!(profile[i,:,j,k,],heights[i,:,j,k], legend=:false))
                end
            end
        end
        display(p1)
=#

       #=
        # gEt canopy profiles
        targ_loc        = Canopy_profile_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:2:length(0.0:0.5:ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))];
        for i=1:length(targ_loc)
            targ_loc[i] = targ_loc[i] + (rand(Uniform(-1,1)) .* (0.2367/2))
        end
        targ_ref        = Canopy_profiles[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:2:length(0.0:0.5:ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))];
        #targ_ref = transpose(targ_ref ./ maximum(targ_ref))

        =#
        # Define user parameters
        params = UserParameters.inputParameters(
            mode = 1,
            processing_mode = 1,
            #pos_n   = [-16 -12 -8 -4 0 4 8 12 16]*1e3 ,
            #pos_n   = [-70 -65 -60 -55 -50 -45 -40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70]*1e3 ,
            #pos_n   = [ -20 -10  0  10]*1e3 ,
            #pos_n   = [-15 -10 -5  0  5  10  15]*1e3 ,
            #pos_n   = [-15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10 12.5 15]*1e3 ,
            #pos_n   = [-30 -27.5 -25 -22.5 -20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30]*1e3 ,
            #pos_n   = collect(-30:1:30)' .*1e3,#[-30 -27.5 -25 -22.5 -20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30]*1e3 ,
            pos_n   = collect(-15:1:15)' .*1e3,#[-30 -27.5 -25 -22.5 -20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30]*1e3 ,

            #pos_n   = [-30 -28 -26 -24 -22 -20 -18 -16 -14 -12 -10 -8 -6 -4 -2 0 2 4 6 8 10 12 14 16.0 18 20 22 24 26 28 30 ]*1e3 ,


            fp = 40,#40,# 1510, #40,
            s_loc_1 = -32:2:32,#-160:1:160,#-186:2:186,#-20:1:20,
            s_loc_2 = -60:5:60,#-90:1:90,#-132:5:132,#-20:1:20,
            s_loc_3 = -25:1:25,
            SAR_duration = 5,
            SAR_start_time=-2.5,

            look_angle = 30,#30,# lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1],

            #target_pos_mode = "layered-grid-GEDIL2",
            #t_loc_3 = targ_loc,
            #t_ref   = targ_ref, 
            target_pos_mode = "CR",
            t_loc_1 = t_loc_S',
            t_loc_2 = t_loc_C',
            t_loc_3 = t_loc_H',
            t_ref   = t_ref_val', 

            #p_t0_LLH=[0.0;0.0;697e3],
            #pulse_length=40e-6,
            #dt_orbits=0.1,

        )
        # Check consistency of input parameters
        paramsIsValid = UserParameters.validateInputParams(params)

        filt_len = 5 # Covariance matrix box filter size

        # theoretical resolution factor
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

        # Compute orbits time, position, and velocity
        #global orbit_time = orbit_time_all[close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        #global orbit_pos = orbit_pos_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        #global orbit_vel = orbit_vel_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]

        #global SAR_start_time = orbit_time_all[close_val_lat_lon[2]] - (params.SAR_duration / 2)

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

       #try
 
        # Generate range spread function (matched filter output)
        min_range, max_range    = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
        Trx                     = 2*(max_range-min_range)/c + 5*params.pulse_length # s duration of RX window
        Srx, MF, ft, t_rx       = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
        # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

        Srx2 = zeros(length(Srx))
        C,D = findmax(Srx)
        Srx2[D-12:D+12] .= C
        # Generate TomoSAR raw data
        ref_range               = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
        #rawdata2                 = Generate_Raw_Data.main_RSF_slowtime_perf_opt(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)


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
        Srx_shifted2 = zeros(Nft)
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
                     temp_sum .= temp_sum .+  (targets_ref[j].*exp(-im*2*pi/(0.23793052)*(range_tx+range_rx)).*Srx_shifted)

             end
         end

            rawdata[s,i,:].=temp_sum

        end
    end
    

    @everywhere include("../modules/process_raw_data.jl")

            # Process raw data to generate image
            if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
                image_3D        = Process_Raw_Data.main_SAR_tomo_3D_new(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
            elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
                if params.processing_mode == 1
                    SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
                elseif params.processing_mode == 2
                    # for co-flyer
                    SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata[:,2:size(p_xyz)[2],:], s_xyz_3xN, p_xyz, t_rx, ref_range, params)
                end
             image_3D           = Process_Raw_Data.tomo_processing_afterSAR_full(SAR_images_3D)
            end
 

        image_3D_pow = abs.(image_3D .* conj(image_3D))
        image_3D_ph = angle.(image_3D)

        display(heatmap(params.s_loc_2,params.s_loc_1,10 .* log10.(image_3D_pow[:,:,26]),xlabel="C [m]",ylabel="S [m]",
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))

        display(heatmap(params.s_loc_2,params.s_loc_1,(image_3D_ph[:,:,26]),xlabel="C [m]",ylabel="S [m]",
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))

        ref_ip = reshape(t_ref_val,13,17)
        display(heatmap(t_loc_C_range,t_loc_S_range,ref_ip',xlabel="C [m]",ylabel="S [m]",
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))


        display(heatmap(params.s_loc_2,params.s_loc_1,(image_3D_pow[:,:,26])./maximum(image_3D_pow[:,:,26]),xlabel="C [m]",ylabel="S [m]",
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))




        (plot(params.s_loc_2,(image_3D_pow[17,:,26])./maximum((image_3D_pow[17,:,26])),label="Simulation",))
        display(plot!(t_loc_C_range,t_ref_all[1,:],xlabel="C [m]",ylabel="Normalized Intensity", label="Input",
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))

        (plot(params.s_loc_1,(image_3D_pow[:,13,26])./maximum((image_3D_pow[:,13,26])),label="Simulation",))
        display(plot!(t_loc_S_range,t_ref_all[1,:],xlabel="S [m]",ylabel="Normalized Intensity", label="Input",
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))


        A = (image_3D_pow[17,:,26])./maximum((image_3D_pow[17,:,26])) 
        B = t_ref_all[1,:]
        display(plot(t_loc_C_range,10 .* log10.(A[7:19]./B),xlabel="C [m]",ylabel="Log Ratio", label="",
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))


        A = (image_3D_pow[:,13,26])./maximum((image_3D_pow[:,13,26])) 
        B = t_ref_all[1,:]
        display(plot(t_loc_S_range,10 .* log10.(A[9:25]./B),xlabel="C [m]",ylabel="Log Ratio", label="",
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))

        norm_flag  					= 1
        #plot_idx 					= [93,26,1]
        plot_idx 					= [17,13,26]

        #plot_idx 					= [21,21,21]

        #plot_idx 					= [Int64(ceil(length(params.s_loc_1)/2)),Int64(ceil(length(params.s_loc_2)/2)),Int64(ceil(length(params.s_loc_3)/2))]

        figsavepath                 = "/Users/joshil/Documents/Outputs/temp/"
        if ~ispath(figsavepath)
            mkdir(figsavepath)
        end
        Data_Plotting.plot_tomo_output(image_3D_pow, params,norm_flag, plot_idx, figsavepath, 10, "Back-projection output \n Normalized Intensity [dB]", "_orig")


        figsavepath                 = "/Users/joshil/Documents/Outputs/temp/"
        if ~ispath(figsavepath)
            mkdir(figsavepath)
        end
        Data_Plotting.plot_tomo_output_lin(image_3D_pow, params,norm_flag, plot_idx, figsavepath, "Back-projection output \n Normalized Intensity [dB]", "_orig")


        figsavepath                 = "/Users/joshil/Documents/Outputs/temp2/"
        if ~ispath(figsavepath)
            mkdir(figsavepath)
        end
        Data_Plotting.plot_tomo_output_phase(image_3D_ph, params, plot_idx, figsavepath, "Back-projection output \n Normalized Intensity [dB]", "_orig")


        #plot_var_op_bpa1        = image_3D_pow[plot_idx[1],plot_idx[2],:]
        plot_var_op_bpa1        = image_3D_pow[plot_idx[1]-8:plot_idx[1]+8,plot_idx[2]-6:plot_idx[2]+6,:]

        #plot_var_op_bpa1        = image_3D_pow[plot_idx[1]-5:plot_idx[1]+5,plot_idx[2]-5:plot_idx[2]+5,:]
        #plot_var_op_bpa1        = image_3D_pow[plot_idx[1]-2:plot_idx[1]+2,plot_idx[2]-2:plot_idx[2]+2,:]

        plot_var_op_bpa         = mean(mean(plot_var_op_bpa1,dims=1),dims=2)[:]

        temp_op[:,jjj] = plot_var_op_bpa


    end

    p1=plot()
    for i=1:10
        stat_var[:,i] = temp_op[:,i]./maximum(temp_op[:,i])
        p1 = (plot!((stat_var[:,i]),params.s_loc_3,ylabel="Height (H axis) [m]",xlabel="Output level",title ="Back-projection output \n Normalized Intensity", label="Sim", lw=2, legend=:false,
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360),ylim=(-20,25) ))
    end
    display(p1)

    p2=plot()
    for i=1:10
        stat_var[:,i] = temp_op[:,i]./maximum(temp_op[:,i])
        p2 = (plot!(10 .* log10.(stat_var[:,i]),params.s_loc_3,ylabel="Height (H axis) [m]",xlabel="Output level [dB]",title ="Back-projection output \n Normalized Intensity [dB]", label="Sim", lw=2, legend=:false,
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360),ylim=(-20,25) ))
    end
    display(p2)


        p1 = (plot((plot_var_op_bpa./maximum(plot_var_op_bpa)),params.s_loc_3,ylabel="Height (H axis) [m]",xlabel="Output level",title ="Back-projection output \n Normalized Intensity", label="Sim", lw=2, legend=:false,
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360),ylim=(-20,20) ))
        display(p1)


        p1 = (plot(10 .* log10.(plot_var_op_bpa./maximum(plot_var_op_bpa)),params.s_loc_3,ylabel="Height (H axis) [m]",xlabel="Output level [dB]",title ="Back-projection output \n Normalized Intensity [dB]", legend=:false, lw=2,
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
        display(p1)




        display(scatter(t_loc_S,t_loc_C,t_loc_H,xlabel="S [m]",ylabel="C [m]",zlabel="H [m]",markersize=0.2))

        p1 = (plot!((t_ref_all'),collect(t_loc_H_range),ylabel="Height (H axis) [m]",xlabel="Reflectivity input",title ="Back-projection output \n Normalized Intensity", label="", lw=2,
        tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
        display(p1)

        [rmprocs(p) for p in workers()]
