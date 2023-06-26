using Distributed
addprocs(24) 

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


c = 299792458 #TODO does not work without redefining c here
earth_radius = 6378.137e3 # Earth semi-major axis at equator

@everywhere const to = TimerOutput()

@timeit to "Reading L2 and L3 file for canopy heights and profile " begin

#Read canopy heights
#filepath_GEDIL3 = "/Users/joshil/Documents/GEDI_Data/GEDI_L3_LandSurface_Metrics_V2_1952/data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
filepath_GEDIL3 = "/u/epstein-z0/wblr/joshil/DARTS/GEDI_Data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"

grid_res        = 100;
Canopy_heights_L3, Geo_location, size_row, size_col = Global_Scale_Support.read_GEDI_L3_data(filepath_GEDIL3, grid_res)

#@load "/Users/joshil/Documents/GEDI_Data/Outputs_L2/Output_GEDIL2_1year_2022.jld2"
@load "/u/epstein-z0/darts/joshil/GEDI/GEDI_L2/Output_GEDIL2_1year_2022_res100m.jld2"

global Canopy_heights           = Canopy_heights_L2[:,:];

global Canopy_profiles          = Canopy_profiles_L2_itp[:,:,:]
global Canopy_profile_heights   = Canopy_heights_profile_L2_itp[:,:,:]

GC.gc()

end

@timeit to "Initialization " begin
# Define outout variables
global Output_stat_bpa          = SharedArray(zeros(size_row,size_col,4))
global Output_stat_beamforming  = SharedArray(zeros(size_row,size_col,4))
global Output_stat_capon        = SharedArray(zeros(size_row,size_col,4))

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

global res_theory_n    = SharedArray(zeros(size_row,size_col,1))
global res_theory_s    = SharedArray(zeros(size_row,size_col,1))

global amb_H                = SharedArray(zeros(size_row,size_col,1))
global amb_N                = SharedArray(zeros(size_row,size_col,1))
global slnt_range           = SharedArray(zeros(size_row,size_col,1))

end

@timeit to "Read orbit file and get orbits " begin
#region_xlims = [50,110]
#region_ylims = [15,55]
region_xlims = 1:347#90:95
region_ylims = 1:146#30:30
#region_xlims        = 50:110
#region_ylims        = 15:55

lat_lon_idx         = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

#orbit_dataset       = Dataset("./inputs/orbit_output_04052023.nc")
#orbit_dataset       = Dataset("/Users/joshil/Documents/Orbits/Outputs/06152023/3/orbit_output_06152023_3.nc") # "orbit_output_04052023.nc") # Read orbits data in NetCDF format
orbit_dataset       = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/orbit_output_06152023_3.nc")

global mast_plat            = 1
flag_plat           = 1 #descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, mast_plat, flag_plat)

end


lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
global data2= data

global lat2 = lat' .* ones(length(lon))
global lon2  = ones(length(lat))' .* lon

global Lats_p = Geo_location.A[1,:,1]
global Lons_p = Geo_location.A[:,1,2]
global mask = SharedArray(zeros(length(Lons_p),length(Lats_p)) )

@timeit to "Processing loop over all pixels " begin

@sync @distributed for i1 = 1:size(lat_lon_idx,1)   

    close_val_lat_lon_m   = findmin(abs.(lat2.-Lats_p[lat_lon_idx[i1,2]]) + abs.(lon2.-Lons_p[lat_lon_idx[i1,1]]))
    if data2[close_val_lat_lon_m[2]]>0
        global mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] = 1
    end
    
    if mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] !=1
        #Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =0
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN
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

        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        
        global Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN;
    
    else
        if isnan(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
            global Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =0.0
            global Canopy_profiles[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:] .=0.0
            global Canopy_profile_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:] = Canopy_profile_heights[94,30,:]
        end
                
        global close_val_lat_lon   = Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = close_val_lat_lon[2]

        slrng_temp2 = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]]; 0]))
        #global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = acosd(orbit_pos_geo_all[:,close_val_lat_lon[2]][3] / slrng_temp2) 
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = Scene.slantrange_to_lookangle(earth_radius,slrng_temp2,orbit_pos_geo_all[:,close_val_lat_lon[2]][3],0.0)[2]


        #ref_profile_height, ref_profile_value, NoPeaks  = Global_Scale_Support.constrct_reflectivity_profile_exp(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        #interp_fn = LinearInterpolation(ref_profile_height, ref_profile_value, extrapolation_bc = Flat());
        #height_range = 0:1:40

        if Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] > 49
            Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = 49.0
        end


        targ_loc = Canopy_profile_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:2:length(0.0:0.5:ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))];
        targ_ref = Canopy_profiles[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:2:length(0.0:0.5:ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))];
        #targ_ref = transpose(targ_ref ./ maximum(targ_ref))

        # Define user parameters
        params = UserParameters.inputParameters(
            look_angle = lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1],
            target_pos_mode = "layered-grid-GEDIL2",
            t_loc_3 = targ_loc,
            t_ref   = targ_ref, 
        )
       #params = UserParameters.inputParameters(look_angle = lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
       #params = UserParameters.inputParameters()
        filt_len = 1

        # theoretical resolution
        if params.mode == 1 # SAR
            global p_mode = 2
        elseif params.mode == 2 # SIMO
            global p_mode = 1
        elseif params.mode == 3 # MIMO
            global p_mode = 1.38
        end

        if close_val_lat_lon[2] < 11
            close_val_lat_lon = close_val_lat_lon .+ 10
        elseif close_val_lat_lon[2] > (length(orbit_time_all) - 10)        
            close_val_lat_lon = close_val_lat_lon .- 10
        end

        # Check consistency of input parameters
        paramsIsValid = UserParameters.validateInputParams(params)

        # Compute orbits time, position, and velocity
        global orbit_time = orbit_time_all[close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_pos = orbit_pos_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_vel = orbit_vel_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]

        global SAR_start_time = orbit_time_all[close_val_lat_lon[2]] - (params.SAR_duration / 2)

        ref_plat = 1 #incicate the reference platform
        bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,1:end,:], orbit_vel[:,1:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, params.left_right_look,  ref_plat)

        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bnorm) ./ 1e3
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bnorm)) ./ 1e3
        global Norm_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,bnorm)) ./ 1e3
        
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bperp) ./ 1e3
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bperp)) ./ 1e3
        global Perp_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,bperp)) ./ 1e3
        
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(b_at) ./ 1e3
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,b_at)) ./ 1e3
        global Par_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,b_at)) ./ 1e3

        # Read number of platforms (todo: move into a struct)
        Np  = size(orbit_pos)[2] # number of platforms
        avg_sep = maximum(bperp)/(Np - 1) # Change this

        # theoretical resolution along-n
        range_s, range_g = Scene.lookangle_to_range(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/p_mode/ (maximum(bperp)) 

        # theoretical resolution along-track
        mu = 3.986004418e14
        sc_speed = sqrt(mu./(mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:])+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
        Lsa = sc_speed*params.SAR_duration + sc_speed/params.fp
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/2/Lsa

        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/p_mode/avg_sep
        
        # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
        p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, SAR_start_time, params)

        # Create target/scene location
        targets_loc, targets_ref, Nt = Scene.construct_targets_str(params); # Nt: number of targets, targets: structure array containing target locations and reflectivities
        s_loc_3xN  = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops
        #t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions
        #t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, orbit_vel, params) ## calculate avg heading from platform positions/velocities
        
        #For co-flyer cofiguration
        t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, reshape(orbit_pos[:,1,:],(size(orbit_pos)[1],1,size(orbit_pos)[3])), params) ## calculate avg heading from platform positions

        # Apply antenna pattern
        if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
            Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
        end

       try
            
        # Generate range spread function (matched filter output)
        min_range, max_range = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
        Trx = 2*(max_range-min_range)/c + 5*params.pulse_length # s duration of RX window
        Srx, MF, ft, t_rx = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
        # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

        # Generate TomoSAR raw data
        ref_range = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
        rawdata = Generate_Raw_Data.main_RSF_slowtime(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
        if params.enable_thermal_noise # adding random noise based on SNR after range (fast-time) processing
            rawdata = Error_Sources.random_noise(rawdata, params)
        end

        # Add phase error
        sync_osc_coeffs = repeat(params.sync_a_coeff_dB, Np)
        if params.enable_sync_phase_error
            rawdata = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, sync_osc_coeffs, params)
        end

        @timeit to "BPA using new function " begin
            # Process raw data to generate image
            if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
                image_3D = Process_Raw_Data.main_SAR_tomo_3D_new(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
            elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
             SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
             # for co-flyer
             #SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata[:,2:size(p_xyz)[2],:], s_xyz_3xN, p_xyz, t_rx, ref_range, params)
             image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
            end
        end


        # Define parametsr for signal processing - Beamforming , CAPON etc
        Master_platform             = ref_plat # Master platform
        #filt_len                    = 1 # Covariance matrix box filter size 3->5x5 box filter
        azimuth_lim 				= [1, Int64(ceil(length(params.s_loc_1)))]
        srange_lim 					= [1, Int64(ceil(length(params.s_loc_2)))]
        ref_hl                      = 6#Int64(ceil(length(params.s_loc_3)/2)) 
        heights_t 					= -80:params.s_loc_3[2]-params.s_loc_3[1]:80 #params.s_loc_3  #2*minimum(params.s_loc_3):params.s_loc_3[2]-params.s_loc_3[1]:2*maximum(params.s_loc_3);#params.s_loc_3

        # Modify input matrix to desired shape
        input_SP_temp               = SAR_images_3D[:,:,:,ref_hl];
        input_SP                    = permutedims(input_SP_temp,[2,1,3])

        Ns1                         = size(input_SP)[1] # Nst
        Ns2                         = size(input_SP)[2] # Np
        Ns3                         = size(input_SP)[3] # Nr

        Cov_mat                     = zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Ns2, Ns2);
        Corr_mat                    = zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Ns2, Ns2);
        Pbf                         = zeros(Ns1,Ns3,length(heights_t));
        PC                          = zeros(Ns1,Ns3,length(heights_t));

        s_xyz_3xN_2D_all            = zeros(3,Int64(ceil(length(params.s_loc_1))),Int64(ceil(length(params.s_loc_2))),Int64(ceil(length(params.s_loc_3))))
        s_xyz_3xN_2D                = zeros(3,Int64(ceil(length(params.s_loc_1))),Int64(ceil(length(params.s_loc_2))))
        l                           = 1
        for i=1:length(params.s_loc_1)
            for j=1:length(params.s_loc_2)
                for k=1:length(params.s_loc_3)
                    s_xyz_3xN_2D_all[:,i,j,k] = s_xyz_3xN[:,k] 
                    l               = l + 1
                end
            end
        end
        s_xyz_3xN_2D                = s_xyz_3xN_2D_all[:,:,:,ref_hl]

        Cov_mat3, Corr_mat3           = Data_Processing.get_covariance_correlation_matrices_new(input_SP,  Ns2, filt_len);

        steering_mat                = Data_Processing.get_steering_matrix(p_xyz, s_xyz_3xN_2D, azimuth_lim, srange_lim, heights_t, Ns2, Master_platform, params.λ, params.mode);

        Pbf                         = Data_Processing.tomo_beamforming(Cov_mat3, steering_mat, azimuth_lim, srange_lim, [size(Cov_mat)[1] size(Cov_mat)[2] size(heights_t)[1]])  
        Pbf = Pbf[:,:,end:-1:1];
        #Pbf = Pbf[:,end:-1:1,:];
        Pbf2                        = Data_Processing.tomocoordinates_to_scenecoordinates(Pbf, heights_t, params.s_loc_2, params.s_loc_3, lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] , params.left_right_look)

 
        PC                          = Data_Processing.tomo_CAPON(Cov_mat3, steering_mat, azimuth_lim, srange_lim, [size(Cov_mat)[1] size(Cov_mat)[2] size(heights_t)[1]])  
        PC = PC[:,:,end:-1:1];
        #PC = PC[:,end:-1:1,:];
        PC2                         = Data_Processing.tomocoordinates_to_scenecoordinates(PC, heights_t, params.s_loc_2, params.s_loc_3, lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] , params.left_right_look)

        plot_idx 					= [Int64(ceil(length(params.s_loc_1)/2)),61,81] #Int64(ceil(length(params.s_loc_3)/2))]   Int64(ceil(length(params.s_loc_2)/2))

	temp3 = 6
        if length(params.t_loc_3)!=1
            #=
            val_max,ind_max         = findmax(abs.(image_3D))
		    norm_BPA_data           = abs.(image_3D) ./ val_max
		    plot_var_op_bpa         = (norm_BPA_data[plot_idx[1],plot_idx[2],temp3:temp3+50])
            plot_var_op_bpa         = reshape(plot_var_op_bpa,length(plot_var_op_bpa),1)
            plot_var_ip             = targets_ref ./ maximum(targets_ref)
            plot_var_op_bpa         = plot_var_op_bpa[1:length(plot_var_ip)]
            plot_var_ip             = plot_var_ip[1:length(plot_var_op_bpa)]
            val_t                   = (findmin(abs.(targets_loc[3,:] .- ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))))[2]
		    pks_ip, vals_ip         = findmaxima(plot_var_ip[1:val_t])
            pks_op_bpa, vals_op_bpa = findmaxima(plot_var_op_bpa[1:val_t])

            val_max,ind_max         = findmax(abs.(Pbf2))
		    norm_Beamforming_data           = abs.(Pbf2) ./ val_max
		    plot_var_op_beamforming         = (norm_Beamforming_data[plot_idx[1],plot_idx[2],plot_idx[3]:plot_idx[3]+50])
            plot_var_op_beamforming         = reshape(plot_var_op_beamforming,length(plot_var_op_beamforming),1)
            plot_var_op_beamforming         = plot_var_op_beamforming[1:length(plot_var_ip)]
            pks_op_beamforming, vals_op_beamforming = findmaxima(plot_var_op_beamforming[1:val_t])

            val_max,ind_max         = findmax(abs.(PC2))
		    norm_Capon_data         = abs.(PC2) ./ val_max
		    plot_var_op_capon         = (norm_Capon_data[plot_idx[1],plot_idx[2],plot_idx[3]:plot_idx[3]+50])
            plot_var_op_capon         = reshape(plot_var_op_capon,length(plot_var_op_capon),1)
            plot_var_op_capon         = plot_var_op_capon[1:length(plot_var_ip)]
            pks_op_capon, vals_op_capon = findmaxima(plot_var_op_capon[1:val_t])
            =#

            plot_var_ip             = ( targets_ref ./ maximum(targets_ref) )[:]
            pks_ip, vals_ip         = findmaxima(plot_var_ip)

		    plot_var_op_bpa         = image_3D[plot_idx[1],plot_idx[2],temp3:temp3+50]
            plot_var_op_bpa         = plot_var_op_bpa[1:length(plot_var_ip)]
            plot_var_op_bpa         = plot_var_op_bpa ./ maximum(plot_var_op_bpa)
            pks_op_bpa, vals_op_bpa = findmaxima(plot_var_op_bpa)

            plot_var_op_beamforming = Pbf2[plot_idx[1],plot_idx[2],plot_idx[3]:plot_idx[3]+50]
            plot_var_op_beamforming = abs.(plot_var_op_beamforming[1:length(plot_var_ip)])
            plot_var_op_beamforming = plot_var_op_beamforming ./ maximum(plot_var_op_beamforming) 
            pks_op_beamforming, vals_op_beamforming = findmaxima(plot_var_op_beamforming)

            plot_var_op_capon       = PC2[plot_idx[1],plot_idx[2],plot_idx[3]:plot_idx[3]+50]
            plot_var_op_capon       = abs.(plot_var_op_capon[1:length(plot_var_ip)])
            plot_var_op_capon       = plot_var_op_capon ./ maximum(plot_var_op_capon) 
            pks_op_capon, vals_op_capon = findmaxima(plot_var_op_capon)

            global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = cor((plot_var_ip),plot_var_op_bpa)[1] #Correlation
		    #global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =  sqrt( (sum( ((plot_var_ip)-plot_var_op_bpa).^2 )) / length(targets_loc[3,:] ) ) #RMSE
            global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =  Global_Scale_Support.compute_nrmse(plot_var_op_bpa, plot_var_ip, "none")  #RMSE
            global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = length(pks_ip) #Total output peaks
		    global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = length(pks_op_bpa) #Total output peaks

            global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = cor((plot_var_ip),plot_var_op_beamforming)[1] #Correlation
		    global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = Global_Scale_Support.compute_nrmse(plot_var_op_beamforming, plot_var_ip, "none") #RMSE
            global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = length(pks_ip) #Total output peaks
		    global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = length(pks_op_beamforming) #Total output peaks

            global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = cor((plot_var_ip),plot_var_op_capon)[1] #Correlation
		    global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = Global_Scale_Support.compute_nrmse(plot_var_op_capon, plot_var_ip, "none") #RMSE 
            global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = length(pks_ip) #Total output peaks
		    global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = length(pks_op_capon) #Total output peaks

#=
            display(plot((plot_var_ip),targets_loc[3,:] ,xlabel="Vertical profile normalized radar intensity",
            ylabel="Height (m)",title = "Input profile" , xlim=(0,1),linewidth=2,legendfont=font(15),legend=:topleft,
			xtickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, label="Input profile")) 

             display(plot!(plot_var_op_bpa, targets_loc[3,:] ,xlabel="Vertical profile normalized radar intensity",
            ylabel="Height (m)",title = "BPA output profile",xlim=(0,1),linewidth=2,legendfont=font(15), 
            xtickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, label="Output profile - Back projection"))

            display(plot!(plot_var_op_beamforming, targets_loc[3,:] ,xlabel="Vertical profile normalized radar intensity",
            ylabel="Height (m)",title = "BPA output profile",xlim=(0,1),linewidth=2,legendfont=font(15),  
            xtickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, label="Output profile - Beamforming"))

            display(plot!(plot_var_op_capon, targets_loc[3,:] ,xlabel="Vertical profile normalized radar intensity",
            ylabel="Height (m)",title = "BPA output profile",xlim=(0,1),linewidth=2,legendfont=font(10), 
            xtickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, label="Output profile - Capon"))
=#
		end
		
        GC.gc()
        catch
            continue
        end
    end
end

end

to

@save "../Outputs/output_gs_study_res_run_062023_100m_5f_514_5plat_5proc_profiles.jld" Geo_location Output_stat_bpa Output_stat_beamforming Output_stat_capon Canopy_heights orbit_time_all orbit_pos_all orbit_vel_all lookang_all Orbit_index Norm_baseline_max Norm_baseline_min Norm_baseline_mean Perp_baseline_max Perp_baseline_min Perp_baseline_mean Par_baseline_max Par_baseline_min Par_baseline_mean res_theory_n res_theory_s amb_H amb_N slnt_range to 

[rmprocs(p) for p in workers()]



to

