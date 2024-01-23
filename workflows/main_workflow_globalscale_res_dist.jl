using Distributed
addprocs(24) 

@everywhere DEPOT_PATH[1]="/u/epstein-z0/wblr/joshil/Julia/.julia" # for Epstein

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
@everywhere include("../modules/data_processing.jl")
@everywhere include("../modules/global_scale_support.jl")
@everywhere include("../modules/scattering.jl")


@everywhere using NCDatasets
@everywhere using Statistics
@everywhere using Parameters
@everywhere using Dates
@everywhere using StaticArrays
@everywhere using Plots
@everywhere using LinearAlgebra
@everywhere using .UserParameters
@everywhere using TimerOutputs
@everywhere using Interpolations
@everywhere using SharedArrays
@everywhere using Peaks
@everywhere using JLD2
@everywhere using GeoDatasets


c = 299792458 
earth_radius = 6378.137e3 # Earth semi-major axis at equator

const to = TimerOutput()

@timeit to "Reading L3 file for canopy heights " begin

#Read canopy heights from GEDI dataset
filepath_GEDIL3       = "/u/epstein-z0/wblr/joshil/DARTS/GEDI_Data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
#filepath_GEDIL3         = "/Users/joshil/Documents/GEDI_Data/GEDI_L3_LandSurface_Metrics_V2_1952/data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
grid_res                = 100; # resolution for simulations
Canopy_heights, Geo_location, size_row, size_col = Global_Scale_Support.read_GEDI_L3_data(filepath_GEDIL3, grid_res)

GC.gc()

end #timeit

@timeit to "Initialization " begin
# Define output variables
global Output_stat_bpa          = SharedArray(zeros(size_row,size_col,12)) #17
global Output_stat_beamforming  = SharedArray(zeros(size_row,size_col,12)) #17
global Output_stat_capon        = SharedArray(zeros(size_row,size_col,12)) #17
global Orbit_index              = SharedArray(zeros(size_row,size_col,1))
global lookang_all              = SharedArray(zeros(size_row,size_col,1))
#Norm baselines
global Norm_baseline_max        = SharedArray(zeros(size_row,size_col,1))
global Norm_baseline_min        = SharedArray(zeros(size_row,size_col,1))
global Norm_baseline_mean       = SharedArray(zeros(size_row,size_col,1))
# perp baselines
global Perp_baseline_max        = SharedArray(zeros(size_row,size_col,1))
global Perp_baseline_min        = SharedArray(zeros(size_row,size_col,1))
global Perp_baseline_mean       = SharedArray(zeros(size_row,size_col,1))
# AT baselines
global Par_baseline_max         = SharedArray(zeros(size_row,size_col,1))
global Par_baseline_min         = SharedArray(zeros(size_row,size_col,1))
global Par_baseline_mean        = SharedArray(zeros(size_row,size_col,1))

global res_theory_n             = SharedArray(zeros(size_row,size_col,1))
global res_theory_s             = SharedArray(zeros(size_row,size_col,1))

global amb_H                    = SharedArray(zeros(size_row,size_col,1))
global amb_N                    = SharedArray(zeros(size_row,size_col,1))
global slnt_range               = SharedArray(zeros(size_row,size_col,1))

end #timeit

@timeit to "Read orbit file and get orbits " begin

region_xlims            = 1:347 # 100km res 
region_ylims            = 1:146 # 100km res 

lat_lon_idx             = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

# Read orbits data in NetCDF format
orbit_dataset           = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/orbit_output_ROSEL_12032023_1.nc")
#orbit_dataset          = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/ROSE_L_Cartwheel_4sat.nc")
#orbit_dataset          = Dataset("/Users/joshil/Documents/Code/darts-simtool/inputs/orbit_output_ROSEL_12032023_1.nc") # "orbit_output_04052023.nc") 

global mast_plat        = 1 # master platform
flag_plat               = 1 # descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, mast_plat, flag_plat)

# Make sure orbit time starts at 0
if orbit_time_all[1] > 0.0
    orbit_time_all      = orbit_time_all .- orbit_time_all[1]
end

end #timeit

# Generate mask based on land/sea, simulation runs only for points on land
global mask             = SharedArray(zeros(size(Geo_location)[1],size(Geo_location)[2]) )
lon,lat,data            = GeoDatasets.landseamask(;resolution='c',grid=5)
itp                     = LinearInterpolation((lon, lat), data)
for i=1:size(Geo_location)[1]
    for j=1:size(Geo_location)[2]
        mask[i,j]       = itp(Geo_location[i,j,2], Geo_location[i,j,1]) 

    end
end

#i1=12003 # test
#i1=28100 # test

@timeit to "Processing loop over all pixels " begin

@sync @distributed for i1 = 1:size(lat_lon_idx,1)   

    if mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] <1

        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:]           .= NaN
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:]   .= NaN
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:]         .= NaN
        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]               = NaN
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]               = NaN
        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]         = NaN
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]         = NaN
        global Norm_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]        = NaN
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]         = NaN
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]         = NaN
        global Perp_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]        = NaN
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]          = NaN
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]          = NaN
        global Par_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]         = NaN
    
        #global max_baseline_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]           = NaN
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]              = NaN
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]              = NaN
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                     = NaN
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                     = NaN
        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                = NaN    

    else
        if isnan(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
            global Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]        = 0.0
        end

        # Get the orbit point
        global close_val_lat_lon   = Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = close_val_lat_lon[2]

        # Compute slant range and look angle
        slrng_temp2         = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]]; 0]))
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = Scene.slantrange_to_lookangle(earth_radius,slrng_temp2,orbit_pos_geo_all[:,close_val_lat_lon[2]][3],0.0)[2]

        global params       = UserParameters.inputParameters(look_angle = lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        # Check consistency of input parameters
        paramsIsValid       = UserParameters.validateInputParams(params)

        filt_len            = 5 # Covariance matrix box filter size

        # theoretical resolution factor
        if params.mode == 1 # SAR
            global p_mode   = 2
        elseif params.mode == 2 # SIMO
            global p_mode   = 1
        elseif params.mode == 3 # MIMO
            global p_mode   = 1.38
        end
        
        # Check boundary consitions for orbit 
	    if close_val_lat_lon[2] < 11
		    close_val_lat_lon = close_val_lat_lon .+ 10
	    elseif close_val_lat_lon[2] > (length(orbit_time_all) - 10)        
		    close_val_lat_lon = close_val_lat_lon .- 10
	    end

        # Compute orbits time, position, and velocity
        global orbit_time   = orbit_time_all[close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_pos    = orbit_pos_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_vel    = orbit_vel_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]

        global SAR_start_time = orbit_time_all[close_val_lat_lon[2]] - (params.SAR_duration / 2)

        # Read number of platforms 
        Np                  = size(orbit_pos)[2] # number of platforms

        if params.processing_mode == 1
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,1:end,:], orbit_vel[:,1:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, params.left_right_look,  1)
            avg_sep = maximum(bperp)/(Np - 1) 
            ref_plat = 1 #incicate the reference platform
        elseif params.processing_mode == 2
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,2:end,:], orbit_vel[:,2:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, params.left_right_look,  1)
            avg_sep = maximum(bperp)/(Np - 2) 
            ref_plat = 2 #incicate the reference platform
        end

        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = maximum(bnorm) ./ 1e3
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = minimum(filter(!iszero,bnorm)) ./ 1e3
        #global Norm_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = mean(filter(!iszero,bnorm)) ./ 1e3
        
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = maximum(bperp) ./ 1e3
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = minimum(filter(!iszero,bperp)) ./ 1e3
        #global Perp_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = mean(filter(!iszero,bperp)) ./ 1e3
        
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = maximum(b_at) ./ 1e3
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = minimum(filter(!iszero,b_at)) ./ 1e3
        #global Par_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]    = mean(filter(!iszero,b_at)) ./ 1e3

        # theoretical resolution along-n
        range_s, range_g = Scene.lookangle_to_range(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]          = (c/params.fc)*range_s/p_mode/  (maximum(bperp) + avg_sep) 

        # theoretical resolution along-track
        mu              = 3.986004418e14
        sc_speed        = sqrt(mu./(mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:])+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
        Lsa             = sc_speed*params.SAR_duration + sc_speed/params.fp
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]          = (c/params.fc)*range_s/2/Lsa

        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                 = (c/params.fc)*range_s/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                 = (c/params.fc)*range_s/p_mode/avg_sep

        # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
         p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, SAR_start_time, params)

        # Create target/scene location
        targets_loc, targets_ref, Nt = Scene.construct_targets_str(params) # Nt: number of targets, targets: structure array containing target locations and reflectivities
         s_loc_3xN  = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops
        
        # For closed form configuration
         #t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions

        # Target location based only on the reference platform
        t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, reshape(orbit_pos[:,ref_plat,:],(size(orbit_pos)[1],1,size(orbit_pos)[3])), params) ## calculate avg heading from platform positions

        # Apply antenna pattern
        if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
            Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
        end
 
	try 

        # Generate range spread function (matched filter output)
         min_range, max_range   = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
         Trx                    = 2*(max_range-min_range)/c + 50*params.pulse_length # s duration of RX window #100 for co-flyer
         Srx, MF, ft, t_rx      = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
        # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

        # Generate TomoSAR raw data
         ref_range              = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
         rawdata                = Generate_Raw_Data.main_RSF_slowtime_surfaceBRCS(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
        if params.enable_thermal_noise # adding random noise based on SNR after range (fast-time) processing
             rawdata            = Error_Sources.random_noise(rawdata, params)
        end

        # Add phase error
         sync_osc_coeffs        = repeat(params.sync_a_coeff_dB, Np)
        if params.enable_sync_phase_error
             rawdata            = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, t_xyz_3xN, sync_osc_coeffs, params)
        end

        # Process raw data to generate image
        if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
             image_3D           = Process_Raw_Data.main_SAR_tomo_3D_new(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
        elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
            if params.processing_mode == 1
                SAR_images_3D   = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
            elseif params.processing_mode == 2
                # for co-flyer
                SAR_images_3D   = Process_Raw_Data.SAR_processing(rawdata[:,2:size(p_xyz)[2],:], s_xyz_3xN, p_xyz, t_rx, ref_range, params)
            end
            image_3D            = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
        end

        image_3D                = Data_Processing.average_2D_data(image_3D, filt_len)
	

        # Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
        if params.left_right_look == "left";PSFcutdir=-1;elseif params.left_right_look == "right";PSFcutdir=1;end
        params.PSF_direction[3] = PSFcutdir*params.PSF_direction[3]
        scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)

        # Calculate point target performance metrics
        bpa_resolutions, bpa_PSLRs, bpa_ISLRs, bpa_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
	
        # Define parameters for signal processing - Beamforming , CAPON 
        #filt_len                    = 1  #3->5x5 box filter
        azimuth_lim 				= [1, Int64(ceil(length(params.s_loc_1)))]
        srange_lim 					= [1, Int64(ceil(length(params.s_loc_2)))]
        ref_hl                      = Int64(ceil(length(params.s_loc_3)/2))
        heights_t 					= params.s_loc_3 #2*minimum(params.s_loc_3):params.s_loc_3[2]-params.s_loc_3[1]:2*maximum(params.s_loc_3);#params.s_loc_3

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

        # Get covariance matrix
        Cov_mat3, Corr_mat3         = Data_Processing.get_covariance_correlation_matrices_new(input_SP,  Ns2, filt_len);

        # Get steering matrix
        steering_mat                = Data_Processing.get_steering_matrix(p_xyz, s_xyz_3xN_2D, azimuth_lim, srange_lim, heights_t, Ns2, ref_plat, params.λ, params.mode, params.processing_mode);

        Pbf                         = Data_Processing.tomo_beamforming(Cov_mat3, steering_mat, azimuth_lim, srange_lim, [size(Cov_mat)[1] size(Cov_mat)[2] size(heights_t)[1]])
        #Pbf = Pbf[:,:,end:-1:1];
        ##Pbf = Pbf[:,end:-1:1,:];
	    Pbf2                        = Data_Processing.tomocoordinates_to_scenecoordinates(Pbf, heights_t, params.s_loc_2, params.s_loc_3, params.look_angle, params.left_right_look, mean(orbit_pos_geo_all[3,:]))

 
        PC                          = Data_Processing.tomo_CAPON(Cov_mat3, steering_mat, azimuth_lim, srange_lim, [size(Cov_mat)[1] size(Cov_mat)[2] size(heights_t)[1]])
        #PC = PC[:,:,end:-1:1];
        ##PC = PC[:,end:-1:1,:];
	    PC2                         = Data_Processing.tomocoordinates_to_scenecoordinates(PC, heights_t, params.s_loc_2, params.s_loc_3, params.look_angle, params.left_right_look, mean(orbit_pos_geo_all[3,:]))

        # Calculate point target performance metrics - beamforming
        scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(sqrt.(Pbf2), params)
        Beamforming_resolutions, Beamforming_PSLRs, Beamforming_ISLRs, Beamforming_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
   
        # Calculate point target performance metrics - CAPON
        scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(sqrt.(PC2), params)
        CAPON_resolutions, CAPON_PSLRs, CAPON_ISLRs, CAPON_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)

        # Calculate point target performance metrics along tomo axis
        params = UserParameters.inputParameters(
            look_angle = lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1],
            PSF_cuts = 2, 
            PSF_direction = [0 1 -tand(params.inc_angle)]
        )

        # Along tomography axis
        if params.left_right_look == "left";PSFcutdir=-1;elseif params.left_right_look == "right";PSFcutdir=1;end
        params.PSF_direction[3]=PSFcutdir*params.PSF_direction[3]

        image_1D_1, scene_axis11, scene_axis22, scene_axis33 = Scene.obtain_1D_slice_tilted(image_3D, params.s_loc_1, params.s_loc_2, params.s_loc_3, params.PSF_direction)
        scene_res=((scene_axis11[2]-scene_axis11[1])^2+(scene_axis22[2]-scene_axis22[1])^2+(scene_axis33[2]-scene_axis33[1])^2)^0.5 # scene resolution along the PSF direction
        BPA_tomo_resolutions, BPA_tomo_PSLRs, BPA_tomo_ISLRs, BPA_tomo_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)

        image_1D_1, scene_axis11, scene_axis22, scene_axis33 = Scene.obtain_1D_slice_tilted(sqrt.(Pbf2), params.s_loc_1, params.s_loc_2, params.s_loc_3, params.PSF_direction)
        scene_res=((scene_axis11[2]-scene_axis11[1])^2+(scene_axis22[2]-scene_axis22[1])^2+(scene_axis33[2]-scene_axis33[1])^2)^0.5 # scene resolution along the PSF direction
        Beamforming_tomo_resolutions, Beamforming_tomo_PSLRs, Beamforming_tomo_ISLRs, Beamforming_tomo_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
  
        image_1D_1, scene_axis11, scene_axis22, scene_axis33 = Scene.obtain_1D_slice_tilted(sqrt.(PC2), params.s_loc_1, params.s_loc_2, params.s_loc_3, params.PSF_direction)
        scene_res=((scene_axis11[2]-scene_axis11[1])^2+(scene_axis22[2]-scene_axis22[1])^2+(scene_axis33[2]-scene_axis33[1])^2)^0.5 # scene resolution along the PSF direction
        CAPON_tomo_resolutions, CAPON_tomo_PSLRs, CAPON_tomo_ISLRs, CAPON_tomo_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)


	    global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                    = ref_range

        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]               = bpa_resolutions[1]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2]               = bpa_resolutions[2]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3]               = bpa_resolutions[3]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4]               = BPA_tomo_resolutions
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],5]               = bpa_PSLRs[1]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],6]               = bpa_PSLRs[2]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],7]               = bpa_PSLRs[3]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],8]               = BPA_tomo_PSLRs
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],9]               = bpa_ISLRs[1]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],10]              = bpa_ISLRs[2]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],11]              = bpa_ISLRs[3]
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],12]              = BPA_tomo_ISLRs
        #global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],13]              = bpa_loc_errors[1]
        #global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],14]              = bpa_loc_errors[2]
        #global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],15]              = bpa_loc_errors[3]
        #global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],16]              = BPA_tomo_loc_errors
        #global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],17]              = round(maximum(20*log10.(image_3D)),digits=2)


        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]       = Beamforming_resolutions[1]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2]       = Beamforming_resolutions[2]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3]       = Beamforming_resolutions[3]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4]       = Beamforming_tomo_resolutions
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],5]       = Beamforming_PSLRs[1]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],6]       = Beamforming_PSLRs[2]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],7]       = Beamforming_PSLRs[3]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],8]       = Beamforming_tomo_PSLRs
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],9]       = Beamforming_ISLRs[1]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],10]      = Beamforming_ISLRs[2]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],11]      = Beamforming_ISLRs[3]
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],12]      = Beamforming_tomo_ISLRs
        #global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],13]      = Beamforming_loc_errors[1]
        #global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],14]      = Beamforming_loc_errors[2]
        #global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],15]      = Beamforming_loc_errors[3]
        #global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],16]      = Beamforming_tomo_loc_errors
        #global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],17]      = round(maximum(10*log10.(Pbf2)),digits=2)


        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]             = CAPON_resolutions[1]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2]             = CAPON_resolutions[2]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3]             = CAPON_resolutions[3]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4]             = CAPON_tomo_resolutions
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],5]             = CAPON_PSLRs[1]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],6]             = CAPON_PSLRs[2]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],7]             = CAPON_PSLRs[3]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],8]             = CAPON_tomo_PSLRs
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],9]             = CAPON_ISLRs[1]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],10]            = CAPON_ISLRs[2]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],11]            = CAPON_ISLRs[3]
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],12]            = CAPON_tomo_ISLRs
        #global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],13]            = CAPON_loc_errors[1]
        #global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],14]            = CAPON_loc_errors[2]
        #global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],15]            = CAPON_loc_errors[3]
        #global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],16]            = CAPON_tomo_loc_errors
        #global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],17]            = round(maximum(10*log10.(PC2)),digits=2)
	println("Finished: ", i1 )

	catch
		continue
	end
    end
end

end

print(to)

@save "../Outputs/output_gs_study_res_run_122023_100m_5f_102_4plat_3proc_test.jld" Geo_location Output_stat_bpa Output_stat_beamforming Output_stat_capon Canopy_heights orbit_time_all orbit_pos_all orbit_vel_all lookang_all Orbit_index Norm_baseline_max Norm_baseline_min Norm_baseline_mean Perp_baseline_max Perp_baseline_min Perp_baseline_mean Par_baseline_max Par_baseline_min Par_baseline_mean res_theory_n res_theory_s amb_H amb_N slnt_range to

[rmprocs(p) for p in workers()]

to

#=
include("../modules/data_plotting.jl")
norm_flag  					= 1
plot_idx 					= [Int64(ceil(length(params.s_loc_1)/2)),Int64(ceil(length(params.s_loc_2)/2)),Int64(ceil(length(params.s_loc_3)/2))]

figsavepath = "/Users/joshil/Documents/Code/analysis_loop/temp16/"

Data_Plotting.plot_tomo_output(image_3D, params,norm_flag, plot_idx, figsavepath, 20, "Back-projection output \n Normalized Intensity [dB]", "_avg")

Data_Plotting.plot_tomo_output(Pbf, params,norm_flag, plot_idx,figsavepath, 10,  "Beamforming output \n Normalized Intensity [dB]", "_orig", heights_t);

Data_Plotting.plot_tomo_output(Pbf2, params,norm_flag, plot_idx,figsavepath,10,  "Beamforming output \n Normalized Intensity [dB]", "_2", heights_t);

Data_Plotting.plot_tomo_output(PC, params,norm_flag, plot_idx,figsavepath,10,  "CAPON output \n Normalized Intensity [dB] \n", "_orig", heights_t);

Data_Plotting.plot_tomo_output(PC2, params,norm_flag, plot_idx,figsavepath,10, "CAPON output  \n Normalized Intensity [dB]\n", "_2", heights_t);
=#