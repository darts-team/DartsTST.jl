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
@everywhere using Distributions
@everywhere using Measures

c = 299792458 
earth_radius = 6378.137e3 # Earth semi-major axis at equator

@everywhere const to = TimerOutput()

@timeit to "Reading L2 and L3 file for canopy heights and profile " begin

#Read canopy heights
if Sys.islinux()
    filepath_GEDIL3 = "/u/epstein-z0/wblr/joshil/DARTS/GEDI_Data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
else
    filepath_GEDIL3 = "/Users/joshil/Documents/GEDI_Data/GEDI_L3_LandSurface_Metrics_V2_1952/data/GEDI03_rh100_mean_2019108_2021104_002_02.tif" 
end

grid_res        = 100; # resolution for simulations
Canopy_heights_L3, Geo_location, size_row, size_col = Global_Scale_Support.read_GEDI_L3_data(filepath_GEDIL3, grid_res)

if Sys.islinux()
    @load "/u/epstein-z0/darts/joshil/GEDI/GEDI_L2/Output_GEDIL2_1year_2022_res100m.jld2"
else
    @load "/Users/joshil/Documents/GEDI_Data/Outputs_L2/Output_GEDIL2_1year_2022.jld2"
end
global Canopy_heights           = Canopy_heights_L2[:,:];

global Canopy_profiles          = Canopy_profiles_L2_itp[:,:,:]
global Canopy_profile_heights   = Canopy_heights_profile_L2_itp[:,:,:]

GC.gc()

end #timeit

@timeit to "Initialization " begin
# Define output variables
global Output_stat_bpa          = SharedArray(zeros(size_row,size_col,4))
global Output_stat_beamforming  = SharedArray(zeros(size_row,size_col,4))
global Output_stat_capon        = SharedArray(zeros(size_row,size_col,4))

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

region_xlims        = 1:347# 100km res #90:95
region_ylims        = 1:146# 100km res #30:30


lat_lon_idx         = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

# Read orbits data in NetCDF format
if Sys.islinux()
    orbit_dataset       = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/orbit_output_ROSEL_12032023_1.nc")
    #orbit_dataset       = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/ROSEL_orbit_coflier_p4_lag3_theta15_12032023_1.nc")
else
    #orbit_dataset       = Dataset("/Users/joshil/Documents/Orbits/Outputs/06152023/3/orbit_output_06152023_3.nc") # "orbit_output_04052023.nc") # Read orbits data in NetCDF format
    orbit_dataset       = Dataset("/Users/joshil/Documents/Orbits/Outputs/ROSEL-L/orbit_output_ROSEL_12032023_1.nc") 
end

global mast_plat    = 1 # master platform
flag_plat           = 1 # descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, "ECI",mast_plat, flag_plat)

end #timeit

# Generate mask based on land/sea, simulation runs only for points on land
global mask         = SharedArray(zeros(size(Geo_location)[1],size(Geo_location)[2]) )

lon,lat,data        = GeoDatasets.landseamask(;resolution='c',grid=5)
itp = LinearInterpolation((lon, lat), data)
for i=1:size(Geo_location)[1]
    for j=1:size(Geo_location)[2]
        mask[i,j]   = itp(Geo_location[i,j,2], Geo_location[i,j,1])

    end
end

i1=13608

@timeit to "Processing loop over all pixels " begin

@sync @distributed for i1 = 1:size(lat_lon_idx,1)   

    if mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] <1
        
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4]         .= NaN
        global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN
        global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4]       .= NaN
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

        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]              = NaN
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]              = NaN
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                     = NaN
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                     = NaN
        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                = NaN
        
        global Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN;
    
    else

        if isnan(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
            global Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =0.0
            global Canopy_profiles[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:] .=0.0
            global Canopy_profile_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:] = Canopy_profile_heights[94,30,:]
        end
                
        # Get the orbit point
        global close_val_lat_lon   = Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = close_val_lat_lon[2]

        # Compute slant range and look angle
        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]]; 0]))
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = Scene.slantrange_to_lookangle(earth_radius,slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1],orbit_pos_geo_all[:,close_val_lat_lon[2]][3],0.0)[2]

        # Custom reflectivity profiles
        #ref_profile_height, ref_profile_value, NoPeaks  = Global_Scale_Support.constrct_reflectivity_profile_exp(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        #interp_fn = LinearInterpolation(ref_profile_height, ref_profile_value, extrapolation_bc = Flat());
        #height_range = 0:1:40

        # Boundary on canopy heights
        if Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] > 49
            Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = 49.0
        end


        targ_loc        = Canopy_profile_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:2:length(0.0:0.5:ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))];
        targ_ref        = Canopy_profiles[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:2:length(0.0:0.5:ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))];
        targ_ref = targ_ref./maximum(targ_ref)
        #targ_ref = [1;0.3;0.4;0.7;1;0.7;0.4;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1]

        if length(targ_loc) ==1
            itp_ht = LinearInterpolation( [targ_loc[1],targ_loc[1]+1.0],  [targ_ref[1],targ_ref[1]+1.0], extrapolation_bc=Flat())
        else
            itp_ht = LinearInterpolation(targ_loc, targ_ref, extrapolation_bc=Flat())
        end


        num_targ_vol            = 1
        t_loc_S_range           = 0#-2:1:2
        t_loc_C_range           = 0#-2:1:2
        t_loc_H_range           = Canopy_profile_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:2:length(0.0:0.5:ceil(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]))];

        tgrid_midS = 1
        tgrid_midC = 1

        t_loc_S                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_C                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_loc_H                 =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
        t_ref_val               =  zeros(length(t_loc_S_range) * length(t_loc_C_range) * length(t_loc_H_range) * num_targ_vol)
    
        m=1
        for i=1:length(t_loc_S_range)
            for j= 1:length(t_loc_C_range)
                for k=1:length(t_loc_H_range)
                    for l=1:num_targ_vol
                        global t_loc_S[m] = t_loc_S_range[i] #+ (rand(Uniform(-1,1)) .* 0.5 )#(t_loc_S_range[2]-t_loc_S_range[1]/2) )
                        global t_loc_C[m] = t_loc_C_range[j] #+ (rand(Uniform(-1,1)) .* 0.5)#(t_loc_C_range[2]-t_loc_C_range[1]/2) )
                        global t_loc_H[m] = t_loc_H_range[k] #+ (rand(Uniform(-1,1)) .* 0.5)#(0.24/2))
                        global t_ref_val[m] = itp_ht(t_loc_H[m])
                        m=m+1;
                    end
                end
            end
        end

        
        
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
            look_angle = lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1],
            #target_pos_mode = "layered-grid-GEDIL2",
            #t_loc_3 = targ_loc,
            #t_ref   = targ_ref, 
            target_pos_mode = "CR",
            t_loc_1 = t_loc_S',
            t_loc_2 = t_loc_C',
            t_loc_3 = t_loc_H',
            t_ref   = t_ref_val', 
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

        # Check boundary consitions for orbit 
        if close_val_lat_lon[2] < 11
            close_val_lat_lon = close_val_lat_lon .+ 10
        elseif close_val_lat_lon[2] > (length(orbit_time_all) - 10)        
            close_val_lat_lon = close_val_lat_lon .- 10
        end

        # Compute orbits time, position, and velocity
        global orbit_time = orbit_time_all[close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_pos = orbit_pos_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_vel = orbit_vel_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]

        global SAR_start_time = orbit_time_all[close_val_lat_lon[2]] - (params.SAR_duration / 2)

        # Read number of platforms 
        Np  = size(orbit_pos)[2] # number of platforms

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
        global Norm_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]    = mean(filter(!iszero,bnorm)) ./ 1e3
        
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = maximum(bperp) ./ 1e3
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = minimum(filter(!iszero,bperp)) ./ 1e3
        global Perp_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]    = mean(filter(!iszero,bperp)) ./ 1e3
        
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = maximum(b_at) ./ 1e3
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = minimum(filter(!iszero,b_at)) ./ 1e3
        global Par_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = mean(filter(!iszero,b_at)) ./ 1e3


        # theoretical resolution along-n
        range_s, range_g = Scene.lookangle_to_range(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/p_mode/ (maximum(bperp) + avg_sep) 

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
        min_range, max_range    = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
        Trx                     = 2*(max_range-min_range)/c + 5*params.pulse_length # s duration of RX window
        Srx, MF, ft, t_rx       = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
        # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

        # Generate TomoSAR raw data
        ref_range               = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
        rawdata                 = Generate_Raw_Data.main_RSF_slowtime_perf_opt(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
        #rawdata                 = Generate_Raw_Data.main_RSF_slowtime(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
        if params.enable_thermal_noise # adding random noise based on SNR after range (fast-time) processing
            rawdata             = Error_Sources.random_noise(rawdata, params)
        end

        # Add phase error
        sync_osc_coeffs         = repeat(params.sync_a_coeff_dB, Np)
        if params.enable_sync_phase_error
            rawdata             = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, sync_osc_coeffs, params)
        end

        @timeit to "BPA using new function " begin
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
             image_3D           = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
            end
        end


        # Define parameters for signal processing - Beamforming , CAPON 
        height_res                  = params.s_loc_3[2]-params.s_loc_3[1]
        #filt_len                    = 1 # Covariance matrix box filter size 3->5x5 box filter
        azimuth_lim 				= [1, Int64(ceil(length(params.s_loc_1)))]
        srange_lim 					= [1, Int64(ceil(length(params.s_loc_2)))]
        ref_hl                      = 6 # this is hard coded change this #Int64(ceil(length(params.s_loc_3)/2)) 
        heights_t 					= -80:height_res:80 #params.s_loc_3  #2*minimum(params.s_loc_3):params.s_loc_3[2]-params.s_loc_3[1]:2*maximum(params.s_loc_3);#params.s_loc_3

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
        Cov_mat3, Corr_mat3           = Data_Processing.get_covariance_correlation_matrices_new(input_SP,  Ns2, filt_len);

         # Get steering matrix
        steering_mat                = Data_Processing.get_steering_matrix(p_xyz, s_xyz_3xN_2D, azimuth_lim, srange_lim, heights_t, Ns2, ref_plat, params.λ, params.mode, params.processing_mode);

        Pbf                         = Data_Processing.tomo_beamforming(Cov_mat3, steering_mat, azimuth_lim, srange_lim, [size(Cov_mat)[1] size(Cov_mat)[2] size(heights_t)[1]])  
        #Pbf = Pbf[:,:,end:-1:1];
        ##Pbf = Pbf[:,end:-1:1,:];
        Pbf2                        = Data_Processing.tomocoordinates_to_scenecoordinates(Pbf, heights_t, params.s_loc_2, params.s_loc_3, lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] , params.left_right_look, mean(orbit_pos_geo_all[3,:]))

 
        PC                          = Data_Processing.tomo_CAPON(Cov_mat3, steering_mat, azimuth_lim, srange_lim, [size(Cov_mat)[1] size(Cov_mat)[2] size(heights_t)[1]])  
        #PC = PC[:,:,end:-1:1];
        ##PC = PC[:,end:-1:1,:];
        PC2                         = Data_Processing.tomocoordinates_to_scenecoordinates(PC, heights_t, params.s_loc_2, params.s_loc_3, lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] , params.left_right_look, mean(orbit_pos_geo_all[3,:]))

        plot_idx 					= [Int64(ceil(length(params.s_loc_1)/2)),61,Int64(ceil(length(heights_t)/2))] #61 #Int64(ceil(length(params.s_loc_3)/2))]   Int64(ceil(length(params.s_loc_2)/2))

	    temp3 = ref_hl

        if length(params.t_loc_3)!=1

            targets_ref_ip_grid     = reshape(targets_ref,num_targ_vol,length(t_loc_H_range),length(t_loc_C_range),length(t_loc_S_range))
            targets_ref_ip_grid_mid = targets_ref_ip_grid[1,:,tgrid_midC,tgrid_midS]
            plot_var_ip             = ( targets_ref_ip_grid_mid ./ maximum(targets_ref_ip_grid_mid) )[:]
            pks_ip, vals_ip         = findmaxima(plot_var_ip)

            if (t_loc_H_range[2]-t_loc_H_range[1]) != (heights_t[2]-heights_t[1])
                itp                 = linear_interpolation(targ_loc, plot_var_ip,extrapolation_bc=0) 
                targ_loc_new        = t_loc_H_range[1]:height_res:t_loc_H_range[end]
                plot_var_ip         = itp.(targ_loc_new)
            else
                targ_loc_new        = collect(t_loc_H_range[1]:height_res:t_loc_H_range[end])
            end

            slim = 0
            clim = 2

            plot_var_op_bpa1        = image_3D[plot_idx[1]-slim:plot_idx[1]-slim,plot_idx[2]-clim:plot_idx[2]+clim,temp3:Int64(ceil(length(params.s_loc_3)))]
            plot_var_op_bpa         = mean(mean(plot_var_op_bpa1,dims=1),dims=2)[:]
            plot_var_op_bpa         = plot_var_op_bpa[1:length(t_loc_H_range)]
            plot_var_op_bpa         = plot_var_op_bpa ./ maximum(plot_var_op_bpa)

            #plot_var_op_bpa         = image_3D[plot_idx[1],plot_idx[2],temp3:Int64(ceil(length(params.s_loc_3)))]
            #plot_var_op_bpa         = plot_var_op_bpa[1:length(t_loc_H_range)]
            #plot_var_op_bpa         = plot_var_op_bpa ./ maximum(plot_var_op_bpa)
            pks_op_bpa, vals_op_bpa = findmaxima(plot_var_op_bpa)

            #plot_var_op_beamforming = Pbf2[plot_idx[1],plot_idx[2],plot_idx[3]:Int64(ceil(length(heights_t)))]
            plot_var_op_beamforming1        = Pbf2[plot_idx[1]-slim:plot_idx[1]+slim,plot_idx[2]-clim:plot_idx[2]+clim,plot_idx[3]:Int64(ceil(length(heights_t)))]
            plot_var_op_beamforming         = mean(mean(plot_var_op_beamforming1,dims=1),dims=2)[:]
            plot_var_op_beamforming = sqrt.(abs.(plot_var_op_beamforming[1:length(t_loc_H_range)]))
            plot_var_op_beamforming = plot_var_op_beamforming ./ maximum(plot_var_op_beamforming) 
            pks_op_beamforming, vals_op_beamforming = findmaxima(plot_var_op_beamforming)

            #plot_var_op_capon       = PC2[plot_idx[1],plot_idx[2],plot_idx[3]:Int64(ceil(length(heights_t)))]
            plot_var_op_capon1        = PC2[plot_idx[1]-slim:plot_idx[1]-slim,plot_idx[2]-clim:plot_idx[2]+clim,plot_idx[3]:Int64(ceil(length(heights_t)))]
            plot_var_op_capon         = mean(mean(plot_var_op_capon1,dims=1),dims=2)[:]
            plot_var_op_capon       = sqrt.(abs.(plot_var_op_capon[1:length(t_loc_H_range)]))
            plot_var_op_capon       = plot_var_op_capon ./ maximum(plot_var_op_capon) 
            pks_op_capon, vals_op_capon = findmaxima(plot_var_op_capon)

            global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2]           = cor((plot_var_ip),plot_var_op_bpa)[1] #Correlation
		    #global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]          =  sqrt( (sum( ((plot_var_ip)-plot_var_op_bpa).^2 )) / length(targets_loc[3,:] ) ) #RMSE
            global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]           =  Global_Scale_Support.compute_nrmse(plot_var_op_bpa, plot_var_ip, "none")  #RMSE
            global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3]           = length(pks_ip) #Total output peaks
		    global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4]           = length(pks_op_bpa) #Total output peaks

            global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2]   = cor((plot_var_ip),plot_var_op_beamforming)[1] #Correlation
		    global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = Global_Scale_Support.compute_nrmse(plot_var_op_beamforming, plot_var_ip, "none") #RMSE
            global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3]   = length(pks_ip) #Total output peaks
		    global Output_stat_beamforming[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4]   = length(pks_op_beamforming) #Total output peaks

            global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2]         = cor((plot_var_ip),plot_var_op_capon)[1] #Correlation
		    global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]         = Global_Scale_Support.compute_nrmse(plot_var_op_capon, plot_var_ip, "none") #RMSE 
            global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3]         = length(pks_ip) #Total output peaks
		    global Output_stat_capon[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4]         = length(pks_op_capon) #Total output peaks

 # DEBUG
            p2=(plot((plot_var_ip),targ_loc_new,xlabel="Vertical profile normalized radar intensity",
            ylabel="Height (m)",title = "Input profile",linewidth=2 ,legend=:topleft, label="Input profile")) 

            p2= (plot!(plot_var_op_bpa, targ_loc_new,xlabel="Vertical profile normalized radar intensity",
            ylabel="Height (m)",title = "BPA output profile",  label="Output profile - Back projection"))

            p2=(plot!(plot_var_op_beamforming,  targ_loc_new ,xlabel="Vertical profile normalized radar intensity",
            ylabel="Height (m)",title = "BPA output profile", label="Output profile - Beamforming"))

            p2=(plot!(plot_var_op_capon, targ_loc_new ,xlabel="Vertical profile normalized radar intensity",
            ylabel="Height (m)",title = "Profile",xlim=(0,1),legendfont=font(10), rightmargin=5mm, leftmargin=5mm,
            xtickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, label="Output profile - Capon"))

            display(p2)
            #img_path = "../../Outputs/Plots_test/"
            #savefig(p2, img_path*"Profile_"*string(Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]][1])*"_"*string(Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]][2])*".png")

		end

        GC.gc()
        catch
            #continue
        end
    end
end

end

print(to)

@save "../../Outputs/output_gs_study_res_run_122023_100m_5f_106_4plat_3proc_profiles_test.jld" Geo_location Output_stat_bpa Output_stat_beamforming Output_stat_capon Canopy_heights orbit_time_all orbit_pos_all orbit_vel_all lookang_all Orbit_index Norm_baseline_max Norm_baseline_min Norm_baseline_mean Perp_baseline_max Perp_baseline_min Perp_baseline_mean Par_baseline_max Par_baseline_min Par_baseline_mean res_theory_n res_theory_s amb_H amb_N slnt_range to 

[rmprocs(p) for p in workers()]

to

