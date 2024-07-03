using Distributed
if Sys.islinux()
    addprocs(10) 
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
@everywhere include("../modules/scattering.jl")

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

#TODO add separate roughness parameters for within swath and outside swath
function surface_brcs(surface_case,λ,tx_incidence,tx_azimuth,rx_incidence,rx_azimuth,polarization::Int64=3,patch_area_eff::Float64 = 1.0)
    k = 2*pi/λ
    θᵥ = 0.2
    if surface_case == 1 #high roughness
        # σ = λ*10
        # l = 1
        m = 0.5   # this is about 7.5 lambda
        l = 5
        σ = m *l / sqrt(2)
        σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = Scattering.BRCS_KA.(λ,tx_incidence,tx_azimuth,rx_incidence,rx_azimuth,σ,l,θᵥ)
    elseif surface_case == 2 #moderate roughness
        m = 0.25
        l = 1
        σ = m *l / sqrt(2)
        # σ = λ
        # l = 1
        σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = Scattering.BRCS_KA.(λ,tx_incidence,tx_azimuth,rx_incidence,rx_azimuth,σ,l,θᵥ)
    elseif surface_case == 3 #low roughness
        # kσ = 0.1
        # kl = 1
        kσ = 0.3
        kl = 3
        
        l = kl /k
        σ = kσ /k
        σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = Scattering.BRCS_SPM_tsang.(λ,tx_incidence,tx_azimuth,rx_incidence,rx_azimuth,l,σ,θᵥ)
        #
        if abs(tx_incidence) < 1.0 # if near nadir
            sigma_coherent = Scattering.RCS_coherent.(σ,θᵥ,λ,tx_incidence,tx_azimuth,rx_incidence,rx_azimuth)
            # tx_range = height
            # rx_range = height# currently assuming all platforms at same height. TODO fix assumption based on inputs
            # fresnel_minor_axis = sqrt.(λ*tx_range.*rx_range./(tx_range.+rx_range))
            # fresnel_major_axis = fresnel_minor_axis./cosd.(rx_incidence)
            # fresnel_area = fresnel_minor_axis .* fresnel_major_axis .* π
            fresnel_area = 1;
            σʳ_vh = σʳ_vh + sigma_coherent[1] ./ fresnel_area
            σʳ_hv = σʳ_hv + sigma_coherent[2] ./ fresnel_area
            σʳ_vv = σʳ_vv + sigma_coherent[3] ./ fresnel_area
            σʳ_hh = σʳ_hh + sigma_coherent[4] ./ fresnel_area
        end
    elseif surface_case == 4 #use Ulaby values - note, is backscatter direction only. used for comparions with monostatic codes
        fc = c/λ
        
        terrain = "soil" #fixed for now
        σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = Scattering.Ulaby_book_terrain_backscatter_values(rx_incidence,fc,terrain)
    
    end#if case

    if polarization == 1 #vh
        brcs = σʳ_vh .* patch_area_eff
    elseif polarization == 2 # hv
        brcs = σʳ_hv .* patch_area_eff
    elseif polarization == 3 # vv
        brcs = σʳ_vv .* patch_area_eff
    elseif polarization == 4 # hh
        brcs = σʳ_hh .* patch_area_eff
    end    

    return brcs #using the effective patch area = 1 (the default) returns the nbrcs

end#function

# Function to convert units of difference in latitude, degrees to meters
function deg_to_m_lat(deg_value)
    delta_lat   = deg_value*pi*earth_radius/180
    return delta_lat
end

# Function to convert units of difference in longitude, degrees to meters
# Depends on the latitude 
function deg_to_m_lon(deg_value, lat=0)
    delta_lon   = (deg_value*(pi*earth_radius * cos(pi*lat/180))/180)
    return delta_lon
end


for monte_carlo_idx = 1:1 # Just running one simulation currently

    @timeit to "Overall time" begin

	platform_ref_point      = [34.207;-121.83;697.5e3]
	platform_heading_ip     = 0.0
	look_angle_ip 		= 30.0  
	lat_res                	= 0.000033
	lon_res                 = 0.000033
	lat_extent              = 0.000333
	lon_extent              = 0.05
        target_mode             = 2     # 1: target fixed in center, 2: Distributed target, 3: Distributed target with 1 dominant scatterer
        num_targ_vol            = 3     # number of targets in each voxel

	ground_range_ini        = Scene.lookangle_to_range(look_angle_ip,platform_ref_point[3],0.0, earth_radius)[2]
	displacement_N          = ground_range_ini .* cosd(platform_heading_ip + 90) 
	displacement_E          = ground_range_ini .* sind(platform_heading_ip + 90)
	trg_ref_lat             = platform_ref_point[1] + (displacement_N/earth_radius*180/pi)
	trg_ref_lon             = platform_ref_point[2] + (displacement_E/(earth_radius+cosd(platform_ref_point[1]))*180/pi) 

        t_loc_S_range           = trg_ref_lat:lat_res:(trg_ref_lat+lat_extent)-lat_res	# S-grid range
        t_loc_C_range           = trg_ref_lon:lon_res:(trg_ref_lon+lon_extent)-lon_res	# C-grid range
        t_loc_H_range           = 0 	# H-grid range

        # Read Copernicus 30m resolution DEM
        dataset                 = AG.read("/home/mlavalle/dat/nisar-dem-copernicus/EPSG4326.vrt")
        DEM_orig                = AG.getband(dataset,1)

        gt                      = AG.getgeotransform(dataset)
        ulx, uly                = gt[1], gt[4]
        x_pixel_size, y_pixel_size = gt[2], gt[6]
        num_cols                = AG.width(DEM_orig)
        num_rows                = AG.height(DEM_orig)

        Lat_vals_full           = zeros(num_rows,1);
        Lon_vals_full           = zeros(num_cols,1);
        for i=1:num_rows
            Lat_vals_full[i]    = uly + i * y_pixel_size;
        end
        for i=1:num_cols
            Lon_vals_full[i]    = ulx + i * x_pixel_size;
        end
	
	# Get DEM for the region of interest by interpolation
	DEM_region_o            = ArchGDAL.read("/home/mlavalle/dat/nisar-dem-copernicus/EPSG4326.vrt") do source
           ArchGDAL.gdalwarp([source], [
                                      "-te", "$(trg_ref_lon-lon_res)", "$(trg_ref_lat-lat_res)",
                                      "$(trg_ref_lon+lon_extent+lon_res)", "$(trg_ref_lat+lat_extent+lat_res)",
                                      "-tr", "$(lon_res)", "$(lat_res)",
                                      "-r", "bilinear"]) do warped
               DEM_data         = ArchGDAL.getband(warped, 1)

               ArchGDAL.read(DEM_data)
           end
       end
       DEM_region               = DEM_region_o[2:end-1,2:end-1]'

	DEM                     = DEM_region_o

	Geo_location_lon        = repeat(collect(trg_ref_lon-lon_res:lon_res:trg_ref_lon+lon_extent),1,length(t_loc_S_range)+2)
	Geo_location_lat        = repeat(collect(trg_ref_lat-lat_res:lat_res:trg_ref_lat+lat_extent)',length(t_loc_C_range)+2,1)


	# Compute East-west slope along longitide
	DEM_temp1                   = zeros(size(DEM)) 
	DEM_temp1[1:end-1,:,1]      = DEM[2:end,:,1] # Shift DEM 1 pixel left
	DEM_temp2                   = zeros(size(DEM)) 
	DEM_temp2[2:end,:,1]        = DEM[1:end-1,:,1] # Shift DEM 1 pixel right
	# Compute the height difference
	ht_diff_mat_lon             = DEM_temp2 .- DEM_temp1

	Geo_location_lon_lshift                 = copy(Geo_location_lon)
	Geo_location_lon_lshift[1:end-1,:,1]    = Geo_location_lon[2:end,:,1] # Shift Geo_location_lon 1 pixel left
	Geo_location_lon_lshift[end,:,1]        = 0 .* Geo_location_lon_lshift[end,:,1]

	Geo_location_lon_rshift                 = copy(Geo_location_lon)
	Geo_location_lon_rshift[2:end,:,1]      = Geo_location_lon[1:end-1,:,1] # Shift Geo_location_lon 1 pixel right
	Geo_location_lon_rshift[1,:,1]          = 0 .* Geo_location_lon_rshift[1,:,1]

	# Compute the longitude difference
	dist_diff_mat_lon           = (Geo_location_lon_rshift - Geo_location_lon_lshift)
	dist_lon                    = deg_to_m_lon.(dist_diff_mat_lon, mean(Geo_location_lat))

	# Compute slope along longitude (X)
	slope_lon                   = ht_diff_mat_lon ./ dist_lon		
	slope_lon_region  	    = slope_lon[2:end-1,2:end-1]'

	# Compute North-south slope along latitude
	DEM_temp1                   = zeros(size(DEM)) 
	DEM_temp1[:,1:end-1,1]      = DEM[:,2:end,1] # Shift DEM 1 pixel left
	DEM_temp2                   = zeros(size(DEM)) 
	DEM_temp2[:,2:end,1]        = DEM[:,1:end-1,1] # Shift DEM 1 pixel right
	# Compute the height difference
	ht_diff_mat_lat             = DEM_temp2 .- DEM_temp1

	Geo_location_lat_lshift                 = copy(Geo_location_lat) 
	Geo_location_lat_lshift[:,1:end-1,1]    = Geo_location_lat[:,2:end,1] # Shift Geo_location_lon 1 pixel left
	Geo_location_lat_lshift[:,end,1]        = 0 .* Geo_location_lat_lshift[:,end,1]

	Geo_location_lat_rshift                 = copy(Geo_location_lat)
	Geo_location_lat_rshift[:,2:end,1]      = Geo_location_lat[:,1:end-1,1] # Shift Geo_location_lon 1 pixel right
	Geo_location_lat_rshift[:,1,1]          = 0 .* Geo_location_lat_rshift[:,1,1]

	# Compute the latitude difference
	dist_diff_mat_lat           = (Geo_location_lat_rshift - Geo_location_lat_lshift)
	dist_lat                    = deg_to_m_lat.(dist_diff_mat_lat)

	# Compute slope along latitude (Y)
	slope_lat                   = ht_diff_mat_lat ./ dist_lat
	slope_lat_region 	    = slope_lat[2:end-1,2:end-1]' 

	# Compute terrain normal vector 
	#N = [-dh/dx, -dh/dy, 1] where dh/dx is the west-east slope and dh/dy is the south-north slope wrt to the geogrid. 
	#dh, dx, and dy are computed in meters using ENU coordinates (xyz).
	N                           =  zeros(size(DEM_region)[1],size(DEM_region)[2],3)
        N_all                       = zeros(size(DEM_region)[1]*size(DEM_region)[2],3)
        k = 1
	for i=1:size(DEM_region)[1]
    		for j=1:size(DEM_region)[2]
        		N[i,j,:]            = [-1*slope_lon_region[i,j,1], -1*slope_lat_region[i,j,1], 1]
                        N_all[k,:]          = [-1*slope_lon_region[i,j,1], -1*slope_lat_region[i,j,1], 1]
    		k = k + 1
                end
	end

	### Get the slope angle along lat and lon
	slope_lon_region2          = atand.(slope_lon_region)
	slope_lat_region2          = atand.(slope_lat_region)



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

            s_loc_1             = trg_ref_lat:lat_res:(trg_ref_lat+lat_extent)-lat_res, 
            s_loc_2             = trg_ref_lon:lon_res:(trg_ref_lon+lon_extent)-lon_res,
            s_loc_3             = 480.0,

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
            fp                  = 1550,
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


        #geometry computations
        rng_slope                   = 0
        slant_range_all             = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        look_angle_all              = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        incidence_angle_all         = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        Critical_baseline_all       = zeros(size(s_xyz_3xN,2),size(p_xyz,2))
        Correlation_theo_all        = zeros(size(s_xyz_3xN,2),1)
        Perp_baseline_all           = zeros(size(s_xyz_3xN,2),1)
        Vert_wavnum_all             = zeros(size(s_xyz_3xN,2),1)
	local_incidence_angle       = zeros(size(s_xyz_3xN,2),1)
	range_slope_angle           = zeros(size(s_xyz_3xN,2),1)
        targets_ref_corr            = zeros(size(s_xyz_3xN,2)*num_targ_vol,1)

        pixel_area = (lat_res*pi/180*earth_radius) * (lon_res*pi/180*(earth_radius+cosd(params.p_heading)))

        for oi = 1:size(p_xyz,2)
            mean_plats_pos              = mean(p_xyz[:,oi,:], dims=2)
            for ti = 1:size(s_xyz_3xN,2)
                    slant_range_all[ti,oi]       = Geometry.distance( mean_plats_pos  , s_xyz_3xN[:,ti] )
                    look_angle_all[ti,oi]        = Scene.slantrange_to_lookangle(earth_radius,slant_range_all[ti,oi],Geometry.xyz_to_geo(mean_plats_pos)[3],Geometry.xyz_to_geo(t_xyz_3xN[:,ti])[3])[2]
                    incidence_angle_all[ti,oi]   = Scene.lookangle_to_incangle(look_angle_all[ti,oi],Geometry.xyz_to_geo(mean_plats_pos)[3],Geometry.xyz_to_geo(t_xyz_3xN[:,ti])[3],earth_radius)
                    Critical_baseline_all[ti,oi] = params.λ * ((2*params.bandwidth)/c) * slant_range_all[ti,oi] * tand(incidence_angle_all[ti,oi] - rng_slope) / p_mode
            end
        end

        ti2 = 1
        for ti=1:size(s_xyz_3xN,2)
            bs_perp2, bs_at2, bs_norm2          = Orbits.get_perp_baselines_new(mean(p_xyz[:,:,:],dims=3), mean(v_xyz[:,:,:],dims=3), look_angle_all[ti,1], 0.0, params.left_right_look, 1)
            Perp_baseline_all[ti]               = bs_perp2[1,2,1]
            Correlation_theo_all[ti]            = 1 - (Perp_baseline_all[ti] ./ (Critical_baseline_all[ti,1]))

            Va                                  = mean(p_xyz[:, 1, :], dims=2) - s_xyz_3xN[:, ti]
            Vb                                  = mean(p_xyz[:, 2, :], dims=2) - s_xyz_3xN[:, ti]

            angle_ip                            = Data_Processing.angle_2vec(Va, Vb) * 1

            if params.mode == 1
                Vert_wavnum_all[ti]             = (4 * pi * (angle_ip * pi / 180) ) / (params.λ * sind(look_angle_all[ti,1]))
            elseif params.mode == 2
                Vert_wavnum_all[ti]             = (2 * pi * (angle_ip * pi / 180) ) / (params.λ * sind(look_angle_all[ti,1]))
            end

            plat_pt_xyz                         =  mean(mean(p_xyz,dims=2),dims=3)[:]
            look_vec_xyz                        = (plat_pt_xyz - s_xyz_3xN[:,ti])
            look_vec_xyz_norm                   = (plat_pt_xyz - s_xyz_3xN[:,ti]) / Geometry.distance( plat_pt_xyz,s_xyz_3xN[:,ti])

            Geo_location                        = Geometry.xyz_to_geo(s_xyz_3xN[:,ti])
            pegθ                                = Geo_location[1]*π/180
            pegϕ                                = Geo_location[2]*π/180
            #ENU to XYZ transformation matrix
            Menu_xyz                = [-sin(pegϕ) -sin(pegθ)*cos(pegϕ) cos(pegθ)*cos(pegϕ);
                                    cos(pegϕ) -sin(pegθ)*sin(pegϕ) cos(pegθ)*sin(pegϕ);
                                    0            cos(pegθ)             sin(pegθ)]
            #XYZ to ENU transformation matrix
            Mxyz_enu                = [-sin(pegϕ)           cos(pegϕ)             0;
                                    -sin(pegθ)*cos(pegϕ) -sin(pegθ)*sin(pegϕ)  cos(pegθ)  ;
                                    cos(pegθ)*cos(pegϕ)   cos(pegθ)*sin(pegϕ)  sin(pegθ)]

            Nxyz                                = Menu_xyz * N_all[ti,:];
            look_vec_enu                        = Mxyz_enu * look_vec_xyz

            local_incidence_angle[ti]           = Data_Processing.angle_2vec(look_vec_xyz_norm, Nxyz)

            look_direction_norm                 = sqrt( look_vec_enu[1] * look_vec_enu[1]  + look_vec_enu[2] * look_vec_enu[2])
            range_slope_angle[ti]= atand((N_all[ti,1] * (look_vec_enu[1] / look_direction_norm)) + (N_all[ti,2] * (look_vec_enu[2] / look_direction_norm)))     

            targets_ref_corr[ti2:ti2+2] .= pixel_area .*  surface_brcs(2,params.λ,local_incidence_angle[ti],params.p_heading,local_incidence_angle[ti],params.p_heading,3,1.0)
 	    ti2 = ti2+3
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
                    if targets_ref_corr[j]!=0
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
                            temp_sum .= temp_sum .+  (targets_ref_corr[j].*exp(-im*2*pi/(0.23793052222222222)*(range_tx+range_rx)).*Srx_shifted)
    
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

	#=
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
	=#

        savepath                 = "/u/epstein-z0/darts/joshil/Outputs/SC_grid/"
        if ~ispath(savepath)
            mkdir(savepath)
        end
        
        @save savepath*"Output58_10km_30la_testgeo8.jld" image_3D SAR_images_3D stat_var_all_1p stat_var_all_2p rawdata params slant_range_all look_angle_all incidence_angle_all Critical_baseline_all Correlation_theo_all Perp_baseline_all Vert_wavnum_all rng_slope s_xyz_3xN p_xyz t_rx ref_range t_loc_H t_xyz_3xN slope_lon_region slope_lat_region DEM_region N_all local_incidence_angle range_slope_angle targets_ref_corr

    end
    
    end

    [rmprocs(p) for p in workers()]

println(to)
            

