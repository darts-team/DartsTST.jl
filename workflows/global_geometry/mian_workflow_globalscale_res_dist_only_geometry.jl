using Distributed
addprocs(8) 

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


c = 299792458 #TODO does not work without redefining c here
earth_radius = 6378.137e3 # Earth semi-major axis at equator

const to = TimerOutput()

@timeit to "Reading L3 file for canopy heights " begin

#Read canopy heights
#filepath_GEDIL3 = "/u/epstein-z0/wblr/joshil/DARTS/GEDI_Data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
filepath_GEDIL3 = "/Users/joshil/Documents/GEDI_Data/GEDI_L3_LandSurface_Metrics_V2_1952/data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
grid_res        = 100;
Canopy_heights, Geo_location, size_row, size_col = Global_Scale_Support.read_GEDI_L3_data(filepath_GEDIL3, grid_res)

GC.gc()

end

@timeit to "Initialization " begin
# Define outout variables
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

global res_theory_n    = SharedArray(zeros(size_row,size_col,4))
global res_theory_s    = SharedArray(zeros(size_row,size_col,4))

global amb_H                = SharedArray(zeros(size_row,size_col,4))
global amb_N                = SharedArray(zeros(size_row,size_col,4))
global slnt_range           = SharedArray(zeros(size_row,size_col,4))

end

@timeit to "Read orbit file and get orbits " begin
#region_xlims = [50,110]
#region_ylims = [15,55]
#region_xlims = 59:61
#region_ylims = 17:19
region_xlims = 1:347 #265:506# 1:347#1:868 
region_ylims = 1:146 #85:271#1:146#1:366
#region_xlims = 530:1150 #530:1150
#region_ylims = 170:540 #170:600
#region_xlims        = 80:81
#region_ylims        = 37:37

lat_lon_idx         = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

#orbit_dataset       = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/NISAR_orbit_coflier_p5_lag4_theta15_06152023_3.nc") # orbit_output_06132023_1.nc#
#orbit_dataset       = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/orbit_output_06152023_3.nc")
#orbit_dataset       = Dataset("/Users/joshil/Documents/Orbits/Outputs/06152023/3/orbit_output_06152023_3.nc") # "orbit_output_04052023.nc") # Read orbits data in NetCDF format

orbit_dataset       = Dataset("/Users/joshil/Documents/Orbits/Outputs/08282023/1/NISAR_orbit_coflier_p4_lag3_theta15_08282023_1.nc") # Read orbits data in NetCDF format
#orbit_dataset       = Dataset("/Users/joshil/Documents/Orbits/Outputs/08282023/1/orbit_output_08282023_1.nc") # Read orbits data in NetCDF format


global mast_plat            = 1
flag_plat           = 1 #descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, mast_plat, flag_plat);

end

global mask = SharedArray(zeros(size(Geo_location)[1],size(Geo_location)[2]) )

lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
itp = LinearInterpolation((lon, lat), data)
for i=1:size(Geo_location)[1]
    for j=1:size(Geo_location)[2]
        mask[i,j] = itp(Geo_location[i,j,2], Geo_location[i,j,1]) 

    end
end

#=lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
global data2= data

global lat2 = lat' .* ones(length(lon))
global lon2  = ones(length(lat))' .* lon

global Lats_p = Geo_location.A[1,:,1]
global Lons_p = Geo_location.A[:,1,2]
global mask = SharedArray(zeros(length(Lons_p),length(Lats_p)) )
=#
#=
@sync @distributed for i1 =1:size(lat_lon_idx,1)   
    
        close_val_lat_lon   = findmin(abs.(lat2.-Lats_p[lat_lon_idx[i1,2]]) + abs.(lon2.-Lons_p[lat_lon_idx[i1,1]]))
        if data2[close_val_lat_lon[2]]>0
            global mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] = 1
        end
end
=#

Canopy_heights_orig = Canopy_heights


global Geo_location2 = zeros(size_row,size_col,2);
global Geo_location2[:,:,1] = Geo_location[:,:,1];
global Geo_location2[:,:,2] = Geo_location[:,:,2];



@timeit to "Processing loop over all pixels1 " begin

@sync @distributed for i1 = 1:size(lat_lon_idx,1)   
    
    if mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] <1
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
    
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN
        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1:4] .= NaN    
    else

        global close_val_lat_lon     = Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

        if close_val_lat_lon[2] < 11
		    close_val_lat_lon       = close_val_lat_lon .+ 10
	    elseif close_val_lat_lon[2] > (length(orbit_time_all) - 10)        
		    close_val_lat_lon       = close_val_lat_lon .- 10
	    end

        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = close_val_lat_lon[2]

        slrng_temp                  = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]]; 0]))
        lookang_tx                  = Scene.slantrange_to_lookangle(earth_radius,slrng_temp,orbit_pos_geo_all[:,close_val_lat_lon[2]][3],0.0)[2]

        global params               = UserParameters.inputParameters(look_angle = lookang_tx)
        # Check consistency of input parameters
        paramsIsValid               = UserParameters.validateInputParams(params)

        # Compute orbits time, position, and velocity
        global orbit_time           = orbit_time_all[close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_pos            = orbit_pos_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_vel            = orbit_vel_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]

        global SAR_start_time       = orbit_time_all[close_val_lat_lon[2]] - (params.SAR_duration / 2)

        # Read number of platforms (todo: move into a struct)
        Np                          = size(orbit_pos)[2] # number of platforms

        # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
        p_xyz, Nst, slow_time      = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, SAR_start_time, params)

        # Create target/scene location
        targets_loc, targets_ref, Nt = Scene.construct_targets_str(params) # Nt: number of targets, targets: structure array containing target locations and reflectivities
         s_loc_3xN                   = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops
        #t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions

        #For co-flyer cofiguration
        t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, reshape(orbit_pos[:,1,:],(size(orbit_pos)[1],1,size(orbit_pos)[3])), params) ## calculate avg heading from platform positions

        # Generate range spread function (matched filter output)
         min_range, max_range = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)

        ## Therorectical computations 
        # theoretical resolution
        if params.mode == 1 # SAR
            global p_mode           = 2
        elseif params.mode == 2 # SIMO
            global p_mode           = 1
        elseif params.mode == 3 # MIMO
             global p_mode           = 1.38
        end 

        # Generate TomoSAR raw data
        ref_range1 = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) 
	    global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = ref_range1

        ref_range2 = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(p_xyz[:,1,:], dims=3)) 
	    global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2]  = ref_range2

        ref_range3 = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz[:,2:end,:],dims=2), dims=3)) 
	    global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3]  = ref_range3

        ref_range4 = ( ref_range2 + ref_range3 ) / 2
	    global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4]  = ref_range4


        if params.processing_mode == 1
            ht_la = Geometry.xyz_to_geo(mean(mean(p_xyz[:,:,:],dims=2),dims=3))[3]
            global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = Scene.slantrange_to_lookangle(earth_radius,ref_range1,ht_la,0.0)[2]    
        elseif params.processing_mode == 2
            ht_la = Geometry.xyz_to_geo(mean(mean(p_xyz[:,2:end,:],dims=2),dims=3))[3]
            global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = Scene.slantrange_to_lookangle(earth_radius,ref_range3,ht_la,0.0)[2]
        end

        if params.processing_mode == 1
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,1:end,:], orbit_vel[:,1:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, params.left_right_look,  1)
            avg_sep = maximum(bperp)/(Np - 1) # Change this
            ref_plat = 1 #incicate the reference platform
        elseif params.processing_mode == 2
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,2:end,:], orbit_vel[:,2:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, params.left_right_look,  1)
            avg_sep = maximum(bperp)/(Np - 2) # Change this
            ref_plat = 2 #incicate the reference platform
        end

        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bnorm) ./ 1e3
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bnorm)) ./ 1e3
        
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bperp) ./ 1e3
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bperp)) ./ 1e3
        
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(b_at) ./ 1e3
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,b_at)) ./ 1e3

        # theoretical resolution along-n
        #range_s, range_g = Scene.lookangle_to_range(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
        
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*ref_range1/p_mode/  (maximum(bperp) + avg_sep) 
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = (c/params.fc)*ref_range2/p_mode/  (maximum(bperp) + avg_sep) 
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = (c/params.fc)*ref_range3/p_mode/  (maximum(bperp) + avg_sep) 
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = (c/params.fc)*ref_range4/p_mode/  (maximum(bperp) + avg_sep) 

        # theoretical resolution along-track
        mu = 3.986004418e14
        sc_speed = sqrt(mu./(mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:])+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
        Lsa = sc_speed*params.SAR_duration + sc_speed/params.fp
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*ref_range1/2/Lsa
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = (c/params.fc)*ref_range2/2/Lsa
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = (c/params.fc)*ref_range3/2/Lsa
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = (c/params.fc)*ref_range4/2/Lsa


        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*ref_range1/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = (c/params.fc)*ref_range2/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = (c/params.fc)*ref_range3/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = (c/params.fc)*ref_range4/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])

        
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*ref_range1/p_mode/avg_sep
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = (c/params.fc)*ref_range2/p_mode/avg_sep
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = (c/params.fc)*ref_range3/p_mode/avg_sep
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = (c/params.fc)*ref_range4/p_mode/avg_sep

    end

end

end

print(to)

@save "../Outputs/outputs_geometry_single_target_configuration_4.jld" Geo_location orbit_time_all orbit_pos_all orbit_vel_all lookang_all Orbit_index Norm_baseline_max Norm_baseline_min Norm_baseline_mean Perp_baseline_max Perp_baseline_min Perp_baseline_mean Par_baseline_max Par_baseline_min Par_baseline_mean res_theory_n res_theory_s amb_H amb_N slnt_range to

[rmprocs(p) for p in workers()]
