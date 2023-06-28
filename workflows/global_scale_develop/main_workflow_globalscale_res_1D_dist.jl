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
@everywhere include("../modules/tomographic_ISLR.jl")

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
mu = 3.986004418e14

const to = TimerOutput()

@timeit to "Reading L3 file for canopy heights " begin

#Read canopy heights
#ßfilepath_GEDIL3 = "/u/epstein-z0/wblr/joshil/DARTS/GEDI_Data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
filepath_GEDIL3 = "/Users/joshil/Documents/GEDI_Data/GEDI_L3_LandSurface_Metrics_V2_1952/data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
grid_res        = 100;
Canopy_heights, Geo_location, size_row, size_col = Global_Scale_Support.read_GEDI_L3_data(filepath_GEDIL3, grid_res)

GC.gc()

end

@timeit to "Initialization " begin
# Define outout variables
global Output_stat_bpa          = SharedArray(zeros(size_row,size_col,4))
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

#global max_baseline_n    = SharedArray(zeros(size_row,size_col,1))
global res_theory_n    = SharedArray(zeros(size_row,size_col,1))
global res_theory_s    = SharedArray(zeros(size_row,size_col,1))

end

@timeit to "Read orbit file and get orbits " begin
#region_xlims = [50,110]
#region_ylims = [15,55]
#region_xlims = 59:61
#region_ylims = 17:19
region_xlims = 1:347#93:96#1:347
region_ylims = 1:146#30:31#1:146
#region_xlims = 530:1150
#region_ylims = 170:600
#region_xlims        = 80:81
#region_ylims        = 37:37

lat_lon_idx         = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

orbit_dataset       = Dataset("inputs/NISAR_orbit_coflier_lag3_theta_15_new.nc") # "orbit_output_04052023.nc") # Read orbits data in NetCDF format
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

#=
@sync @distributed for i1 =1:size(lat_lon_idx,1)   
    
        close_val_lat_lon   = findmin(abs.(lat2.-Lats_p[lat_lon_idx[i1,2]]) + abs.(lon2.-Lons_p[lat_lon_idx[i1,1]]))
        if data2[close_val_lat_lon[2]]>0
            global mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] = 1
        end
end
=#


@timeit to "Processing loop over all pixels " begin

@sync @distributed for i1 = 1:size(lat_lon_idx,1)   
    
    close_val_lat_lon_m   = findmin(abs.(lat2.-Lats_p[lat_lon_idx[i1,2]]) + abs.(lon2.-Lons_p[lat_lon_idx[i1,1]]))
    if data2[close_val_lat_lon_m[2]]>0
        global mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] = 1
    end

    if mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] !=1
    #if isnan(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        #Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =0
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = NaN
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = NaN
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = NaN
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
    
        #global max_baseline_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
    
    else
        if isnan(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
            global Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =0.0
        end

        global close_val_lat_lon   = Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = close_val_lat_lon[2]

        slrng_temp2 = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]]; 0]))
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = acosd(orbit_pos_geo_all[:,close_val_lat_lon[2]][3] / slrng_temp2) 

        global params = UserParameters.inputParameters(look_angle = lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        # Check consistency of input parameters
        paramsIsValid = UserParameters.validateInputParams(params)
        # theoretical resolution
        if params.mode == 1 # SAR
            global p_mode = 2
        elseif params.mode == 2 # SIMO
            global p_mode = 1
        elseif params.mode == 3 # MIMO
            global p_mode = 1.38
        end
                

        # Compute orbits time, position, and velocity
        global orbit_time = orbit_time_all[close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_pos = orbit_pos_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_vel = orbit_vel_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]

        global SAR_start_time = orbit_time_all[close_val_lat_lon[2]] - (params.SAR_duration / 2)

        ref_plat = 2 #incicate the reference platform
        bperp, b_at, bnorm = Orbits.get_perp_baselines(orbit_pos[:,2:end,:], orbit_vel[:,2:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], ref_plat)

        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bnorm) ./ 1e3
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bnorm)) ./ 1e3
        #global Norm_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,bnorm)) ./ 1e3
        
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bperp) ./ 1e3
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bperp)) ./ 1e3
        #global Perp_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,bperp)) ./ 1e3
        
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(b_at) ./ 1e3
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,b_at)) ./ 1e3
        #global Par_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,b_at)) ./ 1e3

        # theoretical resolution along-n
        range_s, range_g = Scene.lookangle_to_range(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/p_mode/ maximum(bperp) 

        # theoretical resolution along-track
        sc_speed = sqrt(mu./(mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:])+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
        Lsa = sc_speed*params.SAR_duration + sc_speed/params.fp
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/2/Lsa


        scene_res       = 0.1
        freq            = c/params.λ
        max_veg_H       =  100
        altitude        = mean(Geometry.xyz_to_geo(orbit_pos[:,1,:])[3,:])
        local_slope     = 0
        look_angle      = lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]
        plat_distr      = bperp[:,:,1][2,:]
        mode            = params.mode
        try
        PSF,ISLR,scene_axis                         = tomographic_ISLR.main(mode,plat_distr,look_angle,local_slope,altitude,max_veg_H,freq,scene_res)

        image_1D_itp,scene_res_itp,scene_axis_itp   = Performance_Metrics.upsample_PSFcut(PSF[:],scene_axis[2]-scene_axis[1],100)
        res_1,res_ind_1,res_ind_2                   = Performance_Metrics.resolution_1D(image_1D_itp,scene_res_itp,4)
        PSLR_1,ISLR_1                               = Performance_Metrics.sidelobe_1D(image_1D_itp,1,res_ind_1,res_ind_2)
        loc_error_1                                 = Performance_Metrics.location_error(image_1D_itp,0,scene_axis_itp)

        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = res_1
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = PSLR_1
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = ISLR
        global Output_stat_bpa[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = loc_error_1
        catch
            continue
        end
        #println("Finished: ", i1 )
    end
end

end

to



@save "../Outputs/output_gs_study_1D_run_062023_2.jld" Geo_location Output_stat_bpa  Canopy_heights orbit_time_all orbit_pos_all orbit_vel_all lookang_all Orbit_index Norm_baseline_max Norm_baseline_min Norm_baseline_mean Perp_baseline_max Perp_baseline_min Perp_baseline_mean Par_baseline_max Par_baseline_min Par_baseline_mean res_theory_n res_theory_s to

[rmprocs(p) for p in workers()]



to
