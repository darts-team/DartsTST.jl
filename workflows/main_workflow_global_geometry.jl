using Distributed
addprocs(24) 

@everywhere include("../modules/geometry.jl")
@everywhere include("../modules/scene.jl")
@everywhere include("../modules/orbits.jl")
@everywhere include("../modules/antenna.jl")
@everywhere include("../modules/simsetup.jl")
@everywhere include("../modules/global_scale_support.jl")

@everywhere using StaticArrays
@everywhere using LinearAlgebra
@everywhere using TimerOutputs
@everywhere using SharedArrays
@everywhere using JLD2
@everywhere using GeoDatasets
@everywhere using GeoArrays
@everywhere using NCDatasets
@everywhere using Statistics


c = 299792458 #TODO does not work without redefining c here
earth_radius = 6378.137e3 # Earth semi-major axis at equator

const to = TimerOutput()

@timeit to "Reading L3 file for canopy heights " begin

#Read canopy heights
filepath_GEDIL3 = "/u/epstein-z0/wblr/joshil/DARTS/GEDI_Data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
#filepath_GEDIL3 = "/Users/joshil/Documents/GEDI_Data/GEDI_L3_LandSurface_Metrics_V2_1952/data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
grid_res        = 10;
Canopy_heights, Geo_location, size_row, size_col = Global_Scale_Support.read_GEDI_L3_data(filepath_GEDIL3, grid_res)

GC.gc()

end

@timeit to "Initialization " begin

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

region_xlims = 1:3470 
region_ylims = 1:1462 

lat_lon_idx         = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

#orbit_dataset       = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/NISAR_orbit_coflier_p5_lag4_theta15_06152023_3.nc") # orbit_output_06132023_1.nc#
orbit_dataset       = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/orbit_output_06152023_3.nc")
#orbit_dataset       = Dataset("/Users/joshil/Documents/Orbits/Outputs/06152023/3/orbit_output_06152023_3.nc") # "orbit_output_04052023.nc") # Read orbits data in NetCDF format

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

Canopy_heights_orig = Canopy_heights


@timeit to "Processing loop over all pixels " begin

@sync @distributed for i1 = 1:size(lat_lon_idx,1)   
    
    close_val_lat_lon_m   = findmin(abs.(lat2.-Lats_p[lat_lon_idx[i1,2]]) + abs.(lon2.-Lons_p[lat_lon_idx[i1,1]]))
    if data2[close_val_lat_lon_m[2]]>0
        global mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] = 1
    end

    if mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] !=1
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
    else
        if isnan(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
            global Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =0.0
        end

        global close_val_lat_lon   = Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = close_val_lat_lon[2]

        slrng_temp2 = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]]; 0]))
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = Scene.slantrange_to_lookangle(earth_radius,slrng_temp2,orbit_pos_geo_all[:,close_val_lat_lon[2]][3],0.0)[2]

        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = slrng_temp2

        p_mode = 2
        processing_mode = 1
        SAR_duration = 5
        fc = 1.25e9
        fp = 10


        
	    if close_val_lat_lon[2] < 11
		    close_val_lat_lon = close_val_lat_lon .+ 10
	    elseif close_val_lat_lon[2] > (length(orbit_time_all) - 10)        
		    close_val_lat_lon = close_val_lat_lon .- 10
	    end

        # Compute orbits time, position, and velocity
        global orbit_time = orbit_time_all[close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_pos = orbit_pos_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_vel = orbit_vel_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]

        global SAR_start_time = orbit_time_all[close_val_lat_lon[2]] - (SAR_duration / 2) #SAR duartion

        # Read number of platforms (todo: move into a struct)
        Np  = size(orbit_pos)[2] # number of platforms

        ref_plat = 1 #incicate the reference platform
        if processing_mode == 1
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,1:end,:], orbit_vel[:,1:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, "left",  ref_plat)
            avg_sep = maximum(bperp)/(Np - 1) # Change this
        elseif processing_mode == 2
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos[:,2:end,:], orbit_vel[:,2:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, "left",  ref_plat)
            avg_sep = maximum(bperp)/(Np - 2) # Change this
        end

        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bnorm) ./ 1e3
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bnorm)) ./ 1e3
        
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bperp) ./ 1e3
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bperp)) ./ 1e3
        
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(b_at) ./ 1e3
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,b_at)) ./ 1e3

        # theoretical resolution along-n
        range_s, range_g = Scene.lookangle_to_range(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/fc)*range_s/p_mode/ maximum(bperp) 

        # theoretical resolution along-track
        mu = 3.986004418e14
        sc_speed = sqrt(mu./(mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:])+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
        Lsa = sc_speed*SAR_duration + sc_speed/fp
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/fc)*range_s/2/Lsa

        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/fc)*range_s/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/fc)*range_s/p_mode/avg_sep

    end
end

end

to


temp = GeoArray(Canopy_heights_orig[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Canopy_heights_GEDIL3.tif", temp)

temp = GeoArray(Norm_baseline_max[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Baseline_max.tif", temp)

temp = GeoArray(Norm_baseline_min[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Baseline_min.tif", temp)

temp = GeoArray(Perp_baseline_max[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Baseline_perp_max.tif", temp)

temp = GeoArray(Perp_baseline_min[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Baseline_perp_min.tif", temp)

temp = GeoArray(Par_baseline_max[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Baseline_AT_max.tif", temp)

temp = GeoArray(Par_baseline_min[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Baseline_AT_min.tif", temp)

temp = GeoArray(amb_H[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Theo_ambiguity_height.tif", temp)

temp = GeoArray(res_theory_n[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Theo_resolution_n.tif", temp)

temp = GeoArray(lookang_all[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/Look_angle.tif", temp)

temp = GeoArray(Geo_location[:,:,1])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/temp_geo_location_lat.tif", temp)

temp = GeoArray(Geo_location[:,:,2])
bbox!(temp, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))  # roughly the Netherlands
epsg!(temp, 6933)  # in WGS84
GeoArrays.write!("../Outputs/global_geometry_10m/temp_geo_location_lon.tif", temp)

@save "../Outputs/global_geometry_10m/global_geometry_10km_1.jld" Geo_location Canopy_heights Canopy_heights_orig orbit_time_all orbit_pos_all orbit_vel_all lookang_all Orbit_index Norm_baseline_max Norm_baseline_min Norm_baseline_mean Perp_baseline_max Perp_baseline_min Perp_baseline_mean Par_baseline_max Par_baseline_min Par_baseline_mean res_theory_n res_theory_s amb_H amb_N slnt_range to

[rmprocs(p) for p in workers()]

