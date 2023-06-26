using Distributed
addprocs(24) 


@everywhere include("../modules/geometry.jl")
@everywhere include("../modules/orbits.jl")
@everywhere include("../modules/scene.jl")
@everywhere include("../modules/global_scale_support.jl")
@everywhere using NCDatasets
@everywhere using Statistics
@everywhere using Parameters
@everywhere using Dates
@everywhere using StaticArrays
@everywhere using SharedArrays
@everywhere using Interpolations
@everywhere using Plots
@everywhere using TimerOutputs
#@everywhere using PyPlot
@everywhere using GeoDatasets
@everywhere using JLD2

c = 299792458 #TODO does not work without redefining c here
earth_radius = 6378.137e3 # Earth semi-major axis at equator

const to = TimerOutput()

@timeit to "Reading L3 file for canopy heights " begin

#Read canopy heights
filepath_GEDIL3 = "/u/epstein-z0/wblr/joshil/DARTS/GEDI_Data/GEDI03_rh100_mean_2019108_2021104_002_02.tif"
grid_res        = 40;
Canopy_heights, Geo_location, size_row, size_col = Global_Scale_Support.read_GEDI_L3_data(filepath_GEDIL3, grid_res)


end

@timeit to "Initialization " begin
# Define outout variables
global Output_stat          = SharedArray(zeros(size_row,size_col,4))
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

global amb_H                = SharedArray(zeros(size_row,size_col,1))
global amb_N                = SharedArray(zeros(size_row,size_col,1))
global slnt_range           = SharedArray(zeros(size_row,size_col,1))

end

@timeit to "Read orbit file and get orbits " begin
#region_xlims = [50,110]
#region_ylims = [15,55]
#region_xlims = 1:347
#region_ylims = 1:146
region_xlims = 1:868
region_ylims = 1:366
#region_xlims = 1:3470
#region_ylims = 1:1462

#region_xlims        = 80:81
#region_ylims        = 37:37

lat_lon_idx         = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

orbit_dataset       = Dataset("/u/epstein-z0/darts/joshil/Orbits/darts-orbits-latest/DARTS/orbit_output_06152023_3.nc") # Read orbits data in NetCDF format
global mast_plat            = 1
flag_plat           = 1 #descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, mast_plat, flag_plat)

end

lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
global data2= data

lat2 = lat' .* ones(length(lon))
lon2  = ones(length(lat))' .* lon

Lats_p = Geo_location.A[1,:,1]
Lons_p = Geo_location.A[:,1,2]

global mask = SharedArray(zeros(length(Lats_p),length(Lons_p)) )

@sync @distributed for i1 =1:size(lat_lon_idx,1)   
    
        close_val_lat_lon   = findmin(abs.(lat2.-Lats_p[lat_lon_idx[i1,2]]) + abs.(lon2.-Lons_p[lat_lon_idx[i1,1]]))
        if data2[close_val_lat_lon[2]]>0
            global mask[lat_lon_idx[i1,2],lat_lon_idx[i1,1]] = 1
        end
end


@timeit to "Processing loop over all pixels " begin

@sync @distributed for i1 = 1:size(lat_lon_idx,1)   

    try
        if mask[lat_lon_idx[i1,2],lat_lon_idx[i1,1]] !=1
        #if isnan(Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
            #Canopy_heights[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] =0
            global Output_stat[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
            global Output_stat[lat_lon_idx[i1,1],lat_lon_idx[i1,2],2] = NaN
            global Output_stat[lat_lon_idx[i1,1],lat_lon_idx[i1,2],3] = NaN
            global Output_stat[lat_lon_idx[i1,1],lat_lon_idx[i1,2],4] = NaN
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

            global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
            global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
            global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        
        else

                    
            global close_val_lat_lon   = Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

            global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = close_val_lat_lon[2]

            slrng_temp2 = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2]]; 0]))
            global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = acosd(orbit_pos_geo_all[:,close_val_lat_lon[2]][3] / slrng_temp2) 

            # Compute orbits time, position, and velocity
            if close_val_lat_lon[2] > 10
                min_lim = close_val_lat_lon[2]-10
                max_lim = close_val_lat_lon[2]+10
            else
                min_lim = 1
                max_lim = close_val_lat_lon[2]+10
            end
            global orbit_time = orbit_time_all[min_lim:max_lim]
            global orbit_pos = orbit_pos_all[:,:,min_lim:max_lim]
            global orbit_vel = orbit_vel_all[:,:,min_lim:max_lim]

            ref_plat = 1 #incicate the reference platform
            bperp, b_at, bnorm = Orbits.get_perp_baselines_new(orbit_pos, orbit_vel, lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, "right", 1)
            #bperp, b_at, bnorm = Orbits.get_perp_baselines(orbit_pos, orbit_vel, lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], ref_plat)

            global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bnorm) ./ 1e3
            global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bnorm)) ./ 1e3
            global Norm_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,bnorm)) ./ 1e3
            
            global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(bperp) ./ 1e3
            global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,bperp)) ./ 1e3
            global Perp_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,bperp)) ./ 1e3
            
            global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(b_at) ./ 1e3
            global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = minimum(filter(!iszero,b_at)) ./ 1e3
            global Par_baseline_mean[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = mean(filter(!iszero,b_at)) ./ 1e3

            # theoretical resolution along-n
            range_s, range_g = Scene.lookangle_to_range(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
            global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/1.26e9)*range_s/1/ maximum(bperp) 

            # Read number of platforms (todo: move into a struct)
            Np  = size(orbit_pos)[2] # number of platforms

            avg_sep = maximum(bperp)/(Np - 1)
            global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
            global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = (c/params.fc)*range_s/p_mode/avg_sep
                    
            global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = slrng_temp2

        end 

    catch
        continue
    end

end


end

@save "/u/epstein-z0/darts/joshil/Outputs/calc_baselines/output_baselines_06152023_3.jld" Geo_location Output_stat Canopy_heights orbit_time_all orbit_pos_all orbit_vel_all lookang_all Orbit_index Norm_baseline_max Norm_baseline_min Norm_baseline_mean Perp_baseline_max Perp_baseline_min Perp_baseline_mean Par_baseline_max Par_baseline_min Par_baseline_mean res_theory_n res_theory_s to


[rmprocs(p) for p in workers()]

