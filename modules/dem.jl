module DEM

using Statistics
using ArchGDAL
using GeoArrays

import ArchGDAL as AG


c               = 299792458
earth_radius    = 6378.137e3 # Earth semi-major axis at equator

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

function custom_DEM(ref_lat, ref_lon, lat_extent, lon_extent, lat_res, lon_res)

    lat_range               = ref_lat:lat_res:(ref_lat+lat_extent)-lat_res   
    lon_range               = ref_lon:lon_res:(ref_lon+lon_extent)-lon_res 

    if iseven(length(lat_range))
        lat_range               = ref_lat:lat_res:(ref_lat+lat_extent+lat_res)-lat_res 
    end  
    if iseven(length(lon_range))
        lon_range               = ref_lon:lon_res:(ref_lon+lon_extent+lon_res)-lon_res  
    end 

    lon_res_dist            = deg_to_m_lon(lon_res, ref_lat)
    lon_vert_dist_up        = lon_res_dist .* tand(3.9)
    lon_vert_dist_down      = lon_res_dist .* tand(-1.9)
    mid_idx                 = Int64(ceil(length(lon_range)/2))

    DEM_lat                 = zeros(length(lat_range))
    DEM_lon_up              = collect((0:mid_idx-1)*lon_vert_dist_up)
    DEM_lon_down            = collect(DEM_lon_up[end] .+ ((1:mid_idx-1)*lon_vert_dist_down))
    DEM_lon                 = vcat(DEM_lon_up,DEM_lon_down)

    DEM_full_lat            = repeat(DEM_lat',length(DEM_lon),1)
    DEM_full_lon            = repeat(DEM_lon,1,length(DEM_lat))

    DEM_full                = DEM_full_lat + DEM_full_lon

    Geo_location_lat_mat    = repeat(collect(lat_range)',length(lon_range),1)
    Geo_location_lon_mat    = repeat(collect(lon_range),1,length(lat_range))

    ag_geotransform = [ref_lon-(lon_res)/2;lon_res;0.0;ref_lat-(lat_res/2);0.0;lat_res]
    ag_ref = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]"

    return DEM_full, Geo_location_lat_mat, Geo_location_lon_mat, ag_geotransform, ag_ref 

end

function custom_DEM2(ref_lat, ref_lon, lat_extent, lon_extent, lat_res, lon_res)

    lat_range               = ref_lat:lat_res:(ref_lat+lat_extent)-lat_res   
    lon_range               = ref_lon:lon_res:(ref_lon+lon_extent)-lon_res 

    if iseven(length(lat_range))
        lat_range               = ref_lat:lat_res:(ref_lat+lat_extent+lat_res)-lat_res 
    end  
    if iseven(length(lon_range))
        lon_range               = ref_lon:lon_res:(ref_lon+lon_extent+lon_res)-lon_res  
    end 

    lon_res_dist            = deg_to_m_lon(lon_res, ref_lat)
    lon_vert_dist_up        = lon_res_dist .* tand(3.9)
    lon_vert_dist_down      = lon_res_dist .* tand(-1.9)
    mid_idx                 = Int64(ceil(length(lon_range)/2))
    strat_idx               = Int64(ceil(length(lon_range)/3))

    DEM_lat                 = zeros(length(lat_range))
    DEM_lon                 = zeros(length(lon_range))

    DEM_lon[strat_idx:mid_idx-1]              = collect((strat_idx:mid_idx-1)*lon_vert_dist_up)
    DEM_lon[mid_idx+1:strat_idx*2]            = collect(DEM_lon[mid_idx-1] .+ ((mid_idx+1:strat_idx*2)*lon_vert_dist_down))

    DEM_full_lat            = repeat(DEM_lat',length(DEM_lon),1)
    DEM_full_lon            = repeat(DEM_lon,1,length(DEM_lat))

    DEM_full                = DEM_full_lat + DEM_full_lon

    Geo_location_lat_mat    = repeat(collect(lat_range)',length(lon_range),1)
    Geo_location_lon_mat    = repeat(collect(lon_range),1,length(lat_range))

    ag_geotransform = [ref_lon-(lon_res)/2;lon_res;0.0;ref_lat-(lat_res/2);0.0;lat_res]
    ag_ref = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]"

    return DEM_full, Geo_location_lat_mat, Geo_location_lon_mat, ag_geotransform, ag_ref 

end


function create_slope_DEM(ref_lat, ref_lon, lat_extent, lon_extent, lat_res, lon_res, slope_lat, slope_lon) 

    lat_range               = ref_lat:lat_res:(ref_lat+lat_extent)-lat_res   
    lon_range               = ref_lon:lon_res:(ref_lon+lon_extent)-lon_res  

    if iseven(length(lat_range))
        lat_range               = ref_lat:lat_res:(ref_lat+lat_extent+lat_res)-lat_res 
    end  
    if iseven(length(lon_range))
        lon_range               = ref_lon:lon_res:(ref_lon+lon_extent+lon_res)-lon_res  
    end  
    
    lat_res_dist            = deg_to_m_lat(lat_res)
    lat_vert_dist           = lat_res_dist .* tand(slope_lat)
    lon_res_dist            = deg_to_m_lon(lon_res, ref_lat)
    lon_vert_dist           = lon_res_dist .* tand(slope_lon)

    if lat_vert_dist == 0
        DEM_lat             = zeros(length(lat_range))
    else
        DEM_lat             = collect(0:lat_vert_dist:(lat_vert_dist*length(lat_range))-lat_vert_dist)
    end

    if lon_vert_dist == 0
        DEM_lon             = zeros(length(lon_range))
    else
        DEM_lon             = collect(0:lon_vert_dist:(lon_vert_dist*length(lon_range))-lon_vert_dist)
    end

    DEM_full_lat            = repeat(DEM_lat',length(DEM_lon),1)
    DEM_full_lon            = repeat(DEM_lon,1,length(DEM_lat))

    DEM_full                = DEM_full_lat + DEM_full_lon

    Geo_location_lat_mat    = repeat(collect(lat_range)',length(lon_range),1)
    Geo_location_lon_mat    = repeat(collect(lon_range),1,length(lat_range))

    ag_geotransform = [ref_lon-(lon_res)/2;lon_res;0.0;ref_lat-(lat_res/2);0.0;lat_res]
    ag_ref = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]"

    return DEM_full, Geo_location_lat_mat, Geo_location_lon_mat, ag_geotransform, ag_ref 
end


##
function read_DEM_from_source(dem_source)

    dataset                     = AG.read(dem_source)
    DEM                         = AG.getband(dataset,1)
    ref                         = AG.getproj(dataset)
    gt                          = AG.getgeotransform(dataset)
    ulx, uly                    = gt[1], gt[4]
    x_pixel_size, y_pixel_size  = gt[2], gt[6]
    num_cols                    = AG.width(DEM)
    num_rows                    = AG.height(DEM)
    
    Geo_location_lat            = zeros(num_rows,1);
    Geo_location_lon            = zeros(num_cols,1);
    for i=1:num_rows
        Geo_location_lat[i]     = uly + i * y_pixel_size;
    end
    for i=1:num_cols
        Geo_location_lon[i]     = ulx + i * x_pixel_size;
    end

    return DEM, Geo_location_lat, Geo_location_lon, gt, ref

end

##
function read_interp_DEM_from_source(dem_source, trg_ref_lat, trg_ref_lon, lat_extent, lon_extent, lat_res, lon_res)
        
    # Get DEM for the region of interest by interpolation
    DEM, gt, ref            = AG.read(dem_source) do source
        AG.gdalwarp([source], [
                                "-te", "$(trg_ref_lon-lon_res)", "$(trg_ref_lat-lat_res)",
                                "$(trg_ref_lon+lon_extent+lon_res)", "$(trg_ref_lat+lat_extent+lat_res)",
                                "-tr", "$(lon_res)", "$(lat_res)",
                                "-r", "bilinear"]) do warped
            DEM_data         = AG.getband(warped, 1)
            gt_data          = AG.getgeotransform(warped)
            ref_data         = AG.getproj(warped)

            AG.read(DEM_data), gt_data, ref_data
        end
    end

    ulx, uly                    = gt[1], gt[4]
    x_pixel_size, y_pixel_size  = gt[2], gt[6]
    num_cols                    = size(DEM)[1]
    num_rows                    = size(DEM)[2]
    
    Geo_location_lat            = zeros(num_rows,1);
    Geo_location_lon            = zeros(num_cols,1);
    for i=1:num_rows
        Geo_location_lat[i]     = uly + i * y_pixel_size;
    end
    for i=1:num_cols
        Geo_location_lon[i]     = ulx + i * x_pixel_size;
    end

    Geo_location_lat_mat = repeat(Geo_location_lat', num_cols)
    Geo_location_lon_mat = repeat(Geo_location_lon, 1,num_rows)

    return DEM, Geo_location_lat_mat, Geo_location_lon_mat, gt, ref

end

function read_interp_DEM_from_source(dem_source, trg_ref_lat, trg_ref_lon, lat_extent, lon_extent, lat_res, lon_res, check_size_flag)

    DEM, Geo_location_lat_mat, Geo_location_lon_mat, gt, ref = read_interp_DEM_from_source(dem_source, trg_ref_lat, trg_ref_lon, lat_extent, lon_extent, lat_res, lon_res)

    if check_size_flag == 1
        if iseven(size(DEM)[2]) || iseven(size(DEM)[1])
            lat_extent_new = lat_extent
            lon_extent_new = lon_extent
            if iseven(size(DEM)[2])
                lat_extent_new = lat_extent + lat_res
            end
            if iseven(size(DEM)[1])
                lon_extent_new = lon_extent + lon_res
            end
            DEM, Geo_location_lat_mat, Geo_location_lon_mat, gt, ref = read_interp_DEM_from_source(dem_source, trg_ref_lat, trg_ref_lon, lat_extent_new, lon_extent_new, lat_res, lon_res)
        end
    end

    return DEM, Geo_location_lat_mat, Geo_location_lon_mat, gt, ref

end


##
function get_slopes_from_DEM(DEM, Geo_location_lon, Geo_location_lat)

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


    ### Get the slope in angles
    #slope_lon          = atand.(slope_lon)
    #slope_lat          = atand.(slope_lat)

    # Return slopes
    return slope_lat, slope_lon

end

##
function get_terrain_norm(DEM, slope_lat, slope_lon)

    # Compute terrain normal vector 
    #N = [-dh/dx, -dh/dy, 1] where dh/dx is the west-east slope and dh/dy is the south-north slope wrt to the geogrid. 
    #dh, dx, and dy are computed in meters using ENU coordinates (xyz).
    N_all                       = zeros(size(DEM)[1]*size(DEM)[2],3)
    k = 1
    for i=1:size(DEM)[1]
        for j=1:size(DEM)[2]
            #N_all[k,:]          = [-1*slope_lat[i,j,1], -1*slope_lon[i,j,1], 1]
            N_all[k,:]          = [-1*slope_lon[i,j,1], -1*slope_lat[i,j,1], 1]
            k = k + 1
        end
    end

    return N_all

end


end