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

    #t_loc_lat_range           = trg_ref_lat:lat_res:(trg_ref_lat+lat_extent)-lat_res	# S-grid range
    #t_loc_lon_range           = trg_ref_lon:lon_res:(trg_ref_lon+lon_extent)-lon_res	# C-grid range
    #Geo_location_lat_mat        = repeat(collect(trg_ref_lat-lat_res:lat_res:trg_ref_lat+lat_extent)',length(t_loc_lon_range)+2,1)
    #Geo_location_lon_mat        = repeat(collect(trg_ref_lon-lon_res:lon_res:trg_ref_lon+lon_extent),1,length(t_loc_lat_range)+2)

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
            N_all[k,:]          = [-1*slope_lat[i,j,1], -1*slope_lon[i,j,1], 1]
            k = k + 1
        end
    end

    return N_all

end


end