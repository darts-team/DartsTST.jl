using Distributed
addprocs(4) 

#@everywhere DEPOT_PATH[1]="/u/epstein-z0/wblr/joshil/Julia/.julia" # for Epstein

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
@everywhere using Interpolations
@everywhere using Setfield

c               = 299792458 # Speed of light
earth_radius    = 6378.137e3 # Earth semi-major axis at equator

const to        = TimerOutput()

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

# Function to compute angel between two vectors
@everywhere function angle_2vec(a, b)
    return acosd(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
end

res_sim                     = 10 #resolution grid in km

file_path                   = "/Users/joshil/Documents/Code/Inputs/EPSG4326_v1.0_resamp_cubic_1.03km.tif" # DEM file path

DEM_orig                    = GeoArrays.read(file_path) # Read DEM
DEM                         = DEM_orig # Make a copy of the DEM
Geo_location                = GeoArrays.coords(DEM) # Get coordinates corresponding to DEM

DEM                         = DEM[1:res_sim:end,1:res_sim:end,:] # Skip 10 samples for reducing the resolution 1km -> 10km
Geo_location                = Geo_location[1:res_sim:end,1:res_sim:end,:] # Skip 10 samples for reducing the resolution 1km -> 10km

#Himalayan region
#DEM                         = DEM[27001:29000,5351:6850,:] # Remove the first few and last few samples, near the poles
#Geo_location                = Geo_location[27001:29000,5351:6850,:] # Remove the first few and last few samples, near the poles

#DEM                         = DEM[:,100:end-100,:] # Remove the first few and last few samples, near the poles
#Geo_location                = Geo_location[:,100:end-100,:] # Remove the first few and last few samples, near the poles

size_row                    = size(DEM)[1]
size_col                    = size(DEM)[2]

# Seperate lat and lon from Geo_location for computations
Geo_location_lon            = zeros(size(Geo_location)[1], size(Geo_location)[2],1)
Geo_location_lat            = zeros(size(Geo_location)[1], size(Geo_location)[2],1)

for i=1:size(Geo_location)[1]
    for j=1:size(Geo_location)[2]
        Geo_location_lon[i,j,1] = Geo_location[i,j,1][1]
        Geo_location_lat[i,j,1] = Geo_location[i,j,1][2]
    end
end

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
dist_lon                    = deg_to_m_lon.(dist_diff_mat_lon, Geo_location_lat)

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

# Compute terrain normal vector 
#N = [-dh/dx, -dh/dy, 1] where dh/dx is the west-east slope and dh/dy is the south-north slope wrt to the geogrid. 
#dh, dx, and dy are computed in meters using ENU coordinates (xyz).
N                           =  zeros(size(Geo_location)[1], size(Geo_location)[2],3)
for i=1:size(Geo_location)[1]
    for j=1:size(Geo_location)[2]
        N[i,j,:]            = [-1*slope_lon[i,j,1], -1*slope_lat[i,j,1], 1]
    end
end

### Get the slope angle along lat and lon
slope_lon2          = atand.(slope_lon)
slope_lat2          = atand.(slope_lat)


@timeit to "Initialization " begin

global Orbit_index          = SharedArray(zeros(size_row,size_col,1)) # Orbit time sample index corresponding to each geolocation pixel
global lookang_all          = SharedArray(zeros(size_row,size_col,1)) # Look angle 
#Norm baselines
global Norm_baseline_max    = SharedArray(zeros(size_row,size_col,1)) # Maximum baseline
global Norm_baseline_min    = SharedArray(zeros(size_row,size_col,1)) # Minumum baseline
# perp baselines
global Perp_baseline_max    = SharedArray(zeros(size_row,size_col,1)) # Maximum prependicular baseline
global Perp_baseline_min    = SharedArray(zeros(size_row,size_col,1)) # Minumum perpendicular baseline
# AT baselines
global Par_baseline_max     = SharedArray(zeros(size_row,size_col,1)) # Maximum along track baseline
global Par_baseline_min     = SharedArray(zeros(size_row,size_col,1)) # Minumum along track baseline

global res_theory_n         = SharedArray(zeros(size_row,size_col,1)) # Theoretical resolution - tomo axis
global res_theory_s         = SharedArray(zeros(size_row,size_col,1)) # Theoretical resolution - azimuth axis

global amb_H                = SharedArray(zeros(size_row,size_col,1)) # Ambiguity along height
global amb_N                = SharedArray(zeros(size_row,size_col,1)) # Ambiguity along tomo axis
global slnt_range           = SharedArray(zeros(size_row,size_col,1)) # Slant range
global altitude             = SharedArray(zeros(size_row,size_col,1)) # Slant range

global local_incidence_angle= SharedArray(zeros(size_row,size_col,1)) # Local incidence angle
global range_slope_angle    = SharedArray(zeros(size_row,size_col,1)) # Range slope angle
global incidence_angle      = SharedArray(zeros(size_row,size_col,1)) # incidence angle
global slope_angle          = SharedArray(zeros(size_row,size_col,1)) # slope angle


end

@timeit to "Read orbit file and get orbits " begin
# Read teh orbit file and obtain orbit information
region_xlims                = 1:size_row 
region_ylims                = 1:size_col 

lat_lon_idx                 = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

#orbit_dataset               = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/orbit_output_06152023_3.nc")
#orbit_dataset               = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/ROSE_L_Helical_4sat.nc") 

#orbit_dataset               = Dataset("/Users/joshil/Documents/Orbits/From_Eric/cartwheel/ROSE_L_Cartwheel_4sat.nc")
orbit_dataset               = Dataset("/Users/joshil/Documents/Orbits/From_Eric/helix/ROSE_L_Helical_4sat.nc")

#orbit_dataset               = Dataset("/Users/joshil/Documents/Orbits/Outputs/06152023/3/orbit_output_06152023_3.nc") # "orbit_output_04052023.nc") # Read orbits data in NetCDF format

global mast_plat            = 1
flag_plat                   = 1 #descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset,"ECI", mast_plat, flag_plat)

end

@timeit to "Generate mask " begin
# Create a global mask based on land-ocean boundary to remove data over ocean
global mask                 = SharedArray(zeros(size(Geo_location)[1],size(Geo_location)[2]) )

lon,lat,data                = GeoDatasets.landseamask(;resolution='c',grid=5)
itp                         = LinearInterpolation((lon, lat), data, extrapolation_bc=0)
for i=1:size(Geo_location)[1]
    for j=1:size(Geo_location)[2]
        mask[i,j]           = itp(Geo_location[i,j,1][1], Geo_location[i,j,1][2]) 

    end
end

end


@timeit to "Processing loop over all pixels " begin

@sync @distributed for i1 = 1:size(lat_lon_idx,1)   
     
    if mask[lat_lon_idx[i1,1],lat_lon_idx[i1,2]] <1
        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]       = NaN
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]       = NaN
        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = NaN
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = NaN
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = NaN
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = NaN
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]             = NaN
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]             = NaN
        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]        = NaN    
        global local_incidence_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]        = NaN
        global range_slope_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global incidence_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = NaN
        global slope_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]       = NaN
        global altitude[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]       = NaN

    else

	try

        # Obtain closest orbit point
        global close_val_lat_lon                                    =  Global_Scale_Support.find_close_val_lat_lon_test(Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1][end:-1:1], orbit_pos_all, orbit_pos_geo_all)
        #global close_val_lat_lon                                    =  Global_Scale_Support.find_close_val_lat_lon(Geo_location, lat_lon_idx[i1,:], orbit_pos_all, orbit_pos_geo_all)

        # Obtain orbit index corresponding to closest orbit point
        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = close_val_lat_lon[2]
        # Compute slant range and look angle
        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]    = Geometry.distance(orbit_pos_all[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1][2];Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1][1]; Float64(DEM[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])]))
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = Scene.slantrange_to_lookangle(earth_radius,slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1],orbit_pos_geo_all[:,close_val_lat_lon[2]][3],0.0)[2]

        global altitude[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = orbit_pos_geo_all[:,close_val_lat_lon[2]][3]

        # Define sim parameters
        p_mode                  = 1 # 2: SISO, 1: SIMO
        processing_mode         = 1 # Processing configuration 1: use all platforms
        SAR_duration            = 5 # SAR duration
        fc                      = 1.25e9 # Carrier frequency
        fp                      = 10 # Pulse repitition time

        ref_plat                = 1 #incicate the reference platform

        # Compute orbits time, position, and velocity
        global orbit_time       = orbit_time_all[close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_pos        = orbit_pos_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        global orbit_vel        = orbit_vel_all[:,:,close_val_lat_lon[2]-10:close_val_lat_lon[2]+10]
        
        # Compute SAR start time
        global SAR_start_time   = orbit_time_all[close_val_lat_lon[2]] - (SAR_duration / 2) #SAR duartion

        # Read number of platforms (todo: move into a struct)
        Np                      = size(orbit_pos)[2] # number of platforms

        if processing_mode == 1
            bperp, b_at, bnorm  = Orbits.get_perp_baselines_new(orbit_pos[:,1:end,:], orbit_vel[:,1:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, "left",  ref_plat)
            avg_sep             = maximum(bperp)/(Np - 1)
        elseif processing_mode == 2
            bperp, b_at, bnorm  = Orbits.get_perp_baselines_new(orbit_pos[:,2:end,:], orbit_vel[:,2:end,:], lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], 0.0, "left",  ref_plat)
            avg_sep             = maximum(bperp)/(Np - 2) 
        end

        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = maximum(bnorm) ./ 1e3
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = minimum(filter(!iszero,bnorm)) ./ 1e3
        
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = maximum(bperp) ./ 1e3
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = minimum(filter(!iszero,bperp)) ./ 1e3
        
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]  = maximum(b_at) ./ 1e3
        global Par_baseline_min[lat_lon_idx[i1,1],la t_lon_idx[i1,2],1]  = minimum(filter(!iszero,b_at)) ./ 1e3

        # theoretical resolution along-n
        range_s, range_g = Scene.lookangle_to_range(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1], mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:]), 0, earth_radius)
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = (c/fc)*range_s/p_mode/ maximum(bperp) 

        # theoretical resolution along-track
        mu                      = 3.986004418e14
        sc_speed                = sqrt(mu./(mean(Geometry.xyz_to_geo(orbit_pos[:,mast_plat,:])[3,:])+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
        Lsa                     = sc_speed*SAR_duration + sc_speed/fp
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = (c/fc)*range_s/2/Lsa

        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]             = (c/fc)*range_s/p_mode/avg_sep*sind(lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]             = (c/fc)*range_s/p_mode/avg_sep

        plat_pt_xyz             = orbit_pos_all[:,1,close_val_lat_lon[2]]
        geo_pt_xyz              = Geometry.geo_to_xyz([Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1][2], Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1][1], Float64(DEM[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1])])
        look_vec_xyz            = (plat_pt_xyz - geo_pt_xyz) 
        look_vec_xyz_norm       = (plat_pt_xyz - geo_pt_xyz) / Geometry.distance(geo_pt_xyz, plat_pt_xyz)

        pegθ                    = Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1][2]*π/180
        pegϕ                    = Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1][1]*π/180
        #ENU to XYZ transformation matrix
        Menu_xyz                = [-sin(pegϕ) -sin(pegθ)*cos(pegϕ) cos(pegθ)*cos(pegϕ);
                                    cos(pegϕ) -sin(pegθ)*sin(pegϕ) cos(pegθ)*sin(pegϕ);
                                    0            cos(pegθ)             sin(pegθ)]
        #XYZ to ENU transformation matrix
        Mxyz_enu                = [-sin(pegϕ)           cos(pegϕ)             0;
                                    -sin(pegθ)*cos(pegϕ) -sin(pegθ)*sin(pegϕ)  cos(pegθ)  ;
                                    cos(pegθ)*cos(pegϕ)   cos(pegθ)*sin(pegϕ)  sin(pegθ)]

        Nxyz                    = Menu_xyz * N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:];
        look_vec_enu            = Mxyz_enu * look_vec_xyz


        global local_incidence_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]         = angle_2vec(look_vec_xyz_norm, Nxyz)

        look_direction_norm                                                         = sqrt( look_vec_enu[1] * look_vec_enu[1]  + look_vec_enu[2] * look_vec_enu[2])
        global range_slope_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]             = atand((N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:][1] * (look_vec_enu[1] / look_direction_norm)) + (N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:][2] * (look_vec_enu[2] / look_direction_norm)))
       #global range_slope_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]             = atand((N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:][1] * (look_vec_enu[1] / norm(look_vec_enu))) + (N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],:][2] * (look_vec_enu[2] / norm(look_vec_enu)) )) # Previous implementation incorrect


        #testing code
        #temp                = (geo_pt_xyz) / norm((geo_pt_xyz))
        temp                    = Menu_xyz * [0;0;1];
        global slope_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]       = angle_2vec(temp, Nxyz) 
        global incidence_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = angle_2vec(temp, look_vec_xyz_norm) 

    catch
        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]           = NaN
        global lookang_all[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]           = NaN
        global Norm_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = NaN
        global Norm_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = NaN
        global Perp_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = NaN
        global Perp_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = NaN
        global Par_baseline_max[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = NaN
        global Par_baseline_min[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]      = NaN
        global res_theory_n[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]          = NaN
        global res_theory_s[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]          = NaN
        global amb_H[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                 = NaN
        global amb_N[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]                 = NaN
        global slnt_range[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]            = NaN    
        global local_incidence_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1] = NaN
        global range_slope_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]     = NaN
        global incidence_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]       = NaN
        global slope_angle[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]           = NaN
        global altitude[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]           = NaN

    end

    end
   
end

end

to

#Save the required outputs to TIF files

#save_path = "/Users/joshil/Documents/Outputs/geometry/global_geometry_10km_Geogrid/Outputs/"
save_path = "../Outputs/global_geometry_1km_helix/"

temp = GeoArray(DEM[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"DEM_1km.tif", temp)

temp = GeoArray(Norm_baseline_max[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Baseline_max.tif", temp)

temp = GeoArray(Norm_baseline_min[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Baseline_min.tif", temp)

temp = GeoArray(Perp_baseline_max[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Baseline_perp_max.tif", temp)

temp = GeoArray(Perp_baseline_min[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Baseline_perp_min.tif", temp)

temp = GeoArray(Par_baseline_max[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Baseline_AT_max.tif", temp)

temp = GeoArray(Par_baseline_min[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Baseline_AT_min.tif", temp)

temp = GeoArray(amb_H[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Theo_ambiguity_height.tif", temp)

temp = GeoArray(res_theory_n[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Theo_resolution_n.tif", temp)

temp = GeoArray(lookang_all[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Look_angle.tif", temp)

temp = GeoArray(local_incidence_angle[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"local_incidence_angle.tif", temp)

temp = GeoArray(range_slope_angle[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"range_slope_angle.tif", temp)

temp = GeoArray(incidence_angle[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"incidence_angle.tif", temp)

temp = GeoArray(slope_angle[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"slope_angle.tif", temp)

temp = GeoArray(slope_lon[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Slope_lon.tif", temp)

temp = GeoArray(slope_lat[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Slope_lat.tif", temp)

temp = GeoArray(slope_lon[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Slope_lon_angle.tif", temp)

temp = GeoArray(slope_lat[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Slope_lat_angle.tif", temp)

temp = GeoArray(altitude[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Altitude.tif", temp)

temp = GeoArray(slnt_range[:,:,1])
temp.f = DEM_orig.f
temp.crs = DEM_orig.crs
temp = @set temp.f.linear = (temp.f.linear .* res_sim)
temp = @set temp.f.translation = Geo_location[1,1]
GeoArrays.write!(save_path*"Slant_range.tif", temp)

# Write variables to output file
@save save_path*"global_geometry_1km_Helix_4plat.jld"  lookang_all Orbit_index Norm_baseline_max Norm_baseline_min Perp_baseline_max Perp_baseline_min Par_baseline_max Par_baseline_min res_theory_n res_theory_s amb_H amb_N slnt_range local_incidence_angle range_slope_angle incidence_angle slope_angle slope_lon slope_lat altitude to

# Remove workers initiated for distributed processing
[rmprocs(p) for p in workers()]

#END
