using Distributed
addprocs(8) 

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
@everywhere using Interpolations


c               = 299792458 # Speed of light
earth_radius    = 6378.137e3 # Earth semi-major axis at equator

const to        = TimerOutput()


res_sim                     = 10 #resolution grid in km

file_path                   = "inputs/EPSG4326_v1.0_resamp_cubic_1.03km.tif" # DEM file path

DEM_orig                    = GeoArrays.read(file_path) # Read DEM
DEM                         = DEM_orig # Make a copy of the DEM
Geo_location                = GeoArrays.coords(DEM) # Get coordinates corresponding to DEM

DEM                         = DEM[1:res_sim:end,1:res_sim:end,:] # Skip 10 samples for reducing the resolution 1km -> 10km
Geo_location                = Geo_location[1:res_sim:end,1:res_sim:end,:] # Skip 10 samples for reducing the resolution 1km -> 10km



size_row                    = size(DEM)[1]
size_col                    = size(DEM)[2]

global Orbit_index          = SharedArray(zeros(size_row,size_col,1)) # Orbit time sample index corresponding to each geolocation pixel


@timeit to "Read orbit file and get orbits " begin
# Read teh orbit file and obtain orbit information
region_xlims                = 1:size_row 
region_ylims                = 1:size_col 

lat_lon_idx                 = Global_Scale_Support.get_lat_lon_idx(region_xlims, region_ylims)

#orbit_dataset               = Dataset("/u/epstein-z0/darts/joshil/code/darts-simtool/inputs/orbit_output_06152023_3.nc")
orbit_dataset               = Dataset("/Users/joshil/Documents/Orbits/Outputs/06152023/3/orbit_output_06152023_3.nc") # "orbit_output_04052023.nc") # Read orbits data in NetCDF format

global mast_plat            = 1
flag_plat                   = 1 #descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, mast_plat, flag_plat)

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
    else
        # Obtain closest orbit point
        global close_val_lat_lon                                    =  Global_Scale_Support.find_close_val_lat_lon_test(Geo_location[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1][end:-1:1], orbit_pos_all, orbit_pos_geo_all)
        # Obtain orbit index corresponding to closest orbit point
        global Orbit_index[lat_lon_idx[i1,1],lat_lon_idx[i1,2],1]   = close_val_lat_lon[2]
    end

end

end

to


#save_path = "/Users/joshil/Documents/Outputs/geometry/global_geometry_10km_Geogrid/Outputs/"
save_path = "../Outputs/global_geometry_10m_new/"

# Write variables to output file
@save save_path*"orbit_index_10km2.jld"  Orbit_index Geo_location

# Remove workers initiated for distributed processing
[rmprocs(p) for p in workers()]

#END