include("modules/geometry.jl")
include("modules/scene.jl")
include("modules/raw_data.jl")
using Plots
pyplot()

c=299792458 # speed of light (m/s)
#planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)
# radar parameters
mode=1 #1: SAR, 2:SIMO, 3:MIMO
tx_element=1 # which element transmits for SIMO (max value N)

fc=1e9 # center frequency (Hz)
h=500e3 # altitude (m)
#θ_l=60 # look angle wrt nadir (deg)
α_b=30 # baseline tilt wrt horizontal (deg)

# platform locations
pb=-10e3:1000:10e3 # platform locations along baseline

# target volume grid on surface defined in geo (θϕh)
t_θ=30:1:40 # deg
t_ϕ=0:2:20 # deg
t_h=0:100:1000 # m target heights
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,a,e) # convert target volume from geo to xyz
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],camera=(30,50),markersize=0.1,xlim=(0,a+1e5),ylim=(0,a+1e5),zlim=(0,a+1e5)) #display grid in 3D

# generate raw data
rawdata=Raw_Data.main(t_xyz_grid,pb,fc,h,α_b,mode,tx_element)

# plot raw data

# save raw data
