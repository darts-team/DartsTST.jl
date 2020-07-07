include("modules/geometry.jl")
include("modules/scene.jl")
using Plots

## planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)


# target volume grid on surface defined in geo (θϕh)
t_θ=30:1:60 # deg
t_ϕ=0:2:60 # deg
t_h=0:100:3000 # m
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops (took 50 sec)
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing (took 33 sec)
# convert target volume from geo to xyz
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,a,e)
#display grid in 3D
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],markersize=1)

# target volume grid on surface defined in azimuth-elevation-height (az-el is for the look angle)
θ_l=60 # deg look vector look angle
ϕ_l=0 # deg look vector azimuth angle
t_h=0 # m target height
p_h=500e3 # m platform altitude
p_θ=0 # deg platform latitude
p_ϕ=0 # deg platform longitude
γ=90 # deg track angle (heading), 0 deg is north, 90 deg is east
peg=[p_θ,p_ϕ,γ] # deg peg point is the nadir point of platform at the center of SAR aperture
p_geo=[p_θ,p_ϕ,p_h]
t_azelh_grid=Scene.form3Dgrid_for(ϕ_l,θ_l,t_h) # using 3 nested for loops
#t_azelh_grid=Scene.form3Dgrid_array(ϕ_l,θ_l,t_h) # using array processing
t_xyz_grid=Scene.azelh_to_xyz(t_azelh_grid,p_geo,peg,a,e)
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],markersize=1)
