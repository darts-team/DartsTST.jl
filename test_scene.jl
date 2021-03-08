include("modules/geometry.jl")
include("modules/scene.jl")
using Plots
pyplot()

## planetary shape constants
<<<<<<< Updated upstream
a=6378.137e3
e=sqrt(0.00669437999015)

#planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)
=======
earth_radius=6378.137e3
earth_eccentricity=sqrt(0.00669437999015)
>>>>>>> Stashed changes

# target volume grid on surface defined in geo (θϕh)
t_θ=-90:2:90 # deg
t_ϕ=-180:4:180 # deg
t_h=0:100:9000 # m target heights
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops (took 50 sec)
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing (took 33 sec)
<<<<<<< Updated upstream
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,a,e) # convert target volume from geo to xyz
=======
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,earth_radius,earth_eccentricity) # convert target volume from geo to xyz
>>>>>>> Stashed changes
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],camera=(0,0),markersize=0.1,xlim=(-a-1e5,a+1e5),ylim=(-a-1e5,a+1e5),zlim=(-a-1e5,a+1e5)) #display grid in 3D

# target volume grid on surface defined in geo (θϕh)
t_θ=30:1:60 # deg
t_ϕ=0:2:60 # deg
t_h=0:100:3000 # m
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops (took 50 sec)
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing (took 33 sec)
# convert target volume from geo to xyz
<<<<<<< Updated upstream
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,a,e)
=======
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,earth_radius,earth_eccentricity)
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
t_xyz_grid=Scene.azelh_to_xyz(t_azelh_grid,p_geo,peg,a,e)
=======
t_xyz_grid=Scene.azelh_to_xyz(t_azelh_grid,p_geo,peg,earth_radius,earth_eccentricity)
>>>>>>> Stashed changes
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],markersize=1)

t_h=0:100:3000 # m target heights
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing
<<<<<<< Updated upstream
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,a,e) # convert target volume from geo to xyz
=======
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,earth_radius,earth_eccentricity) # convert target volume from geo to xyz
>>>>>>> Stashed changes
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],camera=(30,50),markersize=0.1,xlim=(0,a+1e5),ylim=(0,a+1e5),zlim=(0,a+1e5)) #display grid in 3D

# target volume grid on surface defined in azimuth-look angle-height (azimuth and look angle are for the look angle, target height)
θ_l=0:1:30#-60:1:60 # deg look vector look angle
ϕ_l=-180:10:180 # deg look vector azimuth angle
t_h=0:2e6:4e6 # m target heights
p_θϕh=[0,50,6000e3] # deg platform latitude, longitude, altitude
γ=0 # deg track angle (heading), 0 deg is north, 90 deg is east
peg=[p_θϕh[1],p_θϕh[2],γ] # deg peg point is the nadir point of platform at the center of SAR aperture
t_lookh_grid=Scene.form3Dgrid_for(ϕ_l,θ_l,t_h) # using 3 nested for loops
#t_azelh_grid=Scene.form3Dgrid_array(ϕ_l,θ_l,t_h) # using array processing
<<<<<<< Updated upstream
t_xyz_grid=Scene.lookh_to_xyz(t_lookh_grid,p_θϕh,peg,a,e) # convert target volume from lookh to xyz
=======
t_xyz_grid=Scene.lookh_to_xyz(t_lookh_grid,p_θϕh,peg,earth_radius,earth_eccentricity) # convert target volume from lookh to xyz
>>>>>>> Stashed changes
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],camera=(20,20),markersize=0.5,xlim=(-a-1e5,a+1e5),ylim=(-a-1e5,a+1e5),zlim=(-a-1e5,a+1e5)) #display grid in 3D

# target volume grid on surface defined in ch (of sch) and azimuth rotation angles (rotated around platform position vector)
t_c=1e6:5e5:5e6 # m target positions in c of sch
t_h=0:1e5:5e5 # m target heights
rot_P=-120:5:120#-60:10:60 # deg rotation angles around platform position vector
p_θϕh=[60,-30,500e3]
γ=0 # deg track angle (heading), 0 deg is north, 90 deg is east
peg=[p_θϕh[1],p_θϕh[2],γ] # deg peg point is the nadir point of platform at the center of SAR aperture
<<<<<<< Updated upstream
t_xyz_grid_rot=Scene.chP_to_xyz_grid(t_c,t_h,rot_P,p_θϕh,peg,a,e)# convert target volume from geo to xyz
=======
t_xyz_grid_rot=Scene.chP_to_xyz_grid(t_c,t_h,rot_P,p_θϕh,peg,earth_radius,earth_eccentricity)# convert target volume from geo to xyz
>>>>>>> Stashed changes
scatter(t_xyz_grid_rot[1,:],t_xyz_grid_rot[2,:],t_xyz_grid_rot[3,:],camera=(30,30),markersize=1
,xlim=(-a-1e5,a+1e5),ylim=(-a-1e5,a+1e5),zlim=(-a-1e5,a+1e5))#display grid in 3D

include("modules/geometry.jl")
include("modules/scene.jl")
t_sch=[0,100e3,0]
p_θϕh=[0,0,500e3]
rot_P=90
γ=0 # deg track angle (heading), 0 deg is north, 90 deg is east
peg=[p_θϕh[1],p_θϕh[2],γ] # deg peg point is the nadir point of platform at the center of SAR aperture
<<<<<<< Updated upstream
t_xyz=Scene.chP_to_xyz_single(t_sch,rot_P,p_θϕh,peg,a,e)
=======
t_xyz=Scene.chP_to_xyz_single(t_sch,rot_P,p_θϕh,peg,earth_radius,earth_eccentricity
>>>>>>> Stashed changes
println(t_xyz)
