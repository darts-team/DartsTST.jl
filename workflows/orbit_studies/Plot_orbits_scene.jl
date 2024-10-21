


include("../modules/geometry.jl")
include("../modules/scene.jl")
include("../modules/orbits.jl")
include("../modules/antenna.jl")
include("../modules/simsetup.jl")
include("../modules/global_scale_support.jl")

using Plots
using NCDatasets


orbit_dataset               = Dataset("/Users/joshil/Documents/Orbits/From_Eric/cartwheel/ROSE_L_Cartwheel_4sat.nc")
#orbit_dataset               = Dataset("/Users/joshil/Documents/Orbits/From_Eric/helix/ROSE_L_Helical_4sat.nc")

global mast_plat            = 1
flag_plat                   = 1 #descending orbit
orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, mast_plat, flag_plat)


orbit_pos_geo = zeros(size(orbit_pos_all)) 
orbit_dist_12 = zeros(size(orbit_pos_all)[2],size(orbit_pos_all)[3]) 

for i=1:size(orbit_pos_all)[2]
    for j=1:size(orbit_pos_all)[3]
        orbit_pos_geo[:,i,j] = Geometry.xyz_to_geo(orbit_pos_all[:,i,j])
        orbit_dist_12[i,j] = Geometry.distance(orbit_pos_all[:,1,j],orbit_pos_all[:,i,j])
    end
end

plotly()


#display(histogram(orbit_dist_12[:,:], xlabel="Baseline [km]"))


start_idx = 1900
stop_idx = 2005

p1 = plot()
for i=1:size(orbit_pos_all)[2]
    p1 = (plot!(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],orbit_pos_geo[3,i,start_idx:stop_idx],leg=false,camera=(20,40),markersize=0.2,title="Platform Positions")) #display  position of each platform at each pulse in 3D
    #p1 = (plot!(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],zeros(length(orbit_pos_geo[3,i,start_idx:stop_idx])),lc=:grey,leg=false,camera=(20,40),markersize=3,title="Platform Positions",background_color = :white)) #display  position of each platform at each pulse in 3D

    p1 = (scatter!(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],orbit_pos_geo[3,i,start_idx:stop_idx],leg=false,camera=(70,20),markersize=0.2,title="Platform Positions")) #display  position of each platform at each pulse in 3D
end
display(p1)


master_p = 1
start_idx = 1
stop_idx = 4400
p1 = plot()
for i=1:size(orbit_pos_all)[2]
    p1 = (scatter!(orbit_pos_geo[1,i,start_idx:stop_idx].-orbit_pos_geo[1,master_p,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx].-orbit_pos_geo[2,master_p,start_idx:stop_idx],orbit_pos_geo[3,i,start_idx:stop_idx].-orbit_pos_geo[3,master_p,start_idx:stop_idx],leg=false,camera=(20,40),markersize=0.2,title="Platform Positions")) #display  position of each platform at each pulse in 3D
    #p1 = (plot!(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],zeros(length(orbit_pos_geo[3,i,start_idx:stop_idx])),lc=:grey,leg=false,camera=(20,40),markersize=3,title="Platform Positions",background_color = :white)) #display  position of each platform at each pulse in 3D

    #p1 = (scatter!(orbit_pos_geo[1,i,start_idx:stop_idx].-orbit_pos_geo[1,1,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx].-orbit_pos_geo[2,1,start_idx:stop_idx],leg=false,camera=(70,20),markersize=3,title="Platform Positions")) #display  position of each platform at each pulse in 3D
end
display(p1)



start_idx = 1
stop_idx = 4400
p1 = plot()
for i=1:size(orbit_pos_all)[2]
    p1 = (scatter!(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],orbit_pos_geo[3,i,start_idx:stop_idx],leg=false,camera=(20,40),markersize=0.2,title="Platform Positions")) #display  position of each platform at each pulse in 3D
    #p1 = (plot!(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],zeros(length(orbit_pos_geo[3,i,start_idx:stop_idx])),lc=:grey,leg=false,camera=(20,40),markersize=3,title="Platform Positions",background_color = :white)) #display  position of each platform at each pulse in 3D

    #p1 = (scatter!(orbit_pos_geo[1,i,start_idx:stop_idx].-orbit_pos_geo[1,1,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx].-orbit_pos_geo[2,1,start_idx:stop_idx],leg=false,camera=(70,20),markersize=3,title="Platform Positions")) #display  position of each platform at each pulse in 3D
end
display(p1)



master_p = 1
start_idx = 1
stop_idx = 4400
p1 = plot()
for i=4#size(orbit_pos_all)[2]
    p1 = (scatter!(orbit_pos_all[1,i,start_idx:stop_idx].-orbit_pos_all[1,master_p,start_idx:stop_idx],orbit_pos_all[2,i,start_idx:stop_idx].-orbit_pos_all[2,master_p,start_idx:stop_idx],orbit_pos_all[3,i,start_idx:stop_idx].-orbit_pos_all[3,master_p,start_idx:stop_idx],leg=false,camera=(20,40),markersize=0.2,title="Platform Positions")) #display  position of each platform at each pulse in 3D
    #p1 = (plot!(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],zeros(length(orbit_pos_geo[3,i,start_idx:stop_idx])),lc=:grey,leg=false,camera=(20,40),markersize=3,title="Platform Positions",background_color = :white)) #display  position of each platform at each pulse in 3D

    #p1 = (scatter!(orbit_pos_geo[1,i,start_idx:stop_idx].-orbit_pos_geo[1,1,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx].-orbit_pos_geo[2,1,start_idx:stop_idx],leg=false,camera=(70,20),markersize=3,title="Platform Positions")) #display  position of each platform at each pulse in 3D
end
display(p1)
p1= plot!(xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]")





start_idx = 2000
stop_idx = 2005
p1 = plot()
for i=1:size(orbit_pos_all)[2]
    p1 = (plot!(orbit_pos_all[1,i,start_idx:stop_idx],orbit_pos_all[2,i,start_idx:stop_idx],orbit_pos_all[3,i,start_idx:stop_idx],leg=false,camera=(20,40),markersize=3,title="Platform Positions")) #display  position of each platform at each pulse in 3D
    #p1 = (plot!(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],zeros(length(orbit_pos_geo[3,i,start_idx:stop_idx])),lc=:grey,leg=false,camera=(20,40),markersize=3,title="Platform Positions",background_color = :white)) #display  position of each platform at each pulse in 3D

    #p1 = (scatter(orbit_pos_geo[1,i,start_idx:stop_idx],orbit_pos_geo[2,i,start_idx:stop_idx],orbit_pos_geo[3,i,start_idx:stop_idx],leg=false,camera=(70,20),markersize=3,title="Platform Positions")) #display  position of each platform at each pulse in 3D
end
display(p1)
p1= plot!(xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]")







start_idx = 1
stop_idx = 4440
p1 = plot()

earth_radius    = 6378.137e3 # Earth semi-major axis at equator
n = 100
u = range(-π, π; length = n)
v = range(0, π; length = n)
x = earth_radius * cos.(u) * sin.(v)'
y = earth_radius * sin.(u) * sin.(v)'
z = earth_radius * ones(n) * cos.(v)'

my_cg = cgrad([:blue,:blue])
p1 = surface!(x, y, z,c=my_cg,colorbar=false)

for i=1:size(orbit_pos_all)[2]
    p1 = (scatter!(orbit_pos_all[1,i,start_idx:stop_idx],orbit_pos_all[2,i,start_idx:stop_idx],orbit_pos_all[3,i,start_idx:stop_idx],leg=false,camera=(20,40),markersize=0.2,title="Platform Positions")) #display  position of each platform at each pulse in 3D
end
display(p1)

p1= plot!(xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]")
