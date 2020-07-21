include("modules/generate_raw_data.jl")
include("modules/process_raw_data.jl")
include("modules/geometry.jl")
include("modules/scene.jl")
using Plots
pyplot()

c=299792458 # speed of light (m/s)
#planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)
# radar parameters
mode=1 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)
fc=1e9 # center frequency (Hz)

## PLATFORM AND TARGET LOCATIONS
# platform locations (volumetric grid) defined in geo (θϕh)
p_θ=0 # deg latitude
p_ϕ=-0.1:0.01:0.1 # deg longitude
p_h=500e3 # m  heights
# target locations (volumetric grid) defined in geo (θϕh)
t_θ=0.1-0.01:0.004:0.1+0.01 # deg latitude
t_ϕ=-0.01:0.004:0.01 # deg longitude
t_h=0:400:2000 # m  heights

t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops
p_geo_grid=Scene.form3Dgrid_for(p_θ,p_ϕ,p_h) # using 3 nested for loops
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing
#p_geo_grid=Scene.form3Dgrid_array(p_θ,p_ϕ,p_h) # using array processing
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,a,e)
p_xyz_grid=Geometry.geo_to_xyz(p_geo_grid,a,e)
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.1) #display grid in 3D
display(scatter!(p_xyz_grid[1,:],p_xyz_grid[2,:],p_xyz_grid[3,:],leg=false,camera=(20,40),markersize=1,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Platforms and Targets")) #display grid in 3D
#savefig("platforms_and_targets.png")
display(scatter(p_xyz_grid[1,:],p_xyz_grid[2,:],p_xyz_grid[3,:],leg=false,camera=(20,40),markersize=3,xlim=(6878132-12000,6878132+12000),ylim=(-12000,12000),zlim=(-12000,12000),xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Platforms")) #display grid in 3D
#savefig("platforms.png")
display(scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlim=(6377100,6380100),ylim=(-1500,1500),zlim=(9000,12000),xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Targets")) #display grid in 3D
#savefig("targets.png")
## GENERATE RAW DATA
rawdata=Generate_Raw_Data.main(t_xyz_grid,p_xyz_grid,mode,tx_el,fc,a,e)
# plot raw data
if mode==3
    display(heatmap(abs.(rawdata),c=cgrad([:black,:white]),xlabel="RX platform",ylabel="TX platform",title="raw data amplitude"))
    savefig("rawdata_magnitude.png")
    display(heatmap(angle.(rawdata)*180/pi,c=cgrad([:black,:white]),xlabel="RX platform",ylabel="TX platform",title="raw data phase (deg)"))
    savefig("rawdata_phase.png")
else
    display(plot(abs.(rawdata),leg=false,xlabel="TX-RX platform pairs",ylabel="relative amplitude",title="raw data amplitude"))
    savefig("rawdata_magnitude.png")
    display(plot(angle.(rawdata)*180/pi,leg=false,xlabel="TX-RX platform pairs",ylabel="phase (deg)",title="raw data phase"))
    savefig("rawdata_phase.png")
end

## DEFINE IMAGE SCENE

s_θ=0.1-0.012:0.001:0.1+0.012 # deg latitude
s_ϕ=-0.012:0.001:0.012 # deg longitude
s_h=-200:100:2200 # m  heights

s_geo_grid=Scene.form3Dgrid_for(s_θ,s_ϕ,s_h) # using 3 nested for loops
#s_geo_grid=Scene.form3Dgrid_array(s_θ,s_ϕ,s_h) # using array processing
s_xyz_grid=Geometry.geo_to_xyz(s_geo_grid,a,e)
display(scatter(s_xyz_grid[1,:],s_xyz_grid[2,:],s_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlim=(6377100,6380100),ylim=(-1500,1500),zlim=(9000,12000),xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Scene Pixel Locations")) #display grid in 3D
## PROCESS RAW DATA TO GENERATE IMAGE
processed_image=Process_Raw_Data.main(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc,a,e)
# display and save image
#display(heatmap(abs.(processed_image),c=cgrad([:black,:white]),xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Processed 3D Image"))
display(scatter(s_xyz_grid[1,:],s_xyz_grid[2,:],s_xyz_grid[3,:],marker_z=abs.(processed_image),leg=false,camera=(20,40),markersize=0.3,xlim=(6377100,6380100),ylim=(-1500,1500),zlim=(9000,12000),xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Processed Image")) #display grid in 3D
#savefig("image.png")
