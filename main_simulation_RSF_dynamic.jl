include("modules/generate_raw_data.jl")
include("modules/process_raw_data.jl")
include("modules/geometry.jl")
include("modules/scene.jl")
include("inputs/input_parameters_RSF_orbits_nadirlooking.jl")
include("modules/range_spread_function.jl") # as RSF
include("modules/orbits.jl")
using NCDatasets
using Plots
pyplot()
## RANGE SPREAD FUNCTION (matched filter output)
Srx,MF,ft,t_rx=RSF.ideal_RSF(τ,Δt,B,Trx) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
# Srx,MF,ft,t_rx=RSF.non_ideal_RSF(τ,Δt,B,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing
display(plot(ft*1e6,20*log10.(abs.(MF)),ylims=(-100+20*log10(B*τ),20*log10(B*τ)),leg=false,xlabel="fast time (μs)",ylabel="amplitude (dB)",title="Matched Filter Output (Range Spread Function)"))
display(plot(t_rx*1e6,20*log10.(abs.(Srx')),ylims=(-100+20*log10(B*τ),20*log10(B*τ)),leg=false,xlabel="fast time (μs)",ylabel="amplitude (dB)",title="Receive Window/Signal"))
## PLATFORM LOCATIONS
orbit_dataset=Dataset("inputs/orbitOutput_082020.nc") # Read orbits data in NetCDF format
t12_orbits=orbit_dataset["time"][1:2] # first two time samples
dt_orbits=t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
orbit_time_index=collect(Int(round(SAR_start_time/dt_orbits))+1:1:Int(round((SAR_start_time+SAR_duration)/dt_orbits))+1) # index range for orbit times for time interval of interest
orbit_time=orbit_dataset["time"][orbit_time_index] # read in time data
orbit_pos=orbit_dataset["position"][:,:,orbit_time_index] # read in position data
slow_time=collect(SAR_start_time:1/fp:SAR_start_time+SAR_duration) # create slow time axis
orbit_pos_interp=Orbits.interp_orbit(orbit_time,orbit_pos,slow_time) # interpolate orbit to slow time
#plot(orbit_time,orbit_pos[:,1,:]',xaxis=("time (sec)"),ylabel=("ECEF position (km)")) # plot the ECI orbit in the limited time range
#plot(slow_time, orbit_pos_interp[:,1,:]', xaxis=("time (sec)"), ylabel=("ECEF position (km)")) # plot orbits in slow-time over SAR duration with PRI sampling
# TODO modify  Generate_Raw_Data and Process_Raw_Data to accept 3D array with slow-time dimension
#p_geo_grid=Scene.form3Dgrid_for(p_θ,p_ϕ,p_h) # using 3 nested for loops
#p_geo_grid=Scene.form3Dgrid_array(p_θ,p_ϕ,p_h) # using array processing
#p_xyz_grid=Geometry.geo_to_xyz(p_geo_grid,a,e)
display(scatter(p_xyz_grid[1,:],p_xyz_grid[2,:],p_xyz_grid[3,:],leg=false,camera=(20,40),markersize=3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Platforms")) #display grid in 3D
#savefig("platforms.png")
## TARGET LOCATIONS
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,a,e)
display(scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Targets")) #display grid in 3D
#savefig("targets.png")
## DISPLAY PLATFORM AND TARGET LOCATIONS ON THE SAME PLOT
scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.1) #display grid in 3D
display(scatter!(p_xyz_grid[1,:],p_xyz_grid[2,:],p_xyz_grid[3,:],leg=false,camera=(20,40),markersize=1,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Platforms and Targets")) #display grid in 3D
#savefig("platforms_and_targets.png")
## GENERATE RAW DATA
#rawdata=Generate_Raw_Data.main(t_xyz_grid,p_xyz_grid,mode,tx_el,fc,a,e) # without RSF
ref_range=520e3 # reference range (ideally average of min and max ranges: TODO calculate from all target vs platform positions)
rawdata=Generate_Raw_Data.main_RSF(t_xyz_grid,p_xyz_grid,mode,tx_el,fc,a,e,Srx,t_rx,ref_range) # with RSF
#= plot raw data (no RSF)
if mode==3
    display(heatmap(abs.(rawdata),c=cgrad([:black,:white]),xlabel="RX platform",ylabel="TX platform",title="raw data amplitude"))
    #savefig("rawdata_magnitude.png")
    display(heatmap(angle.(rawdata)*180/pi,c=cgrad([:black,:white]),xlabel="RX platform",ylabel="TX platform",title="raw data phase (deg)"))
    #savefig("rawdata_phase.png")
else
    display(plot(abs.(rawdata),leg=false,xlabel="TX-RX platform pairs",ylabel="relative amplitude",title="raw data amplitude"))
    #savefig("rawdata_magnitude.png")
    display(plot(angle.(rawdata)*180/pi,leg=false,xlabel="TX-RX platform pairs",ylabel="phase (deg)",title="raw data phase"))
    #savefig("rawdata_phase.png")
end=#
# plot raw data (RSF)
display(heatmap(t_rx,1:size(p_xyz_grid)[2],20*log10.(abs.(rawdata)),c=cgrad([:black,:white]),xlabel="fast-time (s)",ylabel="TX/RX platform pairs",title="raw data amplitude (dB)"))
display(heatmap(t_rx,1:size(p_xyz_grid)[2],angle.(rawdata)*180/pi,c=cgrad([:black,:white]),xlabel="fast-time (s)",ylabel="TX/RX platform pairs",title="raw data phase (deg)"))
## IMAGE SCENE
Ns_θ=length(s_θ)
Ns_ϕ=length(s_ϕ)
Ns_h=length(s_h)
s_geo_grid=Scene.form3Dgrid_for(s_θ,s_ϕ,s_h) # using 3 nested for loops
#s_geo_grid=Scene.form3Dgrid_array(s_θ,s_ϕ,s_h) # using array processing
s_xyz_grid=Geometry.geo_to_xyz(s_geo_grid,a,e)
display(scatter(s_geo_grid[1,:],s_geo_grid[2,:],s_geo_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="latitude (deg)",ylabel="longitude (deg)",zlabel="height (m)",title="Scene Pixel Locations in GEO")) #display grid in 3D
display(scatter(s_xyz_grid[1,:],s_xyz_grid[2,:],s_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Scene Pixel Locations in XYZ")) #display grid in 3D
## PROCESS RAW DATA TO GENERATE IMAGE
#image_3xN=Process_Raw_Data.main(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc,a,e) # without RSF
image_3xN=Process_Raw_Data.main_RSF(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc,a,e,t_rx,ref_range)  # with RSF
image_3D=Scene.convert_image_3xN_to_3D(image_3xN,Ns_θ,Ns_ϕ,Ns_h)
## DISPLAY AND SAVE IMAGE
display(scatter(s_geo_grid[1,:],s_geo_grid[2,:],s_geo_grid[3,:],marker_z=image_3xN/maximum(image_3xN),leg=false,camera=(20,40),markersize=1,markerstrokewidth=0,xlabel="latitude (deg)",ylabel="longitude (deg)",zlabel="height (m)",title="3D Image in GEO")) #display grid in 3D
display(scatter(s_xyz_grid[1,:],s_xyz_grid[2,:],s_xyz_grid[3,:],marker_z=image_3xN/maximum(image_3xN),leg=false,camera=(20,40),markersize=1,markerstrokewidth=0,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="3D Image in XYZ")) #display grid in 3D
for k=1:Ns_h # height slices from the scene
    display(heatmap(s_ϕ,s_θ,image_3D[:,:,k],ylabel="latitude (deg)",xlabel="longitude (deg)",title="Lat/Lon 2D Image at Height="*string(s_h[k])*"m",c=cgrad([:black,:white])))
end
for k=1:Ns_θ # latitude slices from the scene
    display(heatmap(s_h,s_ϕ,image_3D[k,:,:],ylabel="longitude (deg)",xlabel="heights (m)",title="Lon/Height 2D Image at Lat="*string(s_θ[k])*"deg",c=cgrad([:black,:white])))
end
for k=1:Ns_ϕ # longitude slices from the scene
    display(heatmap(s_h,s_θ,image_3D[:,k,:],ylabel="latitude (deg)",xlabel="heights (m)",title="Lat/Height 2D Image at Lon="*string(s_ϕ[k])*"deg",c=cgrad([:black,:white])))
end
#savefig("image1.png")
