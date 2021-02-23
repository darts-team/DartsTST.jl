include("modules/generate_raw_data.jl")
include("modules/process_raw_data.jl")
include("modules/geometry.jl")
include("modules/scene.jl")
include("inputs/input_parameters_RSF_orbits_nadirlooking.jl")
include("modules/range_spread_function.jl") # as RSF
include("modules/orbits.jl")
include("modules/error_sources.jl")
include("modules/Performance_Metrics.jl")
using NCDatasets
using Plots
using Statistics
#pyplot()
plotly()
#gr()
## RANGE SPREAD FUNCTION (matched filter output)
if enable_fast_time # matched filter gain is included in Srx
    Srx,MF,ft,t_rx=RSF.ideal_RSF(pulse_length,Δt,bandwidth,Trx) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
    # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing
    #display(plot(ft*1e6,20*log10.(abs.(MF)),ylims=(-100+20*log10(bandwidth*pulse_length),20*log10(bandwidth*pulse_length)),leg=false,xlabel="fast time (μs)",ylabel="amplitude (dB)",title="Matched Filter Output (Range Spread Function)",size=(1600,900)))
    #display(plot(t_rx*1e6,20*log10.(abs.(Srx)),ylims=(-100+20*log10(bandwidth*pulse_length),20*log10(bandwidth*pulse_length)),leg=false,xlabel="fast time (μs)",ylabel="amplitude (dB)",title="Receive Window/Signal",size=(1600,900)))
end
## PLATFORM LOCATIONS
orbit_dataset=Dataset("inputs/orbitOutput_082020.nc") # Read orbits data in NetCDF format
t12_orbits=orbit_dataset["time"][1:2] # first two time samples
dt_orbits=t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
orbit_time_index=(Int(round(SAR_start_time/dt_orbits))+1:1:Int(round((SAR_start_time+SAR_duration)/dt_orbits))+1) # index range for orbit times for time interval of interest
orbit_time=orbit_dataset["time"][orbit_time_index] # read in time data
orbit_pos=orbit_dataset["position"][:,:,orbit_time_index] # read in position data
#TODO convert ECI to ECEF
slow_time=(SAR_start_time:1/fp:SAR_start_time+SAR_duration) # create slow time axis
orbit_pos_interp=Orbits.interp_orbit(orbit_time,orbit_pos,slow_time) # interpolate orbit to slow time
p_xyz=1e3*orbit_pos_interp # convert km to m
#display(plot(orbit_time,orbit_pos[1,:,:]',xaxis=("time (sec)"),ylabel=("ECI X position (km)"),size=(1600,900))) # plot the ECI orbit in the limited time range
#display(plot(orbit_time,orbit_pos[2,:,:]',xaxis=("time (sec)"),ylabel=("ECI Y position (km)"),size=(1600,900))) # plot the ECI orbit in the limited time range
#display(plot(orbit_time,orbit_pos[3,:,:]',xaxis=("time (sec)"),ylabel=("ECI Z position (km)"),size=(1600,900))) # plot the ECI orbit in the limited time range
Np=size(orbit_pos)[2] # number of platforms
Nst=size(slow_time)[1] # number of slow-time samples (pulses processed)
orbit_pos_all=reshape(p_xyz,3,Np*Nst) # platform positions in xyz; for each platform, its position at each PRF treated as a different platform; output loops over platforms first, then slow-time
#display(scatter(orbit_pos_all[1,:],orbit_pos_all[2,:],orbit_pos_all[3,:],leg=false,camera=(20,40),markersize=3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Platform Positions at Each Pulse",size=(1600,900))) #display  position of each platform at each pulse in 3D
#savefig("platforms.png")
## TARGET LOCATIONS
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,earth_radius,earth_eccentricity)
#display(scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Targets",size=(1600,900))) #display grid in 3D
#savefig("targets.png")
## DISPLAY PLATFORM AND TARGET LOCATIONS ON THE SAME PLOT
#display(scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],leg=false,camera=(20,40),markersize=1,size=(1600,900))) #display grid in 3D
#display(scatter!(orbit_pos_all[1,:],orbit_pos_all[2,:],orbit_pos_all[3,:],leg=false,camera=(20,40),markersize=1,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Platforms and Targets")) #display grid in 3D
#savefig("platforms_and_targets.png")
## GENERATE RAW DATA
#rawdata=Generate_Raw_Data.main(t_xyz_grid,p_xyz_grid,mode,tx_el,fc) # without RSF
ref_range=Generate_Raw_Data.distance(mean(t_xyz_grid,dims=2),mean(mean(p_xyz,dims=2),dims=3)) # reference range
#rawdata=Generate_Raw_Data.main_RSF(t_xyz_grid,p_xyz,mode,tx_el,fc,Srx,t_rx,ref_range) # with fasttime, without slowtime #TODO use structure as input
if enable_fast_time # with fastime and slowtime; matched filter gain is included in Srx
    rawdata=Generate_Raw_Data.main_RSF_slowtime(t_xyz_grid,p_xyz,mode,tx_el,fc,Srx,t_rx,ref_range)
else # without fastime, with slowtime; matched filter gain is included inside the function
    rawdata=Generate_Raw_Data.main_noRSF_slowtime(t_xyz_grid,p_xyz,mode,tx_el,fc)
end
if !enable_fast_time
    SNR=SNR*pulse_length*bandwidth # SNR increases after matched filter
end
if enable_thermal_noise
    rawdata=Error_Sources.random_noise(rawdata,SNR,enable_fast_time,mode)
end
# plot raw data (RSF)
if mode==1 || mode==2
#    if enable_fast_time;display(heatmap(t_rx,1:size(p_xyz)[2]*size(p_xyz)[3],20*log10.(abs.(reshape(rawdata,size(p_xyz)[2]*size(p_xyz)[3],size(t_rx)[1]))),c=cgrad([:black,:white]),xlabel="fast-time (s)",ylabel="TX/RX platform pairs",title="raw data amplitude (dB)",size=(1600,900)))
#    else;display(heatmap(1:size(p_xyz)[2],1:size(p_xyz)[3],20*log10.(abs.(rawdata)),c=cgrad([:black,:white]),xlabel="platforms",ylabel="pulse number",title="raw data amplitude (dB)",size=(1600,900)));end
elseif mode==3 #TODO
end
#display(heatmap(t_rx,1:size(p_xyz)[2],angle.(rawdata)*180/pi,c=cgrad([:black,:white]),xlabel="fast-time (s)",ylabel="TX/RX platform pairs",title="raw data phase (deg)",size=(1600,900)))
## IMAGE SCENE
Ns_θ=length(s_θ)
Ns_ϕ=length(s_ϕ)
Ns_h=length(s_h)
s_geo_grid=Scene.form3Dgrid_for(s_θ,s_ϕ,s_h) # using 3 nested for loops
#s_geo_grid=Scene.form3Dgrid_array(s_θ,s_ϕ,s_h) # using array processing
s_xyz_grid=Geometry.geo_to_xyz(s_geo_grid,earth_radius,earth_eccentricity)
#display(scatter(s_geo_grid[1,:],s_geo_grid[2,:],s_geo_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="latitude (deg)",ylabel="longitude (deg)",zlabel="height (m)",title="Scene Pixel Locations in GEO",size=(1600,900))) #display grid in 3D
#display(scatter(s_xyz_grid[1,:],s_xyz_grid[2,:],s_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Scene Pixel Locations in XYZ",size=(1600,900),xlim=(minimum(s_xyz_grid[1,:]),maximum(s_xyz_grid[1,:])))) #display grid in 3D
## PROCESS RAW DATA TO GENERATE IMAGE
#image_3xN=Process_Raw_Data.main(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc) # without fastime, without slowtime
#image_3xN=Process_Raw_Data.main_RSF(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc,t_rx,ref_range)  # with fastime, without slowtime
if enable_fast_time # with fastime, with slowtime
    image_3xN=Process_Raw_Data.main_RSF_slowtime(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc,t_rx,ref_range)
else # without fastime, with slowtime
    image_3xN=Process_Raw_Data.main_noRSF_slowtime(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc)
end
image_3D=Scene.convert_image_3xN_to_3D(image_3xN,Ns_θ,Ns_ϕ,Ns_h)
## DISPLAY AND SAVE IMAGE
#display(scatter(s_geo_grid[1,:],s_geo_grid[2,:],s_geo_grid[3,:],marker_z=image_3xN/maximum(image_3xN),leg=false,camera=(20,40),markersize=1,markerstrokewidth=0,xlabel="latitude (deg)",ylabel="longitude (deg)",zlabel="height (m)",title="3D Image in GEO",size=(1600,900))) #display grid in 3D
#display(scatter(s_xyz_grid[1,:],s_xyz_grid[2,:],s_xyz_grid[3,:],marker_z=image_3xN/maximum(image_3xN),leg=false,camera=(20,40),markersize=1,markerstrokewidth=0,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="3D Image in XYZ",size=(1600,900))) #display grid in 3D
for k=1:Ns_h # height slices from the scene
    display(heatmap(s_ϕ,s_θ,image_3D[:,:,k],ylabel="latitude (deg)",xlabel="longitude (deg)",title="Lat/Lon 2D Image at Height="*string(s_h[k])*"m",c=cgrad([:black,:white]),aspect_ratio=:equal,size=(1600,900)))
end
#=for k=1:Ns_θ # latitude slices from the scene
    display(heatmap(s_h,s_ϕ,image_3D[k,:,:],ylabel="longitude (deg)",xlabel="heights (m)",title="Lon/Height 2D Image at Lat="*string(s_θ[k])*"deg",c=cgrad([:black,:white]),aspect_ratio=:equal,size=(1600,900)))
end
for k=1:Ns_ϕ # longitude slices from the scene
    display(heatmap(s_h,s_θ,image_3D[:,k,:],ylabel="latitude (deg)",xlabel="heights (m)",title="Lat/Height 2D Image at Lon="*string(s_ϕ[k])*"deg",c=cgrad([:black,:white]),aspect_ratio=:equal,size=(1600,900)))
end=#
#savefig("image1.png")
## PERFORMANCE METRICS
# Resolution
include("modules/Performance_Metrics.jl")
target_location=[t_θ t_ϕ t_h] # point target location for resolution calculation
if size(target_location)[1]==1 # resolution is calculated when there is only one point target
    resolutions=Performance_Metrics.resolution(image_3D,res_dB,target_location,s_θ,s_ϕ,s_h) # resolutions in each of the 3 axes
    print(resolutions)
end
# Sidelobes
