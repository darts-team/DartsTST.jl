module Plotting

using Plots

function plot_RSF_rawdata(enable_fast_time,mode,ft,t_rx,MF,Srx,Np,Nst,rawdata)
    gr()
    if enable_fast_time
        display(plot(ft*1e6,20*log10.(abs.(MF)),ylims=(-100+20*log10(maximum(MF)),20*log10(maximum(MF))),leg=false,xlabel="fast time (μs)",ylabel="amplitude (dB)",title="Matched Filter Output (Range Spread Function)",size=(1600,900)))
        display(plot(t_rx*1e6,20*log10.(abs.(Srx)),ylims=(-100+20*log10(maximum(Srx)),20*log10(maximum(Srx))),leg=false,xlabel="fast time (μs)",ylabel="amplitude (dB)",title="Receive Window/Signal",size=(1600,900)))
    end
    if mode==1 || mode==2
        if enable_fast_time
            display(heatmap(t_rx,1:Np*Nst,20*log10.(abs.(reshape(rawdata,Np*Nst,size(t_rx)[1]))),c=cgrad([:black,:white]),xlabel="fast-time (s)",ylabel="TX/RX platform pairs",title="raw data amplitude (dB)",size=(1600,900)))
        else
            display(heatmap(1:Np,1:Nst,20*log10.(abs.(rawdata)),c=cgrad([:black,:white]),xlabel="platforms",ylabel="pulse number",title="raw data amplitude (dB)",size=(1600,900)))
        end
    elseif mode==3 #TODO
    end
end

function plot_geometry(orbit_time,orbit_pos,p_xyz,t_xyz_grid,s_geo_grid,s_xyz_grid)
    orbit_pos_all=reshape(p_xyz,3,size(p_xyz)[2]*size(p_xyz)[3]) # platform positions in xyz; for each platform, its position at each pulse (PRI) is plotted; output loops over platforms first, then slow-time
    gr()
    display(plot(orbit_time,orbit_pos[1,:,:]',xaxis=("time (sec)"),ylabel=("ECI X position (km)"),size=(1600,900))) # plot the ECI orbit in the limited time range
    display(plot(orbit_time,orbit_pos[2,:,:]',xaxis=("time (sec)"),ylabel=("ECI Y position (km)"),size=(1600,900))) # plot the ECI orbit in the limited time range
    display(plot(orbit_time,orbit_pos[3,:,:]',xaxis=("time (sec)"),ylabel=("ECI Z position (km)"),size=(1600,900))) # plot the ECI orbit in the limited time range
    plotly()
    display(scatter(orbit_pos_all[1,:],orbit_pos_all[2,:],orbit_pos_all[3,:],leg=false,camera=(20,40),markersize=3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Platform Positions at Each Pulse",size=(1600,900))) #display  position of each platform at each pulse in 3D
    display(scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Targets",size=(1600,900))) #display grid in 3D
    #DISPLAY PLATFORM AND TARGET LOCATIONS ON THE SAME PLOT
    display(scatter(t_xyz_grid[1,:],t_xyz_grid[2,:],t_xyz_grid[3,:],leg=false,camera=(20,40),markersize=1,size=(1600,900))) #display grid in 3D
    display(scatter!(orbit_pos_all[1,:],orbit_pos_all[2,:],orbit_pos_all[3,:],leg=false,camera=(20,40),markersize=1,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Platforms and Targets")) #display grid in 3D
    display(scatter(s_geo_grid[1,:],s_geo_grid[2,:],s_geo_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="latitude (deg)",ylabel="longitude (deg)",zlabel="height (m)",title="Scene Pixel Locations in GEO",size=(1600,900))) #display grid in 3D
    display(scatter(s_xyz_grid[1,:],s_xyz_grid[2,:],s_xyz_grid[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Scene Pixel Locations in XYZ",size=(1600,900),xlim=(minimum(s_xyz_grid[1,:]),maximum(s_xyz_grid[1,:])))) #display grid in 3D
end

function plot_tomogram(display_tomograms,image_3D,s_θ,s_ϕ,s_h,s_geo_grid,s_xyz_grid)
    brightest=maximum(image_3D)
    faintest=minimum(image_3D)
    Ns_θ=length(s_θ)
    Ns_ϕ=length(s_ϕ)
    Ns_h=length(s_h)
    if display_tomograms==1
        gr()
        k=Int(round(Ns_h/2));display(heatmap(s_ϕ,s_θ,image_3D[:,:,k],ylabel="latitude (deg)",xlabel="longitude (deg)",title="Lat/Lon 2D Image at Height="*string(s_h[k])*"m",c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,900))) #aspect_ratio=:equal
        k=Int(round(Ns_θ/2));display(heatmap(s_h,s_ϕ,image_3D[k,:,:],ylabel="longitude (deg)",xlabel="heights (m)",title="Lon/Height 2D Image at Lat="*string(s_θ[k])*"deg",c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,900))) #aspect_ratio=:equal
        k=Int(round(Ns_ϕ/2));display(heatmap(s_h,s_θ,image_3D[:,k,:],ylabel="latitude (deg)",xlabel="heights (m)",title="Lat/Height 2D Image at Lon="*string(s_ϕ[k])*"deg",c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,900))) #aspect_ratio=:equal
    elseif display_tomograms==2
        gr()
        for k=1:Ns_h # height slices from the scene
            display(heatmap(s_ϕ,s_θ,image_3D[:,:,k],ylabel="latitude (deg)",xlabel="longitude (deg)",title="Lat/Lon 2D Image at Height="*string(s_h[k])*"m",c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,900))) #aspect_ratio=:equal
        end
        for k=1:Ns_θ # latitude slices from the scene
            display(heatmap(s_h,s_ϕ,image_3D[k,:,:],ylabel="longitude (deg)",xlabel="heights (m)",title="Lon/Height 2D Image at Lat="*string(s_θ[k])*"deg",c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,900))) #aspect_ratio=:equal
        end
        for k=1:Ns_ϕ # longitude slices from the scene
            display(heatmap(s_h,s_θ,image_3D[:,k,:],ylabel="latitude (deg)",xlabel="heights (m)",title="Lat/Height 2D Image at Lon="*string(s_ϕ[k])*"deg",c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,900))) #aspect_ratio=:equal
        end
    elseif display_tomograms==3
        plotly()
        display(scatter(s_geo_grid[1,:],s_geo_grid[2,:],s_geo_grid[3,:],marker_z=image_1xN/maximum(image_1xN),leg=false,camera=(20,40),markersize=1,markerstrokewidth=0,xlabel="latitude (deg)",ylabel="longitude (deg)",zlabel="height (m)",title="3D Image in GEO",size=(1600,900))) #display grid in 3D
        display(scatter(s_xyz_grid[1,:],s_xyz_grid[2,:],s_xyz_grid[3,:],marker_z=image_1xN/maximum(image_1xN),leg=false,camera=(20,40),markersize=1,markerstrokewidth=0,xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="3D Image in XYZ",size=(1600,900))) #display grid in 3D
    end
    #savefig("tomogram.png")
end

end
