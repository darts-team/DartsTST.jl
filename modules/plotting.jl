module Plotting

using Plots

function plot_RSF_rawdata(enable_fast_time,mode,ft,t_rx,MF,Srx,Np,Nst,rawdata) #TODO add descriptions
    gr()
    if enable_fast_time # plot range-compressed pulse (matched filter output)
        display(plot(ft*1e6,20*log10.(abs.(MF)),ylims=(-100+20*log10(maximum(MF)),20*log10(maximum(MF))),leg=false,xlabel="fast time (μs)",ylabel="amplitude (dB)",title="Matched Filter Output (Range Spread Function)",size=(1600,900)))
        display(plot(t_rx*1e6,20*log10.(abs.(Srx)),ylims=(-100+20*log10(maximum(Srx)),20*log10(maximum(Srx))),leg=false,xlabel="fast time (μs)",ylabel="amplitude (dB)",title="Receive Window/Signal",size=(1600,900)))
    end
    # plot rawdata
    if mode==1 || mode==2
        if enable_fast_time #rawdata is displayed as a 2D image of size Nft x Np*Nst
            display(heatmap(t_rx*1e6,1:Np*Nst,20*log10.(abs.(reshape(rawdata,Np*Nst,size(t_rx)[1]))),c=cgrad([:black,:white]),xlabel="fast-time (us)",ylabel="TX/RX platform pairs",title="raw data amplitude (dB)",size=(1600,900)))
        else #rawdata is displayed as a 2D image of size Np x Nst
            display(heatmap(1:Np,1:Nst,20*log10.(abs.(rawdata)),c=cgrad([:black,:white]),xlabel="platforms",ylabel="pulse number",title="raw data amplitude (dB)",size=(1600,900)))
        end
    elseif mode==3 #TODO how to display it?
        if enable_fast_time # rawdata is a 4D array of size Nst x Np(RX) x Np(TX) x Nft
        else # rawdata is a 3D array of size Nst x Np(RX) x Np(TX)
        end
    end
end

function coordinates(coord_sys)
    if coord_sys=="LLH"
        coords=["latitude","longitude","height"]
    elseif coord_sys=="SCH"
        coords=["along-track","cross-track","height"]
    elseif coord_sys=="XYZ"
        coords=["ECEF-X","ECEF-Y","ECEF-Z"]
    end
end

function plot_geometry(orbit_time,orbit_pos,p_loc,t_loc,s_loc,coords) #TODO smart plotting for limited set (only platforms, only targets, only scene)
    orbit_pos_all=reshape(p_loc,3,size(p_loc)[2]*size(p_loc)[3]) # platform positions in xyz; for each platform, its position at each pulse (PRI) is plotted; output loops over platforms first, then slow-time
    gr()
    platform_labels=Array{String}(undef,1,size(orbit_pos)[2])
    for i=1:size(orbit_pos)[2];platform_labels[i]=string("platform-",i);end
    display(plot(orbit_time,orbit_pos[1,:,:]'/1e3,xaxis=("time (sec)"),ylabel=("ECI X position (km)"),size=(1600,900),labels=platform_labels)) # plot the ECI orbit in the limited time range
    display(plot(orbit_time,orbit_pos[2,:,:]'/1e3,xaxis=("time (sec)"),ylabel=("ECI Y position (km)"),size=(1600,900),labels=platform_labels)) # plot the ECI orbit in the limited time range
    display(plot(orbit_time,orbit_pos[3,:,:]'/1e3,xaxis=("time (sec)"),ylabel=("ECI Z position (km)"),size=(1600,900),labels=platform_labels)) # plot the ECI orbit in the limited time range
    plotly()
    display(scatter(orbit_pos_all[1,:],orbit_pos_all[2,:],orbit_pos_all[3,:],leg=false,camera=(20,40),markersize=3,xlabel=coords[1],ylabel=coords[2],zlabel=coords[3],title="Platform Positions at Each Pulse",size=(1600,900))) #display  position of each platform at each pulse in 3D
    display(scatter(t_loc[1,:],t_loc[2,:],t_loc[3,:],leg=false,camera=(20,40),markersize=2,xlabel=coords[1],ylabel=coords[2],zlabel=coords[3],title="Target Locations",size=(1600,900))) #display grid in 3D
    #DISPLAY PLATFORM AND TARGET LOCATIONS ON THE SAME PLOT
    scatter(t_loc[1,:],t_loc[2,:],t_loc[3,:],leg=false,camera=(20,40),markersize=1,size=(1600,900)) #display grid in 3D
    #scatter!([avg_peg.pegLat],[avg_peg.pegLon],[0],markersize=3)
    display(scatter!(orbit_pos_all[1,:],orbit_pos_all[2,:],orbit_pos_all[3,:],leg=false,camera=(20,40),markersize=1,xlabel=coords[1],ylabel=coords[2],zlabel=coords[3],title="Platforms and Targets")) #display grid in 3D
    # DISPLAY SCENE
    display(scatter(s_loc[1,:],s_loc[2,:],s_loc[3,:],leg=false,camera=(20,40),markersize=0.3,xlabel=coords[1],ylabel=coords[2],zlabel=coords[3],title="Scene Pixel Locations",size=(1600,900))) #display grid in 3D
    #DISPLAY PLATFORM AND TARGET AND SCENE ON THE SAME PLOT
    scatter(t_loc[1,:],t_loc[2,:],t_loc[3,:],leg=false,camera=(20,40),markersize=1,size=(1600,900)) #display grid in 3D
    scatter!(orbit_pos_all[1,:],orbit_pos_all[2,:],orbit_pos_all[3,:],markersize=1) #display grid in 3D
    display(scatter!(s_loc[1,:],s_loc[2,:],s_loc[3,:],markersize=0.3,xlabel=coords[1],ylabel=coords[2],zlabel=coords[3],title="Platforms and Targets and Scene")) #display grid in 3D
end

function plot_tomogram(PSF_image_point,display_tomograms,image_1xN,image_3D,s_loc_1,s_loc_2,s_loc_3,s_loc_3xN,t_loc_1,t_loc_2,t_loc_3,coords,mode,scene_axis11,scene_axis22,scene_axis33,PSF_cuts,PSF_metrics)
    image_3D=image_3D/maximum(image_3D)
    brightest=maximum(image_3D)
    faintest=minimum(image_3D)
    Ns_1=length(s_loc_1)
    Ns_2=length(s_loc_2)
    Ns_3=length(s_loc_3)
    if mode==1;mode_txt="SAR-SISO";elseif mode==2;mode_txt="SIMO";elseif mode==3;mode_txt="MIMO";end
    if display_tomograms==1
        gr()
        if PSF_image_point==3 # center of scene
            k1=Int(ceil(Ns_1/2))
            k2=Int(ceil(Ns_2/2))
            k3=Int(ceil(Ns_3/2))
        elseif PSF_image_point==1 # peak location
            max_ind=findall(image_3D .==maximum(image_3D))
            k1=max_ind[1][1]
            k2=max_ind[1][2]
            k3=max_ind[1][3]
        elseif PSF_image_point==2 # target location
            k1=findall(t_loc_1 .==s_loc_1);k1=k1[1]
            k2=findall(t_loc_2 .==s_loc_2);k2=k2[1]
            k3=findall(t_loc_3 .==s_loc_3);k3=k3[1]
        end
        col=[:black,:yellow,:red]
        psize=(900,900)
        if Ns_2>1 && Ns_3>1
            display(heatmap(s_loc_2,s_loc_3,image_3D[k1,:,:]',xguidefontsize=20,yguidefontsize=20,xtickfontsize=18,ytickfontsize=18,titlefont=24;ylabel=coords[3]*" (m)",xlabel=coords[2]*" (m)",title="2D Tomogram ("*mode_txt*" mode)",
            c=:thermal,xticks=s_loc_2[1]:10:s_loc_2[end],yticks=s_loc_3[1]:10:s_loc_3[end],clims=(faintest,brightest),size=psize));savefig("tomograms/tomogram_"*mode_txt*".png") #title="2D Image at Loc-1="*string(s_loc_1[k1])
            if PSF_metrics
                heatmap(s_loc_2,s_loc_3,image_3D[k1,:,:]',xguidefontsize=20,yguidefontsize=20,xtickfontsize=18,ytickfontsize=18,titlefont=24;ylabel=coords[3]*" (m)",xlabel=coords[2]*" (m)",title="2D Tomogram ("*mode_txt*" mode)",
                c=:thermal,xticks=s_loc_2[1]:10:s_loc_2[end],yticks=s_loc_3[1]:10:s_loc_3[end],clims=(faintest,brightest),size=psize)
                if PSF_cuts==1
                    plot!(s_loc_2[k2]*ones(Ns_3),s_loc_3,lc=[:white],leg=false)
                    display(plot!(s_loc_2,s_loc_3[k3]*ones(Ns_2),lc=[:white],leg=false))
                elseif PSF_cuts==2;display(plot!(scene_axis22,scene_axis33,lc=[:white],leg=false));end
            end
        end
        if Ns_1>1 && Ns_3>1
            display(heatmap(s_loc_1,s_loc_3,image_3D[:,k2,:]',xguidefontsize=20,yguidefontsize=20,xtickfontsize=18,ytickfontsize=18,titlefont=24;ylabel=coords[3]*" (m)",xlabel=coords[1]*" (m)",title="2D Tomogram ("*mode_txt*" mode)",
            c=:thermal,xticks=s_loc_1[1]:10:s_loc_1[end],yticks=s_loc_3[1]:10:s_loc_3[end],clims=(faintest,brightest),size=psize));savefig("tomograms/tomogram_"*mode_txt*".png") #title="2D Image at Loc-2="*string(s_loc_2[k2])
            if PSF_metrics
                heatmap(s_loc_1,s_loc_3,image_3D[:,k2,:]',xguidefontsize=20,yguidefontsize=20,xtickfontsize=18,ytickfontsize=18,titlefont=24;ylabel=coords[3]*" (m)",xlabel=coords[1]*" (m)",title="2D Tomogram ("*mode_txt*" mode)",
                c=:thermal,xticks=s_loc_1[1]:10:s_loc_1[end],yticks=s_loc_3[1]:10:s_loc_3[end],clims=(faintest,brightest),size=psize)
                if PSF_cuts==1
                    plot!(s_loc_1[k1]*ones(Ns_3),s_loc_3,lc=[:white],leg=false)
                    display(plot!(s_loc_1,s_loc_3[k3]*ones(Ns_1),lc=[:white],leg=false))
                elseif PSF_cuts==2;display(plot!(scene_axis11,scene_axis33,lc=[:white],leg=false));end
            end
        end
        if Ns_1>1 && Ns_2>1
            display(heatmap(s_loc_2,s_loc_1,image_3D[:,:,k3],xguidefontsize=20,yguidefontsize=20,xtickfontsize=18,ytickfontsize=18,titlefont=24;ylabel=coords[1]*" (m)",xlabel=coords[2]*" (m)",title="2D Tomogram ("*mode_txt*" mode)",
            c=:thermal,xticks=s_loc_2[1]:10:s_loc_2[end],yticks=s_loc_1[1]:10:s_loc_1[end],clims=(faintest,brightest),size=psize));savefig("tomograms/tomogram_"*mode_txt*".png") #title="2D Image at Loc-3="*string(s_loc_3[k3])
            if PSF_metrics
                heatmap(s_loc_2,s_loc_1,image_3D[:,:,k3],xguidefontsize=20,yguidefontsize=20,xtickfontsize=18,ytickfontsize=18,titlefont=24;ylabel=coords[1]*" (m)",xlabel=coords[2]*" (m)",title="2D Tomogram ("*mode_txt*" mode)",
                c=:thermal,xticks=s_loc_2[1]:10:s_loc_2[end],yticks=s_loc_1[1]:10:s_loc_1[end],clims=(faintest,brightest),size=psize)
                if PSF_cuts==1
                    plot!(s_loc_2[k2]*ones(Ns_1),s_loc_1,lc=[:white],leg=false)
                    display(plot!(s_loc_2,s_loc_1[k1]*ones(Ns_2),lc=[:white],leg=false))
                elseif PSF_cuts==2;display(plot!(scene_axis22,scene_axis11,lc=[:white],leg=false));end
            end
        end
    elseif display_tomograms==2
        gr()
        for k=1:Ns_3 # height slices from the scene
            display(heatmap(s_loc_2,s_loc_1,image_3D[:,:,k],ylabel=coords[1],xlabel=coords[2],title="2D Image at Loc-3="*string(s_loc_3[k]),c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,1200))) #aspect_ratio=:equal
        end
        for k=1:Ns_1 # latitude slices from the scene
            display(heatmap(s_loc_3,s_loc_2,image_3D[k,:,:],ylabel=coords[2],xlabel=coords[3],title="2D Image at Loc-1="*string(s_loc_1[k]),c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,1200))) #aspect_ratio=:equal
        end
        for k=1:Ns_2 # longitude slices from the scene
            display(heatmap(s_loc_3,s_loc_1,image_3D[:,k,:],ylabel=coords[1],xlabel=coords[3],title="2D Image at Loc-2="*string(s_loc_2[k]),c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,1200))) #aspect_ratio=:equal
        end
    elseif display_tomograms==3
        gr()
        #display(scatter(s_loc_3xN[1,:],s_loc_3xN[2,:],s_loc_3xN[3,:],marker_z=image_1xN/maximum(image_1xN),leg=false,camera=(20,40),markersize=1,markerstrokewidth=0,xlabel=coords[1],ylabel=coords[2],zlabel=coords[3],title="3D Image",size=(1600,900))) #display grid in 3D
        #display(scatter(image_3D,size=(1600,1200)))
        for i=30:5:60;display(scatter(image_3D,camera=(i,30),size=(1600,1200)));end;for i=30:5:60;display(scatter(image_3D,camera=(60,i),size=(1600,1200)));end;for i=60:-5:30;display(scatter(image_3D,camera=(i,60),size=(1600,1200)));end
    end
end

function plot_input_scene(inputscene_3D,s_loc_1,s_loc_2,s_loc_3,coords)
    brightest=maximum(inputscene_3D)
    faintest=minimum(inputscene_3D)
    Ns_1=length(s_loc_1)
    Ns_2=length(s_loc_2)
    Ns_3=length(s_loc_3)
    k1=Int(ceil(Ns_1/2))
    k2=Int(ceil(Ns_2/2))
    k3=Int(ceil(Ns_3/2))
    gr()
    #display(scatter(inputscene_3D,size=(1600,1200)))
    if Ns_2>1 && Ns_3>1;display(heatmap(s_loc_2,s_loc_3,inputscene_3D[k1,:,:]',ylabel=coords[3],xlabel=coords[2],title="2D Input Scene at Loc-1="*string(s_loc_1[k1]),c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,1200)));end #aspect_ratio=:equal
    if Ns_1>1 && Ns_3>1;display(heatmap(s_loc_1,s_loc_3,inputscene_3D[:,k2,:]',ylabel=coords[3],xlabel=coords[1],title="2D Input Scene at Loc-2="*string(s_loc_2[k2]),c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,1200)));end #aspect_ratio=:equal
    if Ns_1>1 && Ns_2>1;display(heatmap(s_loc_2,s_loc_1,inputscene_3D[:,:,k3],ylabel=coords[1],xlabel=coords[2],title="2D Input Scene at Loc-3="*string(s_loc_3[k3]),c=cgrad([:black,:white]),clims=(faintest,brightest),size=(1600,1200)));end #aspect_ratio=:equal
end

end
