module Data_Plotting

using Plots
using Parameters
##
c = 299792458
plotly()

##
"""
Plotting the Doppler frequencies and Slant ranges

"""
function plot_DoppFreq(fDopp, Np::Int64=1)
    for i=1:Np
        if i==Np
            if Np==1
                p1 = plot(fDopp[i,:],ylabel="Doppler frequency (Hz)",xlabel="Azimuth/Along-track (samples)",title = "Doppler frequency versus Azimuth/Along-track")
            else
                p1 = plot!(fDopp[i,:],ylabel="Doppler frequency (Hz)",xlabel="Azimuth/Along-track (samples)",title = "Doppler frequency versus Azimuth/Along-track")
            end
            display(p1)
        elseif i==1
            p1 = plot(fDopp[i,:],ylabel="Doppler frequency (Hz)",xlabel="Azimuth/Along-track (samples)",title = "Doppler frequency versus Azimuth/Along-track")
        else
            p1 = plot!(fDopp[i,:],ylabel="Doppler frequency (Hz)",xlabel="Azimuth/Along-track (samples)",title = "Doppler frequency versus Azimuth/Along-track")
        end
        #png("../Plots/Dopp_freq.png")
        #savefig(p1,"/Users/joshil/Documents/Presentations/plots/plot_tests.png")
    end   
end

function plot_Slrng(Slrng,Np)
    for i=1:Np
        if i==Np
            if Np==1
                p1 = plot(Slrng[i,:],ylabel="Target slant range (m)",xlabel="Azimuth/Along-track(samples)",title = "Target slant range versus Azimuth/Along-track")
            else
                p1 = plot!(Slrng[i,:],ylabel="Target slant range (m)",xlabel="Azimuth/Along-track (samples)",title = "Target slant range versus Azimuth/Along-track")
            end
            display(p1)
        elseif i==1
            p1 = plot(Slrng[i,:],ylabel="Target slant range (m)",xlabel="Azimuth/Along-track (samples)",title = "Target slant range versus Azimuth/Along-track")
        else
            p1 = plot!(Slrng[i,:],ylabel="Target slant range (m)",xlabel="Azimuth/Along-track (samples)",title = "Target slant range versus Azimuth/Along-track")
        end
    end
    #display(plot(Slrng[4,:],ylabel="Slant range (m)",xlabel="along-track (samples)",title = "" ))
end

##
"""
Plotting SAR RAW received signal

"""
function plot_raw_received_data(rawdata, plot_tx_platform_idx, plot_rx_platform_idx, Ncr, Veff, slow_time, params_densesim, norm_flag, fasttime_lim)
    @unpack mode, fs, SAR_duration, SAR_start_time   = params_densesim

    x_axis_fasttime 				= ((1:Ncr).*(c/(fs)) )/1000
    y_axis_slowtime 				= ((slow_time.-(SAR_start_time+SAR_duration/2)) .* Veff[plot_rx_platform_idx,:])/1000

    if(norm_flag == 1)
        val_max,ind_max = findmax(abs.(rawdata))
        plot_data = abs.(rawdata) ./ val_max
        c_lim=(-50,0)
    else
        plot_data = abs.(rawdata)
        c_lim=(-250,-150)
    end

    if mode == 3
        display( heatmap(x_axis_fasttime,y_axis_slowtime, 20*log10.(plot_data[:,plot_tx_platform_idx,plot_rx_platform_idx,:]),ylabel="Azimuth/Along-track (km)",xlabel="Slant range/Cross-track (km)",title = "Raw output signal" )); #, clim=(-80,40),aspect_ratio=:equal
    else
    	display( heatmap(x_axis_fasttime[fasttime_lim[1]:fasttime_lim[2]],y_axis_slowtime, 20*log10.(plot_data[:,plot_rx_platform_idx,fasttime_lim[1]:fasttime_lim[2]]),ylabel="Azimuth/Along-track (km)",xlabel="Slant range/Cross-track (km)",title = "Raw output signal", clim=c_lim )); #,aspect_ratio=:equal
    end

end

##
"""
Plotting range compresed received signal

"""
function plot_rangecompressed_received_data(RC_signal, plot_tx_platform_idx, plot_rx_platform_idx, Ncr, Veff, slow_time, params_densesim, norm_flag, fasttime_lim)
    @unpack mode, fs, SAR_duration, SAR_start_time   = params_densesim

    x_axis_fasttime 				= ((1:Ncr).*(c/(fs)) )/1000
    y_axis_slowtime 				= ((slow_time.-(SAR_start_time+SAR_duration/2)) .* Veff[plot_rx_platform_idx,:])/1000

    if(norm_flag == 1)
        val_max,ind_max = findmax(abs.(RC_signal))
        plot_data = abs.(RC_signal) ./ val_max
        c_lim=(-50,0)
    else
        plot_data = abs.(RC_signal)
        c_lim=(-200,-100)
    end

    if mode ==3
        display( heatmap(x_axis_fasttime,y_axis_slowtime, 20*log10.(plot_data[:,plot_tx_platform_idx,plot_rx_platform_idx,:]),ylabel="Azimuth/Along-track (km)",xlabel="Slant range/Cross-track (km)",title = "Range compressed output signal" )); #, clim=(-80,40),aspect_ratio=:equal
    else
    	display( heatmap(x_axis_fasttime[fasttime_lim[1]:fasttime_lim[2]],y_axis_slowtime,20*log10.(plot_data[:,plot_rx_platform_idx,fasttime_lim[1]:fasttime_lim[2]]),ylabel="Azimuth/Along-track (km)",xlabel="Slant range/Cross-track (km)",title = "Range compressed output signal", clim=c_lim));#,,xlim=(9,10))); #, clim=(-80,40),aspect_ratio=:equal
    end

end

##
"""
Plotting azimuth compresed received signal

"""
function plot_azimuthcompressed_received_data(AC_signal, plot_tx_platform_idx, plot_rx_platform_idx, Ncr, Veff, slow_time, params_densesim, norm_flag, fasttime_lim)
    @unpack mode, fs, SAR_duration, SAR_start_time   = params_densesim

    x_axis_fasttime 				= ((1:Ncr).*(c/(fs)) )/1000
    y_axis_slowtime 				= ((slow_time.-(SAR_start_time+SAR_duration/2)) .* Veff[plot_rx_platform_idx,:])/1000

    if(norm_flag == 1)
        val_max,ind_max = findmax(abs.(AC_signal))
        plot_data = abs.(AC_signal) ./ val_max
        c_lim=(-100,0)
    else
        plot_data = abs.(AC_signal)
        c_lim=(-135,-35)
    end

    if mode ==3
        display( heatmap(x_axis_fasttime,y_axis_slowtime, 20*log10.(plot_data[:,plot_tx_platform_idx,plot_rx_platform_idx,:]),ylabel="Azimuth/Along-track (km)",xlabel="Slant range/Cross-track (km)",title = "Azimuth compressed output signal" )); #, clim=(-80,40),aspect_ratio=:equal
    else
    	display( heatmap(x_axis_fasttime[fasttime_lim[1]:fasttime_lim[2]],y_axis_slowtime,20*log10.(plot_data[:,plot_rx_platform_idx,fasttime_lim[1]:fasttime_lim[2]]),ylabel="Azimuth/Along-track (km)",xlabel="Slant range/Cross-track (km)",title = "Azimuth compressed output signal",clim=c_lim));#,,xlim=(9,10))); #, clim=(-80,40),aspect_ratio=:equal
    end

end

##
"""
Plotting Back projection algorithm output signal

"""
function plot_tomo_output(image_3D, params_densesim, norm_flag, plot_idx, savepath, dB_level = 10, title_string=" ", save_name_ext=" ", heights_scene =0)
    @unpack s_loc_1, s_loc_2 = params_densesim
    xlimp = [minimum(s_loc_1) maximum(s_loc_1)]
    ylimp = [minimum(s_loc_2) maximum(s_loc_2)]
    if heights_scene == 0
        @unpack  s_loc_3 = params_densesim
        zlimp = [minimum(s_loc_3) maximum(s_loc_3)]
    else
        s_loc_3 = heights_scene
        zlimp = [minimum(params_densesim.s_loc_3) maximum(params_densesim.s_loc_3)]
    end

    xlimp = [-20 20]
    ylimp = [-20 20]
    zlimp = [-20 20]

    if(norm_flag == 1)
        val_max,ind_max = findmax(abs.(image_3D))
        plot_data = abs.(image_3D) ./ val_max
        c_lim=(-60,0)
    else
        plot_data = abs.(image_3D)
        c_lim=(-200,0)
    end
    gr()
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p1 = (heatmap(s_loc_1,s_loc_2,dB_level .* log10.(transpose(plot_data[:,:,plot_idx[3]])),ylabel="Ground range (C axis) [m]",xlabel="Along-track distance (S axis) [m]",title =title_string, clim=c_lim,
    #c=:imola, xlim=(xlimp[1],xlimp[2]),ylim=(ylimp[1],ylimp[2]),tickfont=font(18), xtickfont=font(18), ytickfont=font(18), guidefont=font(18), titlefontsize=20, size=(650,500) ))
    xlim=(xlimp[1],xlimp[2]),ylim=(ylimp[1],ylimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p11 = plot(p1,p0,layout=l)
    savefig(p11, savepath*title_string[1:12]*save_name_ext*"_1.png")
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p2 = (heatmap(s_loc_1,s_loc_3,dB_level .* log10.(transpose(plot_data[:,plot_idx[2],:])),ylabel="Height (H axis) [m]",xlabel="Along-track distance (S axis) [m]",title = title_string, clim=c_lim,
    xlim=(xlimp[1],xlimp[2]),ylim=(zlimp[1],zlimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p22 = plot(p2,p0,layout=l)
    savefig(p22, savepath*title_string[1:12]*save_name_ext*"_2.png")
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p3 = (heatmap(s_loc_2,s_loc_3,dB_level .* log10.(transpose(plot_data[plot_idx[1],:,:])),ylabel="Height (H axis) [m]",xlabel="Ground range (C axis) [m]",title = title_string, clim=c_lim,
    xlim=(ylimp[1],ylimp[2]),ylim=(zlimp[1],zlimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p33 = plot(p3,p0,layout=l)
    savefig(p33, savepath*title_string[1:12]*save_name_ext*"_3.png")
    p4 = (plot(s_loc_2,dB_level .* log10.(plot_data[plot_idx[1],:,plot_idx[3]]),xlabel="Ground range (C axis) [m]",ylabel="Output level [dB]",title = title_string,ylim=c_lim,
    xlim=(ylimp[1],ylimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false, linewidth=3))
    savefig(p4, savepath*title_string[1:12]*save_name_ext*"_4.png")
    p5 = (plot(s_loc_1,dB_level .* log10.(plot_data[:,plot_idx[2],plot_idx[3]]),xlabel="Along-track distance (S axis) [m]",ylabel="Output level [dB]",title = title_string,ylim=c_lim,
    xlim=(xlimp[1],xlimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15,legend=false, linewidth=3))
    savefig(p5, savepath*title_string[1:12]*save_name_ext*"_5.png")
    p6 = (plot(dB_level .* log10.(plot_data[plot_idx[1],plot_idx[2],:]),s_loc_3,xlabel="Output level [dB]",ylabel="Height (H axis) [m]",title = title_string,xlim=c_lim,
    ylim=(zlimp[1],zlimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false, linewidth=3))
    savefig(p6, savepath*title_string[1:12]*save_name_ext*"_6.png")

end

function plot_tomo_output_lin(image_3D, params_densesim, norm_flag, plot_idx, savepath, title_string=" ", save_name_ext=" ", heights_scene =0)
    @unpack s_loc_1, s_loc_2 = params_densesim
    xlimp = [minimum(s_loc_1) maximum(s_loc_1)]
    ylimp = [minimum(s_loc_2) maximum(s_loc_2)]
    if heights_scene == 0
        @unpack  s_loc_3 = params_densesim
        zlimp = [minimum(s_loc_3) maximum(s_loc_3)]
    else
        s_loc_3 = heights_scene
        zlimp = [minimum(params_densesim.s_loc_3) maximum(params_densesim.s_loc_3)]
    end

    xlimp = [-30 30]
    ylimp = [-30 30]
    zlimp = [-20 20]

    if(norm_flag == 1)
        val_max,ind_max = findmax(abs.(image_3D))
        plot_data = abs.(image_3D) ./ val_max
        c_lim=(0,1)
    else
        plot_data = abs.(image_3D)
        c_lim=(0,1)
    end
    gr()
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p1 = (heatmap(s_loc_1,s_loc_2,(transpose(plot_data[:,:,plot_idx[3]])),ylabel="Ground range (C axis) [m]",xlabel="Along-track distance (S axis) [m]",title =title_string, clim=c_lim,
    #c=:imola, xlim=(xlimp[1],xlimp[2]),ylim=(ylimp[1],ylimp[2]),tickfont=font(18), xtickfont=font(18), ytickfont=font(18), guidefont=font(18), titlefontsize=20, size=(650,500) ))
    xlim=(xlimp[1],xlimp[2]),ylim=(ylimp[1],ylimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p11 = plot(p1,p0,layout=l)
    savefig(p11, savepath*title_string[1:12]*save_name_ext*"_91.png")
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p2 = (heatmap(s_loc_1,s_loc_3,(transpose(plot_data[:,plot_idx[2],:])),ylabel="Height (H axis) [m]",xlabel="Along-track distance (S axis) [m]",title = title_string, clim=c_lim,
    xlim=(xlimp[1],xlimp[2]),ylim=(zlimp[1],zlimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p22 = plot(p2,p0,layout=l)
    savefig(p22, savepath*title_string[1:12]*save_name_ext*"_92.png")
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p3 = (heatmap(s_loc_2,s_loc_3,(transpose(plot_data[plot_idx[1],:,:])),ylabel="Height (H axis) [m]",xlabel="Ground range (C axis) [m]",title = title_string, clim=c_lim,
    xlim=(ylimp[1],ylimp[2]),ylim=(zlimp[1],zlimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p33 = plot(p3,p0,layout=l)
    savefig(p33, savepath*title_string[1:12]*save_name_ext*"_93.png")
    p4 = (plot(s_loc_2,(plot_data[plot_idx[1],:,plot_idx[3]]),xlabel="Ground range (C axis) [m]",ylabel="Output level [dB]",title = title_string,ylim=c_lim,
    xlim=(ylimp[1],ylimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false, linewidth=3))
    savefig(p4, savepath*title_string[1:12]*save_name_ext*"_94.png")
    p5 = (plot(s_loc_1,(plot_data[:,plot_idx[2],plot_idx[3]]),xlabel="Along-track distance (S axis) [m]",ylabel="Output level [dB]",title = title_string,ylim=c_lim,
    xlim=(xlimp[1],xlimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15,legend=false, linewidth=3))
    savefig(p5, savepath*title_string[1:12]*save_name_ext*"_95.png")
    p6 = (plot((plot_data[plot_idx[1],plot_idx[2],:]),s_loc_3,xlabel="Output level [dB]",ylabel="Height (H axis) [m]",title = title_string,xlim=c_lim,
    ylim=(zlimp[1],zlimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false, linewidth=3))
    savefig(p6, savepath*title_string[1:12]*save_name_ext*"_96.png")

end


function plot_tomo_output_phase(image_3D_ph, params_densesim, plot_idx, savepath, title_string=" ", save_name_ext=" ", heights_scene =0)
    @unpack s_loc_1, s_loc_2 = params_densesim
    xlimp = [minimum(s_loc_1) maximum(s_loc_1)]
    ylimp = [minimum(s_loc_2) maximum(s_loc_2)]
    if heights_scene == 0
        @unpack  s_loc_3 = params_densesim
        zlimp = [minimum(s_loc_3) maximum(s_loc_3)]
    else
        s_loc_3 = heights_scene
        zlimp = [minimum(params_densesim.s_loc_3) maximum(params_densesim.s_loc_3)]
    end

    xlimp = [-20 20]
    ylimp = [-20 20]
    zlimp = [-20 20]

    c_lim=(-3.2,3.2)

    plot_data = image_3D_ph

    gr()
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p1 = (heatmap(s_loc_1,s_loc_2,(transpose(plot_data[:,:,plot_idx[3]])),ylabel="Ground range (C axis) [m]",xlabel="Along-track distance (S axis) [m]",title =title_string, clim=c_lim,
    #c=:imola, xlim=(xlimp[1],xlimp[2]),ylim=(ylimp[1],ylimp[2]),tickfont=font(18), xtickfont=font(18), ytickfont=font(18), guidefont=font(18), titlefontsize=20, size=(650,500) ))
    xlim=(xlimp[1],xlimp[2]),ylim=(ylimp[1],ylimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p11 = plot(p1,p0,layout=l)
    savefig(p11, savepath*title_string[1:12]*save_name_ext*"_1.png")
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p2 = (heatmap(s_loc_1,s_loc_3,(transpose(plot_data[:,plot_idx[2],:])),ylabel="Height (H axis) [m]",xlabel="Along-track distance (S axis) [m]",title = title_string, clim=c_lim,
    xlim=(xlimp[1],xlimp[2]),ylim=(zlimp[1],zlimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p22 = plot(p2,p0,layout=l)
    savefig(p22, savepath*title_string[1:12]*save_name_ext*"_2.png")
    l = @layout[grid(1,1) a{0.0001w}]
    p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
    p3 = (heatmap(s_loc_2,s_loc_3,(transpose(plot_data[plot_idx[1],:,:])),ylabel="Height (H axis) [m]",xlabel="Ground range (C axis) [m]",title = title_string, clim=c_lim,
    xlim=(ylimp[1],ylimp[2]),ylim=(zlimp[1],zlimp[2]),tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
    p33 = plot(p3,p0,layout=l)
    savefig(p33, savepath*title_string[1:12]*save_name_ext*"_3.png")
    p4 = (plot(s_loc_2,(plot_data[plot_idx[1],:,plot_idx[3]]),xlabel="Ground range (C axis) [m]",ylabel="Output level [dB]",title = title_string,ylim=c_lim,
    xlim=(ylimp[1],ylimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false, linewidth=3))
    savefig(p4, savepath*title_string[1:12]*save_name_ext*"_4.png")
    p5 = (plot(s_loc_1,(plot_data[:,plot_idx[2],plot_idx[3]]),xlabel="Along-track distance (S axis) [m]",ylabel="Output level [dB]",title = title_string,ylim=c_lim,
    xlim=(xlimp[1],xlimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15,legend=false, linewidth=3))
    savefig(p5, savepath*title_string[1:12]*save_name_ext*"_5.png")
    p6 = (plot((plot_data[plot_idx[1],plot_idx[2],:]),s_loc_3,xlabel="Output level [dB]",ylabel="Height (H axis) [m]",title = title_string,xlim=c_lim,
    ylim=(zlimp[1],zlimp[2]), tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false, linewidth=3))
    savefig(p6, savepath*title_string[1:12]*save_name_ext*"_6.png")

end

##
end
