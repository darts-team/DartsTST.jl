
using Plots
using LinearAlgebra
using Statistics
using Interpolations
using Measures
using StatsPlots
using NaNStatistics
using SpecialFunctions
using Distributions
using CurveFit
using SpecialFunctions
using HypergeometricFunctions
using ColorSchemes
using GeoArrays
using ArchGDAL
using DelimitedFiles
using DSP

function plot_image(x_axis,y_axis,Data,unit_flag, savepath, savename, figure_title)

    if ~ispath(savepath)
        mkdir(savepath)
    end    

    if unit_flag == "log"
        p1=(heatmap(x_axis,y_axis,10 .* log10.(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]",title=figure_title,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #500,360
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "lin"
        p1=(heatmap(x_axis,y_axis,(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "phase"
        p1=(heatmap(x_axis,y_axis,(Data),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200
        savefig(p1, savepath*savename*".png")
    end
end

function plot_image2(x_axis,y_axis,Data,unit_flag, savepath, savename, figure_title, colorbar)

    if ~ispath(savepath)
        mkdir(savepath)
    end    

    if unit_flag == "log"
        p1=(heatmap(x_axis,y_axis,10 .* log10.(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]",title=figure_title,c=colorbar,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #500,360
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "lin"
        p1=(heatmap(x_axis,y_axis,(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title,c=colorbar,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "phase"
        p1=(heatmap(x_axis,y_axis,(Data),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title, #c=:twilight,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200
        savefig(p1, savepath*savename*".png")
    end
end

function plot_profile(x_axis, Data1, Data2, savepath, savename, xlabel, ylabel, figure_title, label1, label2)

    if Data2==""
        p1=(plot(x_axis, Data1,xlabel=xlabel,ylabel=ylabel,title=figure_title,legend=:topleft, lc=:black, label=label1,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #500,360
        savefig(p1, savepath*savename*".png")
    else
        p1=(plot(x_axis, Data1,xlabel=xlabel,ylabel=ylabel,title=figure_title, legend=:topleft, lc=:black, label=label1,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #500,360
        p1=(plot!(x_axis, Data2,xlabel=xlabel,ylabel=ylabel,title=figure_title, legend=:topleft,lc=:blue,label=label2,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #500,360
        savefig(p1, savepath*savename*".png")
    end

end

function plot_histogram_only(plot_var_hist,xlabel_text, savepath, savename)

    bins_range              = range(minimum(plot_var_hist),maximum(plot_var_hist), length=51)

    # histogram plot
    p1= (histogram(plot_var_hist,xlabel=xlabel_text,ylabel="Count",title = "Histogram", colorbar=true, bins= bins_range,
        tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
    savefig(p1, savepath*savename*"_Hist_plot"*".png")
    
    return plot_var_hist

end


function plot_histogram_statistics_amplitude(plot_var_hist, savepath, savename)

    if ~ispath(savepath)
        mkdir(savepath)
    end    

    mean_p                 = mean(plot_var_hist)
    plot_var_hist          = plot_var_hist ./ mean_p
    bins_range             = range(0,maximum(plot_var_hist), length=51)

    # histogram plot
    p1= (histogram(plot_var_hist,xlabel="Amplitude",ylabel="Count",title = "Histogram", colorbar=true, bins= bins_range,
        tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
    savefig(p1, savepath*savename*"_Hist_plot1"*".png")
    
    N, bin                 = histcountindices(plot_var_hist,bins_range)
    pdf1                   = (N / sum(N)) / (bins_range[2]-bins_range[1])
    #axis_plot              = range(0+(maximum(plot_var_hist)/50),maximum(plot_var_hist)-(maximum(plot_var_hist)/50), length=50)
    axis_plot              = range(0+(maximum(plot_var_hist)/(length(bins_range)-1)),maximum(plot_var_hist)-(maximum(plot_var_hist)/(length(bins_range)-1)), length=(length(bins_range)-1))

    theory_val             = zeros(length(axis_plot),1)
    sigma                  = 1 * sqrt(2/pi)
    for i =1:length(axis_plot)
        #theory_val[i]      = 2 .* pdf(Chisq(2),axis_plot[i])
        theory_val[i]      = (axis_plot[i] * exp((-axis_plot[i]^2)/(2*(sigma^2)))) / (sigma^2)
    end

    (plot(axis_plot./sigma , sigma .*pdf1, ylim=(0,1.1), xlim=(0,maximum(plot_var_hist)),xlabel="Amplitude / sigma",ylabel="sigma * PDF",title = ""
        ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
    p3 = (plot!(axis_plot./sigma, sigma .*theory_val, ylim=(0,1.1), xlim=(0,maximum(plot_var_hist)),xlabel="Amplitude / sigma",ylabel="sigma * PDF",title = ""
        ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
    savefig(p3, savepath*savename*"_Hist_plot2"*".png")

    return plot_var_hist 

end


function plot_histogram_statistics_power(plot_var_hist, savepath, savename)

    mean_p                  = mean(plot_var_hist)
    plot_var_hist           = plot_var_hist ./ mean_p
    bins_range              = range(0,maximum(plot_var_hist), length=51)

    # histogram plot
    p1= (histogram(plot_var_hist,xlabel="Amplitude",ylabel="Count",title = "Histogram", colorbar=true, bins= bins_range,
        tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
    savefig(p1, savepath*savename*"_Hist_plot1"*".png")
    
    N, bin                  = histcountindices(plot_var_hist,bins_range)
    pdf1                    = (N / sum(N)) / (bins_range[2]-bins_range[1])
    #axis_plot               = range(0+(maximum(plot_var_hist)/50),maximum(plot_var_hist)-(maximum(plot_var_hist)/50), length=50)
    axis_plot              = range(0+(maximum(plot_var_hist)/(length(bins_range)-1)),maximum(plot_var_hist)-(maximum(plot_var_hist)/(length(bins_range)-1)), length=(length(bins_range)-1))

    theory_val              = zeros(length(axis_plot),1)
    N                       = 1
    for i =1:length(axis_plot)
        #theory_val[i]       = 2 .* pdf(Chisq(2),axis_plot[i])
        theory_val[i]       = ((axis_plot[i]^(N-1)) * exp((-axis_plot[i])/(1))) / ((1^N) * gamma(1))
    end
    
    (plot(axis_plot, pdf1, ylim=(0,1.1), xlim=(0,maximum(plot_var_hist)),xlabel="Power",ylabel="PDF",title = ""
        ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
    p3=(plot!(axis_plot, theory_val, ylim=(0,1.1), xlim=(0,maximum(plot_var_hist)),xlabel="Power",ylabel="PDF",title = ""
        ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
    savefig(p3, savepath*savename*"_Hist_plot2"*".png")
    
    return plot_var_hist

end


function plot_histogram_statistics_phase(plot_var_hist, savepath, savename)

    bins_range              = range(-pi,pi, length=51)

    #Phase histogram plot
    p2 = (histogram(plot_var_hist,xlabel="Phase",ylabel="Count",title = "Histogram", colorbar=true, bins=bins_range,
    tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
    savefig(p2, savepath*savename*"_Hist_plot1"*".png")

    N, bin                  = histcountindices(plot_var_hist,bins_range)
    pdf1                    = (N / sum(N)) / (bins_range[2]-bins_range[1])
    #axis_plot               = range(-pi+(pi/50),pi-(pi/50), length=50)
    axis_plot               = range(-pi+(pi/(length(bins_range)-1)),pi-(pi/(length(bins_range)-1)), length=(length(bins_range)-1))
    theory_val              = repeat([1/(2*pi)], size(axis_plot)[1],1) 

    (plot(axis_plot, pdf1, ylim=(0,1/pi), xlim=(-1.5*pi,1.5*pi),xlabel="Phase",ylabel="PDF",title = ""
        ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
    p3  =(plot!(axis_plot, theory_val, ylim=(0,1/pi), xlim=(-1.5*pi,1.5*pi),xlabel="Phase",ylabel="PDF",title = ""
        ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
    savefig(p3, savepath*savename*"_Hist_plot2"*".png")

    return plot_var_hist

end


function plot_interferogram_histogram_statistics_phase_han1(plot_var_hist, gamma_ip, L, phase_ip_ref, xlims, savepath, savename)

    bins_range              = range(-pi,pi, length=501)

    #Phase histogram plot
    p2 = (histogram(plot_var_hist,xlabel="Phase",ylabel="Count",title = "Histogram", colorbar=true, bins=bins_range, xlim=xlims,
    tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
    savefig(p2, savepath*savename*"_Hist_plot1_Phase"*".png")

    N, bin                  = histcountindices(plot_var_hist,bins_range)
    pdf1                    = (N / sum(N)) / (bins_range[2]-bins_range[1])
    axis_plot               = range(-pi+(pi/(length(bins_range)-1)),pi-(pi/(length(bins_range)-1)), length=(length(bins_range)-1))

        #Hanssen - phase
        theory_val       = zeros( length(bins_range))

        for phase_ip_idx=1:length(bins_range)
            beta_ip     = gamma_ip * cos(bins_range[phase_ip_idx] - phase_ip_ref)
            beta_sqrt_t = (1-(beta_ip^2))

            pdf_ph_t1 = ((gamma(2*L-1)) / ((gamma(L)^2) * (2^(2*(L-1))))) * ( ((((2*L-1) * beta_ip) / (beta_sqrt_t^(L+0.5))) * ( acos(-beta_ip))) + (1/beta_sqrt_t^L))

            if L == 1
                pdf_ph_t2 = 0# (1 / (2*(L-1))) 
            else
                sum_term = 0
                for i=0:L-2
                    sum_term = sum_term + ( ( (gamma(L-0.5)) * (gamma(L-i-1)) * (1+(((2*i)+1)*(beta_ip^2)))) / ( (gamma(L-0.5-i)) * (gamma(L-1)) * (beta_sqrt_t^(i+2))))
                end    
                pdf_ph_t2 = (1 / (2*(L-1))) * sum_term
            end
            theory_val[phase_ip_idx] = (((1 - gamma_ip^2)^L)/ (2*pi)) .* (pdf_ph_t1 + pdf_ph_t2)

        end

    (plot(axis_plot, pdf1, xlim=xlims,xlabel="Phase",ylabel="PDF",title = ""
        ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
    p3  =(plot!(bins_range, theory_val, xlim=xlims,xlabel="Phase",ylabel="PDF",title = ""
        ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
    savefig(p3, savepath*savename*"_PDF_comp_IntPhase"*".png")

    return plot_var_hist

end

function plot_interferogram_histogram_statistics_phase(plot_var_hist, gamma_ip, L, phase_ip_ref, xlims, savepath, savename)

    bins_range              = range(-pi,pi, length=501)

    #Phase histogram plot
    p2 = (histogram(plot_var_hist,xlabel="Phase",ylabel="Count",title = "Histogram", colorbar=true, bins=bins_range, xlim=xlims,
    tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
    savefig(p2, savepath*savename*"_Hist_plot1_Phase"*".png")

    N, bin                  = histcountindices(plot_var_hist,bins_range)
    pdf1                    = (N / sum(N)) / (bins_range[2]-bins_range[1])
    axis_plot               = range(-pi+(pi/(length(bins_range)-1)),pi-(pi/(length(bins_range)-1)), length=(length(bins_range)-1))

        #Hanssen - phase
        theory_val       = zeros( length(bins_range))

        for phase_ip_idx=1:length(bins_range)
            beta_ip     = gamma_ip * cos(bins_range[phase_ip_idx] - phase_ip_ref)
            beta_sqrt_t = (1-(beta_ip^2))
    
            Term1 = (gamma(L + 0.5) * ((1-(gamma_ip^2))^L) * abs.(gamma_ip) *  cos(bins_range[phase_ip_idx] - phase_ip_ref)) / ( 2 * sqrt(pi) * gamma(L) * (beta_sqrt_t^(L+0.5)) )
            Term2 = ((1-(gamma_ip^2))^L) / (2*pi) * _₂F₁(L,1,0.5,beta_ip^2)
    
            theory_val[phase_ip_idx] = Term1 + Term2

        end

    (plot(axis_plot, pdf1, xlim=xlims,xlabel="Phase",ylabel="PDF",title = ""
        ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
    p3  =(plot!(bins_range, theory_val, xlim=xlims,xlabel="Phase",ylabel="PDF",title = ""
        ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
    savefig(p3, savepath*savename*"_PDF_comp_IntPhase"*".png")

    return plot_var_hist

end

function plot_interferogram_histogram_statistics_magnitude(plot_var_hist, gamma_ip, L, xlims, savepath, savename)

    bins_range              = range(0,1, length=251)

    #Phase histogram plot
    p2 = (histogram(plot_var_hist,xlabel="Magnitude",ylabel="Count",title = "Histogram", colorbar=true, bins=bins_range, xlim=xlims,
    tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
    savefig(p2, savepath*savename*"_Hist_plot1_Mag"*".png")

    N, bin                  = histcountindices(plot_var_hist,bins_range)
    pdf1                    = (N / sum(N)) / (bins_range[2]-bins_range[1])
    axis_plot               = range(0+(1/(length(bins_range)-1)),1-(1/(length(bins_range)-1)), length=(length(bins_range)-1))

        #Hanssen - magnitude
        theory_val       = zeros( length(bins_range))

        for gamm_idx=1:length(bins_range)
            
            theory_val[gamm_idx]         = 2 .* (L-1) .* ((1-(gamma_ip^2))^L) .* bins_range[gamm_idx] .* ((1-(bins_range[gamm_idx]^2))^(L-2)) * _₂F₁(L,L,1,(bins_range[gamm_idx]^2) * (gamma_ip^2) )
            
        end

    (plot(axis_plot, pdf1, xlim=xlims,xlabel="Magnitude",ylabel="PDF",title = ""
        ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
    p3  =(plot!(bins_range, theory_val, xlim=xlims,xlabel="Magnitude",ylabel="PDF",title = ""
        ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm,legend=:topleft))
    savefig(p3, savepath*savename*"_PDF_comp_IntMag"*".png")

    return plot_var_hist

end

function plot_geometry_variables(ref_data, Lon_vals, Lat_vals, maxind_val_lon, maxind_val_lat, profile_flag, savepath)


    #profile_flag = 1, mean of profiles along-track
    #profile_flag = 0, 1st sample along-track

    #DEM
    DEM_data                            = ref_data[:,:,14]
    if profile_flag == 1
        DEM_region_rangeprofile         = mean(DEM_data[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        DEM_region_rangeprofile         = DEM_data[1:maxind_val_lon,1]
    end
    plot_image(Lon_vals,Lat_vals,DEM_data[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "DEM", "")
    plot_profile(Lon_vals, DEM_region_rangeprofile[:], "", savepath, "DEM_profile", "Lon [deg]", "DEM height [m]","", "", "", )

    #BRCS
    targets_ref_data                    = ref_data[:,:,13]
    if profile_flag == 1
        targets_ref_corr_rangeprofile   = mean(targets_ref_data[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        targets_ref_corr_rangeprofile   = targets_ref_data[1:maxind_val_lon,1]
    end
    plot_image(Lon_vals,Lat_vals,targets_ref_data[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "RCS", "")
    plot_profile(Lon_vals, targets_ref_corr_rangeprofile[:], "", savepath, "RCS_profile", "Lon [deg]", "Reflectivity from RCS","", "", "", )

    #Look angle
    Lookangle_P1                        = ref_data[:,:,3]
    Lookangle_P2                        = ref_data[:,:,4]
    Lookangle_Diff                      = Lookangle_P1 - Lookangle_P2

    if profile_flag == 1
        Lookangle_P1_rangeprofile       = mean(Lookangle_P1[1:maxind_val_lon,:],dims=2)
        Lookangle_P2_rangeprofile       = mean(Lookangle_P2[1:maxind_val_lon,:],dims=2)
        Lookangle_Diff_rangeprofile     = mean(Lookangle_Diff[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        Lookangle_P1_rangeprofile       =  Lookangle_P1[1:maxind_val_lon,1] 
        Lookangle_P2_rangeprofile       =  Lookangle_P2[1:maxind_val_lon,1]
        Lookangle_Diff_rangeprofile     =  Lookangle_Diff[1:maxind_val_lon,1]
    end

    plot_image(Lon_vals,Lat_vals,Lookangle_P1[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Lookangle_P1", "")
    plot_image(Lon_vals,Lat_vals,Lookangle_P2[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Lookangle_P2", "")
    plot_image(Lon_vals,Lat_vals,Lookangle_Diff[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Lookangle_Diff", "")

    plot_profile(Lon_vals, Lookangle_P1_rangeprofile[:], Lookangle_P2_rangeprofile[:], savepath, "Lookangle_profile1", "Lon [deg]", "Look Angle [deg]", "", "Platform 1", "Platform 2")
    plot_profile(Lon_vals, Lookangle_Diff_rangeprofile[:], "", savepath, "Lookangle_profile", "Lon [deg]", "Look Angle Diff [deg]", "", "Difference", "")

    #Incidence angle
    Incangle_P1                         = ref_data[:,:,5]
    Incangle_P2                         = ref_data[:,:,6]
    Incangle_Diff                       = Incangle_P1 - Incangle_P2

    if profile_flag == 1
        Incangle_P1_rangeprofile        = mean(Incangle_P1[1:maxind_val_lon,:],dims=2)
        Incangle_P2_rangeprofile        = mean(Incangle_P2[1:maxind_val_lon,:],dims=2)
        Incangle_Diff_rangeprofile      = mean(Incangle_Diff[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        Incangle_P1_rangeprofile        =  Incangle_P1[1:maxind_val_lon,1]
        Incangle_P2_rangeprofile        =  Incangle_P2[1:maxind_val_lon,1] 
        Incangle_Diff_rangeprofile      =  Incangle_Diff[1:maxind_val_lon,1] 
    end

    plot_image(Lon_vals,Lat_vals,Incangle_P1[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Incangle_P1", "")
    plot_image(Lon_vals,Lat_vals,Incangle_P2[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Incangle_P2", "")
    plot_image(Lon_vals,Lat_vals,Incangle_Diff[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Incangle_Diff", "")

    plot_profile(Lon_vals, Incangle_P1_rangeprofile[:], Incangle_P2_rangeprofile[:], savepath, "Incangle_profile1", "Lon [deg]", "Incidence Angle [deg]", "", "Platform 1", "Platform 2")
    plot_profile(Lon_vals, Incangle_Diff_rangeprofile[:], "", savepath, "Incangle_profile", "Lon [deg]", "Incidence Angle Diff [deg]", "", "Difference", "")

    #Local Incidence angle
    LocalIncangle                       = ref_data[:,:,9]
    if profile_flag == 1
        LocalIncangle_rangeprofile      = mean(LocalIncangle[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        LocalIncangle_rangeprofile      =  LocalIncangle[1:maxind_val_lon,1] 
    end
    plot_image(Lon_vals,Lat_vals,LocalIncangle[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "LocalIncangle", "")
    plot_profile(Lon_vals, LocalIncangle_rangeprofile[:], "", savepath, "LocalIncangle_profile", "Lon [deg]", "Local incidence Angle [deg]","", "", "", )

    #Range slope angle
    range_slope                         = ref_data[:,:,10]
    if profile_flag == 1
        range_slope_rangeprofile        = mean(range_slope[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        range_slope_rangeprofile        =  range_slope[1:maxind_val_lon,1]
    end
    plot_image(Lon_vals,Lat_vals,range_slope[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "RangeSlopeangle", "")
    plot_profile(Lon_vals, range_slope_rangeprofile[:], "", savepath, "RangeSlopeangle_profile", "Lon [deg]", "Range Slope Angle [deg]","", "", "", )

    #Slant range
    Slantrange_P1                       = ref_data[:,:,1]
    Slantrange_P2                       = ref_data[:,:,2]
    Slantrange_Diff                     = Slantrange_P1 - Slantrange_P2

    if profile_flag == 1
        Slantrange_P1_rangeprofile      = mean(Slantrange_P1[1:maxind_val_lon,:],dims=2)
        Slantrange_P2_rangeprofile      = mean(Slantrange_P2[1:maxind_val_lon,:],dims=2)
        Slantrange_Diff_rangeprofile    = mean(Slantrange_Diff[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        Slantrange_P1_rangeprofile      =  Slantrange_P1[1:maxind_val_lon,1] 
        Slantrange_P2_rangeprofile      =  Slantrange_P2[1:maxind_val_lon,1] 
        Slantrange_Diff_rangeprofile    =  Slantrange_Diff[1:maxind_val_lon,1] 
    end

    plot_image(Lon_vals,Lat_vals,Slantrange_P1[1:maxind_val_lon,1:maxind_val_lat]'./1e3,"lin", savepath, "Slantrange_P1", "")
    plot_image(Lon_vals,Lat_vals,Slantrange_P2[1:maxind_val_lon,1:maxind_val_lat]'./1e3,"lin", savepath, "Slantrange_P2", "")
    plot_image(Lon_vals,Lat_vals,Slantrange_Diff[1:maxind_val_lon,1:maxind_val_lat]',"phase", savepath, "Slantrange_Diff", "")

    plot_profile(Lon_vals, Slantrange_P1_rangeprofile[:], Slantrange_P2_rangeprofile[:], savepath, "Slantrange_profile1", "Lon [deg]", "Slant range [m]", "", "Platform 1", "Platform 2")
    plot_profile(Lon_vals, Slantrange_Diff_rangeprofile[:], "", savepath, "Slantrange_profile", "Lon [deg]", "Slant range [m]", "", "Difference", "")

    #Synthetic fringes
    Unwrapped_data                      = ( (4 * pi )/0.23793052222222222).*Slantrange_Diff[1:maxind_val_lon,1:maxind_val_lat]
    Wrapped_data                        = mod.(Unwrapped_data[1:maxind_val_lon,1:maxind_val_lat],2*pi) #.- pi

    plot_image(Lon_vals,Lat_vals,Unwrapped_data',"phase", savepath, "Slantrange_Diff_unwrapped", "")
    plot_image(Lon_vals,Lat_vals,Wrapped_data',"phase", savepath, "Slantrange_Diff_wrapped", "")

    if profile_flag == 1
        Wrapped_data_rangeprofile       = mean(Wrapped_data[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        Wrapped_data_rangeprofile       =  Wrapped_data[1:maxind_val_lon,1] 
    end
    plot_profile(Lon_vals, Wrapped_data_rangeprofile[:], "", savepath, "Slantrange_profile_wrapped", "Lon [deg]", "", "", "", "")

    # Critical baseline
    Critical_baseline                   = ref_data[:,:,11]
    if profile_flag == 1
        Critical_baseline_rangeprofile  = mean(Critical_baseline[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        Critical_baseline_rangeprofile  =  Critical_baseline[1:maxind_val_lon,1] 
    end

    plot_image(Lon_vals,Lat_vals,Critical_baseline[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Critical_baseline_P1", "")
    plot_profile(Lon_vals, Critical_baseline_rangeprofile[:], "", savepath, "Critical_baseline_profile", "Lon [deg]", "Critical baseline  [m]", "", "", "")

    p1=(heatmap(Lon_vals,Lat_vals,Critical_baseline[1:maxind_val_lon,1:maxind_val_lat]',xlabel="Lon [deg]",ylabel="Lat [deg]", title="", clim=(0e3,40e3),
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200
    savefig(p1, savepath*"Critical_baseline_P1_2"*".png")


    # Perp baseline
    Perp_baseline                       = ref_data[:,:,7]
    if profile_flag == 1
        Perp_baseline_rangeprofile      = mean(Perp_baseline[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        Perp_baseline_rangeprofile      =  Perp_baseline[1:maxind_val_lon,1]   
    end

    plot_image(Lon_vals,Lat_vals,Perp_baseline[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Perp_baseline_P1", "")
    plot_profile(Lon_vals, Perp_baseline_rangeprofile[:], "", savepath, "Perp_baseline_profile","Lon [deg]", "Perp baseline  [m]", "", "", "")

    # vertical wavenumber
    Vert_wavnum                         = ref_data[:,:,8]
    if profile_flag == 1
        Vert_wavnum_rangeprofile        = mean(Vert_wavnum[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        Vert_wavnum_rangeprofile        =  Vert_wavnum[1:maxind_val_lon,1] 
    end

    plot_image(Lon_vals,Lat_vals,Vert_wavnum[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Vert_wavnum_P1", "")
    plot_profile(Lon_vals, Vert_wavnum_rangeprofile[:], "", savepath, "Vert_wavnum_profile", "Lon [deg]", "Vert wave num  ", "", "", "")

    # Theoretical correlation 
    Correlation_theo                    =  ref_data[:,:,12]
    if profile_flag == 1
        Correlation_theo_rangeprofile   = mean(Correlation_theo[1:maxind_val_lon,:],dims=2)
    elseif profile_flag == 0
        Correlation_theo_rangeprofile   =  Correlation_theo[1:maxind_val_lon,1] 
    end

    plot_image(Lon_vals,Lat_vals,Correlation_theo[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Correlation_theo_P1", "")
    plot_profile(Lon_vals, Correlation_theo_rangeprofile[:], "", savepath, "Correlation_theo_profile", "Lon [deg]", "Correlation theory  ", "", "", "")

    p1=(heatmap(Lon_vals,Lat_vals,Correlation_theo[1:maxind_val_lon,1:maxind_val_lat]',xlabel="Lon [deg]",ylabel="Lat [deg]", title="", clim=(0,1),
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200
    savefig(p1, savepath*"Correlation_theo_P1_2"*".png")

end


function get_max_ind_vec(data_len,looks)
    for i=1:100
        data_len = data_len -1
        if mod(data_len,looks)==0
            return data_len
        end
    end
    return 0
end
#------------------------------------------------------------------------------------------

#Load the output file

s_geom_filepath         = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/50/50_scene_geometry_1.tif" 
t_geom_filepath         = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/50/50_target_geometry_1.tif" 

opdata_filepath_s       = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/50/50_sim_output_main_scene_1.tif" 
opdata_filepath_t       = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/50/50_sim_output_main_target_1.tif" 

savepath                = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/50/"
#=
s_geom_filepath         = "/Users/joshil/Documents/Outputs/InSAR Outputs/Intrepid_Geo_outputs_set1/1/1_scene_geometry_1.tif" 
t_geom_filepath         = "/Users/joshil/Documents/Outputs/InSAR Outputs/Intrepid_Geo_outputs_set1/1/1_target_geometry_1.tif" 

opdata_filepath_s       = "/Users/joshil/Documents/Outputs/InSAR Outputs/Intrepid_Geo_outputs_set1/1/1_sim_output_main_1.tif" 
opdata_filepath_t       = "/Users/joshil/Documents/Outputs/InSAR Outputs/Intrepid_Geo_outputs_set1/1/1_sim_output_main_1.tif" 

savepath                = "/Users/joshil/Documents/Outputs/InSAR Outputs/Intrepid_Geo_outputs_set1/1/"
=#

process_plot_output_flag      = "S"

profile_flag            = 1

Looks_along_Lat         = 4#2#2
Looks_along_Lon         = 12#4 #15

oversampling_factor_looks = 1#2.58*2.18


if process_plot_output_flag == "S"

    savepath                 = savepath*"Plots_s/"
    if ~ispath(savepath)
        mkdir(savepath)
    end    

    s_geom_data                     = GeoArrays.read(s_geom_filepath)
    t_geom_data                     = GeoArrays.read(t_geom_filepath)

    coords_op                       = collect(GeoArrays.coords(s_geom_data))

    Lon_vals_all                    = [x[1] for x in coords_op]
    Lon_vals                        = (Lon_vals_all[:,1])

    Lat_vals_all                    = [x[2] for x in coords_op]
    Lat_vals                        = reverse(Lat_vals_all[1,:])

    maxind_val_lon                  = get_max_ind_vec(length(Lon_vals),Looks_along_Lon)
    maxind_val_lat                  = get_max_ind_vec(length(Lat_vals),Looks_along_Lat)

    Lon_vals                        = Lon_vals[1:maxind_val_lon]
    Lat_vals                        = Lat_vals[1:maxind_val_lat]

    geom_savepath                   = savepath*"Geom_plots_s/"
    plot_geometry_variables(s_geom_data, Lon_vals, Lat_vals, maxind_val_lon, maxind_val_lat, profile_flag, geom_savepath)

    data_op                         = ArchGDAL.read(opdata_filepath_s)

elseif process_plot_output_flag == "T"

    savepath                 = savepath*"Plots_t/"
    if ~ispath(savepath)
        mkdir(savepath)
    end  

    t_geom_data                     = GeoArrays.read(t_geom_filepath)

    coords_op                       = collect(GeoArrays.coords(t_geom_data))

    Lon_vals_all                    = [x[1] for x in coords_op]
    Lon_vals                        = Lon_vals_all[:,1]

    Lat_vals_all                    = [x[2] for x in coords_op]
    Lat_vals                        = reverse(Lat_vals_all[1,:])

    maxind_val_lon                  = get_max_ind_vec(length(Lon_vals),Looks_along_Lon)
    maxind_val_lat                  = get_max_ind_vec(length(Lat_vals),Looks_along_Lat)

    Lon_vals                        = Lon_vals[1:maxind_val_lon]
    Lat_vals                        = Lat_vals[1:maxind_val_lat]

    geom_savepath                   = savepath*"Geom_plots_t/"
    plot_geometry_variables(t_geom_data, Lon_vals, Lat_vals, maxind_val_lon, maxind_val_lat, profile_flag, geom_savepath)

    data_op                         = ArchGDAL.read(opdata_filepath_t)

end


Data_1                          = ArchGDAL.getband(data_op, 1)
Data_2                          = ArchGDAL.getband(data_op, 2)

Data_1                          = Data_1[1:maxind_val_lon,1:maxind_val_lat]
Data_2                          = Data_2[1:maxind_val_lon,1:maxind_val_lat]


# Amplitude, phase and power plots and statistics
Amp_1                           = abs.(Data_1)
Amp_2                           = abs.(Data_2)

Pow_1                           = abs.(Amp_1 .* conj(Amp_1))
Pow_2                           = abs.(Amp_2 .* conj(Amp_2))

Phase_1                         = angle.(Data_1)
Phase_2                         = angle.(Data_2)

gr()

plot_image(Lon_vals,Lat_vals,Pow_1',"log", savepath*"Power_plots/", "Power_P1_log", "")
plot_image(Lon_vals,Lat_vals,Pow_2',"log", savepath*"Power_plots/", "Power_P2_log", "")

plot_image(Lon_vals,Lat_vals,Pow_1',"lin", savepath*"Power_plots/", "Power_P1_lin", "")
plot_image(Lon_vals,Lat_vals,Pow_2',"lin", savepath*"Power_plots/", "Power_P2_lin", "")

plot_image(Lon_vals,Lat_vals,(Pow_1./mean(Pow_1))',"lin", savepath*"Power_plots/", "Power_P1_lin_norm", "")
plot_image(Lon_vals,Lat_vals,(Pow_2./mean(Pow_2))',"lin", savepath*"Power_plots/", "Power_P2_lin_norm", "")

if profile_flag == 1
    Pow_1_rangeprofile          = mean(Pow_1[1:maxind_val_lon,:],dims=2)
    Pow_2_rangeprofile          = mean(Pow_2[1:maxind_val_lon,:],dims=2)
elseif profile_flag == 0
    Pow_1_rangeprofile          =  Pow_1[:,1] 
    Pow_2_rangeprofile          =  Pow_2[:,1] 
end 

plot_profile(Lon_vals, 10 .* log10.(Pow_1_rangeprofile[:]), "", savepath*"Power_plots/", "Pow_1_profile", "Lon [deg]", "Power [dB] ", "", "", "")
plot_profile(Lon_vals, 10 .* log10.(Pow_2_rangeprofile[:]), "", savepath*"Power_plots/", "Pow_2_profile", "Lon [deg]", "Power [dB]", "", "", "")


plot_image(Lon_vals,Lat_vals,Amp_1',"log", savepath*"Amplitude_plots/", "Amp_P1_log", "")
plot_image(Lon_vals,Lat_vals,Amp_2',"log", savepath*"Amplitude_plots/", "Amp_P2_log", "")

plot_image(Lon_vals,Lat_vals,Amp_1',"lin", savepath*"Amplitude_plots/", "Amp_P1_lin", "")
plot_image(Lon_vals,Lat_vals,Amp_2',"lin", savepath*"Amplitude_plots/", "Amp_P2_lin", "")

plot_image(Lon_vals,Lat_vals,(Amp_1./mean(Amp_1))',"lin", savepath*"Amplitude_plots/", "Amp_P1_lin_norm", "")
plot_image(Lon_vals,Lat_vals,(Amp_2./mean(Amp_2))',"lin", savepath*"Amplitude_plots/", "Amp_P2_lin_norm", "")

if profile_flag == 1
    Amp_1_rangeprofile          = mean(Amp_1[1:maxind_val_lon,:],dims=2)
    Amp_2_rangeprofile          = mean(Amp_2[1:maxind_val_lon,:],dims=2)
elseif profile_flag == 0
    Amp_1_rangeprofile          =  Amp_1[:,1]
    Amp_2_rangeprofile          =  Amp_2[:,1]
end

plot_profile(Lon_vals, 10 .* log10.(Amp_1_rangeprofile[:]), "", savepath*"Amplitude_plots/", "Amp_1_profile","Lon [deg]", "Amplitude [dB] ", "", "", "")
plot_profile(Lon_vals, 10 .* log10.(Amp_2_rangeprofile[:]), "", savepath*"Amplitude_plots/", "Amp_2_profile", "Lon [deg]", "Amplitude [dB]", "", "", "")

plot_image(Lon_vals,Lat_vals,Phase_1',"phase", savepath*"Phase_plots/", "Phase_P1_lin", "")
plot_image(Lon_vals,Lat_vals,Phase_2',"phase", savepath*"Phase_plots/", "Phase_P2_lin", "")

# Amplitude statistics plots 
plot_var                        = Amp_1[:]
A = plot_histogram_statistics_amplitude(plot_var[:], savepath*"Amplitude_plots/", "P1_amp")
plot_var                        = Amp_2[:]
A = plot_histogram_statistics_amplitude(plot_var[:], savepath*"Amplitude_plots/", "P2_amp")

# Power statistics plots 
plot_var                        = Pow_1[:]
A = plot_histogram_statistics_power(plot_var[:], savepath*"Power_plots/", "P1_pow")
plot_var                        = Pow_2[:]
A = plot_histogram_statistics_power(plot_var[:], savepath*"Power_plots/", "P2_pow")

# Phase statistics plots 
plot_var                        = Phase_1[:]
A = plot_histogram_statistics_phase(plot_var[:], savepath*"Phase_plots/", "P1_phase")
plot_var                        = Phase_2[:]
A = plot_histogram_statistics_phase(plot_var[:], savepath*"Phase_plots/", "P2_phase")


# Multi looking power stat and plot
Lat_length_ml                   = Int(maxind_val_lat/Looks_along_Lat)
Lon_length_ml                   = Int(maxind_val_lon/Looks_along_Lon)

Pow_1_multilooked               = zeros(Lon_length_ml,Lat_length_ml)
Pow_2_multilooked               = zeros(Lon_length_ml,Lat_length_ml)

Lat_vals_multilooked            = zeros(Lat_length_ml)
Lon_vals_multilooked            = zeros(Lon_length_ml)

norm_mean_pow1                  = maximum(Pow_1)
norm_mean_pow2                  = maximum(Pow_1)

k=1
for i=1:Lat_length_ml
    l=1
    for j=1:Lon_length_ml
        Pow_1_multilooked[j,i]  = mean(Pow_1[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1]) ./norm_mean_pow1
        Pow_2_multilooked[j,i]  = mean(Pow_2[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1]) ./norm_mean_pow2
        Lat_vals_multilooked[i] = mean(Lat_vals[k:k+Looks_along_Lat-1])
        Lon_vals_multilooked[j] = mean(Lon_vals[l:l+Looks_along_Lon-1])
        l=l+Looks_along_Lon
    end
    k=k+Looks_along_Lat
end

plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Pow_1_multilooked',"log", savepath*"Power_plots/", "ML_Power_P1_log", "")
plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Pow_2_multilooked',"log", savepath*"Power_plots/", "ML_Power_P2_log", "")
        
plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Pow_1_multilooked',"lin", savepath*"Power_plots/", "ML_Power_P1_lin", "")
plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Pow_2_multilooked',"lin", savepath*"Power_plots/", "ML_Power_P2_lin", "")

# Power statistics plots 
plot_var                        = Pow_1_multilooked[:]
A = plot_histogram_only(plot_var[:],"Power", savepath*"Power_plots/", "ML_P1_pow")
plot_var                        = Pow_2_multilooked[:]
A = plot_histogram_only(plot_var[:],"Power", savepath*"Power_plots/", "ML_P2_pow")



#Interferometry
Data_12 		                = Data_1 .* conj(Data_2)

Pow_12                          = abs.(Data_12)
Phase_12                        = angle.(Data_12)

plot_image(Lon_vals,Lat_vals,Pow_12',"log", savepath*"Int_plots/", "Int_Power_P12_log", "")
plot_image(Lon_vals,Lat_vals,Pow_12',"lin", savepath*"Int_plots/", "Int_Power_P12_lin", "")

plot_image(Lon_vals,Lat_vals,(Pow_12 ./ ((Pow_1.*Pow_2).^0.5))',"log", savepath*"Int_plots/", "Int_Coh_Power_P12_log", "")
plot_image(Lon_vals,Lat_vals,(Pow_12 ./ ((Pow_1.*Pow_2).^0.5))',"lin", savepath*"Int_plots/", "Int_Coh_Power_P12_lin", "")

plot_image(Lon_vals,Lat_vals,Phase_12',"lin", savepath*"Int_plots/", "Int_Phase_P12_lin", "")

if profile_flag == 1
    Pow_12_rangeprofile         =  mean(Pow_12,dims=2)
elseif profile_flag == 0
    Pow_12_rangeprofile         =  Pow_12[:,1]
end
plot_profile(Lon_vals, 10 .* log10.(Pow_12_rangeprofile[:]), "", savepath*"Int_plots/", "Pow_12_profile", "Lon [deg]", "Power [dB] ", "", "", "")


Int_Pow_multilooked             = zeros(Lon_length_ml,Lat_length_ml)
Int_Phase_multilooked           = zeros(Lon_length_ml,Lat_length_ml)
Correlation                     = zeros(Lon_length_ml,Lat_length_ml)

Lat_vals_multilooked            = zeros(Lat_length_ml)
Lon_vals_multilooked            = zeros(Lon_length_ml)

k=1
for i=1:Int((length(Lat_vals))/Looks_along_Lat)
    l=1
    for j=1:Int((length(Lon_vals))/Looks_along_Lon)
        complex_coherence = mean(Data_12[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1])./ ((mean(Pow_1[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1]) .* mean(Pow_2[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1])).^0.5)
        Int_Pow_multilooked[j,i] = abs.(complex_coherence)
        Int_Phase_multilooked[j,i] = angle.(complex_coherence) 

        Lat_vals_multilooked[i] = mean(Lat_vals[k:k+Looks_along_Lat-1])
        Lon_vals_multilooked[j] = mean(Lon_vals[l:l+Looks_along_Lon-1])
        l=l+Looks_along_Lon
    end
    k=k+Looks_along_Lat
end

replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v) 

Int_Pow_multilooked = replace_nan(Int_Pow_multilooked)
Int_Phase_multilooked= replace_nan(Int_Phase_multilooked)



plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Int_Pow_multilooked',"log", savepath*"Int_plots/", "ML_Power_P12_log", "")
plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Int_Pow_multilooked',"lin", savepath*"Int_plots/", "ML_Power_P12_lin", "")
        
plot_image2(Lon_vals_multilooked,Lat_vals_multilooked,Int_Phase_multilooked',"phase", savepath*"Int_plots/", "ML_Phase_P12_lin", "",":twilight")
#plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Correlation,"lin", savepath*"Int_plots/", "ML_Correlation_lin", "")

if profile_flag == 1
    Correlation_rangeprofile    =  mean(Int_Pow_multilooked,dims=2)
    Int_phase_rangeprofile      =  mean(Int_Phase_multilooked,dims=2)
elseif profile_flag == 0
    Correlation_rangeprofile    = Int_Pow_multilooked[:,1] 
    Int_phase_rangeprofile      = Int_Phase_multilooked[:,1] 
end

plot_profile(Lon_vals_multilooked, Correlation_rangeprofile[:], "", savepath*"Int_plots/", "Int_mag_profile", "Lon [deg]", "Magnitude  ", "", "", "")
plot_profile(Lon_vals_multilooked, Int_phase_rangeprofile[:], "", savepath*"Int_plots/", "Int_phase__profile", "Lon [deg]", "Phase [rad]  ", "", "", "")

A = plot_histogram_only(Int_Pow_multilooked[:], "Magnitude ",savepath*"Int_plots/", "ML_Int_mag")
A = plot_histogram_only(Int_Phase_multilooked[:], "Phase",savepath*"Int_plots/", "ML_Int_ph")

if process_plot_output_flag == "S"
#gamma_ip = mean(Int_Pow_multilooked[:])
gamma_ip = mean(s_geom_data[:,:,12])
elseif process_plot_output_flag == "T"
    gamma_ip = mean(t_geom_data[:,:,12])
end


A = plot_interferogram_histogram_statistics_phase(Int_Phase_multilooked[:], gamma_ip, (Looks_along_Lat*Looks_along_Lon)/oversampling_factor_looks,0.0, (-0.5,0.5), savepath*"Int_plots/", "PDF_Phase_Han1")
A = plot_interferogram_histogram_statistics_phase_han1(Int_Phase_multilooked[:], gamma_ip, (Looks_along_Lat*Looks_along_Lon)/oversampling_factor_looks, 0.0, (-0.5,0.5), savepath*"Int_plots/", "PDF_Phase_Han2")
B = plot_interferogram_histogram_statistics_magnitude(Int_Pow_multilooked[:], gamma_ip, (Looks_along_Lat*Looks_along_Lon)/oversampling_factor_looks, (0.5,1), savepath*"Int_plots/", "PDF_Mag_Han1")

Correlation_rangeprofile_window = movmean(Correlation_rangeprofile, 20)
Int_phase_rangeprofile_window   = movmean(Int_phase_rangeprofile, 20)

plot_profile(Lon_vals_multilooked, Correlation_rangeprofile_window[:], "", savepath*"Int_plots/", "Int_mag_profile_window", "Lon [deg]", "Magnitude  ", "", "", "")
plot_profile(Lon_vals_multilooked, Int_phase_rangeprofile_window[:], "", savepath*"Int_plots/", "Int_phase_profile_window", "Lon [deg]", "Phase [rad]  ", "", "", "")


CF_A, CF_B = linear_fit(Lon_vals_multilooked,Correlation_rangeprofile_window)

p1=(plot(Lon_vals_multilooked, Correlation_rangeprofile_window[:],xlabel="Lon [deg]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:black, label="",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
p1=(plot!(Lon_vals_multilooked, (Lon_vals_multilooked .* CF_B).+CF_A,xlabel="Lon [deg]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:red, label="Linear Fit",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p1, savepath*"Int_plots/"*"Int_mag_profile_window2"*".png")


if process_plot_output_flag == "S"
#Estimate height from wrapped interferometric phase

if profile_flag == 1
    int_wrapped_data        =  mean(Int_Phase_multilooked,dims=2)
elseif profile_flag == 0
    int_wrapped_data        = Int_Phase_multilooked[:,1]
end
plot_profile(Lon_vals_multilooked, int_wrapped_data[:], "", savepath*"Int_plots/", "Int_wrapped_data_profile", "Lon [deg]", "  ", "", "", "")

int_unwrapped_data          = DSP.unwrap(int_wrapped_data[:]) .* -1
plot_profile(Lon_vals_multilooked, int_unwrapped_data[:], "", savepath*"Int_plots/", "Int_unwrapped_data_profile", "Lon [deg]", "  ", "", "", "")

if profile_flag == 1
    slant_range_profile         =  mean(s_geom_data[1:Looks_along_Lon:maxind_val_lon,:,1],dims=2)
    look_angle_profile          =  mean(s_geom_data[1:Looks_along_Lon:maxind_val_lon,:,3],dims=2)
    perp_baseline_profile       =  mean(s_geom_data[1:Looks_along_Lon:maxind_val_lon,:,7],dims=2)
elseif profile_flag == 0
    slant_range_profile         = s_geom_data[1:Looks_along_Lon:maxind_val_lon,1,1]
    look_angle_profile          = s_geom_data[1:Looks_along_Lon:maxind_val_lon,1,3]
    perp_baseline_profile       = s_geom_data[1:Looks_along_Lon:maxind_val_lon,1,7]
end

int_unwrapped_height_2      = int_unwrapped_data .* (0.23793052222222222/(4*pi)) .* (slant_range_profile .*sind.(look_angle_profile) ./perp_baseline_profile)
plot_profile(Lon_vals_multilooked, int_unwrapped_height_2[:], "", savepath*"Int_plots/", "Int_unwrapped_height_profile", "Lon [deg]", "Height [m]", "", "", "")


end
#=
p2 = (plot(Lon_vals, Correlation_theo_P1_rangeprofile[:], label="Theory"))
p2 = (plot!(Lon_vals_multilooked, (Lon_vals_multilooked .* CF_B).+CF_A, label="Sim"))

AASR = 10 ^(-26/10)
p2 = (plot!(Lon_vals_multilooked, ((Lon_vals_multilooked .* CF_B).+CF_A) * (1/1+(AASR)), label="Sim - AASR -26 dB", legend=:topleft, xlabel="Lon [deg]", ylabel="Correlation",
topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p2, savepath*"Int_plots/"*"Correlation_comparison_1"*".png")

AASR = 10 ^(-21.9/10)
p2 = (plot!(Lon_vals_multilooked, ((Lon_vals_multilooked .* CF_B).+CF_A) * (1/1+(AASR)), label="Sim - AASR -22 dB", legend=:topleft, xlabel="Lon [deg]", ylabel="Correlation",
topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p2, savepath*"Int_plots/"*"Correlation_comparison_2"*".png")
=#


#=
using JLD2
using FFTW
@load "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/56/Simulation_info_56_1.jld"

F_Srx = fftshift(fft(Srx))
freq_Axis = LinRange(-1/(2*params.Δt),1/(2*params.Δt),length(Srx))


display(plot(freq_Axis./1e6,abs.(F_Srx), ylabel="FFT Magnitude", xlabel="Frequency [MHz]", label=:false))
display(plot!(tickfont=font(12), ytickfont=font(12), legendfontsize=12, guidefont=font(12), titlefontsize=12, legend=:bottomright, leftmargin=2mm, bottommargin=2mm))




F_rawdata = fftshift(fft(rawdata[1,1,:]))
freq_Axis = LinRange(-1/(2*params.Δt),1/(2*params.Δt),length(rawdata[1,1,:]))


display(plot(freq_Axis./1e6,abs.(F_rawdata), ylabel="FFT Magnitude", xlabel="Frequency [MHz]", label=:false, xlim=(-50,50)))
display(plot!(tickfont=font(12), ytickfont=font(12), legendfontsize=12, guidefont=font(12), titlefontsize=12, legend=:bottomright, leftmargin=2mm, bottommargin=2mm))
=#

targets_ref_data                    = t_geom_data[:,:,13]
targets_ref_corr_rangeprofile   = targets_ref_data[1:maxind_val_lon,1]

Amp_1_rangeprofile          =  Amp_1[:,1]

plot_profile(Lon_vals, (Amp_1_rangeprofile[:]), "", savepath, "Amp_profile1_comp", "Lon [deg]", "Amplitude", "", "", "")

plot_profile(Lon_vals, (targets_ref_corr_rangeprofile[:]), "", savepath, "RCS_comp", "Lon [deg]", "RCS", "", "", "")




if process_plot_output_flag == "S"
    #Estimate height from wrapped interferometric phase

    Int_unwrapped_data_new          = DSP.unwrap(Int_Phase_multilooked,dims=1) .* -1
    plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Int_unwrapped_data_new',"lin", savepath*"Int_plots/", "Test_1", "")


        slant_range_profile         = s_geom_data[1:Looks_along_Lon:maxind_val_lon,1:Looks_along_Lat:maxind_val_lat,1]
        look_angle_profile          = s_geom_data[1:Looks_along_Lon:maxind_val_lon,1:Looks_along_Lat:maxind_val_lat,3]
        perp_baseline_profile       = s_geom_data[1:Looks_along_Lon:maxind_val_lon,1:Looks_along_Lat:maxind_val_lat,7]

    int_unwrapped_height_2      = Int_unwrapped_data_new .* (0.23793052222222222/(4*pi)) .* (slant_range_profile[:,:,1] .*sind.(look_angle_profile[:,:,1]) ./perp_baseline_profile[:,:,1])
    plot_image(Lon_vals_multilooked,Lat_vals_multilooked,int_unwrapped_height_2',"lin", savepath*"Int_plots/", "Test_3", "")
    
    plot_image(Lon_vals_multilooked,Lat_vals_multilooked,int_unwrapped_height_2'.+2000,"lin", savepath*"Int_plots/", "Test_4", "")


    ref_DEM = s_geom_data[1:Looks_along_Lon:maxind_val_lon,1:Looks_along_Lat:maxind_val_lat,14]

    plot_image(Lon_vals_multilooked,Lat_vals_multilooked,ref_DEM[:,:,1]' .- (int_unwrapped_height_2'.+2000),"lin", savepath*"Int_plots/", "Test_5", "")

end

#using JLD2
#@save "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/50/test_data.jld" Int_Phase_multilooked  Int_Pow_multilooked

#
#using PyCall
#snaphu = pyimport("snaphu")
