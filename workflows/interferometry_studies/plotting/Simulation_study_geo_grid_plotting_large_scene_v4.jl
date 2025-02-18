
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

function plot_geometry_variables(ref_data, Lon_vals, Lat_vals, maxind_val_lon, maxind_val_lat, savepath)

    #DEM
    DEM_data                            = ref_data[:,:,14]
    DEM_region_rangeprofile             = DEM_data[1:maxind_val_lon,1]
    plot_image(Lon_vals,Lat_vals,DEM_data[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "DEM", "")
    plot_profile(Lon_vals, DEM_region_rangeprofile[:], "", savepath, "DEM_profile", "Lon [deg]", "DEM height [m]","", "", "", )

    #BRCS
    targets_ref_data                    = ref_data[:,:,13]
    targets_ref_corr_rangeprofile       = targets_ref_data[1:maxind_val_lon,1]
    plot_image(Lon_vals,Lat_vals,targets_ref_data[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "RCS", "")
    plot_profile(Lon_vals, targets_ref_corr_rangeprofile[:], "", savepath, "RCS_profile", "Lon [deg]", "Reflectivity from RCS","", "", "", )

    #Look angle
    Lookangle_P1                        = ref_data[:,:,3]
    Lookangle_P2                        = ref_data[:,:,4]
    Lookangle_Diff                      = Lookangle_P1 - Lookangle_P2

    Lookangle_P1_rangeprofile           =  Lookangle_P1[1:maxind_val_lon,1] 
    Lookangle_P2_rangeprofile           =  Lookangle_P2[1:maxind_val_lon,1]
    Lookangle_Diff_rangeprofile         =  Lookangle_Diff[1:maxind_val_lon,1]

    plot_image(Lon_vals,Lat_vals,Lookangle_P1[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Lookangle_P1", "")
    plot_image(Lon_vals,Lat_vals,Lookangle_P2[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Lookangle_P2", "")
    plot_image(Lon_vals,Lat_vals,Lookangle_Diff[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Lookangle_Diff", "")

    plot_profile(Lon_vals, Lookangle_P1_rangeprofile[:], Lookangle_P2_rangeprofile[:], savepath, "Lookangle_profile1", "Lon [deg]", "Look Angle [deg]", "", "Platform 1", "Platform 2")
    plot_profile(Lon_vals, Lookangle_Diff_rangeprofile[:], "", savepath, "Lookangle_profile", "Lon [deg]", "Look Angle Diff [deg]", "", "Difference", "")

    #Incidence angle
    Incangle_P1                         = ref_data[:,:,5]
    Incangle_P2                         = ref_data[:,:,6]
    Incangle_Diff                       = Incangle_P1 - Incangle_P2

    Incangle_P1_rangeprofile            =  Incangle_P1[1:maxind_val_lon,1]
    Incangle_P2_rangeprofile            =  Incangle_P2[1:maxind_val_lon,1] 
    Incangle_Diff_rangeprofile          =  Incangle_Diff[1:maxind_val_lon,1] 

    plot_image(Lon_vals,Lat_vals,Incangle_P1[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Incangle_P1", "")
    plot_image(Lon_vals,Lat_vals,Incangle_P2[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Incangle_P2", "")
    plot_image(Lon_vals,Lat_vals,Incangle_Diff[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Incangle_Diff", "")

    plot_profile(Lon_vals, Incangle_P1_rangeprofile[:], Incangle_P2_rangeprofile[:], savepath, "Incangle_profile1", "Lon [deg]", "Incidence Angle [deg]", "", "Platform 1", "Platform 2")
    plot_profile(Lon_vals, Incangle_Diff_rangeprofile[:], "", savepath, "Incangle_profile", "Lon [deg]", "Incidence Angle Diff [deg]", "", "Difference", "")

    #Local Incidence angle
    LocalIncangle                       = ref_data[:,:,9]
    LocalIncangle_rangeprofile          =  LocalIncangle[1:maxind_val_lon,1] 
    plot_image(Lon_vals,Lat_vals,LocalIncangle[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "LocalIncangle", "")
    plot_profile(Lon_vals, LocalIncangle_rangeprofile[:], "", savepath, "LocalIncangle_profile", "Lon [deg]", "Local incidence Angle [deg]","", "", "", )

    #Range slope angle
    range_slope                         = ref_data[:,:,10]
    range_slope_rangeprofile            =  range_slope[1:maxind_val_lon,1] 
    plot_image(Lon_vals,Lat_vals,range_slope[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "RangeSlopeangle", "")
    plot_profile(Lon_vals, range_slope_rangeprofile[:], "", savepath, "RangeSlopeangle_profile", "Lon [deg]", "Range Slope Angle [deg]","", "", "", )

    #Slant range
    Slantrange_P1                       = ref_data[:,:,1]
    Slantrange_P2                       = ref_data[:,:,2]
    Slantrange_Diff                     = Slantrange_P1 - Slantrange_P2

    Slantrange_P1_rangeprofile          =  Slantrange_P1[1:maxind_val_lon,1] 
    Slantrange_P2_rangeprofile          =  Slantrange_P2[1:maxind_val_lon,1] 
    Slantrange_Diff_rangeprofile        =  Slantrange_Diff[1:maxind_val_lon,1] 

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

    Wrapped_data_rangeprofile           =  Wrapped_data[1:maxind_val_lon,1] 
    plot_profile(Lon_vals, Wrapped_data_rangeprofile[:], "", savepath, "Slantrange_profile_wrapped", "Lon [deg]", "", "", "", "")

    # Critical baseline
    Critical_baseline                   = ref_data[:,:,11]
    Critical_baseline_rangeprofile      =  Critical_baseline[1:maxind_val_lon,1] 

    plot_image(Lon_vals,Lat_vals,Critical_baseline[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Critical_baseline_P1", "")
    plot_profile(Lon_vals, Critical_baseline_rangeprofile[:], "", savepath, "Critical_baseline_profile", "Lon [deg]", "Critical baseline  [m]", "", "", "")

    # Perp baseline
    Perp_baseline                       = ref_data[:,:,7]
    Perp_baseline_rangeprofile          =  Perp_baseline[1:maxind_val_lon,1]  

    plot_image(Lon_vals,Lat_vals,Perp_baseline[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Perp_baseline_P1", "")
    plot_profile(Lon_vals, Perp_baseline_rangeprofile[:], "", savepath, "Perp_baseline_profile","Lon [deg]", "Perp baseline  [m]", "", "", "")

    # vertical wavenumber
    Vert_wavnum                         = ref_data[:,:,8]
    Vert_wavnum_rangeprofile            =  Vert_wavnum[1:maxind_val_lon,1]  

    plot_image(Lon_vals,Lat_vals,Vert_wavnum[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Vert_wavnum_P1", "")
    plot_profile(Lon_vals, Vert_wavnum_rangeprofile[:], "", savepath, "Vert_wavnum_profile", "Lon [deg]", "Vert wave num  ", "", "", "")

    # Theoretical correlation 
    Correlation_theo                    =  ref_data[:,:,12]
    Correlation_theo_rangeprofile       =  Correlation_theo[1:maxind_val_lon,1] 

    plot_image(Lon_vals,Lat_vals,Correlation_theo[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Correlation_theo_P1", "")
    plot_profile(Lon_vals, Correlation_theo_rangeprofile[:], "", savepath, "Correlation_theo_profile", "Lon [deg]", "Correlation theory  ", "", "", "")

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


function combine_tif_files(file_list_path, file_list_name)

    file_list_path = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/11/s_all/"
    file_list_name = "s_all_list.txt"
    
    list =  readdlm(file_list_path*file_list_name)

    display(heatmap())
    for i=1:size(list)[1]

        
        file_t1 = file_list_path*list[i]

        s_geom_data             = GeoArrays.read(file_t1)
        
        coords_op               = collect(GeoArrays.coords(s_geom_data))
        
        Lon_vals_all            = [x[1] for x in coords_op]
        Lon_vals                = Lon_vals_all[:,1]
        
        Lat_vals_all            = [x[2] for x in coords_op]
        Lat_vals                = Lat_vals_all[1,:]

        p1=(heatmap!(Lon_vals,Lat_vals,s_geom_data[:,:,3],xlabel="Lon [deg]",ylabel="Lat [deg]", title="",
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200
    end
    savefig(p1, file_list_path*"LA.png")



end


#------------------------------------------------------------------------------------------

#Load the output file


s_file_list_path = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/11/s_all/"
s_file_list_name = "s_all_list.txt"

s_list =  readdlm(s_file_list_path*s_file_list_name)



file_list_path = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/11/sim_all/"
file_list_name = "sim_all_list.txt"

list =  readdlm(file_list_path*file_list_name)

display(heatmap())
for i=1:size(list)[1]


    s_geom_data             = GeoArrays.read(s_file_list_path*s_list[i])

    coords_op               = collect(GeoArrays.coords(s_geom_data))
    
    Lon_vals_all            = [x[1] for x in coords_op]
    Lon_vals                = Lon_vals_all[:,1]
    
    Lat_vals_all            = [x[2] for x in coords_op]
    Lat_vals                = (Lat_vals_all[1,:])
    


    println(i)
file_t1 = file_list_path*list[i]

Looks_along_Lat         = 4#2
Looks_along_Lon         = 12 #15

maxind_val_lon          = get_max_ind_vec(length(Lon_vals),Looks_along_Lon)
maxind_val_lat          = get_max_ind_vec(length(Lat_vals),Looks_along_Lat)

Lon_vals                = Lon_vals[1:maxind_val_lon]
Lat_vals                = Lat_vals[1:maxind_val_lat]


data_op                 = ArchGDAL.read(file_t1)

Data_1                  = ArchGDAL.getband(data_op, 3)
Data_2                  = ArchGDAL.getband(data_op, 4)

Data_1                  = Data_1[1:maxind_val_lon,1:maxind_val_lat]
Data_2                  = Data_2[1:maxind_val_lon,1:maxind_val_lat]

Amp_1                   = abs.(Data_1)
Amp_2                   = abs.(Data_2)


Pow_1                   = abs.(Amp_1 .* conj(Amp_1))
Pow_2                   = abs.(Amp_2 .* conj(Amp_2))

#Interferometry
Data_12 		        = Data_1 .* conj(Data_2)

Pow_12                  = abs.(Data_12)
Phase_12                = angle.(Data_12)


Lat_length_ml           = Int(maxind_val_lat/Looks_along_Lat)
Lon_length_ml           = Int(maxind_val_lon/Looks_along_Lon)


Pow_1_multilooked       = zeros(Lon_length_ml,Lat_length_ml)
Pow_2_multilooked       = zeros(Lon_length_ml,Lat_length_ml)

Lat_vals_multilooked    = zeros(Lat_length_ml)
Lon_vals_multilooked    = zeros(Lon_length_ml)

norm_mean_pow1 = maximum(Pow_1)
norm_mean_pow2 = maximum(Pow_1)

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

Int_Pow_multilooked         = zeros(Lon_length_ml,Lat_length_ml)
Int_Phase_multilooked       = zeros(Lon_length_ml,Lat_length_ml)
Correlation                 = zeros(Lon_length_ml,Lat_length_ml)

Lat_vals_multilooked        = zeros(Lat_length_ml)
Lon_vals_multilooked        = zeros(Lon_length_ml)

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
        

global p1=(heatmap!(Lon_vals_multilooked,(Lat_vals_multilooked),Int_Phase_multilooked',xlabel="Lon [deg]",ylabel="Lat [deg]", title="",
topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200

global p2=(plot!(Lon_vals_multilooked,Int_Phase_multilooked[:,1][:],xlabel="Lon [deg]",ylabel="", title="", ylim=(-0.25,0.25), legend=:false,
topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,500) )) #1200


end
savefig(p1, file_list_path*"Ph1.png")
savefig(p2, file_list_path*"Ph2.png")



