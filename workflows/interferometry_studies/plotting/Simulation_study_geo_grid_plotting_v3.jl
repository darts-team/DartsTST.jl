
using JLD2
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

function plot_image(x_axis,y_axis,Data,unit_flag, savepath, savename, figure_title)

    if ~ispath(savepath)
        mkdir(savepath)
    end    

    if unit_flag == "log"
        p1=(heatmap(x_axis,y_axis,10 .* log10.(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]",title=figure_title,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #500,360
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "lin"
        p1=(heatmap(x_axis,y_axis,(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #1200
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "phase"
        p1=(heatmap(x_axis,y_axis,(Data),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #1200
        savefig(p1, savepath*savename*".png")
    end
end

function plot_image2(x_axis,y_axis,Data,unit_flag, savepath, savename, figure_title, colorbar)

    if ~ispath(savepath)
        mkdir(savepath)
    end    

    if unit_flag == "log"
        p1=(heatmap(x_axis,y_axis,10 .* log10.(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]",title=figure_title,c=colorbar,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #500,360
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "lin"
        p1=(heatmap(x_axis,y_axis,(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title,c=colorbar,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #1200
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "phase"
        p1=(heatmap(x_axis,y_axis,(Data),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title, #c=:twilight,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #1200
        savefig(p1, savepath*savename*".png")
    end
end

function plot_profile(x_axis, Data1, Data2, savepath, savename, xlabel, ylabel, figure_title, label1, label2)

    if Data2==""
        p1=(plot(x_axis, Data1,xlabel=xlabel,ylabel=ylabel,title=figure_title,legend=:topleft, lc=:black, label=label1,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #500,360
        savefig(p1, savepath*savename*".png")
    else
        p1=(plot(x_axis, Data1,xlabel=xlabel,ylabel=ylabel,title=figure_title, legend=:topleft, lc=:black, label=label1,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #500,360
        p1=(plot!(x_axis, Data2,xlabel=xlabel,ylabel=ylabel,title=figure_title, legend=:topleft,lc=:blue,label=label2,
        topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) )) #500,360
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


#------------------------------------------------------------------------------------------

#Load the output file
@load "/Users/joshil/Documents/Outputs/InSAR Outputs/3t_distributed_DEM/test_la2/Output81_10km_30la_geo994.jld" 

# Mention the directory to save plots
savepath                 = "/Users/joshil/Documents/Outputs/InSAR Outputs/3t_distributed_DEM/test_la2/Plots1/"

if ~ispath(savepath)
    mkdir(savepath)
end

s_loc_S             = params.s_loc_1
s_loc_C             = params.s_loc_2
s_loc_H             = params.s_loc_3

Looks_along_S = 3 #2
Looks_along_C = 5 #15

#Data_1              = trg_stat_var_all_1p
#Data_2              = trg_stat_var_all_2p

Data_1              = stat_var_all_1p
Data_2              = stat_var_all_2p

#s_loc_C = s_loc_C[30:end-31]
#Data_1 =  Data_1[:,30:end-31]
#Data_2 =  Data_2[:,30:end-31]

Amp_1               = abs.(Data_1)
Amp_2               = abs.(Data_2)

Pow_1               = abs.(Amp_1 .* conj(Amp_1))
Pow_2               = abs.(Amp_2 .* conj(Amp_2))

Phase_1             = angle.(Data_1)
Phase_2             = angle.(Data_2)

gr()

plot_image(s_loc_C,s_loc_S,Pow_1,"log", savepath*"Power_plots/", "Power_P1_log", "")
plot_image(s_loc_C,s_loc_S,Pow_2,"log", savepath*"Power_plots/", "Power_P2_log", "")

plot_image(s_loc_C,s_loc_S,Pow_1,"lin", savepath*"Power_plots/", "Power_P1_lin", "")
plot_image(s_loc_C,s_loc_S,Pow_2,"lin", savepath*"Power_plots/", "Power_P2_lin", "")

plot_image(s_loc_C,s_loc_S,Pow_1./mean(Pow_1),"lin", savepath*"Power_plots/", "Power_P1_lin_norm", "")
plot_image(s_loc_C,s_loc_S,Pow_2./mean(Pow_2),"lin", savepath*"Power_plots/", "Power_P2_lin_norm", "")

Pow_1_rangeprofile   =  Pow_1[1,:] #mean(Pow_1,dims=1)
Pow_2_rangeprofile   =  Pow_2[1,:] #mean(Pow_2,dims=1)

plot_profile(s_loc_C, 10 .* log10.(Pow_1_rangeprofile[:]), "", savepath*"Power_plots/", "Pow_1_profile", "Lon [deg]", "Power [dB] ", "", "", "")
plot_profile(s_loc_C, 10 .* log10.(Pow_2_rangeprofile[:]), "", savepath*"Power_plots/", "Pow_2_profile", "Lon [deg]", "Power [dB]", "", "", "")


plot_image(s_loc_C,s_loc_S,Amp_1,"log", savepath*"Amplitude_plots/", "Amp_P1_log", "")
plot_image(s_loc_C,s_loc_S,Amp_2,"log", savepath*"Amplitude_plots/", "Amp_P2_log", "")

plot_image(s_loc_C,s_loc_S,Amp_1,"lin", savepath*"Amplitude_plots/", "Amp_P1_lin", "")
plot_image(s_loc_C,s_loc_S,Amp_2,"lin", savepath*"Amplitude_plots/", "Amp_P2_lin", "")

plot_image(s_loc_C,s_loc_S,Amp_1./mean(Amp_1),"lin", savepath*"Amplitude_plots/", "Amp_P1_lin_norm", "")
plot_image(s_loc_C,s_loc_S,Amp_2./mean(Amp_2),"lin", savepath*"Amplitude_plots/", "Amp_P2_lin_norm", "")

Amp_1_rangeprofile   =  Amp_1[1,:] #mean(Amp_1,dims=1)
Amp_2_rangeprofile   =  Amp_2[1,:] #mean(Amp_2,dims=1)

plot_profile(s_loc_C, 10 .* log10.(Amp_1_rangeprofile[:]), "", savepath*"Amplitude_plots/", "Amp_1_profile","Lon [deg]", "Amplitude [dB] ", "", "", "")
plot_profile(s_loc_C, 10 .* log10.(Amp_2_rangeprofile[:]), "", savepath*"Amplitude_plots/", "Amp_2_profile", "Lon [deg]", "Amplitude [dB]", "", "", "")



plot_image(s_loc_C,s_loc_S,Phase_1,"phase", savepath*"Phase_plots/", "Phase_P1_lin", "")
plot_image(s_loc_C,s_loc_S,Phase_2,"phase", savepath*"Phase_plots/", "Phase_P2_lin", "")

# Amplitude statistics plots 
plot_var                = Amp_1[:]
A = plot_histogram_statistics_amplitude(plot_var[:], savepath*"Amplitude_plots/", "P1_amp")
plot_var                = Amp_2[:]
A = plot_histogram_statistics_amplitude(plot_var[:], savepath*"Amplitude_plots/", "P2_amp")

# Power statistics plots 
plot_var                = Pow_1[:]
A = plot_histogram_statistics_power(plot_var[:], savepath*"Power_plots/", "P1_pow")
plot_var                = Pow_2[:]
A = plot_histogram_statistics_power(plot_var[:], savepath*"Power_plots/", "P2_pow")


# Phase statistics plots 
plot_var                = Phase_1[:]
A = plot_histogram_statistics_phase(plot_var[:], savepath*"Phase_plots/", "P1_phase")
plot_var                = Phase_2[:]
A = plot_histogram_statistics_phase(plot_var[:], savepath*"Phase_plots/", "P2_phase")



# Multi looking
Block_size_S = Looks_along_S
Block_size_C = Looks_along_C


Pow_1_multilooked = zeros(Int(length(s_loc_S)/Block_size_S),Int(length(s_loc_C)/Block_size_C))
Pow_2_multilooked = zeros(Int(length(s_loc_S)/Block_size_S),Int(length(s_loc_C)/Block_size_C))

s_loc_1_multilooked = zeros(Int(length(s_loc_S)/Block_size_S))
s_loc_2_multilooked = zeros(Int(length(s_loc_C)/Block_size_C))

norm_mean_pow1 = maximum(Pow_1)
norm_mean_pow2 = maximum(Pow_1)


        k=1
        for i=1:Int((length(s_loc_S))/Block_size_S)
            l=1
            for j=1:Int((length(s_loc_C))/Block_size_C)
                Pow_1_multilooked[i,j] = mean(Pow_1[k:k+Block_size_S-1,l:l+Block_size_C-1]) ./norm_mean_pow1
                Pow_2_multilooked[i,j] = mean(Pow_2[k:k+Block_size_S-1,l:l+Block_size_C-1]) ./norm_mean_pow2
                s_loc_1_multilooked[i] = mean(s_loc_S[k:k+Block_size_S-1])
                s_loc_2_multilooked[j] = mean(s_loc_C[l:l+Block_size_C-1])
                l=l+Block_size_C
            end
            k=k+Block_size_S
        end

plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Pow_1_multilooked,"log", savepath*"Power_plots/", "ML_Power_P1_log", "")
plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Pow_2_multilooked,"log", savepath*"Power_plots/", "ML_Power_P2_log", "")
        
plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Pow_1_multilooked,"lin", savepath*"Power_plots/", "ML_Power_P1_lin", "")
plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Pow_2_multilooked,"lin", savepath*"Power_plots/", "ML_Power_P2_lin", "")

# Power statistics plots 
plot_var                = Pow_1_multilooked[:]
A = plot_histogram_only(plot_var[:],"Power", savepath*"Power_plots/", "ML_P1_pow")
plot_var                = Pow_2_multilooked[:]
A = plot_histogram_only(plot_var[:],"Power", savepath*"Power_plots/", "ML_P2_pow")



#Interferometry

Data_12 		        = Data_1 .* conj(Data_2)

Pow_12                  = abs.(Data_12)
Phase_12                = angle.(Data_12)

plot_image(s_loc_C,s_loc_S,Pow_12,"log", savepath*"Int_plots/", "Int_Power_P12_log", "")
plot_image(s_loc_C,s_loc_S,Pow_12,"lin", savepath*"Int_plots/", "Int_Power_P12_lin", "")

plot_image(s_loc_C,s_loc_S,Pow_12 ./ ((Pow_1.*Pow_2).^0.5),"log", savepath*"Int_plots/", "Int_Coh_Power_P12_log", "")
plot_image(s_loc_C,s_loc_S,Pow_12 ./ ((Pow_1.*Pow_2).^0.5),"lin", savepath*"Int_plots/", "Int_Coh_Power_P12_lin", "")

plot_image(s_loc_C,s_loc_S,Phase_12,"lin", savepath*"Int_plots/", "Int_Phase_P12_lin", "")


Pow_12_rangeprofile   =  mean(Pow_12,dims=1)

plot_profile(s_loc_C, 10 .* log10.(Pow_12_rangeprofile[:]), "", savepath*"Int_plots/", "Pow_12_profile", "Lon [deg]", "Power [dB] ", "", "", "")



Block_size_S = Looks_along_S
Block_size_C = Looks_along_C

Int_Pow_multilooked = zeros(Int(length(s_loc_S)/Block_size_S),Int(length(s_loc_C)/Block_size_C))
Int_Phase_multilooked = zeros(Int(length(s_loc_S)/Block_size_S),Int(length(s_loc_C)/Block_size_C))
Correlation = zeros(Int(length(s_loc_S)/Block_size_S),Int(length(s_loc_C)/Block_size_C))

s_loc_1_multilooked = zeros(Int(length(s_loc_S)/Block_size_S))
s_loc_2_multilooked = zeros(Int(length(s_loc_C)/Block_size_C))

        k=1
        for i=1:Int((length(s_loc_S))/Block_size_S)
            l=1
            for j=1:Int((length(s_loc_C))/Block_size_C)
                complex_coherence = mean(Data_12[k:k+Block_size_S-1,l:l+Block_size_C-1])./ ((mean(Pow_1[k:k+Block_size_S-1,l:l+Block_size_C-1]) .* mean(Pow_2[k:k+Block_size_S-1,l:l+Block_size_C-1])).^0.5)
                Int_Pow_multilooked[i,j] = abs.(complex_coherence)
                Int_Phase_multilooked[i,j] = angle.(complex_coherence) 

                s_loc_1_multilooked[i] = mean(s_loc_S[k:k+Block_size_S-1])
                s_loc_2_multilooked[j] = mean(s_loc_C[l:l+Block_size_C-1])
                l=l+Block_size_C
            end
            k=k+Block_size_S
        end


plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Int_Pow_multilooked,"log", savepath*"Int_plots/", "ML_Power_P12_log", "")
plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Int_Pow_multilooked,"lin", savepath*"Int_plots/", "ML_Power_P12_lin", "")
        
plot_image2(s_loc_2_multilooked,s_loc_1_multilooked,Int_Phase_multilooked,"phase", savepath*"Int_plots/", "ML_Phase_P12_lin", "",":twilight")
#plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Correlation,"lin", savepath*"Int_plots/", "ML_Correlation_lin", "")

Correlation_rangeprofile   = Int_Pow_multilooked[1,:] # mean(Int_Pow_multilooked,dims=1)
Int_phase_rangeprofile   = Int_Phase_multilooked[1,:] # mean(Int_Phase_multilooked,dims=1)

plot_profile(s_loc_2_multilooked, Correlation_rangeprofile[:], "", savepath*"Int_plots/", "Int_mag_profile", "Lon [deg]", "Magnitude  ", "", "", "")
plot_profile(s_loc_2_multilooked, Int_phase_rangeprofile[:], "", savepath*"Int_plots/", "Int_phase__profile", "Lon [deg]", "Phase [rad]  ", "", "", "")


A = plot_histogram_only(Int_Pow_multilooked[:], "Magnitude ",savepath*"Int_plots/", "ML_Int_mag")
A = plot_histogram_only(Int_Phase_multilooked[:], "Phase",savepath*"Int_plots/", "ML_Int_ph")

#gamma_ip = mean(Int_Pow_multilooked[:])
gamma_ip = mean(Correlation_theo_all[:,1])

A = plot_interferogram_histogram_statistics_phase(Int_Phase_multilooked[:], gamma_ip, Looks_along_S*Looks_along_C,0.0, (-0.5,0.5), savepath*"Int_plots/", "PDF_Phase_Han1")
A = plot_interferogram_histogram_statistics_phase_han1(Int_Phase_multilooked[:], gamma_ip, Looks_along_S*Looks_along_C, 0.0, (-0.5,0.5), savepath*"Int_plots/", "PDF_Phase_Han2")
B = plot_interferogram_histogram_statistics_magnitude(Int_Pow_multilooked[:], gamma_ip, Looks_along_S*Looks_along_C, (0.5,1), savepath*"Int_plots/", "PDF_Mag_Han1")

Correlation_rangeprofile_window = movmean(Correlation_rangeprofile, 40)
Int_phase_rangeprofile_window = movmean(Int_phase_rangeprofile, 40)

plot_profile(s_loc_2_multilooked, Correlation_rangeprofile_window[:], "", savepath*"Int_plots/", "Int_mag_profile_window", "Lon [deg]", "Magnitude  ", "", "", "")
plot_profile(s_loc_2_multilooked, Int_phase_rangeprofile_window[:], "", savepath*"Int_plots/", "Int_phase_profile_window", "Lon [deg]", "Phase [rad]  ", "", "", "")


CF_A, CF_B = linear_fit(s_loc_2_multilooked,Correlation_rangeprofile_window)

p1=(plot(s_loc_2_multilooked, Correlation_rangeprofile_window[:],xlabel="Lon [deg]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:black, label="",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
p1=(plot!(s_loc_2_multilooked, (s_loc_2_multilooked .* CF_B).+CF_A,xlabel="Lon [deg]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:red, label="Linear Fit",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p1, savepath*"Int_plots/"*"Int_mag_profile_window2"*".png")



#=
p2 = (plot(s_loc_C, Correlation_theo_P1_rangeprofile[:], label="Theory"))
p2 = (plot!(s_loc_2_multilooked, (s_loc_2_multilooked .* CF_B).+CF_A, label="Sim"))

AASR = 10 ^(-26/10)
p2 = (plot!(s_loc_2_multilooked, ((s_loc_2_multilooked .* CF_B).+CF_A) * (1/1+(AASR)), label="Sim - AASR -26 dB", legend=:topleft, xlabel="Lon [deg]", ylabel="Correlation",
topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p2, savepath*"Int_plots/"*"Correlation_comparison_1"*".png")

AASR = 10 ^(-21.9/10)
p2 = (plot!(s_loc_2_multilooked, ((s_loc_2_multilooked .* CF_B).+CF_A) * (1/1+(AASR)), label="Sim - AASR -22 dB", legend=:topleft, xlabel="Lon [deg]", ylabel="Correlation",
topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p2, savepath*"Int_plots/"*"Correlation_comparison_2"*".png")
=#


## Geometry plots

#DEM
DEM_region_rangeprofile    = DEM_region[1,:] # mean(DEM_region,dims=1)
plot_image(s_loc_C,s_loc_S,DEM_region,"lin", savepath*"Geom_plots/", "DEM", "")
plot_profile(s_loc_C, DEM_region_rangeprofile[:], "", savepath*"Geom_plots/", "DEM_profile", "Lon [deg]", "DEM height [m]","", "", "", )


#BRCS
targets_ref_corr = reshape(trg_targets_ref_corr[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
targets_ref_corr_rangeprofile    = targets_ref_corr[1,:] # mean(targets_ref_corr,dims=1)
plot_image(s_loc_C,s_loc_S,targets_ref_corr,"lin", savepath*"Geom_plots/", "RCS", "")
plot_profile(s_loc_C, targets_ref_corr_rangeprofile[:], "", savepath*"Geom_plots/", "RCS_profile", "Lon [deg]", "Reflectivity from RCS","", "", "", )

targets_ref_corr = 10 .* log10.(reshape(trg_targets_ref_corr[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]')
targets_ref_corr = targets_ref_corr ./ maximum(targets_ref_corr)
targets_ref_corr_rangeprofile    = targets_ref_corr[1,:] # mean(targets_ref_corr,dims=1)
plot_image(s_loc_C,s_loc_S,targets_ref_corr,"lin", savepath*"Geom_plots/", "RCS_2", "")
plot_profile(s_loc_C, targets_ref_corr_rangeprofile[:], "", savepath*"Geom_plots/", "RCS_profile_2", "Lon [deg]", "Reflectivity from RCS","", "", "", )



#Look angle
Lookangle_P1                = reshape(look_angle_all[:,1],length(s_loc_C),length(s_loc_S))'
Lookangle_P2                = reshape(look_angle_all[:,2],length(s_loc_C),length(s_loc_S))'
Lookangle_Diff              = Lookangle_P1 - Lookangle_P2

Lookangle_P1_rangeprofile   =  Lookangle_P1[1,:] #mean(Lookangle_P1,dims=1)
Lookangle_P2_rangeprofile   =  Lookangle_P2[1,:] #mean(Lookangle_P2,dims=1)
Lookangle_Diff_rangeprofile =  Lookangle_Diff[1,:] #mean(Lookangle_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Lookangle_P1,"lin", savepath*"Geom_plots/", "Lookangle_P1", "")
plot_image(s_loc_C,s_loc_S,Lookangle_P2,"lin", savepath*"Geom_plots/", "Lookangle_P2", "")
plot_image(s_loc_C,s_loc_S,Lookangle_Diff,"lin", savepath*"Geom_plots/", "Lookangle_Diff", "")

plot_profile(s_loc_C, Lookangle_P1_rangeprofile[:], Lookangle_P2_rangeprofile[:], savepath*"Geom_plots/", "Lookangle_profile1", "Lon [deg]", "Look Angle [deg]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Lookangle_Diff_rangeprofile[:], "", savepath*"Geom_plots/", "Lookangle_profile2", "Lon [deg]", "Look Angle Diff [deg]", "", "Difference", "")

#Look angle -trg
Lookangle_P1                = reshape(trg_look_angle_all[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Lookangle_P2                = reshape(trg_look_angle_all[:,2],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Lookangle_Diff              = Lookangle_P1 - Lookangle_P2

Lookangle_P1_rangeprofile   =  Lookangle_P1[1,:] #mean(Lookangle_P1,dims=1)
Lookangle_P2_rangeprofile   =  Lookangle_P2[1,:] #mean(Lookangle_P2,dims=1)
Lookangle_Diff_rangeprofile =  Lookangle_Diff[1,:] #mean(Lookangle_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Lookangle_P1,"lin", savepath*"Geom_plots/", "trg_Lookangle_P1", "")
plot_image(s_loc_C,s_loc_S,Lookangle_P2,"lin", savepath*"Geom_plots/", "trg_Lookangle_P2", "")
plot_image(s_loc_C,s_loc_S,Lookangle_Diff,"lin", savepath*"Geom_plots/", "trg_Lookangle_Diff", "")

plot_profile(s_loc_C, Lookangle_P1_rangeprofile[:], Lookangle_P2_rangeprofile[:], savepath*"Geom_plots/", "trg_Lookangle_profile1", "Lon [deg]", "Look Angle [deg]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Lookangle_Diff_rangeprofile[:], "", savepath*"Geom_plots/", "trg_Lookangle_profile2", "Lon [deg]", "Look Angle Diff [deg]", "", "Difference", "")


#Incidence angle
Incangle_P1                 = reshape(incidence_angle_all[:,1],length(s_loc_C),length(s_loc_S))'
Incangle_P2                 = reshape(incidence_angle_all[:,2],length(s_loc_C),length(s_loc_S))'
Incangle_Diff               = Incangle_P1 - Incangle_P2

Incangle_P1_rangeprofile    =  Incangle_P1[1,:] #mean(Incangle_P1,dims=1)
Incangle_P2_rangeprofile    =  Incangle_P2[1,:] #mean(Incangle_P2,dims=1)
Incangle_Diff_rangeprofile  =  Incangle_Diff[1,:] #mean(Incangle_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Incangle_P1,"lin", savepath*"Geom_plots/", "Incangle_P1", "")
plot_image(s_loc_C,s_loc_S,Incangle_P2,"lin", savepath*"Geom_plots/", "Incangle_P2", "")
plot_image(s_loc_C,s_loc_S,Incangle_Diff,"lin", savepath*"Geom_plots/", "Incangle_Diff", "")

plot_profile(s_loc_C, Incangle_P1_rangeprofile[:], Incangle_P2_rangeprofile[:], savepath*"Geom_plots/", "Incangle_profile1", "Lon [deg]", "Incidence Angle [deg]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Incangle_Diff_rangeprofile[:], "", savepath*"Geom_plots/", "Incangle_profile2", "Lon [deg]", "Incidence Angle Diff [deg]", "", "Difference", "")

#Incidence angle -trg
Incangle_P1                 = reshape(trg_incidence_angle_all[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Incangle_P2                 = reshape(trg_incidence_angle_all[:,2],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Incangle_Diff               = Incangle_P1 - Incangle_P2

Incangle_P1_rangeprofile    =  Incangle_P1[1,:] #mean(Incangle_P1,dims=1)
Incangle_P2_rangeprofile    =  Incangle_P2[1,:] #mean(Incangle_P2,dims=1)
Incangle_Diff_rangeprofile  =  Incangle_Diff[1,:] #mean(Incangle_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Incangle_P1,"lin", savepath*"Geom_plots/", "trg_Incangle_P1", "")
plot_image(s_loc_C,s_loc_S,Incangle_P2,"lin", savepath*"Geom_plots/", "trg_Incangle_P2", "")
plot_image(s_loc_C,s_loc_S,Incangle_Diff,"lin", savepath*"Geom_plots/", "trg_Incangle_Diff", "")

plot_profile(s_loc_C, Incangle_P1_rangeprofile[:], Incangle_P2_rangeprofile[:], savepath*"Geom_plots/", "trg_Incangle_profile1", "Lon [deg]", "Incidence Angle [deg]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Incangle_Diff_rangeprofile[:], "", savepath*"Geom_plots/", "trg_Incangle_profile2", "Lon [deg]", "Incidence Angle Diff [deg]", "", "Difference", "")


#Local Incidence angle
LocalIncangle_P1                 = reshape(local_incidence_angle[:,1],length(s_loc_C),length(s_loc_S))'
LocalIncangle_P1_rangeprofile    =  LocalIncangle_P1[1,:] #mean(LocalIncangle_P1,dims=1)
plot_image(s_loc_C,s_loc_S,LocalIncangle_P1,"lin", savepath*"Geom_plots/", "LocalIncangle", "")
plot_profile(s_loc_C, LocalIncangle_P1_rangeprofile[:], "", savepath*"Geom_plots/", "LocalIncangle_profile", "Lon [deg]", "Local incidence Angle [deg]","", "", "", )

#Local Incidence angle -trg
LocalIncangle_P1                 = reshape(trg_local_incidence_angle[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
LocalIncangle_P1_rangeprofile    =  LocalIncangle_P1[1,:] #mean(LocalIncangle_P1,dims=1)
plot_image(s_loc_C,s_loc_S,LocalIncangle_P1,"lin", savepath*"Geom_plots/", "trg_LocalIncangle", "")
plot_profile(s_loc_C, LocalIncangle_P1_rangeprofile[:], "", savepath*"Geom_plots/", "trg_LocalIncangle_profile", "Lon [deg]", "Local incidence Angle [deg]","", "", "", )

#Range slope angle
range_slope_P1                 = reshape(range_slope_angle[:,1],length(s_loc_C),length(s_loc_S))'
range_slope_P1_rangeprofile    =  range_slope_P1[1,:] #mean(range_slope_P1,dims=1)
plot_image(s_loc_C,s_loc_S,range_slope_P1,"lin", savepath*"Geom_plots/", "RangeSlopeangle", "")
plot_profile(s_loc_C, range_slope_P1_rangeprofile[:], "", savepath*"Geom_plots/", "RangeSlopeangle_profile", "Lon [deg]", "Range Slope Angle [deg]","", "", "", )

#Range slope angle -trg
range_slope_P1                 = reshape(trg_range_slope_angle[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
range_slope_P1_rangeprofile    =  range_slope_P1[1,:] #mean(range_slope_P1,dims=1)
plot_image(s_loc_C,s_loc_S,range_slope_P1,"lin", savepath*"Geom_plots/", "trg_RangeSlopeangle", "")
plot_profile(s_loc_C, range_slope_P1_rangeprofile[:], "", savepath*"Geom_plots/", "trg_RangeSlopeangle_profile", "Lon [deg]", "Range Slope Angle [deg]","", "", "", )



#Slant range
Slantrange_P1 = reshape(slant_range_all[:,1],length(s_loc_C),length(s_loc_S))'
Slantrange_P2 = reshape(slant_range_all[:,2],length(s_loc_C),length(s_loc_S))'
Slantrange_Diff = Slantrange_P1 - Slantrange_P2

Slantrange_P1_rangeprofile =  Slantrange_P1[1,:] #mean(Slantrange_P1,dims=1)
Slantrange_P2_rangeprofile =  Slantrange_P2[1,:] #mean(Slantrange_P2,dims=1)
Slantrange_Diff_rangeprofile =  Slantrange_Diff[1,:] #mean(Slantrange_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Slantrange_P1./1e3,"lin", savepath*"Geom_plots/", "Slantrange_P1", "")
plot_image(s_loc_C,s_loc_S,Slantrange_P2./1e3,"lin", savepath*"Geom_plots/", "Slantrange_P2", "")
plot_image(s_loc_C,s_loc_S,Slantrange_Diff,"phase", savepath*"Geom_plots/", "Slantrange_Diff", "")

plot_profile(s_loc_C, Slantrange_P1_rangeprofile[:], Slantrange_P2_rangeprofile[:], savepath*"Geom_plots/", "Slantrange_profile1", "Lon [deg]", "Slant range [m]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Slantrange_Diff_rangeprofile[:], "", savepath*"Geom_plots/", "Slantrange_profile2", "Lon [deg]", "Slant range [m]", "", "Difference", "")

#Slant range -trg
Slantrange_P12 = reshape(trg_slant_range_all[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Slantrange_P22 = reshape(trg_slant_range_all[:,2],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Slantrange_Diff2 = Slantrange_P1 - Slantrange_P2

Slantrange_P12_rangeprofile =  Slantrange_P12[1,:] #mean(Slantrange_P1,dims=1)
Slantrange_P22_rangeprofile =  Slantrange_P22[1,:] #mean(Slantrange_P2,dims=1)
Slantrange_Diff2_rangeprofile =  Slantrange_Diff2[1,:] #mean(Slantrange_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Slantrange_P12./1e3,"lin", savepath*"Geom_plots/", "trg_Slantrange_P1", "")
plot_image(s_loc_C,s_loc_S,Slantrange_P22./1e3,"lin", savepath*"Geom_plots/", "trg_Slantrange_P2", "")
plot_image(s_loc_C,s_loc_S,Slantrange_Diff2,"phase", savepath*"Geom_plots/", "trg_Slantrange_Diff", "")

plot_profile(s_loc_C, Slantrange_P12_rangeprofile[:], Slantrange_P2_rangeprofile[:], savepath*"Geom_plots/", "trg_Slantrange_profile1", "Lon [deg]", "Slant range [m]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Slantrange_Diff2_rangeprofile[:], "", savepath*"Geom_plots/", "trg_Slantrange_profile2", "Lon [deg]", "Slant range [m]", "", "Difference", "")

#Synthetic fringes
Unwrapped_data = ( (4 * pi )/params.λ).*Slantrange_Diff
Wrapped_data = mod.(Unwrapped_data,2*pi) #.- pi

plot_image(s_loc_C,s_loc_S,Unwrapped_data,"phase", savepath*"Geom_plots/", "Slantrange_Diff_unwrapped", "")
plot_image(s_loc_C,s_loc_S,Wrapped_data,"phase", savepath*"Geom_plots/", "Slantrange_Diff_wrapped", "")

Wrapped_data_rangeprofile =  Wrapped_data[1,:] #mean(Wrapped_data,dims=1)
plot_profile(s_loc_C, Wrapped_data_rangeprofile[:], "", savepath*"Geom_plots/", "Slantrange_profile_wrapped", "Lon [deg]", "", "", "", "")
#plot_profile(s_loc_C[1000:1150], Wrapped_data_rangeprofile[1,1000:1150], "", savepath*"Geom_plots/", "Slantrange_profile_wrapped2", "Lon [deg]", "", "", "", "")

Data_3 = Data_2 .* exp.(1im .* Unwrapped_data)
#Data_3 = Data_2 .* exp.(-1im .* Wrapped_data)

# Critical baseline
Critical_baseline_P1                = reshape(Critical_baseline_all[:,1],length(s_loc_C),length(s_loc_S))'
Critical_baseline_P1_rangeprofile   =  Critical_baseline_P1[1,:] #mean(Critical_baseline_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Critical_baseline_P1,"lin", savepath*"Geom_plots/", "Critical_baseline_P1", "")
plot_profile(s_loc_C, Critical_baseline_P1_rangeprofile[:], "", savepath*"Geom_plots/", "Critical_baseline_profile2", "Lon [deg]", "Critical baseline  [m]", "", "", "")

# Critical baseline - trg
Critical_baseline_P1                = reshape(trg_Critical_baseline_all[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Critical_baseline_P1_rangeprofile   =  Critical_baseline_P1[1,:] #mean(Critical_baseline_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Critical_baseline_P1,"lin", savepath*"Geom_plots/", "trg_Critical_baseline_P1", "")
plot_profile(s_loc_C, Critical_baseline_P1_rangeprofile[:], "", savepath*"Geom_plots/", "trg_Critical_baseline_profile2", "Lon [deg]", "Critical baseline  [m]", "", "", "")

# Perp baseline
Perp_baseline_P1                = reshape(Perp_baseline_all[:,1],length(s_loc_C),length(s_loc_S))'
Perp_baseline_P1_rangeprofile   =  Perp_baseline_P1[1,:] #mean(Perp_baseline_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Perp_baseline_P1,"lin", savepath*"Geom_plots/", "Perp_baseline_P1", "")
plot_profile(s_loc_C, Perp_baseline_P1_rangeprofile[:], "", savepath*"Geom_plots/", "Perp_baseline_profile2","Lon [deg]", "Perp baseline  [m]", "", "", "")

# Perp baseline -trg
Perp_baseline_P1                = reshape(trg_Perp_baseline_all[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Perp_baseline_P1_rangeprofile   =  Perp_baseline_P1[1,:] #mean(Perp_baseline_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Perp_baseline_P1,"lin", savepath*"Geom_plots/", "trg_Perp_baseline_P1", "")
plot_profile(s_loc_C, Perp_baseline_P1_rangeprofile[:], "", savepath*"Geom_plots/", "trg_Perp_baseline_profile2","Lon [deg]", "Perp baseline  [m]", "", "", "")

# vertical wavenumber
Vert_wavnum_P1                = reshape(Vert_wavnum_all[:,1],length(s_loc_C),length(s_loc_S))'
Vert_wavnum_P1_rangeprofile   =  Vert_wavnum_P1[1,:] #mean(Vert_wavnum_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Vert_wavnum_P1,"lin", savepath*"Geom_plots/", "Vert_wavnum_P1", "")
plot_profile(s_loc_C, Vert_wavnum_P1_rangeprofile[:], "", savepath*"Geom_plots/", "Vert_wavnum_profile2", "Lon [deg]", "Vert wave num  ", "", "", "")

# vertical wavenumber - trg
Vert_wavnum_P1                = reshape(trg_Vert_wavnum_all[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Vert_wavnum_P1_rangeprofile   =  Vert_wavnum_P1[1,:] #mean(Vert_wavnum_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Vert_wavnum_P1,"lin", savepath*"Geom_plots/", "trg_Vert_wavnum_P1", "")
plot_profile(s_loc_C, Vert_wavnum_P1_rangeprofile[:], "", savepath*"Geom_plots/", "trg_Vert_wavnum_profile2", "Lon [deg]", "Vert wave num  ", "", "", "")

# Theoretical correlation 
Correlation_theo_P1                = reshape(Correlation_theo_all[:,1],length(s_loc_C),length(s_loc_S))'
Correlation_theo_P1_rangeprofile   =  Correlation_theo_P1[1,:] #mean(Correlation_theo_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Correlation_theo_P1,"lin", savepath*"Geom_plots/", "Correlation_theo_P1", "")
plot_profile(s_loc_C, Correlation_theo_P1_rangeprofile[:], "", savepath*"Geom_plots/", "Correlation_theo_profile2", "Lon [deg]", "Correlation theory  ", "", "", "")

# Theoretical correlation  - trg
Correlation_theo_P1                = reshape(trg_Correlation_theo_all[:,1],3,length(s_loc_C),length(s_loc_S))[1,:,:]'
Correlation_theo_P1_rangeprofile   =  Correlation_theo_P1[1,:] #mean(Correlation_theo_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Correlation_theo_P1,"lin", savepath*"Geom_plots/", "trg_Correlation_theo_P1", "")
plot_profile(s_loc_C, Correlation_theo_P1_rangeprofile[:], "", savepath*"Geom_plots/", "trg_Correlation_theo_profile2", "Lon [deg]", "Correlation theory  ", "", "", "")



# TSlope_lat
Slope_lat                = atand.(-1 .* reshape(N_all[:,1],1,length(s_loc_C),length(s_loc_S))[1,:,:]')
Slope_lat_rangeprofile   =  Slope_lat[:,1] 

plot_image(s_loc_C,s_loc_S,Slope_lat,"lin", savepath*"Geom_plots/", "Slope_lat", "")
plot_profile(s_loc_S, Slope_lat_rangeprofile[:], "", savepath*"Geom_plots/", "Slope_lat_rangeprofile", "Lat [deg]", "Slope Lat  ", "", "", "")


# TSlope_lon
Slope_lon                = atand.(-1 .* reshape(N_all[:,2],1,length(s_loc_C),length(s_loc_S))[1,:,:]')
Slope_lon_rangeprofile   =  Slope_lon[1,:] 

plot_image(s_loc_C,s_loc_S,Slope_lon,"lin", savepath*"Geom_plots/", "Slope_lon", "")
plot_profile(s_loc_C, Slope_lon_rangeprofile[:], "", savepath*"Geom_plots/", "Slope_lon_rangeprofile", "Lon [deg]", "Slope Lat  ", "", "", "")











#Interferometry 2

Data_12 		        = Data_1 .* conj(Data_3)

Pow_12                  = abs.(Data_12)
Phase_12                = angle.(Data_12)

plot_image(s_loc_C,s_loc_S,Pow_12,"log", savepath*"Int_plots/", "2Int_Power_P12_log", "")
plot_image(s_loc_C,s_loc_S,Pow_12,"lin", savepath*"Int_plots/", "2Int_Power_P12_lin", "")

plot_image(s_loc_C,s_loc_S,Pow_12 ./ ((Pow_1.*Pow_2).^0.5),"log", savepath*"Int_plots/", "2Int_Coh_Power_P12_log", "")
plot_image(s_loc_C,s_loc_S,Pow_12 ./ ((Pow_1.*Pow_2).^0.5),"lin", savepath*"Int_plots/", "2Int_Coh_Power_P12_lin", "")

plot_image(s_loc_C,s_loc_S,Phase_12,"lin", savepath*"Int_plots/", "2Int_Phase_P12_lin", "")



Block_size_S = Looks_along_S
Block_size_C = Looks_along_C

Int_Pow_multilooked = zeros(Int(length(s_loc_S)/Block_size_S),Int(length(s_loc_C)/Block_size_C))
Int_Phase_multilooked = zeros(Int(length(s_loc_S)/Block_size_S),Int(length(s_loc_C)/Block_size_C))
Correlation = zeros(Int(length(s_loc_S)/Block_size_S),Int(length(s_loc_C)/Block_size_C))

s_loc_1_multilooked = zeros(Int(length(s_loc_S)/Block_size_S))
s_loc_2_multilooked = zeros(Int(length(s_loc_C)/Block_size_C))

        k=1
        for i=1:Int((length(s_loc_S))/Block_size_S)
            l=1
            for j=1:Int((length(s_loc_C))/Block_size_C)
                complex_coherence = mean(Data_12[k:k+Block_size_S-1,l:l+Block_size_C-1])./ ((mean(Pow_1[k:k+Block_size_S-1,l:l+Block_size_C-1]) .* mean(Pow_2[k:k+Block_size_S-1,l:l+Block_size_C-1])).^0.5)
                Int_Pow_multilooked[i,j] = abs.(complex_coherence)
                Int_Phase_multilooked[i,j] = angle.(complex_coherence) #angle.(mean(Data_12[k:k+Block_size_S-1,l:l+Block_size_C-1]))
                #Correlation[i,j] = abs.(mean(Data_12[k:k+Block_size_S-1,l:l+Block_size_C-1])) ./ ((mean(Pow_1[k:k+Block_size_S-1,l:l+Block_size_C-1]) .* mean(Pow_2[k:k+Block_size_S-1,l:l+Block_size_C-1])).^0.5)

                s_loc_1_multilooked[i] = mean(s_loc_S[k:k+Block_size_S-1])
                s_loc_2_multilooked[j] = mean(s_loc_C[l:l+Block_size_C-1])
                l=l+Block_size_C
            end
            k=k+Block_size_S
        end


plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Int_Pow_multilooked,"log", savepath*"Int_plots/", "2ML_Power_P12_log", "")
plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Int_Pow_multilooked,"lin", savepath*"Int_plots/", "2ML_Power_P12_lin", "")
        
plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Int_Phase_multilooked,"phase", savepath*"Int_plots/", "2ML_Phase_P12_lin", "")
#plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Correlation,"lin", savepath*"Int_plots/", "ML_Correlation_lin", "")

Correlation_rangeprofile   =  mean(Int_Pow_multilooked,dims=1)
Int_phase_rangeprofile   =  mean(Int_Phase_multilooked,dims=1)

plot_profile(s_loc_2_multilooked, Correlation_rangeprofile[:], "", savepath*"Int_plots/", "2Int_mag_profile", "Lon [deg]", "Magnitude  ", "", "", "")
plot_profile(s_loc_2_multilooked, Int_phase_rangeprofile[:], "", savepath*"Int_plots/", "2Int_phase__profile", "Lon [deg]", "Phase [rad]  ", "", "", "")


A = plot_histogram_only(Int_Pow_multilooked[:], "Magnitude ",savepath*"Int_plots/", "2ML_Int_mag")
A = plot_histogram_only(Int_Phase_multilooked[:], "Phase",savepath*"Int_plots/", "2ML_Int_ph")


Correlation_rangeprofile_window = movmean(Correlation_rangeprofile, 40)
Int_phase_rangeprofile_window = movmean(Int_phase_rangeprofile, 40)

plot_profile(s_loc_2_multilooked, Correlation_rangeprofile_window[:], "", savepath*"Int_plots/", "2Int_mag_profile_window", "Lon [deg]", "Magnitude  ", "", "", "")
plot_profile(s_loc_2_multilooked, Int_phase_rangeprofile_window[:], "", savepath*"Int_plots/", "2Int_phase_profile_window", "Lon [deg]", "Phase [rad]  ", "", "", "")


CF_A, CF_B = linear_fit(s_loc_2_multilooked,Correlation_rangeprofile_window)

p1=(plot(s_loc_2_multilooked, Correlation_rangeprofile_window[:],xlabel="Lon [deg]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:black, label="",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
p1=(plot!(s_loc_2_multilooked, (s_loc_2_multilooked .* CF_B).+CF_A,xlabel="Lon [deg]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:red, label="Linear Fit",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p1, savepath*"Int_plots/"*"2Int_mag_profile_window2"*".png")


#gamma_ip = mean(Int_Pow_multilooked[:])
gamma_ip = mean(Correlation_theo_all[:,1])

A = plot_interferogram_histogram_statistics_phase(Int_Phase_multilooked[:], gamma_ip, Looks_along_S*Looks_along_C,0.0, (-0.5,0.5), savepath*"Int_plots/", "2PDF_Phase_Han1")
A = plot_interferogram_histogram_statistics_phase_han1(Int_Phase_multilooked[:], gamma_ip, Looks_along_S*Looks_along_C, 0.0, (-0.5,0.5), savepath*"Int_plots/", "2PDF_Phase_Han2")
B = plot_interferogram_histogram_statistics_magnitude(Int_Pow_multilooked[:], gamma_ip, Looks_along_S*Looks_along_C, (0.5,1), savepath*"Int_plots/", "2PDF_Mag_Han1")

Correlation_rangeprofile_window = movmean(Correlation_rangeprofile, 40)
Int_phase_rangeprofile_window = movmean(Int_phase_rangeprofile, 40)

plot_profile(s_loc_2_multilooked, Correlation_rangeprofile_window[:], "", savepath*"Int_plots/", "2Int_mag_profile_window", "Lon [deg]", "Magnitude  ", "", "", "")
plot_profile(s_loc_2_multilooked, Int_phase_rangeprofile_window[:], "", savepath*"Int_plots/", "2Int_phase_profile_window", "Lon [deg]", "Phase [rad]  ", "", "", "")


CF_A, CF_B = linear_fit(s_loc_2_multilooked,Correlation_rangeprofile_window)

p1=(plot(s_loc_2_multilooked, Correlation_rangeprofile_window[:],xlabel="Lon [deg]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:black, label="",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
p1=(plot!(s_loc_2_multilooked, (s_loc_2_multilooked .* CF_B).+CF_A,xlabel="Lon [deg]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:red, label="Linear Fit",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p1, savepath*"Int_plots/"*"2Int_mag_profile_window2"*".png")



#=
# To do processing again with master as reference
# Processing - 2
include("../modules/process_raw_data.jl")

SAR_images_3D_proc2 = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
image_3D_proc2      = Process_Raw_Data.tomo_processing_afterSAR_full(SAR_images_3D_proc2)


stat_var_all_1p_proc2 = SAR_images_3D_proc2[1,:,:,1]
stat_var_all_2p_proc2 = SAR_images_3D_proc2[2,:,:,1]
=#


