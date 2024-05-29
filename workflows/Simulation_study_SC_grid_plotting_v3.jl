
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

function plot_image(x_axis,y_axis,Data,unit_flag, savepath, savename, figure_title)

    if ~ispath(savepath)
        mkdir(savepath)
    end    

    if unit_flag == "log"
        p1=(heatmap(x_axis,y_axis,10 .* log10.(abs.(Data)),xlabel="C [m]",ylabel="S [m]",title=figure_title,
        topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "lin"
        p1=(heatmap(x_axis,y_axis,(abs.(Data)),xlabel="C [m]",ylabel="S [m]", title=figure_title,
        topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #1200
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "phase"
        p1=(heatmap(x_axis,y_axis,(Data),xlabel="C [m]",ylabel="S [m]", title=figure_title,
        topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #1200
        savefig(p1, savepath*savename*".png")
    end
end

function plot_profile(x_axis, Data1, Data2, savepath, savename, xlabel, ylabel, figure_title, label1, label2)

    if Data2==""
        p1=(plot(x_axis, Data1,xlabel=xlabel,ylabel=ylabel,title=figure_title,legend=:topleft, lc=:black, label=label1,
        topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
        savefig(p1, savepath*savename*".png")
    else
        p1=(plot(x_axis, Data1,xlabel=xlabel,ylabel=ylabel,title=figure_title, legend=:topleft, lc=:black, label=label1,
        topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
        p1=(plot!(x_axis, Data2,xlabel=xlabel,ylabel=ylabel,title=figure_title, legend=:topleft,lc=:blue,label=label2,
        topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
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
    axis_plot              = range(0+(maximum(plot_var_hist)/50),maximum(plot_var_hist)-(maximum(plot_var_hist)/50), length=50)
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
    axis_plot               = range(0+(maximum(plot_var_hist)/50),maximum(plot_var_hist)-(maximum(plot_var_hist)/50), length=50)
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
    axis_plot               = range(-pi+(pi/50),pi-(pi/50), length=50)
    theory_val              = repeat([1/(2*pi)], size(axis_plot)[1],1) 

    (plot(axis_plot, pdf1, ylim=(0,1/pi), xlim=(-1.5*pi,1.5*pi),xlabel="Phase",ylabel="PDF",title = ""
        ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
    p3  =(plot!(axis_plot, theory_val, ylim=(0,1/pi), xlim=(-1.5*pi,1.5*pi),xlabel="Phase",ylabel="PDF",title = ""
        ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
    savefig(p3, savepath*savename*"_Hist_plot2"*".png")

    return plot_var_hist

end


#------------------------------------------------------------------------------------------

#Load the output file
@load "/Users/joshil/Documents/Outputs/epstein_temp/temp_data/Output.jld" 

# Mention the directory to save plots
savepath                 = "/Users/joshil/Documents/Outputs/epstein_temp/temp_data/temp_plot/"
if ~ispath(savepath)
    mkdir(savepath)
end

t_loc_S_range       = -516:8:516
t_loc_C_range       = -520:4:516
t_loc_H_range       = 0

s_loc_S             = params.s_loc_1
s_loc_C             = params.s_loc_2
s_loc_H             = params.s_loc_3

Data_1              = stat_var_all_1p
Data_2              = stat_var_all_2p

#Data_1              = stat_var_all_1p_proc2
#Data_2              = stat_var_all_2p_proc2

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


plot_image(s_loc_C,s_loc_S,Amp_1,"log", savepath*"Amplitude_plots/", "Amp_P1_log", "")
plot_image(s_loc_C,s_loc_S,Amp_2,"log", savepath*"Amplitude_plots/", "Amp_P2_log", "")

plot_image(s_loc_C,s_loc_S,Amp_1,"lin", savepath*"Amplitude_plots/", "Amp_P1_lin", "")
plot_image(s_loc_C,s_loc_S,Amp_2,"lin", savepath*"Amplitude_plots/", "Amp_P2_lin", "")

plot_image(s_loc_C,s_loc_S,Amp_1./mean(Amp_1),"lin", savepath*"Amplitude_plots/", "Amp_P1_lin_norm", "")
plot_image(s_loc_C,s_loc_S,Amp_2./mean(Amp_2),"lin", savepath*"Amplitude_plots/", "Amp_P2_lin_norm", "")

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
Block_size_S = 13
Block_size_C = 10


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



Block_size_S = 13
Block_size_C = 10

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
        
plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Int_Phase_multilooked,"phase", savepath*"Int_plots/", "ML_Phase_P12_lin", "")
#plot_image(s_loc_2_multilooked,s_loc_1_multilooked,Correlation,"lin", savepath*"Int_plots/", "ML_Correlation_lin", "")

Correlation_rangeprofile   =  mean(Int_Pow_multilooked,dims=1)
Int_phase_rangeprofile   =  mean(Int_Phase_multilooked,dims=1)

plot_profile(s_loc_2_multilooked, Correlation_rangeprofile[:], "", savepath*"Int_plots/", "Int_mag_profile", "C [m]", "Magnitude  ", "", "", "")
plot_profile(s_loc_2_multilooked, Int_phase_rangeprofile[:], "", savepath*"Int_plots/", "Int_phase__profile", "C [m]", "Phase [deg]  ", "", "", "")


A = plot_histogram_only(Int_Pow_multilooked[:], "Magnitude ",savepath*"Int_plots/", "ML_Int_mag")
A = plot_histogram_only(Int_Phase_multilooked[:], "Phase",savepath*"Int_plots/", "ML_Int_ph")


Correlation_rangeprofile_window = movmean(Correlation_rangeprofile, 40)
Int_phase_rangeprofile_window = movmean(Int_phase_rangeprofile, 40)

plot_profile(s_loc_2_multilooked, Correlation_rangeprofile_window[:], "", savepath*"Int_plots/", "Int_mag_profile_window", "C [m]", "Magnitude  ", "", "", "")
plot_profile(s_loc_2_multilooked, Int_phase_rangeprofile_window[:], "", savepath*"Int_plots/", "Int_phase_profile_window", "C [m]", "Phase [deg]  ", "", "", "")


CF_A, CF_B = linear_fit(s_loc_2_multilooked,Correlation_rangeprofile_window)

p1=(plot(s_loc_2_multilooked, Correlation_rangeprofile_window[:],xlabel="C [m]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:black, label="",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
p1=(plot!(s_loc_2_multilooked, (s_loc_2_multilooked .* CF_B).+CF_A,xlabel="C [m]",ylabel="Magnitude  ",title="",legend=:topleft, lc=:red, label="Linear Fit",
topmargin=6mm,bottommargin=6mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(1200,360) )) #500,360
savefig(p1, savepath*"Int_plots/"*"Int_mag_profile_window2"*".png")


## Geometry plots

#Look angle
Lookangle_P1                = reshape(look_angle_all[:,1],length(s_loc_C),length(s_loc_S))'
Lookangle_P2                = reshape(look_angle_all[:,2],length(s_loc_C),length(s_loc_S))'
Lookangle_Diff              = Lookangle_P1 - Lookangle_P2

Lookangle_P1_rangeprofile   =  mean(Lookangle_P1,dims=1)
Lookangle_P2_rangeprofile   =  mean(Lookangle_P2,dims=1)
Lookangle_Diff_rangeprofile =  mean(Lookangle_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Lookangle_P1,"lin", savepath*"Geom_plots/", "Lookangle_P1", "")
plot_image(s_loc_C,s_loc_S,Lookangle_P2,"lin", savepath*"Geom_plots/", "Lookangle_P2", "")
plot_image(s_loc_C,s_loc_S,Lookangle_Diff,"lin", savepath*"Geom_plots/", "Lookangle_Diff", "")

plot_profile(s_loc_C, Lookangle_P1_rangeprofile[:], Lookangle_P2_rangeprofile[:], savepath*"Geom_plots/", "Lookangle_profile1", "C [m]", "Look Angle [deg]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Lookangle_Diff_rangeprofile[:], "", savepath*"Geom_plots/", "Lookangle_profile2", "C [m]", "Look Angle Diff [deg]", "", "Difference", "")

#Incidence angle
Incangle_P1                 = reshape(incidence_angle_all[:,1],length(s_loc_C),length(s_loc_S))'
Incangle_P2                 = reshape(incidence_angle_all[:,2],length(s_loc_C),length(s_loc_S))'
Incangle_Diff               = Incangle_P1 - Incangle_P2

Incangle_P1_rangeprofile    =  mean(Incangle_P1,dims=1)
Incangle_P2_rangeprofile    =  mean(Incangle_P2,dims=1)
Incangle_Diff_rangeprofile  =  mean(Incangle_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Incangle_P1,"lin", savepath*"Geom_plots/", "Incangle_P1", "")
plot_image(s_loc_C,s_loc_S,Incangle_P2,"lin", savepath*"Geom_plots/", "Incangle_P2", "")
plot_image(s_loc_C,s_loc_S,Incangle_Diff,"lin", savepath*"Geom_plots/", "Incangle_Diff", "")

plot_profile(s_loc_C, Incangle_P1_rangeprofile[:], Incangle_P2_rangeprofile[:], savepath*"Geom_plots/", "Incangle_profile1", "C [m]", "Incidence Angle [deg]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Incangle_Diff_rangeprofile[:], "", savepath*"Geom_plots/", "Incangle_profile2", "C [m]", "Incidence Angle Diff [deg]", "", "Difference", "")

#Slant range
Slantrange_P1 = reshape(slant_range_all[:,1],length(s_loc_C),length(s_loc_S))'
Slantrange_P2 = reshape(slant_range_all[:,2],length(s_loc_C),length(s_loc_S))'
Slantrange_Diff = Slantrange_P1 - Slantrange_P2

Slantrange_P1_rangeprofile =  mean(Slantrange_P1,dims=1)
Slantrange_P2_rangeprofile =  mean(Slantrange_P2,dims=1)
Slantrange_Diff_rangeprofile =  mean(Slantrange_Diff,dims=1)

plot_image(s_loc_C,s_loc_S,Slantrange_P1./1e3,"lin", savepath*"Geom_plots/", "Slantrange_P1", "")
plot_image(s_loc_C,s_loc_S,Slantrange_P2./1e3,"lin", savepath*"Geom_plots/", "Slantrange_P2", "")
plot_image(s_loc_C,s_loc_S,Slantrange_Diff,"lin", savepath*"Geom_plots/", "Slantrange_Diff", "")

plot_profile(s_loc_C, Slantrange_P1_rangeprofile[:], Slantrange_P2_rangeprofile[:], savepath*"Geom_plots/", "Slantrange_profile1", "C [m]", "Slant range [m]", "", "Platform 1", "Platform 2")
plot_profile(s_loc_C, Slantrange_Diff_rangeprofile[:], "", savepath*"Geom_plots/", "Slantrange_profile2", "C [m]", "Slant range [m]", "", "Difference", "")


Critical_baseline_P1                = reshape(Critical_baseline_all[:,1],length(s_loc_C),length(s_loc_S))'
Critical_baseline_P1_rangeprofile   =  mean(Critical_baseline_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Critical_baseline_P1,"lin", savepath*"Geom_plots/", "Critical_baseline_P1", "")
plot_profile(s_loc_C, Critical_baseline_P1_rangeprofile[:], "", savepath*"Geom_plots/", "Critical_baseline_profile2", "C [m]", "Critical baseline  [m]", "", "", "")


Perp_baseline_P1                = reshape(Perp_baseline_all[:,1],length(s_loc_C),length(s_loc_S))'
Perp_baseline_P1_rangeprofile   =  mean(Perp_baseline_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Perp_baseline_P1,"lin", savepath*"Geom_plots/", "Perp_baseline_P1", "")
plot_profile(s_loc_C, Perp_baseline_P1_rangeprofile[:], "", savepath*"Geom_plots/", "Perp_baseline_profile2", "C [m]", "Perp baseline  [m]", "", "", "")

Vert_wavnum_P1                = reshape(Vert_wavnum_all[:,1],length(s_loc_C),length(s_loc_S))'
Vert_wavnum_P1_rangeprofile   =  mean(Vert_wavnum_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Vert_wavnum_P1,"lin", savepath*"Geom_plots/", "Vert_wavnum_P1", "")
plot_profile(s_loc_C, Vert_wavnum_P1_rangeprofile[:], "", savepath*"Geom_plots/", "Vert_wavnum_profile2", "C [m]", "Vert wave num  ", "", "", "")

Correlation_theo_P1                = reshape(Correlation_theo_all[:,1],length(s_loc_C),length(s_loc_S))'
Correlation_theo_P1_rangeprofile   =  mean(Correlation_theo_P1,dims=1)

plot_image(s_loc_C,s_loc_S,Correlation_theo_P1,"lin", savepath*"Geom_plots/", "Correlation_theo_P1", "")
plot_profile(s_loc_C, Correlation_theo_P1_rangeprofile[:], "", savepath*"Geom_plots/", "Correlation_theo_profile2", "C [m]", "Correlation theory  ", "", "", "")


#Synthetic fringes
Unwrapped_data = ( (4 * pi )/params.λ).*Slantrange_Diff
Wrapped_data = mod.(Unwrapped_data,2*pi) .- pi

plot_image(s_loc_C,s_loc_S,Unwrapped_data,"phase", savepath*"Geom_plots/", "Slantrange_Diff_unwrapped", "")
plot_image(s_loc_C,s_loc_S,Wrapped_data,"phase", savepath*"Geom_plots/", "Slantrange_Diff_wrapped", "")

Wrapped_data_rangeprofile =  mean(Wrapped_data,dims=1)
plot_profile(s_loc_C, Wrapped_data_rangeprofile[:], "", savepath*"Geom_plots/", "Slantrange_profile_wrapped", "C [m]", "", "", "", "")
#plot_profile(s_loc_C[1000:1150], Wrapped_data_rangeprofile[1,1000:1150], "", savepath*"Geom_plots/", "Slantrange_profile_wrapped2", "C [m]", "", "", "", "")




#=
# To do processing again with master as reference
# Processing - 2
include("../modules/process_raw_data.jl")

SAR_images_3D_proc2 = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
image_3D_proc2      = Process_Raw_Data.tomo_processing_afterSAR_full(SAR_images_3D)


stat_var_all_1p_proc2 = SAR_images_3D_proc2[1,:,:,1]
stat_var_all_2p_proc2 = SAR_images_3D_proc2[2,:,:,1]

=#