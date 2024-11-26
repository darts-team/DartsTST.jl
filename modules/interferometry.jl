module Interferometry

include("../modules/geometry.jl")
include("../modules/scene.jl")
include("../modules/orbits.jl")
include("../modules/data_processing.jl")

using Statistics
using Parameters
using Plots
using Measures
using Statistics
using NaNStatistics
using SpecialFunctions
using HypergeometricFunctions

c               = 299792458
earth_radius    = 6378.137e3 # Earth semi-major axis at equator

function get_scene_geometry_values(p_xyz, v_xyz, s_xyz_3xN, N_all, ref_plat, sec_plat, p_mode, params, grid )
    @unpack λ, mode, bandwidth, left_right_look = params

    #geometry computations based on scene
    slant_range_ref                 = zeros(size(s_xyz_3xN,2),1)
    look_angle_ref                  = zeros(size(s_xyz_3xN,2),1)
    incidence_angle_ref             = zeros(size(s_xyz_3xN,2),1)
    slant_range_sec                 = zeros(size(s_xyz_3xN,2),1)
    look_angle_sec                  = zeros(size(s_xyz_3xN,2),1)
    incidence_angle_sec             = zeros(size(s_xyz_3xN,2),1)
    Perp_baseline_ref               = zeros(size(s_xyz_3xN,2),1)
    Vert_wavnum_ref                 = zeros(size(s_xyz_3xN,2),1)
    local_incidence_angle_ref       = zeros(size(s_xyz_3xN,2),1)
    range_slope_angle_ref           = zeros(size(s_xyz_3xN,2),1)
    Critical_baseline_ref           = zeros(size(s_xyz_3xN,2),1)
    Correlation_theo_ref            = zeros(size(s_xyz_3xN,2),1)

    mean_plats_pos_ref              = mean(p_xyz[:,ref_plat,:], dims=2)
    mean_plats_pos_sec              = mean(p_xyz[:,sec_plat,:], dims=2)

    for ti = 1:size(s_xyz_3xN,2)
        slant_range_ref[ti]         = Geometry.distance( mean_plats_pos_ref  , s_xyz_3xN[:,ti] )
        look_angle_ref[ti]          = Scene.slantrange_to_lookangle(earth_radius,slant_range_ref[ti],Geometry.xyz_to_geo(mean_plats_pos_ref)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3])[2]
        incidence_angle_ref[ti]     = Scene.lookangle_to_incangle(look_angle_ref[ti],Geometry.xyz_to_geo(mean_plats_pos_ref)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3],earth_radius)

        slant_range_sec[ti]         = Geometry.distance( mean_plats_pos_sec  , s_xyz_3xN[:,ti] )
        look_angle_sec[ti]          = Scene.slantrange_to_lookangle(earth_radius,slant_range_sec[ti],Geometry.xyz_to_geo(mean_plats_pos_sec)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3])[2]
        incidence_angle_sec[ti]     = Scene.lookangle_to_incangle(look_angle_sec[ti],Geometry.xyz_to_geo(mean_plats_pos_sec)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3],earth_radius)

        bs_perp, bs_at, bs_norm     = Orbits.get_perp_baselines_new(mean(p_xyz[:,:,:],dims=3), mean(v_xyz[:,:,:],dims=3), look_angle_ref[ti], 0.0, left_right_look, 1)
        Perp_baseline_ref[ti]       = bs_perp[1,sec_plat,1]

        Va                          = mean(p_xyz[:, ref_plat, :], dims=2) - s_xyz_3xN[:, ti]
        Vb                          = mean(p_xyz[:, sec_plat, :], dims=2) - s_xyz_3xN[:, ti]
        angle_ip                    = Data_Processing.angle_2vec(Va, Vb) * 1

        if mode == 1
            Vert_wavnum_ref[ti]     = (4 * pi * (angle_ip * pi / 180) ) / (λ * sind(look_angle_ref[ti]))
        elseif mode == 2
            Vert_wavnum_ref[ti]     = (2 * pi * (angle_ip * pi / 180) ) / (λ * sind(look_angle_ref[ti]))
        end

        #plat_pt_xyz                 =  mean(mean(p_xyz,dims=2),dims=3)[:] #??????
        plat_pt_xyz                 =  mean(p_xyz[:,ref_plat,:],dims=2)
        look_vec_xyz                = (plat_pt_xyz - s_xyz_3xN[:,ti])
        look_vec_xyz_norm           = (plat_pt_xyz - s_xyz_3xN[:,ti]) / Geometry.distance( plat_pt_xyz,s_xyz_3xN[:,ti])

        Geo_location                = Geometry.xyz_to_geo(s_xyz_3xN[:,ti])
        pegθ                        = Geo_location[1]*π/180
        pegϕ                        = Geo_location[2]*π/180
        #ENU to XYZ transformation matrix
        Menu_xyz                    = [-sin(pegϕ) -sin(pegθ)*cos(pegϕ) cos(pegθ)*cos(pegϕ);
                                    cos(pegϕ) -sin(pegθ)*sin(pegϕ) cos(pegθ)*sin(pegϕ);
                                    0            cos(pegθ)             sin(pegθ)]
        #XYZ to ENU transformation matrix
        Mxyz_enu                    = [-sin(pegϕ)           cos(pegϕ)             0;
                                    -sin(pegθ)*cos(pegϕ) -sin(pegθ)*sin(pegϕ)  cos(pegθ)  ;
                                    cos(pegθ)*cos(pegϕ)   cos(pegθ)*sin(pegϕ)  sin(pegθ)]
                   
        look_vec_enu                = Mxyz_enu * look_vec_xyz
        look_direction_norm         = sqrt( look_vec_enu[1] * look_vec_enu[1]  + look_vec_enu[2] * look_vec_enu[2])

        if grid == "Flat"
            #N_all[ti,:]             = [0,0,1]
            Nxyz                        = Menu_xyz * [0,0,1];
            local_incidence_angle_ref[ti] = Data_Processing.angle_2vec(look_vec_xyz_norm, Nxyz)
            range_slope_angle_ref[ti]   = atand((0 * (look_vec_enu[1] / look_direction_norm)) + (0 * (look_vec_enu[2] / look_direction_norm)))    
        else
            Nxyz                        = Menu_xyz * N_all[ti,:];
            local_incidence_angle_ref[ti] = Data_Processing.angle_2vec(look_vec_xyz_norm, Nxyz)
            range_slope_angle_ref[ti]   = atand((N_all[ti,1] * (look_vec_enu[1] / look_direction_norm)) + (N_all[ti,2] * (look_vec_enu[2] / look_direction_norm)))          
        end      
        
        #Critical_baseline_ref[ti]   = λ * ((2*bandwidth)/c) * slant_range_ref[ti] * tand(local_incidence_angle_ref[ti] - range_slope_angle_ref[ti]) / p_mode
        Critical_baseline_ref[ti]   = λ * ((2*bandwidth)/c) * slant_range_ref[ti] * tand(local_incidence_angle_ref[ti]) / p_mode

        Correlation_theo_ref[ti]    = 1 - (Perp_baseline_ref[ti] ./ (Critical_baseline_ref[ti]))

    end
    return slant_range_ref, look_angle_ref, incidence_angle_ref, slant_range_sec, look_angle_sec, incidence_angle_sec, Perp_baseline_ref, Vert_wavnum_ref, local_incidence_angle_ref, range_slope_angle_ref, Critical_baseline_ref, Correlation_theo_ref
end

##  Function to get SAR amplitude statistics 
function sar_geogrid_amplitude_statistics(ip_var, savepath, savename, plot_flag)

    # Check if directory is present
    if ~ispath(savepath)
        mkdir(savepath)
    end    

    # Get statistics
    mean_var                = mean(ip_var)
    median_var              = mean(ip_var)
    std_var                 = mean(ip_var)

    # Get histogram
    ip_var                  = ip_var ./ mean_var
    bins_range              = range(0,maximum(ip_var), length=51)

    # Plot histogram
    if plot_flag == 1
        p1 = (histogram(ip_var,xlabel="Amplitude",ylabel="Count",title = "Histogram", colorbar=true, bins= bins_range,
            tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
        savefig(p1, savepath*savename*"_Hist_plot1"*".png")
    end
    
    # get PDF
    N, bin                  = histcountindices(ip_var,bins_range)
    pdf_var                 = (N / sum(N)) / (bins_range[2]-bins_range[1])
    ref_var_vals            = range(0+(maximum(ip_var)/(length(bins_range)-1)),maximum(ip_var)-(maximum(ip_var)/(length(bins_range)-1)), length=(length(bins_range)-1))

    # Compute theoretical values
    theory_val             = zeros(length(ref_var_vals),1)
    sigma                  = 1 * sqrt(2/pi)
    for i =1:length(ref_var_vals)
        #theory_val[i]      = 2 .* pdf(Chisq(2),ref_var_vals[i])
        theory_val[i]      = (ref_var_vals[i] * exp((-ref_var_vals[i]^2)/(2*(sigma^2)))) / (sigma^2)
    end

    # Plot comparison data vs theory
    if plot_flag == 1
        plot(ref_var_vals./sigma , sigma .*pdf_var, ylim=(0,1.1), xlim=(0,maximum(ip_var)),xlabel="Amplitude / sigma",ylabel="sigma * PDF",title = ""
            ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3)
        p3 = (plot!(ref_var_vals./sigma, sigma .*theory_val, ylim=(0,1.1), xlim=(0,maximum(ip_var)),xlabel="Amplitude / sigma",ylabel="sigma * PDF",title = ""
            ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
        savefig(p3, savepath*savename*"_Hist_plot2"*".png")
    end

    data_stat = [mean_var; median_var; std_var]

    return data_stat, ref_var_vals./sigma, sigma.*pdf_var, sigma.*theory_val

end

##  Function to get SAR power statistics 
function sar_geogrid_power_statistics(ip_var, savepath, savename, plot_flag)

    # Check if directory is present
    if ~ispath(savepath)
        mkdir(savepath)
    end    

    # Get statistics
    mean_var                = mean(ip_var)
    median_var              = mean(ip_var)
    std_var                 = mean(ip_var)

    # Get histogram
    ip_var                  = ip_var ./ mean_var
    bins_range              = range(0,maximum(ip_var), length=51)

    # Plot histogram
    if plot_flag == 1
        p1 = (histogram(ip_var,xlabel="Power",ylabel="Count",title = "Histogram", colorbar=true, bins= bins_range,
            tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
        savefig(p1, savepath*savename*"_Hist_plot1"*".png")
    end
    
    # get PDF
    N, bin                  = histcountindices(ip_var,bins_range)
    pdf_var                 = (N / sum(N)) / (bins_range[2]-bins_range[1])
    ref_var_vals            = range(0+(maximum(ip_var)/(length(bins_range)-1)),maximum(ip_var)-(maximum(ip_var)/(length(bins_range)-1)), length=(length(bins_range)-1))

    # Compute theoretical values
    theory_val              = zeros(length(ref_var_vals),1)
    N                       = 1
    for i =1:length(ref_var_vals)
        #theory_val[i]       = 2 .* pdf(Chisq(2),ref_var_vals[i])
        theory_val[i]       = ((ref_var_vals[i]^(N-1)) * exp((-ref_var_vals[i])/(1))) / ((1^N) * gamma(1))
    end

    # Plot comparison data vs theory
    if plot_flag == 1
        plot(ref_var_vals , pdf_var, ylim=(0,1.1), xlim=(0,maximum(ip_var)),xlabel="Power ",ylabel="PDF",title = ""
            ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3)
        p3 = (plot!(ref_var_vals, theory_val, ylim=(0,1.1), xlim=(0,maximum(ip_var)),xlabel="Power ",ylabel="PDF",title = ""
            ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
        savefig(p3, savepath*savename*"_Hist_plot2"*".png")
    end

    data_stat = [mean_var; median_var; std_var]

    return data_stat, ref_var_vals, pdf_var, theory_val

end

##  Function to get SAR phase statistics 
function sar_geogrid_phase_statistics(ip_var, savepath, savename, plot_flag)

    # Check if directory is present
    if ~ispath(savepath)
        mkdir(savepath)
    end    

    # Get statistics
    mean_var                = mean(ip_var)
    median_var              = mean(ip_var)
    std_var                 = mean(ip_var)

    # Get histogram
    bins_range              = range(-pi,pi, length=51)

    # Plot histogram
    if plot_flag == 1
        p1 = (histogram(ip_var,xlabel="Phase",ylabel="Count",title = "Histogram", colorbar=true, bins= bins_range,
            tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
        savefig(p1, savepath*savename*"_Hist_plot1"*".png")
    end

    N, bin                  = histcountindices(ip_var,bins_range)
    pdf_var                 = (N / sum(N)) / (bins_range[2]-bins_range[1])
    ref_var_vals            = range(-pi+(pi/(length(bins_range)-1)),pi-(pi/(length(bins_range)-1)), length=(length(bins_range)-1))
    theory_val              = repeat([1/(2*pi)], size(ref_var_vals)[1],1) 

    # Plot comparison data vs theory
    if plot_flag == 1
        plot(ref_var_vals , pdf_var, ylim=(0,1/pi), xlim=(-1.5*pi,1.5*pi),xlabel="Phase ",ylabel="PDF",title = ""
            ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3)
        p3 = (plot!(ref_var_vals, theory_val, ylim=(0,1/pi), xlim=(-1.5*pi,1.5*pi),xlabel="Phase ",ylabel="PDF",title = ""
            ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
        savefig(p3, savepath*savename*"_Hist_plot2"*".png")
    end

    data_stat = [mean_var; median_var; std_var]

    return data_stat, ref_var_vals, pdf_var, theory_val

end

##  Function to get InSAR magnitude statistics 
function interferogram_statistics_magnitude(ip_var, gamma_ip, L, xlims, savepath, savename, plot_flag)

    # Check if directory is present
    if ~ispath(savepath)
        mkdir(savepath)
    end    
    
    # Get statistics
    mean_var                = mean(ip_var)
    median_var              = mean(ip_var)
    std_var                 = mean(ip_var)
    
    bins_range              = range(0,1, length=251)

    #Phase histogram plot
    if plot_flag == 1
        p2 = (histogram(ip_var,xlabel="Magnitude",ylabel="Count",title = "Histogram", colorbar=true, bins=bins_range, xlim=xlims,
        tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
        savefig(p2, savepath*savename*"_Hist_plot1_Mag"*".png")
    end

    N, bin                  = histcountindices(ip_var,bins_range)
    pdf_var                 = (N / sum(N)) / (bins_range[2]-bins_range[1])
    ref_var_vals            = range(0+(1/(length(bins_range)-1)),1-(1/(length(bins_range)-1)), length=(length(bins_range)-1))

    #Hanssen - magnitude
    theory_val       = zeros( length(bins_range))
    for gamm_idx=1:length(bins_range)  
        theory_val[gamm_idx] = 2 .* (L-1) .* ((1-(gamma_ip^2))^L) .* bins_range[gamm_idx] .* ((1-(bins_range[gamm_idx]^2))^(L-2)) * _₂F₁(L,L,1,(bins_range[gamm_idx]^2) * (gamma_ip^2) )   
    end

    if plot_flag == 1    
        (plot(ref_var_vals, pdf_var, xlim=xlims,xlabel="Magnitude",ylabel="PDF",title = ""
            ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
        p3=(plot!(bins_range, theory_val, xlim=xlims,xlabel="Magnitude",ylabel="PDF",title = ""
            ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm,legend=:topleft))
        savefig(p3, savepath*savename*"_PDF_comp_IntMag"*".png")
    end

    data_stat = [mean_var; median_var; std_var]

    return data_stat, ref_var_vals, pdf_var, theory_val

end

##  Function to get InSAR phase statistics - method 1
function interferogram_statistics_phase(ip_var, gamma_ip, L, phase_ip_ref, xlims, savepath, savename, plot_flag)
    
    # Check if directory is present
    if ~ispath(savepath)
        mkdir(savepath)
    end    
    
    # Get statistics
    mean_var                = mean(ip_var)
    median_var              = mean(ip_var)
    std_var                 = mean(ip_var)
    bins_range              = range(-pi,pi, length=501)

    #Phase histogram plot
    if plot_flag == 1   
        p2 = (histogram(ip_var,xlabel="Phase",ylabel="Count",title = "Histogram", colorbar=true, bins=bins_range, xlim=xlims,
        tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
        savefig(p2, savepath*savename*"_Hist_plot1_Phase"*".png")
    end

    N, bin                  = histcountindices(ip_var,bins_range)
    pdf_var                 = (N / sum(N)) / (bins_range[2]-bins_range[1])
    ref_var_vals            = range(-pi+(pi/(length(bins_range)-1)),pi-(pi/(length(bins_range)-1)), length=(length(bins_range)-1))

    theory_val              = zeros( length(bins_range))

    #Hanssen - phase - method 1
    for phase_ip_idx=1:length(bins_range)
        beta_ip         = gamma_ip * cos(bins_range[phase_ip_idx] - phase_ip_ref)
        beta_sqrt_t     = (1-(beta_ip^2))
        pdf_ph_t1       = ((gamma(2*L-1)) / ((gamma(L)^2) * (2^(2*(L-1))))) * ( ((((2*L-1) * beta_ip) / (beta_sqrt_t^(L+0.5))) * ( acos(-beta_ip))) + (1/beta_sqrt_t^L))
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

    if plot_flag == 1       
        (plot(ref_var_vals, pdf_var, xlim=xlims,xlabel="Phase",ylabel="PDF",title = ""
            ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
        p3  =(plot!(bins_range, theory_val, xlim=xlims,xlabel="Phase",ylabel="PDF",title = ""
            ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
        savefig(p3, savepath*savename*"_PDF_comp_IntPhase"*".png")
    end

    data_stat = [mean_var; median_var; std_var]

    return data_stat, ref_var_vals, pdf_var, theory_val

end


##  Function to get InSAR phase statistics - method 2
function interferogram_statistics_phase_method2(ip_var, gamma_ip, L, phase_ip_ref, xlims, savepath, savename, plot_flag)
    
    # Check if directory is present
    if ~ispath(savepath)
        mkdir(savepath)
    end    
    
    # Get statistics
    mean_var                = mean(ip_var)
    median_var              = mean(ip_var)
    std_var                 = mean(ip_var)
    bins_range              = range(-pi,pi, length=501)

    #Phase histogram plot
    if plot_flag == 1   
        p2 = (histogram(ip_var,xlabel="Phase",ylabel="Count",title = "Histogram", colorbar=true, bins=bins_range, xlim=xlims,
        tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, legend=false,right_margin=10mm))
        savefig(p2, savepath*savename*"_Hist_plot1_Phase"*".png")
    end

    N, bin                  = histcountindices(ip_var,bins_range)
    pdf_var                 = (N / sum(N)) / (bins_range[2]-bins_range[1])
    ref_var_vals            = range(-pi+(pi/(length(bins_range)-1)),pi-(pi/(length(bins_range)-1)), length=(length(bins_range)-1))

    theory_val              = zeros( length(bins_range))

    #Hanssen - phase - method 2
    for phase_ip_idx=1:length(bins_range)
        beta_ip             = gamma_ip * cos(bins_range[phase_ip_idx] - phase_ip_ref)
        beta_sqrt_t         = (1-(beta_ip^2))

        Term1               = (gamma(L + 0.5) * ((1-(gamma_ip^2))^L) * abs.(gamma_ip) *  cos(bins_range[phase_ip_idx] - phase_ip_ref)) / ( 2 * sqrt(pi) * gamma(L) * (beta_sqrt_t^(L+0.5)) )
        Term2               = ((1-(gamma_ip^2))^L) / (2*pi) * _₂F₁(L,1,0.5,beta_ip^2)

        theory_val[phase_ip_idx] = Term1 + Term2

    end

    if plot_flag == 1       
        (plot(ref_var_vals, pdf_var, xlim=xlims,xlabel="Phase",ylabel="PDF",title = ""
            ,label="Data", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3))
        p3  =(plot!(bins_range, theory_val, xlim=xlims,xlabel="Phase",ylabel="PDF",title = ""
            ,label="Theory", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=3,right_margin=10mm))
        savefig(p3, savepath*savename*"_PDF_comp_IntPhase"*".png")
    end

    data_stat = [mean_var; median_var; std_var]

    return data_stat, ref_var_vals, pdf_var, theory_val

end

# Function to unwrap phase based on snaphu python code
function unwrap_phase_snaphu(complex_coherence_mat, Int_Pow_multilooked, NL)

    using PyCall
    shu = pyimport("snaphu")
    unw, conncomp = shu.unwrap(complex_coherence_mat, Int_Pow_multilooked, nlooks=NL, cost="smooth", init="mcf")

    return unw, conncomp 

end

# Function to construct heights from unwrapped phase 
function get_height_from_phase(unwrapped_phase, slant_range, local_inc_angle, perp_baseline, lambda)

    height_est      = (unwrapped_phase .* -1) .* (lambda/(4*pi)) .* (slant_range .*sind.(local_inc_angle) ./ perp_baseline)

    return height_est

end

# Function to multi-look 1D dataset (2nd dimension is the variable list)
function multilook_1D_data(input_data, Looks)

    ip_size                         = size(input_data)
    if length(ip_size) == 2
        data_size2 = ip_size[2]
    else
        data_size2 = 1
    end

    length_ml                       = Int(ip_size[1]/Looks)

    output_data                     = zeros(typeof(input_data[1]),length_ml)

    k=1
    for i=1:length_ml
            for y=1:data_size2
                output_data[i,y] = mean(input_data[k:k+Looks-1,y])
            end

        k=k+Looks
    end

    return output_data
end

# Function to multi-look 2D dataset (3rd dimension is the variable list)
function multilook_2D_data(input_data, Looks_along_X, Looks_along_Y )

    ip_size                         = size(input_data)
    if length(ip_size) == 3
        data_size3 = ip_size[3]
    else
        data_size3 = 1
    end

    X_length_ml                     = Int(ip_size[1]/Looks_along_X)
    Y_length_ml                     = Int(ip_size[2]/Looks_along_Y)

    output_data                     = zeros(typeof(input_data[1]),X_length_ml,Y_length_ml,data_size3)

    k=1
    for i=1:Y_length_ml
        l=1
        for j=1:X_length_ml
            for y=1:data_size3
                output_data[j,i,y] = mean(input_data[l:l+Looks_along_X-1,k:k+Looks_along_Y-1,y])
            end
            l=l+Looks_along_X
        end
        k=k+Looks_along_Y
    end

    if data_size3 == 1
        return output_data[:,:,1]
    else
        return output_data
    end
end



end