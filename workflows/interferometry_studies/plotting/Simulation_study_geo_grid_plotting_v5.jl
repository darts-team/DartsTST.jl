
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

include("../../../modules/interferometry.jl")


function plot_image(x_axis,y_axis,Data,unit_flag, savepath, savename, figure_title)

    fontsize = 18

    if ~ispath(savepath)
        mkdir(savepath)
    end    

    if unit_flag == "log"
        plot_vals = 10 .* log10.(abs.(Data))
        llim = minimum(plot_vals[:])
        hlim = maximum(plot_vals[:])
        if (hlim-llim) > 30
            llim = hlim-30
        end
        p1=(heatmap(x_axis,y_axis,plot_vals,clim=(llim, hlim),xlabel="Lon [deg]",ylabel="Lat [deg]",title=figure_title,xticks=round.(LinRange((minimum(x_axis)),(maximum(x_axis)),4),digits=2), yticks=round.(LinRange((minimum(y_axis)),(maximum(y_axis)),4),digits=2),
        topmargin=6mm,bottommargin=16mm,leftmargin=10mm,rightmargin=15mm,tickfont=font(fontsize), xtickfont=font(fontsize), ytickfont=font(fontsize), guidefont=font(fontsize), titlefontsize=fontsize, size=(1200,500) )) #500,360
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "lin"
        p1=(heatmap(x_axis,y_axis,(abs.(Data)),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title, xticks=round.(LinRange((minimum(x_axis)),(maximum(x_axis)),4),digits=2), yticks=round.(LinRange((minimum(y_axis)),(maximum(y_axis)),4),digits=2),
        topmargin=6mm,bottommargin=16mm,leftmargin=10mm,rightmargin=15mm,tickfont=font(fontsize), xtickfont=font(fontsize), ytickfont=font(fontsize), guidefont=font(fontsize), titlefontsize=fontsize, size=(1200,500) )) #1200
        savefig(p1, savepath*savename*".png")
    elseif unit_flag == "phase"
        p1=(heatmap(x_axis,y_axis,(Data),xlabel="Lon [deg]",ylabel="Lat [deg]", title=figure_title, c=:twilight, xticks=round.(LinRange((minimum(x_axis)),(maximum(x_axis)),4),digits=2), yticks=round.(LinRange((minimum(y_axis)),(maximum(y_axis)),4),digits=2),
        topmargin=6mm,bottommargin=16mm,leftmargin=10mm,rightmargin=15mm,tickfont=font(fontsize), xtickfont=font(fontsize), ytickfont=font(fontsize), guidefont=font(fontsize), titlefontsize=fontsize, size=(1200,500) )) #1200
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
    plot_image(Lon_vals,Lat_vals,Slantrange_Diff[1:maxind_val_lon,1:maxind_val_lat]',"lin", savepath, "Slantrange_Diff", "")

    plot_profile(Lon_vals, Slantrange_P1_rangeprofile[:], Slantrange_P2_rangeprofile[:], savepath, "Slantrange_profile1", "Lon [deg]", "Slant range [m]", "", "Platform 1", "Platform 2")
    plot_profile(Lon_vals, Slantrange_Diff_rangeprofile[:], "", savepath, "Slantrange_profile", "Lon [deg]", "Slant range [m]", "", "Difference", "")

    #Synthetic fringes
    Unwrapped_data                      = ( (4 * pi )/0.23793052222222222).*Slantrange_Diff[1:maxind_val_lon,1:maxind_val_lat]
    Wrapped_data                        = mod.(Unwrapped_data[1:maxind_val_lon,1:maxind_val_lat],2*pi) #.- pi

    plot_image(Lon_vals,Lat_vals,Unwrapped_data',"lin", savepath, "Slantrange_Diff_unwrapped", "")
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


function get_geometry_variables(ref_data)

    #DEM
    DEM_data                            = ref_data[:,:,14]
    #BRCS
    targets_ref_data                    = ref_data[:,:,13]
    #Look angle
    Lookangle_P1                        = ref_data[:,:,3]
    Lookangle_P2                        = ref_data[:,:,4]
    #Incidence angle
    Incangle_P1                         = ref_data[:,:,5]
    Incangle_P2                         = ref_data[:,:,6]
    #Local Incidence angle
    LocalIncangle                       = ref_data[:,:,9]
    #Range slope angle
    range_slope                         = ref_data[:,:,10]
    #Slant range
    Slantrange_P1                       = ref_data[:,:,1]
    Slantrange_P2                       = ref_data[:,:,2]
    # Critical baseline
    Critical_baseline                   = ref_data[:,:,11]
    # Perp baseline
    Perp_baseline                       = ref_data[:,:,7]
    # vertical wavenumber
    Vert_wavnum                         = ref_data[:,:,8]
    # Theoretical correlation 
    Correlation_theo                    = ref_data[:,:,12]

    return DEM_data, targets_ref_data, Lookangle_P1, Lookangle_P2, Incangle_P1, Incangle_P2, LocalIncangle, range_slope, Slantrange_P1, Slantrange_P2, Critical_baseline, Perp_baseline, Vert_wavnum, Correlation_theo

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

folder_path= "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set2/"
folder_index = 202

s_geom_filepath             = folder_path*string(folder_index)*"/"*string(folder_index)*"_scene_geometry_1.tif" 
t_geom_filepath             = folder_path*string(folder_index)*"/"*string(folder_index)*"_target_geometry_1.tif" 

opdata_filepath_s           = folder_path*string(folder_index)*"/"*string(folder_index)*"_sim_output_main_scene_1.tif" 
opdata_filepath_t           = folder_path*string(folder_index)*"/"*string(folder_index)*"_sim_output_main_target_1.tif" 

savepath                    = folder_path*string(folder_index)*"/"


#t_geom_filepath         = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/106/106_target_geometry_1.tif" 

#opdata_filepath_s       = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/106/106_sim_output_main_scene_1.tif" 
#opdata_filepath_t       = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/106/106_sim_output_main_target_1.tif" 

#savepath                = "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/106/"

process_plot_output_flag      = "T"

profile_flag            = 0

Looks_along_Lat         = 4
Looks_along_Lon         = 12

oversampling_factor_looks = 2.78*1.8


if process_plot_output_flag == "S"

    savepath                 = savepath*"Plots_s/"
    if ~ispath(savepath)
        mkdir(savepath)
    end    

    s_geom_data                     = reverse(GeoArrays.read(s_geom_filepath),dims=2)
    t_geom_data                     = reverse(GeoArrays.read(t_geom_filepath),dims=2)

    #s_geom_data                     = (GeoArrays.read(s_geom_filepath))
    #t_geom_data                     = (GeoArrays.read(t_geom_filepath))

    coords_op                       = collect(GeoArrays.coords(s_geom_data))

    Lon_vals_all                    = [x[1] for x in coords_op]
    Lon_vals                        = (Lon_vals_all[:,1])

    Lat_vals_all                    = [x[2] for x in coords_op]
    Lat_vals                        = reverse(Lat_vals_all[1,:])
    #Lat_vals                        = (Lat_vals_all[1,:])


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

    t_geom_data                     = reverse(GeoArrays.read(t_geom_filepath),dims=2)
    #t_geom_data                     = (GeoArrays.read(t_geom_filepath))

    coords_op                       = collect(GeoArrays.coords(t_geom_data))

    Lon_vals_all                    = [x[1] for x in coords_op]
    Lon_vals                        = Lon_vals_all[:,1]

    Lat_vals_all                    = [x[2] for x in coords_op]
    #Lat_vals                        = (Lat_vals_all[1,:])
    Lat_vals                        = reverse(Lat_vals_all[1,:])


    maxind_val_lon                  = get_max_ind_vec(length(Lon_vals),Looks_along_Lon)
    maxind_val_lat                  = get_max_ind_vec(length(Lat_vals),Looks_along_Lat)

    Lon_vals                        = Lon_vals[1:maxind_val_lon]
    Lat_vals                        = Lat_vals[1:maxind_val_lat]

    geom_savepath                   = savepath*"Geom_plots_t/"
    plot_geometry_variables(t_geom_data, Lon_vals, Lat_vals, maxind_val_lon, maxind_val_lat, profile_flag, geom_savepath)

    data_op                         = ArchGDAL.read(opdata_filepath_t)

    DEM_data, targets_ref_data, Lookangle_P1, Lookangle_P2, Incangle_P1, Incangle_P2, LocalIncangle, range_slope, Slantrange_P1, Slantrange_P2, Critical_baseline, Perp_baseline, Vert_wavnum, Correlation_theo = get_geometry_variables(t_geom_data)

end


Data_1                          = ArchGDAL.getband(data_op, 1)
Data_2                          = ArchGDAL.getband(data_op, 2)

Data_1                          = reverse(Data_1[1:maxind_val_lon,1:maxind_val_lat], dims=2)
Data_2                          = reverse(Data_2[1:maxind_val_lon,1:maxind_val_lat], dims=2)


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
A, B, C, D                      = Interferometry.sar_geogrid_amplitude_statistics(Amp_1[:], savepath*"Amplitude_plots/", "P1_amp", 1)
A, B, C, D                      = Interferometry.sar_geogrid_amplitude_statistics(Amp_2[:], savepath*"Amplitude_plots/", "P2_amp", 1)
# Power statistics plots 
A, B, C, D                      = Interferometry.sar_geogrid_power_statistics(Pow_1[:], savepath*"Power_plots/", "P1_pow", 1)
A, B, C, D                      = Interferometry.sar_geogrid_power_statistics(Pow_2[:], savepath*"Power_plots/", "P2_pow", 1)
# Phase statistics plots 
A, B, C, D                      = Interferometry.sar_geogrid_phase_statistics(Phase_1[:], savepath*"Phase_plots/", "P1_phase", 1)
A, B, C, D                      = Interferometry.sar_geogrid_phase_statistics(Phase_2[:], savepath*"Phase_plots/", "P2_phase", 1)

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

complex_coherence_mat = zeros(ComplexF64,Lon_length_ml,Lat_length_ml)

k=1
for i=1:Int((length(Lat_vals))/Looks_along_Lat)
    l=1
    for j=1:Int((length(Lon_vals))/Looks_along_Lon)
        complex_coherence = mean(Data_12[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1])./ ((mean(Pow_1[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1]) .* mean(Pow_2[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1])).^0.5)
        complex_coherence_mat[j,i] = complex_coherence
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
        
plot_image(Lon_vals_multilooked,Lat_vals_multilooked,Int_Phase_multilooked',"phase", savepath*"Int_plots/", "ML_Phase_P12_lin", "")
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


A, B, C, D = Interferometry.interferogram_statistics_magnitude(Int_Pow_multilooked[:], gamma_ip, (Looks_along_Lat*Looks_along_Lon)/oversampling_factor_looks, (0.5,1), savepath*"Int_plots/", "PDF_Mag_Han1", 1)
A, B, C, D = Interferometry.interferogram_statistics_phase(Int_Phase_multilooked[:], gamma_ip, (Looks_along_Lat*Looks_along_Lon)/1, 0.0, (-0.5,0.5), savepath*"Int_plots/", "PDF_Phase_Han1", 1)
A, B, C, D = Interferometry.interferogram_statistics_phase_method2(Int_Phase_multilooked[:], gamma_ip, (Looks_along_Lat*Looks_along_Lon)/1, 0.0, (-0.5,0.5), savepath*"Int_plots/", "PDF_Phase_Han2", 1)


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


# Phase unwrapping and Height estimation
#=
using JLD2
#@save "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/50/test_data.jld" complex_coherence_mat Int_Phase_multilooked Int_Pow_multilooked

#@save "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/50/test_data2.jld" complex_coherence_mat Int_Phase_multilooked Int_Pow_multilooked


@load "/Users/joshil/Documents/Outputs/InSAR Outputs/Geo_outputs_set1/90/test_data.jld" 

plot_image(Lon_vals_multilooked,Lat_vals_multilooked,conncomp',"phase", savepath*"Int_plots/", "SNAPHU_1", "")

plot_image(Lon_vals_multilooked,Lat_vals_multilooked,(unw' .* -1) ,"phase", savepath*"Int_plots/", "SNAPHU_2", "")


int_unwrapped_height_new2      = (unw .* -1) .* (0.23793052222222222/(4*pi)) .* (slant_range_profile[:,:,1] .*sind.(look_angle_profile[:,:,1]) ./perp_baseline_profile[:,:,1])
plot_image(Lon_vals_multilooked,Lat_vals_multilooked,int_unwrapped_height_new2',"phase", savepath*"Int_plots/", "SNAPHU_3", "")

plot_image(Lon_vals_multilooked,Lat_vals_multilooked,int_unwrapped_height_new2'.+2000,"lin", savepath*"Int_plots/", "SNAPHU_4", "")
=#


#=
if process_plot_output_flag == "T"

    include("../../../modules/scattering.jl")

    LocalIncangle_ml             = zeros(Lon_length_ml,Lat_length_ml)
    targets_ref_data_ml             = zeros(Lon_length_ml,Lat_length_ml)

    k=1
    for i=1:Int((length(Lat_vals))/Looks_along_Lat)
        l=1
        for j=1:Int((length(Lon_vals))/Looks_along_Lon)
            LocalIncangle_ml[j,i] = mean(LocalIncangle[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1])
            targets_ref_data_ml[j,i] = mean(targets_ref_data[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1])
            l=l+Looks_along_Lon
        end
        k=k+Looks_along_Lat
    end



    rcs_model_ip = collect(0:60)
    rcs_model = zeros(61)
    for i=1:61
        rcs_model[i] = Scattering.TST_surface_brcs(2,0.23793052222222222,rcs_model_ip[i],0.0,rcs_model_ip[i],180.0-0.0,3,1.0)
    end
    rcs_constant  = 4.5
    rcs_model = (rcs_model .* rcs_constant)

    plot_var = (targets_ref_data[:])
    
    (histogram2d(LocalIncangle[:],plot_var.-maximum(plot_var),colorbar_scale=:log10,xlabel="Local incidence angle",ylabel="RCS", 
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(700,500), )) #1200
    
    display(plot!(rcs_model_ip,rcs_model.-maximum(rcs_model), lw=2, color=:black, label="Model",))
    
    pow_var = 10 .* log10.((Pow_1_multilooked)')
    pow_var = pow_var .- maximum(pow_var[:])
    la_var = LocalIncangle_ml'

    display(histogram2d(la_var[:],(pow_var[:]), colorbar_scale=:log10,xlabel="Local incidence angle",ylabel="Pr",
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(700,500))) 





    Amp_1_multilooked               = zeros(Lon_length_ml,Lat_length_ml)
    Amp_2_multilooked               = zeros(Lon_length_ml,Lat_length_ml)
    Correlation_theo_multilooked    = zeros(Lon_length_ml,Lat_length_ml)
    
    norm_mean_amp1                  = maximum(Amp_1)
    norm_mean_amp2                  = maximum(Amp_1)
    
    k=1
    for i=1:Lat_length_ml
        l=1
        for j=1:Lon_length_ml
            Amp_1_multilooked[j,i]  = mean(Amp_1[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1]) ./norm_mean_amp1
            Amp_2_multilooked[j,i]  = mean(Amp_2[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1]) ./norm_mean_amp2
            Correlation_theo_multilooked[j,i]  = mean(Correlation_theo[l:l+Looks_along_Lon-1,k:k+Looks_along_Lat-1]) 

            l=l+Looks_along_Lon
        end
        k=k+Looks_along_Lat
    end



    pow_var = 10 .* log10.(Amp_1_multilooked')
    pow_var = pow_var .- maximum(pow_var[:])
    la_var = LocalIncangle_ml'

    display(histogram2d(la_var[:],pow_var[:], colorbar_scale=:log10,xlabel="Local incidence angle",ylabel="Amp",
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(700,500))) 



    Correlation_theo[:] = clamp.(Correlation_theo[:],0,1)
    display(histogram2d(LocalIncangle[:],Correlation_theo[:],xlabel="Local incidence angle",ylabel="Theo Corr", 
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(700,500),xlim=(10,40),ylim=(0.8,1) )) #1200



    display(histogram2d(LocalIncangle_ml[:],Int_Pow_multilooked[:],xlabel="Local incidence angle",ylabel="Corr", 
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(700,500),xlim=(10,40),ylim=(0.8,1) )) #1200


    display(histogram2d(Correlation_theo_multilooked[:],Int_Pow_multilooked[:],xlabel="Corr Theory",ylabel="Corr Sim", 
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(700,500), xlim=(0.9,1),ylim=(0.9,1))) #1200

    display(plot!(0:0.01:1,0:0.01:1,lw=2,color=:red, legend=:false))




    A_vals = 10 .* log10.((targets_ref_data[1:maxind_val_lon,1:maxind_val_lat].^2)').+128

    B_vals = 10 .* log10.(Pow_1')

    (histogram2d(A_vals[:],B_vals[:],colorbar_scale=:log10,xlabel="RCS input adjusted [dB]",ylabel="Power [dB]", 
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(700,500), xlim=(40,160),ylim=(40,160) )) #1200
    
    display(plot!(40:160,40:160,lw=2,color=:red, legend=:false))



    A_vals = 10 .* log10.(targets_ref_data[1:maxind_val_lon,1:maxind_val_lat]').+64

    B_vals = 10 .* log10.(Amp_1')

    (histogram2d(A_vals[:],B_vals[:],colorbar_scale=:log10,xlabel="RCS input adjusted [dB]",ylabel="Amplitude [dB]", 
    topmargin=6mm,bottommargin=10mm,leftmargin=6mm,rightmargin=6mm,tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(700,500), xlim=(20,80),ylim=(20,80) )) #1200
    
    display(plot!(20:80,20:80,lw=2,color=:red, legend=:false))

end

=#