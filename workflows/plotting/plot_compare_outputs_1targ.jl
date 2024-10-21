
#Plotting op file
using GeoDatasets
using JLD2
using Plots

@load "../Outputs/outputs_geometry_single_target_configuration_4.jld"


img_path = "/Users/joshil/Documents/Code/Plots_4p/outputs_geometry_single_target_configuration_4/"

#xlim_plot = (-125,-75)
#ylim_plot = (15,50)

xlim_plot = (-180,180)
ylim_plot = (-75,75)

#xlim_plot = (95,140)
#ylim_plot = (-10,15)

Lats_p                  = Geo_location.A[1,:,1]
Lons_p                  = Geo_location.A[:,1,2]

plot_v_la               =  lookang_all[:,:,1]
plot_v_bprpmax          =  Perp_baseline_max[:,:,1]
plot_v_bprpmin          =  Perp_baseline_min[:,:,1]
plot_v_bparmax          =  Par_baseline_max[:,:,1]
plot_v_bparmin          =  Par_baseline_min[:,:,1]
plot_v_bnrmmax          =  Norm_baseline_max[:,:,1]
plot_v_bnrmmin          =  Norm_baseline_min[:,:,1]

theo_res_N_1              =  res_theory_n[:,:,1]
theo_res_S_1              =  res_theory_s[:,:,1]
amb_Ht_1                  =  amb_H[:,:,1]
amb_Nt_1                  =  amb_N[:,:,1]
Slnt_range_1              = slnt_range[:,:,1]./1000

theo_res_N_2              =  res_theory_n[:,:,2]
theo_res_S_2              =  res_theory_s[:,:,2]
amb_Ht_2                  =  amb_H[:,:,2]
amb_Nt_2                  =  amb_N[:,:,2]
Slnt_range_2              = slnt_range[:,:,2]./1000

theo_res_N_3              =  res_theory_n[:,:,3]
theo_res_S_3              =  res_theory_s[:,:,3]
amb_Ht_3                  =  amb_H[:,:,3]
amb_Nt_3                  =  amb_N[:,:,3]
Slnt_range_3              = slnt_range[:,:,3]./1000

theo_res_N_4              =  res_theory_n[:,:,4]
theo_res_S_4              =  res_theory_s[:,:,4]
amb_Ht_4                  =  amb_H[:,:,4]
amb_Nt_4                  =  amb_N[:,:,4]
Slnt_range_4              = slnt_range[:,:,4]./1000

lon_idx1 = 190# 180
lon_idx2 = 210 #223
lat_idx = 74


(plot(Lons_p[lon_idx1:lon_idx2], theo_res_N_1[lon_idx1:lon_idx2,lat_idx],label = "method 1", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], theo_res_N_2[lon_idx1:lon_idx2,lat_idx],label = "method 2", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], theo_res_N_3[lon_idx1:lon_idx2,lat_idx],label = "method 3", linewidth=2))
display(plot!(Lons_p[lon_idx1:lon_idx2], theo_res_N_4[lon_idx1:lon_idx2,lat_idx],xlabel="Longitude [deg]",ylabel="Resolution [m]",title = "Theoretical resolution n"
,label="method 4", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=2))

(plot(Lons_p[lon_idx1:lon_idx2], amb_Ht_1[lon_idx1:lon_idx2,lat_idx],label = "method 1", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], amb_Ht_2[lon_idx1:lon_idx2,lat_idx],label = "method 2", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], amb_Ht_3[lon_idx1:lon_idx2,lat_idx],label = "method 3", linewidth=2))
display(plot!(Lons_p[lon_idx1:lon_idx2], amb_Ht_4[lon_idx1:lon_idx2,lat_idx],xlabel="Longitude [deg]",ylabel="Ambiguity height [m]",title = "Ambiguity along height"
,label="method 4", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=2))

(plot(Lons_p[lon_idx1:lon_idx2], Slnt_range_1[lon_idx1:lon_idx2,lat_idx],label = "method 1", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], Slnt_range_2[lon_idx1:lon_idx2,lat_idx],label = "method 2", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], Slnt_range_3[lon_idx1:lon_idx2,lat_idx],label = "method 3", linewidth=2))
display(plot!(Lons_p[lon_idx1:lon_idx2], Slnt_range_4[lon_idx1:lon_idx2,lat_idx],xlabel="Longitude [deg]",ylabel="Slant range [km]",title = "Slant range"
,label="method 4", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=2))

display(plot(Lons_p[lon_idx1:lon_idx2], plot_v_bprpmax[lon_idx1:lon_idx2,lat_idx],xlabel="Longitude [deg]",ylabel="Max perp. basline [km]",title = "Maximum perpendicular baseline"
,label="method 4", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=2))

display(plot(Lons_p[lon_idx1:lon_idx2], plot_v_la[lon_idx1:lon_idx2,lat_idx],xlabel="Longitude [deg]",ylabel="Look angle [deg]",title = "Look angle"
,label="method 4", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=2))


@load "/Users/joshil/Documents/Code/global_scale_outputs/4p_outputs_1/output_gs_study_res_run_082023_100m_5f_103_4plat_4proc_lag.jld"


plot_v_bpa_s            =  Output_stat_bpa[:,:,1]
plot_v_bpa_c            =  Output_stat_bpa[:,:,2]
plot_v_bpa_h            =  Output_stat_bpa[:,:,3]
plot_v_bpa_n            =  Output_stat_bpa[:,:,4]
plot_v_bpa_n_pslr       =  Output_stat_bpa[:,:,8]
plot_v_bpa_n_islr       =  Output_stat_bpa[:,:,12]
#plot_v_bpa_n_locerr     =  Output_stat_bpa[:,:,16]
#plot_v_bpa_n_pow        =  Output_stat_bpa[:,:,17]

plot_v_beamforming_s    =  Output_stat_beamforming[:,:,1]
plot_v_beamforming_c    =  Output_stat_beamforming[:,:,2]
plot_v_beamforming_h    =  Output_stat_beamforming[:,:,3]
plot_v_beamforming_n    =  Output_stat_beamforming[:,:,4]
plot_v_beamforming_n_pslr       =  Output_stat_beamforming[:,:,8]
plot_v_beamforming_n_islr       =  Output_stat_beamforming[:,:,12]
#plot_v_beamforming_n_locerr     =  Output_stat_beamforming[:,:,16]
#plot_v_beamforming_n_pow        =  Output_stat_beamforming[:,:,17]

plot_v_capon_s          =  Output_stat_capon[:,:,1]
plot_v_capon_c          =  Output_stat_capon[:,:,2]
plot_v_capon_h          =  Output_stat_capon[:,:,3]
plot_v_capon_n          =  Output_stat_capon[:,:,4]
plot_v_capon_n_pslr     =  Output_stat_capon[:,:,8]
plot_v_capon_n_islr     =  Output_stat_capon[:,:,12]
#plot_v_capon_n_locerr   =  Output_stat_capon[:,:,16]
#plot_v_capon_n_pow        =  Output_stat_capon[:,:,17]

(plot(Lons_p[lon_idx1:lon_idx2], theo_res_N_1[lon_idx1:lon_idx2,lat_idx],label = "method 1", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], theo_res_N_2[lon_idx1:lon_idx2,lat_idx],label = "method 2", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], theo_res_N_3[lon_idx1:lon_idx2,lat_idx],label = "method 3", linewidth=2))
(plot!(Lons_p[lon_idx1:lon_idx2], theo_res_N_4[lon_idx1:lon_idx2,lat_idx],xlabel="Longitude [deg]",ylabel="Resolution [m]",title = "Theoretical resolution n"
,label="method 4", legendfontsize=15, tickfont=font(15), ytickfont=font(15), guidefont=font(15), titlefontsize=15, linewidth=2))
display(plot!(Lons_p[lon_idx1:lon_idx2], plot_v_bpa_n[lon_idx1:lon_idx2,lat_idx],label = "BPA", linewidth=2, linestyle=:dash))
#(plot!(Lons_p[lon_idx1:lon_idx2], plot_v_beamforming_n[lon_idx1:lon_idx2,lat_idx],label = "BPA", linewidth=2, linestyle=:dash))

