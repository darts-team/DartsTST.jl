#Plotting op file
using PyPlot
using GeoDatasets
using JLD2

#@load "/Users/joshil/Documents/Code/global_scale_outputs/outputs_5p_profiles/output_gs_study_res_run_062023_100m_5f_513_5plat_5proc_profiles.jld"
@load "/Users/joshil/Documents/Code/global_scale_outputs/outputs_5p_v2/output_gs_study_res_run_062023_100m_5f_901_5plat_5proc_profiles.jld"


img_path = "/Users/joshil/Documents/Code/global_scale_outputs/plots_5p_v2/100m_5f_901_5plat_5proc_profiles/"


@load "/Users/joshil/Documents/GEDI_Data/Outputs_L2/Output_GEDIL2_1year_2022.jld2"

#xlim_plot = (-128,-66)
#ylim_plot = (15,53)

xlim_plot = (-180,180)
ylim_plot = (-52,52)

Lats_p = Geo_location.A[1,:,1]
Lons_p = Geo_location.A[:,1,2]

#plot_v_cr =  Canopy_heights.A[:,:,1]
#plot_v_cr = plot_v_cr[:,:,1]

plot_v_rmse_bpa =  Output_stat_bpa[:,:,1]
plot_v_corr_bpa =  Output_stat_bpa[:,:,2]
plot_v_ipnp_bpa =  Output_stat_bpa[:,:,3]
plot_v_opnp_bpa =  Output_stat_bpa[:,:,4]

plot_v_rmse_beamforming =  Output_stat_beamforming[:,:,1]
plot_v_corr_beamforming =  Output_stat_beamforming[:,:,2]
plot_v_ipnp_beamforming =  Output_stat_beamforming[:,:,3]
plot_v_opnp_beamforming =  Output_stat_beamforming[:,:,4]


plot_v_rmse_capon =  Output_stat_capon[:,:,1]
plot_v_corr_capon =  Output_stat_capon[:,:,2]
plot_v_ipnp_capon =  Output_stat_capon[:,:,3]
plot_v_opnp_capon =  Output_stat_capon[:,:,4]


plot_v_la =  lookang_all[:,:,1]
plot_v_bprpmax =  Perp_baseline_max[:,:,1]
plot_v_bprpmin =  Perp_baseline_min[:,:,1]
plot_v_bparmax =  Par_baseline_max[:,:,1]
plot_v_bparmin =  Par_baseline_min[:,:,1]
plot_v_bnrmmax =  Norm_baseline_max[:,:,1]
plot_v_bnrmmin =  Norm_baseline_min[:,:,1]

theo_res_N =  res_theory_n[:,:,1]
theo_res_S =  res_theory_s[:,:,1]

amb_Ht                  =  amb_H[:,:,1]
amb_Nt                  =  amb_N[:,:,1]
Slnt_range              = slnt_range[:,:,1]./1000






plot_v_rmse_bpa2 = plot_v_rmse_bpa
plot_v_rmse_bpa2[isnan.(plot_v_rmse_bpa2)] .= 0.0
plot_v_rmse_bpa2 = plot_v_rmse_bpa2 ./ maximum(plot_v_rmse_bpa2)
plot_v_rmse_bpa =  Output_stat_bpa[:,:,1]
plot_v_rmse_bpa2[isnan.(plot_v_rmse_bpa)] .= NaN

plot_v_rmse_beamforming2 = plot_v_rmse_beamforming
plot_v_rmse_beamforming2[isnan.(plot_v_rmse_beamforming2)] .= 0.0
plot_v_rmse_beamforming2 = plot_v_rmse_beamforming2 ./ maximum(plot_v_rmse_beamforming2)
plot_v_rmse_beamforming =  Output_stat_beamforming[:,:,1]
plot_v_rmse_beamforming2[isnan.(plot_v_rmse_beamforming)] .= NaN

idx1=203
idx2=101
println(plot_v_rmse_bpa2[idx1,idx2])
println(plot_v_rmse_beamforming2[idx1,idx2])
println(plot_v_corr_bpa[idx1,idx2])
println(plot_v_corr_beamforming[idx1,idx2])
println(plot_v_ipnp_bpa[idx1,idx2])
println(plot_v_ipnp_beamforming[idx1,idx2])
println(plot_v_opnp_bpa[idx1,idx2])
println(plot_v_opnp_beamforming[idx1,idx2])

#=
mask_c_bpa = isnan.(plot_v_corr_bpa)
total_points_corr_bpa = sum(mask_c_bpa.==0)
abv_thre_points_corr_bpa = sum(plot_v_corr_bpa.>0.75)
per_corr_abv_thre_bpa = abv_thre_points_corr_bpa / total_points_corr_bpa

mask_c_beamforming = isnan.(plot_v_corr_beamforming)
total_points_corr_beamforming = sum(mask_c_beamforming.==0)
abv_thre_points_corr_beamforming = sum(plot_v_corr_beamforming.>0.75)
per_corr_abv_thre_beamforming = abv_thre_points_corr_beamforming / total_points_corr_beamforming

=#
#=
plot_v_rmse[isnan.(plot_v_cr)] .= NaN
plot_v_corr[isnan.(plot_v_cr)] .= NaN
plot_v_ipnp[isnan.(plot_v_cr)] .= NaN
plot_v_opnp[isnan.(plot_v_cr)] .= NaN

plot_v_la[isnan.(plot_v_cr)] .= NaN
plot_v_bprpmax[isnan.(plot_v_cr)] .= NaN
plot_v_bprpmin[isnan.(plot_v_cr)] .= NaN
plot_v_bparmax[isnan.(plot_v_cr)] .= NaN
plot_v_bparmin[isnan.(plot_v_cr)] .= NaN
plot_v_bnrmmax[isnan.(plot_v_cr)] .= NaN
plot_v_bnrmmin[isnan.(plot_v_cr)] .= NaN

theo_res_N[isnan.(plot_v_cr)] .= NaN
theo_res_S[isnan.(plot_v_cr)] .= NaN
=#

#display(Plots.heatmap(Lons_p[xlims_p[end]:-1:xlims_p[1]],Lats_p[ylims_p],plot_v_cr[:,:,1]))
#display(Plots.heatmap(rotl90(plot_v_rmse[:,:,1]), clim=(0,1), xlabel="", ylabel="", title="RMSE"))
#display(Plots.heatmap(rotl90(plot_v_corr[:,:,1]), clim=(0,1), xlabel="", ylabel="", title="Corr"))
#display(Plots.heatmap(rotl90(plot_v_ipnp[:,:,1]), clim=(0,5), xlabel="", ylabel="", title="# i/p peaks"))
#display(Plots.heatmap(rotl90(plot_v_opnp[:,:,1]), clim=(0,10), xlabel="", ylabel="", title="# o/p peaks"))

heights_L2 = Canopy_heights_L2[:,:,1]
heights_L2[isnan.(Slnt_range)] .= NaN
#####-- Canopy heights
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],heights_L2'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("Height (m)", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("viridis")
plt.clim(0,30)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Canopy heights", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.clim(0,5)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"L2_H.png")


#RMSE
pygui(true)
fig = plt.figure(figsize=(11, 6))
ii=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_rmse_bpa'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("", size=16)
cbar.ax.tick_params(labelsize=16) 
plt.set_cmap("magma")
plt.clim(0,1)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Back projection - RMSE", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1) # US 1.1
plt.savefig(img_path*"RMSE_BPA.png")

#CORR
pygui(true)
fig = plt.figure(figsize=(11, 6))
tt=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_corr_bpa'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("magma")
plt.clim(0,1)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Back projection - Correlation", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(-50,50)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"CORR_BPA.png")

bounds = [1,2,3]
#Ip_peaks
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_ipnp_bpa'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("Number of peaks", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("cividis")
plt.clim(1,3)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Back projection - Number of Input Peaks", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(-50,50)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"IPP_BPA.png")

#Op_peaks
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_opnp_bpa'))
cbar = colorbar(orientation="vertical", shrink=0.70) #,ticks=bounds
cbar.set_label("Number of peaks", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("cividis")
plt.clim(1,3)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Back projection - Number of Output Peaks", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"OPP_BPA.png")


#RMSE


pygui(true)
fig = plt.figure(figsize=(11, 6))
ii=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_rmse_beamforming'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("", size=16)
cbar.ax.tick_params(labelsize=16) 
plt.set_cmap("magma")
plt.clim(0,1)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Beamforming - RMSE", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1) # US 1.1
plt.savefig(img_path*"RMSE_beamforming.png")

#CORR
pygui(true)
fig = plt.figure(figsize=(11, 6))
tt=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_corr_beamforming'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("magma")
plt.clim(0,1)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Beamforming - Correlation", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(-50,50)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"CORR_beamforming.png")

bounds = [1,2,3]
#Ip_peaks
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_ipnp_beamforming'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("Number of peaks", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("cividis")
plt.clim(1,3)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Beamforming - Number of Input Peaks", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(-50,50)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"IPP_beamforming.png")

#Op_peaks
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_opnp_beamforming'))
cbar = colorbar(orientation="vertical", shrink=0.70) #,ticks=bounds
cbar.set_label("Number of peaks", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("cividis")
plt.clim(1,3)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Beamforming - Number of Output Peaks", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"OPP_beamforming.png")





pygui(true)
fig = plt.figure(figsize=(11, 6))
ii=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_rmse_capon'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("", size=16)
cbar.ax.tick_params(labelsize=16) 
plt.set_cmap("magma")
plt.clim(0,1)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("CAPON - RMSE", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1) # US 1.1
plt.savefig(img_path*"RMSE_capon.png")

#CORR
pygui(true)
fig = plt.figure(figsize=(11, 6))
tt=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_corr_capon'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("magma")
plt.clim(0,1)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("CAPON - Correlation", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(-50,50)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"CORR_capon.png")

bounds = [1,2,3]
#Ip_peaks
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_ipnp_capon'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("Number of peaks", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("cividis")
plt.clim(1,3)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("CAPON - Number of Input Peaks", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(-50,50)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"IPP_capon.png")

#Op_peaks
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_opnp_capon'))
cbar = colorbar(orientation="vertical", shrink=0.70) #,ticks=bounds
cbar.set_label("Number of peaks", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("cividis")
plt.clim(1,3)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("CAPON - Number of Output Peaks", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"OPP_capon.png")



#Look angles
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_la'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("Degrees", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("cividis")
plt.clim(30,44)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Look angle", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"LookAng.png")

#PerpBaselineMax
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_bprpmax'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("km", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(15,40)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Max Perp Baseline", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"PerpBaselineMax.png")

#PerpBaselineMin
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_bprpmin'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("km", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(0.0,15)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Min Perp Baseline", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"PerpBaselineMin.png")

#ParBaselineMax
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_bparmax'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("km", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(0,4)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Max AT Baseline", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"ATBaselineMax.png")

#ParBaselineMin
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_bparmin'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("km", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(0.0,1.5)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Min AT Baseline", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"ATBaselineMin.png")


#NormBaselineMax
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_bnrmmax'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("km", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(20, 45)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Max Baseline", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"NormBaselineMax.png")

#NormBaselineMin
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_bnrmmin'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("km", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(0,15)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Min Baseline", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"NormBaselineMin.png")





#Theo_Res_N
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],theo_res_N'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("m", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(5,12)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Theoretical resolution n", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"Theo_res_n.png")



#NormBaselineMin
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],theo_res_S'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("m", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(2.5, 3.5)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Theoretical resolution along-track", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"Theo_res_s.png")


amb_Ht[amb_Ht.==0.0].=NaN
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],amb_Ht'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("[m]", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(10,25)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Ambiguity along Height", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"Amb_H.png")

amb_Nt[amb_Nt.==0.0].=NaN
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],amb_Nt'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("[m]", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(15,30)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Ambiguity along n", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"Amb_N.png")

#Slnt_range[Slnt_range.==0.0].=NaN
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],Slnt_range'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("[km]", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(850,1050)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Range to target", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"Slnt_Range.png")