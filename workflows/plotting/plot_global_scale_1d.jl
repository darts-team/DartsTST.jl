
#Plotting op file
using PyPlot
using GeoDatasets
using JLD2

@load "/Users/joshil/Documents/Code/Outputs/output_gs_study_1D_run_062023_1.jld"

img_path = "/Users/joshil/Documents/Code/global_scale_outputs/plots/1d1/"

#xlim_plot = (-128,-66)
#ylim_plot = (15,53)

xlim_plot = (-180,180)
ylim_plot = (-75,75)

Lats_p = Geo_location.A[1,:,1]
Lons_p = Geo_location.A[:,1,2]

plot_v_cr =  Canopy_heights.A[:,:,1]
plot_v_cr = plot_v_cr[:,:,1]

plot_v_bpa_res =  Output_stat_bpa[:,:,1]
plot_v_bpa_pslr =  Output_stat_bpa[:,:,2]
plot_v_bpa_islr =  Output_stat_bpa[:,:,3]
plot_v_bpa_loc =  Output_stat_bpa[:,:,4]

plot_v_la       =  lookang_all[:,:,1]
plot_v_bprpmax =  Perp_baseline_max[:,:,1]
plot_v_bprpmin =  Perp_baseline_min[:,:,1]
plot_v_bparmax =  Par_baseline_max[:,:,1]
plot_v_bparmin =  Par_baseline_min[:,:,1]
plot_v_bnrmmax =  Norm_baseline_max[:,:,1]
plot_v_bnrmmin =  Norm_baseline_min[:,:,1]

theo_res_N =  res_theory_n[:,:,1]
theo_res_S =  res_theory_s[:,:,1]

pygui(true)
fig = plt.figure(figsize=(11, 6))
ii=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_bpa_res'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("RMSE", size=16)
cbar.ax.tick_params(labelsize=16) 
plt.set_cmap("magma")
plt.clim(0,30)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Back projection algorithm - Resolution along n", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1) # US 1.1
plt.savefig(img_path*"BPA_res.png")



pygui(true)
fig = plt.figure(figsize=(11, 6))
tt=(pcolormesh(Lons_p[:],Lats_p[:],10 .* log10.(plot_v_bpa_pslr')))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("Correlation", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("magma")
plt.clim(-40,20)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Back projection algorithm - PSLR", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"BPA_pslr.png")

pygui(true)
fig = plt.figure(figsize=(11, 6))
tt=(pcolormesh(Lons_p[:],Lats_p[:],(plot_v_bpa_islr')))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("Correlation", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("magma")
plt.clim(-20,30)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Back projection algorithm - ISLR", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"BPA_islr.png")


pygui(true)
fig = plt.figure(figsize=(11, 6))
tt=(pcolormesh(Lons_p[:],Lats_p[:],(plot_v_bpa_loc')))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("Correlation", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("magma")
plt.clim(0, 1)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
pygui(true)
(PyPlot.contour(lon,lat,data',[0.5],colors=[[0.2,0.2,0.2]],linewidths=1.0))
gca().set_aspect(1)
plt.xlabel("Longitude", fontsize=18)
plt.ylabel("Latitude", fontsize=18)
plt.title("Back projection algorithm - Loc error", fontsize=20)
plt.xlim(xlim_plot)
plt.ylim(ylim_plot)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"BPA_loc.png")





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
plt.clim(10,40)
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
plt.clim(0.0,20)
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
plt.clim(0,10)
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
plt.clim(0.0,5)
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
plt.clim(20, 50)
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
plt.clim(5,25)
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


#####-- Canopy heights
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],plot_v_cr[:,:,1]'))
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
plt.savefig(img_path*"CH_1.png")


pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],Canopy_heights.A[:,:,1]'))
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
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
gca().set_aspect(1.1)
plt.savefig(img_path*"CH_2.png")



#Theo_Res_N
pygui(true)
fig = plt.figure(figsize=(11, 6))
rr=(pcolormesh(Lons_p[:],Lats_p[:],theo_res_N'))
cbar = colorbar(orientation="vertical", shrink=0.70)
cbar.set_label("m", size=16)
cbar.ax.tick_params(labelsize=16) 
gca().set_aspect(1)
plt.set_cmap("jet")
plt.clim(5,20)
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