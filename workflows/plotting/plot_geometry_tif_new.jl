using Plots
using GeoArrays
using Measures

path_dir = "/Users/joshil/Documents/Outputs/geometry/cartwheel/global_geometry_10km_cartwheel/"

img_dir = "/Users/joshil/Documents/Outputs/geometry/cartwheel/global_geometry_10km_cartwheel/plots/"


filepath                            = path_dir*"DEM_10km.tif"
DEM_orig                            = GeoArrays.read(filepath)
filepath                            = path_dir*"Baseline_AT_max.tif"
baseline_AT_max                     = GeoArrays.read(filepath)
filepath                            = path_dir*"Baseline_AT_min.tif"
baseline_AT_min                     = GeoArrays.read(filepath)
filepath                            = path_dir*"Baseline_perp_max.tif"
baseline_perp_max                   = GeoArrays.read(filepath)
filepath                            = path_dir*"Baseline_perp_min.tif"
baseline_perp_min                   = GeoArrays.read(filepath)
filepath                            = path_dir*"Altitude.tif"
altitude                            = GeoArrays.read(filepath)
filepath                            = path_dir*"Slant_range.tif"
slant_range                         = GeoArrays.read(filepath)
filepath                            = path_dir*"Look_angle.tif"
looka                               = GeoArrays.read(filepath)
filepath                            = path_dir*"incidence_angle.tif"
inca                                = GeoArrays.read(filepath)
filepath                            = path_dir*"local_incidence_angle.tif"
l_inca                              = GeoArrays.read(filepath)
filepath                            = path_dir*"range_slope_angle.tif"
ran_slopea                          = GeoArrays.read(filepath)

filepath                            = path_dir*"Slope_lon.tif"
slope_lon                           = GeoArrays.read(filepath)
filepath                            = path_dir*"Slope_lat.tif"
slope_lat                           = GeoArrays.read(filepath)
filepath                            = path_dir*"slope_angle.tif"
slopea                              = GeoArrays.read(filepath)
filepath                            = path_dir*"Slope_lon_angle.tif"
slope_lon2                          = GeoArrays.read(filepath)
filepath                            = path_dir*"Slope_lat_angle.tif"
slope_lat2                          = GeoArrays.read(filepath)

filepath                            = path_dir*"Baseline_max.tif"
baseline_max                     = GeoArrays.read(filepath)
filepath                            = path_dir*"Baseline_min.tif"
baseline_min                     = GeoArrays.read(filepath)


filepath                            = path_dir*"Theo_resolution_n.tif"
theo_res_n                          = GeoArrays.read(filepath)
filepath                            = path_dir*"Theo_ambiguity_height.tif"
hoa                                 = GeoArrays.read(filepath)



DEM_orig[isnan.(looka)]             .= NaN
baseline_AT_max[isnan.(looka)]      .= NaN
baseline_AT_min[isnan.(looka)]      .= NaN
baseline_perp_max[isnan.(looka)]    .= NaN
baseline_perp_min[isnan.(looka)]    .= NaN
altitude[isnan.(looka)]             .= NaN
slant_range[isnan.(looka)]          .= NaN
inca[isnan.(looka)]                 .= NaN
l_inca[isnan.(looka)]               .= NaN
ran_slopea[isnan.(looka)]           .= NaN

baseline_max[isnan.(looka)]      .= NaN
baseline_min[isnan.(looka)]      .= NaN
theo_res_n[isnan.(looka)]      .= NaN
hoa[isnan.(looka)]      .= NaN
#=
slope_lon2 = GeoArray(atand.(slope_lon))
slope_lon2.f = slope_lon.f
slope_lon2.crs = slope_lon.crs
slope_lat2 = GeoArray(atand.(slope_lat))
slope_lat2.f = slope_lat.f
slope_lat2.crs = slope_lat.crs



p = (histogram(slope_lon[:],yaxis=(:log10), xlabel="Slope", ylabel="Frequency", legend=false, xlim=(-0.5,0.5)))
savefig(p, img_dir*"hist_lon.png")

p = (histogram(slope_lat[:],yaxis=(:log10), xlabel="Slope", ylabel="Frequency", legend=false, xlim=(-0.5,0.5)))
savefig(p, img_dir*"hist_lat.png")

p = (histogram(slope_lon2[:],yaxis=(:log10), xlabel="Slope angle lon [deg]", ylabel="Frequency", legend=false, xlim=(-50,50)))
savefig(p, img_dir*"hist_lon_deg.png")

p = (histogram(slope_lat2[:],yaxis=(:log10), xlabel="Slope angle lat [deg]", ylabel="Frequency", legend=false, xlim=(-50,50)))
savefig(p, img_dir*"hist_lat_deg.png")

p = (histogram(ran_slopea[:],yaxis=(:log10), xlabel="Range slope angle [deg]", ylabel="Frequency", legend=false, xlim=(-50,50)))
savefig(p, img_dir*"range_slope_deg.png")
=#

xlim_p = (-180,180)
ylim_p = (-80,80)
img_idx = 1

xlim_p = (70,90)
ylim_p = (35-9,50-9)
img_idx = 2

#1km
xlim_p = (70,88)
ylim_p = (27,40)
img_idx = 2

xlim_p = (-125,-115)
ylim_p = (40-9,55-9)
img_idx = 33

xlim_p = (-80,-60)
ylim_p = (-20-9,-9)
img_idx = 4



p = (heatmap(DEM_orig, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="DEM [m]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_DEM.png")

p = (heatmap(baseline_AT_max, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="Maximum AT baseline [km]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_baseline_AT_max.png")

p = (heatmap(baseline_AT_min, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="Minimum AT baseline [km]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_baseline_AT_min.png")

p = (heatmap(baseline_perp_max, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="Maximum Perp baseline [km]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_baseline_perp_max.png")

p = (heatmap(baseline_perp_min, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="Minimum Perp baseline [km]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_baseline_perp_min.png")

p = (heatmap(altitude, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="Altitude [m]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_altitude.png")

p = (heatmap(slant_range, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="Slant range [m]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_slant_range.png")

p = (heatmap(looka, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(20,60), xlim=xlim_p, ylim=ylim_p, title="Look angle [deg]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_look_angle.png")

p = (heatmap(inca, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(20,60), xlim=xlim_p, ylim=ylim_p, title="Incidence angle [deg]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_incidence_angle.png")

p = (heatmap(l_inca, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(20,60), xlim=xlim_p, ylim=ylim_p, title="Local incidence angle [deg]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_local_incidence_angle.png")

p = (heatmap(ran_slopea, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(-15/3,15/3) , xlim=xlim_p, ylim=ylim_p, title="Range slope angle [deg]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_range_slope_angle.png")


p = (heatmap(baseline_AT_max, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="Maximum baseline [km]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_baseline_max.png")

p = (heatmap(baseline_AT_min, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="Minimum baseline [km]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_baseline_min.png")

p = (heatmap(theo_res_n, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, clim=(10, 40), title="Theoretical resolution n [m]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_theo_res_n.png")

p = (heatmap(hoa, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p,  clim=(20 ,50),title="Height of Ambiguity  [m]",rightmargin=5mm))
savefig(p, img_dir*string(img_idx)*"_hoa.png")

#=
p = (heatmap(slopea, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(-5,5), xlim=xlim_p, ylim=ylim_p, title="Slope angle [deg]",rightmargin=3mm))
savefig(p, img_dir*string(img_idx)*"_slope_angle.png")

p = (heatmap(slope_lon, xlabel="Longitude (deg)", ylabel="Latitude (deg)",clim=(-0.25,0.25), xlim=xlim_p, ylim=ylim_p, title="Slope along longitude (X)",rightmargin=3mm))
savefig(p, img_dir*string(img_idx)*"_Slope_lon.png")

p = (heatmap(slope_lat, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(-0.25,0.25), xlim=xlim_p, ylim=ylim_p, title="Slope along latitude (Y)",rightmargin=3mm))
savefig(p, img_dir*string(img_idx)*"_Slope_lat.png")

p = (heatmap(slope_lon2, xlabel="Longitude (deg)", ylabel="Latitude (deg)",clim=(-15/3,15/3), xlim=xlim_p, ylim=ylim_p, title="Slope along longitude (X)",rightmargin=3mm))
savefig(p, img_dir*string(img_idx)*"_Slope_lon_deg.png")

p = (heatmap(slope_lat2, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(-15/3,15/3), xlim=xlim_p, ylim=ylim_p, title="Slope along latitude (Y)",rightmargin=3mm))
savefig(p, img_dir*string(img_idx)*"_Slope_lat_deg.png")
=#


