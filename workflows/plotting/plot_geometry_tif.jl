using Plots
using GeoArrays

#path_dir = "../Outputs/global_geometry_10m_new/"
#path_dir = "/Users/joshil/Documents/Outputs/geometry/global_geometry_10km_Geo/global_geometry_10m_new_v2/"
#path_dir = "/Users/joshil/Documents/Outputs/geometry/global_geometry_10km_Geo_v2/1km_him_region/global_geometry_1km_new/"
path_dir = "/Users/joshil/Documents/Outputs/geometry/cartwheel_10km/global_geometry_10km_cartwheel/"
#path_dir = "/Users/joshil/Documents/Outputs/geometry/global_geometry_100km_Geogrid/Outputs/"

#img_dir = "/Users/joshil/Documents/Outputs/geometry/global_geometry_10km_Geo/global_geometry_10m_new_v2_plots/"
#img_dir = "/Users/joshil/Documents/Outputs/geometry/global_geometry_10km_Geo_v2/1km_him_region_plots/"
img_dir = "/Users/joshil/Documents/Outputs/geometry/cartwheel_10km/global_geometry_10km_cartwheel/Plots/"
#img_dir = "/Users/joshil/Documents/Outputs/geometry/global_geometry_100km_Geogrid/Plots/"


filepath                            = path_dir*"DEM_10km.tif"
DEM_orig                            = GeoArrays.read(filepath)
filepath                            = path_dir*"Slope_lon.tif"
slope_lon                           = GeoArrays.read(filepath)
filepath                            = path_dir*"Slope_lat.tif"
slope_lat                           = GeoArrays.read(filepath)
filepath                            = path_dir*"Look_angle.tif"
looka                               = GeoArrays.read(filepath)
filepath                            = path_dir*"local_incidence_angle.tif"
l_inca                              = GeoArrays.read(filepath)
filepath                            = path_dir*"incidence_angle.tif"
inca                                = GeoArrays.read(filepath)
filepath                            = path_dir*"range_slope_angle.tif"
ran_slopea                          = GeoArrays.read(filepath)
filepath                            = path_dir*"slope_angle.tif"
slopea                              = GeoArrays.read(filepath)
filepath                            = path_dir*"Slope_lon_angle.tif"
slope_lon2                          = GeoArrays.read(filepath)
filepath                            = path_dir*"Slope_lat_angle.tif"
slope_lat2                          = GeoArrays.read(filepath)


DEM_orig[isnan.(looka)]             .= NaN
slope_lon[isnan.(looka)]            .= NaN
slope_lat[isnan.(looka)]            .= NaN
slope_lon2[isnan.(looka)]            .= NaN
slope_lat2[isnan.(looka)]            .= NaN


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


l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(DEM_orig, xlabel="Longitude (deg)", ylabel="Latitude (deg)", xlim=xlim_p, ylim=ylim_p, title="DEM [m]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_DEM.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(looka, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(20,60), xlim=xlim_p, ylim=ylim_p, title="Look angle [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_look_angle.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (heatmap(l_inca, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(20,60), xlim=xlim_p, ylim=ylim_p, title="Local incidence angle [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_local_incidence_angle.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (heatmap(inca, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(20,60), xlim=xlim_p, ylim=ylim_p, title="Incidence angle [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_incidence_angle.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (heatmap(ran_slopea, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(-15/3,15/3) , xlim=xlim_p, ylim=ylim_p, title="Range slope angle [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_range_slope_angle.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (heatmap(slopea, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(-5,5), xlim=xlim_p, ylim=ylim_p, title="Slope angle [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_slope_angle.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (heatmap(slope_lon, xlabel="Longitude (deg)", ylabel="Latitude (deg)",clim=(-0.25,0.25), xlim=xlim_p, ylim=ylim_p, title="Slope along longitude (X)"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_Slope_lon.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (heatmap(slope_lat, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(-0.25,0.25), xlim=xlim_p, ylim=ylim_p, title="Slope along latitude (Y)"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_Slope_lat.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (heatmap(slope_lon2, xlabel="Longitude (deg)", ylabel="Latitude (deg)",clim=(-15/3,15/3), xlim=xlim_p, ylim=ylim_p, title="Slope along longitude (X)"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_Slope_lon_deg.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (heatmap(slope_lat2, xlabel="Longitude (deg)", ylabel="Latitude (deg)", clim=(-15/3,15/3), xlim=xlim_p, ylim=ylim_p, title="Slope along latitude (Y)"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(img_idx)*"_Slope_lat_deg.png")


lim1 = 925#645#2900#1190#1170
lim2 = 550#500#550#1105#1110

lim1 = 2850
lim2 = 550


lat_tplot = DEM_orig.f.translation[2] + (lim2 * DEM_orig.f.linear[2,2])
lon_tplot = DEM_orig.f.translation[1] + (lim1 * DEM_orig.f.linear[1,1]): DEM_orig.f.linear[1,1]: DEM_orig.f.translation[1] + ((lim1+50) * DEM_orig.f.linear[1,1])

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, DEM_orig[lim1:lim1+50-1,lim2,1], xlabel="Longitude [deg]", ylabel="elevation [m]", label="DEM", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_1.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, looka[lim1:lim1+50-1,lim2,1], label="look angle", linewidth=3))
p = (plot!(lon_tplot, l_inca[lim1:lim1+50-1,lim2,1], xlabel="Longitude [deg]", ylabel="angle [deg]", label="Local incidence angle", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13, ylim=(25,60)))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_2.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, slope_lon2[lim1:lim1+50-1,lim2,1], label="Slope lon", linewidth=3))
p = (plot!(lon_tplot, slope_lat2[lim1:lim1+50-1,lim2,1], xlabel="Longitude [deg]", ylabel="slope", label="Slope lat", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13 ))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_3.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, ran_slopea[lim1:lim1+50-1,lim2,1], xlabel="Longitude [deg]", ylabel="angle [deg]", label="Range slope angle", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_4.png")






lim1 = 1000
lim2 = 1140

length_plot = 500

lat_tplot = DEM_orig.f.translation[2] + (lim2 * DEM_orig.f.linear[2,2])
lon_tplot = DEM_orig.f.translation[1] + (lim1 * DEM_orig.f.linear[1,1]): DEM_orig.f.linear[1,1]: DEM_orig.f.translation[1] + ((lim1+length_plot) * DEM_orig.f.linear[1,1])

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, DEM_orig[lim1:lim1+length_plot-1,lim2,1], xlabel="Longitude [deg]", ylabel="elevation [m]", label="DEM", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_1.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, looka[lim1:lim1+length_plot-1,lim2,1], label="look angle", linewidth=3))
p = (plot!(lon_tplot, l_inca[lim1:lim1+length_plot-1,lim2,1], xlabel="Longitude [deg]", ylabel="angle [deg]", label="Local incidence angle", linewidth=3, legend=:topright, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13, ylim=(10,70)))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_2.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, slope_lon2[lim1:lim1+length_plot-1,lim2,1], label="Slope angle lon", linewidth=3))
p = (plot!(lon_tplot, slope_lat2[lim1:lim1+length_plot-1,lim2,1], xlabel="Longitude [deg]", ylabel="slope", label="Slope angle lat", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13 ))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_3.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, ran_slopea[lim1:lim1+length_plot-1,lim2,1], xlabel="Longitude [deg]", ylabel="angle [deg]", label="Range slope angle", linewidth=3, legend=:topright, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_4.png")

