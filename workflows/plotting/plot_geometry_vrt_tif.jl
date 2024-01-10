using Plots
using GeoArrays

file_name = "/Users/joshil/Documents/Outputs/geometry/test_himalayan_region_240m/Himalayan_region.tif"

img_dir = "/Users/joshil/Documents/Outputs/geometry/test_himalayan_region_240m/plots/"



DEM                             = GeoArrays.read(file_name, band=1)
look_ang                        = GeoArrays.read(file_name, band=2)
local_inca                      = GeoArrays.read(file_name, band=3)
ran_slope                       = GeoArrays.read(file_name, band=4)
slant_range                     = GeoArrays.read(file_name, band=5)
OI                              = GeoArrays.read(file_name, band=6)
test_local_inca                 = GeoArrays.read(file_name, band=7)

Geo_location_full                = GeoArrays.coords(DEM)
# Seperate lat and lon from Geo_location for computations
Geo_location_lon_full            = zeros(size(Geo_location_full)[1], size(Geo_location_full)[2],1)
Geo_location_lat_full             = zeros(size(Geo_location_full)[1], size(Geo_location_full)[2],1)

for i=1:size(Geo_location_full)[1]
    for j=1:size(Geo_location_full)[2]
        Geo_location_lon_full[i,j,1] = Geo_location_full[i,j,1][1]
        Geo_location_lat_full[i,j,1] = Geo_location_full[i,j,1][2]
    end
end

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(DEM,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="DEM [m]",fontsize=16))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"DEM.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(look_ang,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Look angle [deg]", clim=(28,42)))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"look_ang.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(local_inca,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Local incidence angle [deg]", clim=(30,60)))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"local_inca.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(ran_slope,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Range slope angle [deg]", clim=(-30,30)))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"ran_slope.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(slant_range,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Slant range [km]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"slant_range.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(test_local_inca,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Local incidence angle with DEM = 0 [deg]", clim=(30,60)))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"test_local_inca.png")






lim1 = 2800
lim2 = 4400#4850#4450#4000#4400
p_len= 50#1000#50#1000


lat_tplot = DEM.f.translation[2] + (lim2 * DEM.f.linear[2,2])
lon_tplot = DEM.f.translation[1] + (lim1 * DEM.f.linear[1,1]): DEM.f.linear[1,1]: DEM.f.translation[1] + ((lim1+p_len) * DEM.f.linear[1,1])


l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, DEM[lim1:lim1+p_len,lim2,1], xlabel="Longitude [deg]", ylabel="elevation [m]", label="DEM", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_1.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, look_ang[lim1:lim1+p_len,lim2,1], label="look angle", linewidth=3))
p = (plot!(lon_tplot, local_inca[lim1:lim1+p_len,lim2,1], xlabel="Longitude [deg]", ylabel="angle [deg]", label="Local incidence angle", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_2.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, local_inca[lim1:lim1+p_len,lim2,1], label="Local incidence angle", linewidth=3))
p = (plot!(lon_tplot, test_local_inca[lim1:lim1+p_len,lim2,1], xlabel="Longitude [deg]", ylabel="angle [deg]", label="Local incidence angle - DEM=0", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13 ))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_3.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, look_ang[lim1:lim1+p_len,lim2,1], label="look angle", linewidth=3))
p = (plot!(lon_tplot, test_local_inca[lim1:lim1+p_len,lim2,1], xlabel="Longitude [deg]", ylabel="angle [deg]", label="Local incidence angle - DEM=0", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13 ))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_4.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white)
p = (plot(lon_tplot, ran_slope[lim1:lim1+p_len,lim2,1], xlabel="Longitude [deg]", ylabel="angle [deg]", label="Range slope angle", linewidth=3, legend=:topleft, tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), legendfontsize=13))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*string(lim1)*"_"*string(lim2)*"_5.png")






l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(look_ang,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Look angle [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"look_ang2.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(local_inca,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Local incidence angle [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"local_inca2.png")

l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(ran_slope,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Range slope angle [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"ran_slope2.png")



l = @layout[grid(1,1) a{0.0001w}]
p0 = plot(legend=false,grid=false,foreground_color_subplot=:white) 
p = (heatmap(test_local_inca,xlabel="Longitude (deg)", ylabel="Latitude (deg)", title="Local incidence angle with DEM = 0 [deg]"))
p11 = plot(p,p0,layout=l)
savefig(p11, img_dir*"test_local_inca2.png")

