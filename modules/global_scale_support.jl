module Global_Scale_Support

using ..Geometry
using ..Orbits
using ..Scene
using ..Scattering

using PlotlyJS 
using DataFrames 
using CSV
using NCDatasets
using Dates
using Parameters
using GeoArrays
using Proj
using CoordinateTransformations
using LinearAlgebra


"""
Read the GEDI L3 data using the filepath specified,
construct the canopy heights for each pixel and get  
the corresponding geo location.

"""
function read_GEDI_L3_data(filepath_GEDIL3, grid_res)
  #choose the grid resolution based on teh requirements, update teh size row and size col accordingly
  if grid_res == 100 #100km 
    #100km x 100km grid 
    size_row      = 347;
    size_col      = 146;
  elseif grid_res == 40 #40km
    ##40km x 40km grid 
    #size_row     = 868;
    #size_col     = 366;
  else
    throw("Resolution not valid! Change to 100 km or 40 km")
  end

  ga              = GeoArrays.read(filepath_GEDIL3)
  replace!(ga,missing => NaN)
  trans           = Proj.Transformation("EPSG:6933", "WGS84")
  trans2          = Proj.Transformation( "WGS84","EPSG:6933")

  global Canopy_heights = GeoArray(zeros(size_row,size_col)) 
  epsg!(Canopy_heights, 6933)
  bbox!(Canopy_heights, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))

  global Geo_location = GeoArray(zeros(size_row,size_col,2))
  epsg!(Geo_location, 6933)
  bbox!(Geo_location, (min_x=-1.736753044e7, min_y=7.314540830638585e6, max_x=1.736753044e7, max_y=-7.314540830638585e6))

  for i=1:size_row
      for j=1:size_col
        Geo_location[i,j,1],Geo_location[i,j,2] = trans(GeoArrays.coords(Canopy_heights,[i,j]))
        EPSG6933_vals                           = trans2.(Geo_location[i,j,1],Geo_location[i,j,2])
        temp                                    = indices(ga,EPSG6933_vals)
        Canopy_heights[i,j,1]                   = ga[temp[1],temp[2]][1]
      end
  end

  return Canopy_heights, Geo_location, size_row, size_col

end

"""
Construct the lat and lon index
to get a 2D matrix corresponding to lat and lon index

"""
function get_lat_lon_idx(region_xlims, region_ylims)

  lat_lon_idx       = Int64.(zeros(length(region_xlims)*length(region_ylims),2))

  global count_i    = 1;
  for x_lim = region_xlims[1]:region_xlims[end]
    for y_lim = region_ylims[1]:region_ylims[end]
      lat_lon_idx[count_i,:] = Int64.([x_lim,y_lim])
      global count_i         = count_i+1;
    end
  end

  return lat_lon_idx

end


"""
Read the orbit file and get the orbit time, pos and vel
Choose data points corresponding to ascending or descending orbits (flag)
flag 0=ascending orbits, 1=descending orbits

"""
function get_orbit_info_fromfile(orbit_dataset, mast_plat, flag)

  t12_orbits 		        = orbit_dataset["time"][1:2] # first two time samples
  dt_orbits 		        = t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
  orbit_time_index      = Int(1):Int(dt_orbits):length(orbit_dataset["time"])
   # index range for orbit times for time interval of interest
  orbit_time1 		      = orbit_dataset["time"][orbit_time_index] # read in time data
  orbit_pos_ECI 	      = 1e3*orbit_dataset["position"][:,:,orbit_time_index] # read in position data, 3 x Np x Nt
  orbit_vel_ECI         = 1e3*orbit_dataset["velocity"][:,:,orbit_time_index] # read in velocity data, 3 x Np x Nt (used optionally in avg peg and heading calculation)
  dv 				            = orbit_dataset.attrib["epoch"];
  epoch 			          = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
  global dcm 		        = Orbits.eci_dcm(orbit_time1, epoch);
  orbit_pos1,orbit_vel1 = Orbits.ecef_orbitpos(orbit_pos_ECI,orbit_vel_ECI,dcm)
  
  orbit_pos_geo         = Geometry.xyz_to_geo(orbit_pos1[:,mast_plat,:])
  
  temp1 				        = orbit_pos_geo[1,1:end-1]
  temp2 				        = orbit_pos_geo[1,2:end]
  temp3 				        = temp1.>temp2
  temp3 				        = [temp3; temp3[end]]
  orbit_pos_geo 		    = orbit_pos_geo[:,temp3.==flag]
  orbit_time 		        = orbit_time1[temp3.==flag]
  orbit_pos 			      = orbit_pos1[:,:,temp3.==flag]
  orbit_vel 			      = orbit_vel1[:,:,temp3.==flag]

  #DEBUG
  #display(plot(temp3[1:10000]))
  #display(plot(orbit_pos_geo[2,:],orbit_pos_geo[1,:]))

  return orbit_time, orbit_pos, orbit_vel, orbit_pos_geo
end


"""
Find the closest lat and lon location (orbit index) corresponding to the pixel location
look angle lower limit is considered
#TODO, if the look angle goes beyongd the upperlimit?

"""
function find_close_val_lat_lon(Geo_location, lat_lon_idx, orbit_pos, orbit_pos_geo)

  search_lim            = 100 
  look_ang_lower_lim    = 29.9 #look angle lower limit

  #Search for the closest orbit point from the geo point
  close_val_lat_lon     = findmin(abs.(orbit_pos_geo[2,1:end-search_lim].-(Geo_location[lat_lon_idx[1],lat_lon_idx[2],2])) 
  + abs.(orbit_pos_geo[1,1:end-search_lim].-(Geo_location[lat_lon_idx[1],lat_lon_idx[2],1])))

  for i2=1:25 #25 can be changed based on the search required
    #Check whether the closest orbit point is to the left of the geo point
    #if no, get the next closest point left of the current orbit
    #if yes, check for look angle limit and get teh next closest point to left if look angle criterion is not satisfied  
    #To ensure right looking SAR geometry
    if orbit_pos_geo[:,close_val_lat_lon[2]][2] >= (Geo_location[lat_lon_idx[1],lat_lon_idx[2],2])
      close_val_lat_lon   = findmin(abs.(orbit_pos_geo[2,1:end-search_lim].-(Geo_location[lat_lon_idx[1],lat_lon_idx[2],2]-(i2.*0.5))) 
      + abs.(orbit_pos_geo[1,1:end-search_lim].-(Geo_location[lat_lon_idx[1],lat_lon_idx[2],1])))
    else
      slrng_temp        = Geometry.distance(orbit_pos[:,1,close_val_lat_lon[2]], Geometry.geo_to_xyz([Geo_location[lat_lon_idx[1],lat_lon_idx[2]]; 0]))
      look_angle        = acosd(orbit_pos_geo[:,close_val_lat_lon[2]][3] / slrng_temp) 
      if (look_angle < (look_ang_lower_lim)) 
        close_val_lat_lon   = findmin(abs.(orbit_pos_geo[2,1:end-search_lim].-(Geo_location[lat_lon_idx[1],lat_lon_idx[2],2]-(i2.*0.5))) 
        + abs.(orbit_pos_geo[1,1:end-search_lim].-(Geo_location[lat_lon_idx[1],lat_lon_idx[2],1])))
      else
        break
      end
    end
  end

  #Now that we have the approximate orbit location, the next step is to ensure it is approximately in the center of the geometry required
  # search for the +-search_lim points around the orbit location obtained from previous step
  # compute distacnces and choose the pount with the minimum distance as the orbit point for simulations corresponding to the geo location
  dist_idx              = close_val_lat_lon[2]-search_lim:close_val_lat_lon[2]+search_lim
  lookvec_ini           = zeros(3,length(dist_idx))
  slrng_ini             = zeros(length(dist_idx))
  for i3=1:length(dist_idx)
    lookvec_ini[:,i3] 	= Geometry.geo_to_xyz([Geo_location[lat_lon_idx[1],lat_lon_idx[2]]; 0])  .- orbit_pos[:,1,dist_idx[i3]]; # look vector
    slrng_ini[i3] 		  = norm(lookvec_ini[:,i3]); #slant range
  end
  close_val_lat_lon     = (findmin(slrng_ini)[1],Int(dist_idx[findmin(slrng_ini)[2]]))

  return close_val_lat_lon

end


"""
Plotting on map initial version

"""
function maps1(lat_vals, lon_vals, p_var, c1, c2)

    marker = attr(size=[60,60,60,60],
                  color=p_var,
                  cmin=c1,
                  cmax=c2,
                  #symbol="square",
                  colorscale="hot",
                  colorbar=attr(title="Correlation",
                                ticksuffix="  ",
                                showticksuffix="last"),
                  line_color="black")
    trace = scattergeo(;mode="markers", lat=lat_vals, lon=lon_vals,
                        marker=marker, marker_size=3,
                        marker_line_color="black", marker_line_width=2,
                        name="Data")

    layout = Layout(geo_scope="usa", geo_resolution=50, width=900, height=600,
                    margin=attr(l=0, r=0, t=10, b=0))
    plot(trace, layout)
end


end