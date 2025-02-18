module Geometry

#external packages
using ReferenceFrameRotations
using LinearAlgebra
using Statistics
using StaticArrays
using Parameters

"""
Creates a Peg point based on peg coordinates

# Arguments
 - `peg::3xN Float Array`, peg coordinates [peglat, pegLon, pegHeading] in (deg,deg,deg)
 - `earth_radius::Float64`, optional, earth equatorial radius (in meters)
 - `earth_eccentricity::Float64`, optional, earth eccentricity

# Output
 - `peg::PegPoint`, in addition to peg coordinates also contains other parameters necessary for peg calculations
"""
struct PegPoint
    pegLat::Float64 # peg latitude (deg)
    pegLon::Float64 # peg longitude (deg)
    pegHdg::Float64 # peg heading (deg)
    eq_rad::Float64 # earth equatorial Radius
    ecc::Float64 # earth Eccentricity
    Mxyzprime_xyz::Array{Float64,2} #rotation matrix from xyz prime to xyz
    O::Array{Float64,2} # translation vector
    Ra::Float64 # earth radius
    function PegPoint(pegLat, pegLon,pegHdg,eq_rad::Float64=6.378137e6,ecc::Float64=0.08181919084267456)
        e2  = ecc^2 #eccentricity squared
        #break out peg parameters
        pegθ  = pegLat*π/180
        pegϕ  = pegLon*π/180
        peghed  = pegHdg*π/180
        repeg = eq_rad/sqrt(1-e2*sin(pegθ)^2)
        rnpeg = eq_rad*(1-e2)/sqrt((1-e2*sin(pegθ)^2)^3)
        ra = repeg*rnpeg/(repeg*cos(peghed)^2+rnpeg*sin(peghed)^2)

        #ENU to XYZ transformation matrix
        Menu_xyz = [-sin(pegϕ) -sin(pegθ)*cos(pegϕ) cos(pegθ)*cos(pegϕ);
                     cos(pegϕ) -sin(pegθ)*sin(pegϕ) cos(pegθ)*sin(pegϕ);
                      0            cos(pegθ)             sin(pegθ)]

        #X'Y'Z' to ENU transformation matrix
        Mxyzprime_enu = [0 sin(peghed) -cos(peghed);
                         0 cos(peghed) sin(peghed);
                         1    0           0]
        #Up vector in XYZ
        Uxyz = Menu_xyz*[0 0 1]';

        #vector from center of ellipsoid to pegpoint
        P = [repeg*cos(pegθ)*cos(pegϕ), repeg*cos(pegθ)*sin(pegϕ), repeg*(1-e2)*sin(pegθ)]

        #translation vector
        O = P-ra*Uxyz

        #rotation vector X'Y'Z' to XYZ
        Mxyzprime_xyz=Menu_xyz*Mxyzprime_enu
        new(pegLat, pegLon, pegHdg, eq_rad, ecc, Mxyzprime_xyz, O, ra)
    end

end


" Create Quaternion based on rotation angle and axis"
quat(rot_angle, rot_ax) = Quaternion(cosd(rot_angle/2.0), rot_ax*sind(rot_angle/2.0))
" Rotate frame given a rotation quaternion "
rotate_frame(v,q) = convert(Array{Float64,1}, vect(inv(q)*v*q))
" Rotate vector, given a rotation quaternion "
rotate_vec(v,q) = convert(Array{Float64,1}, vect(q*v*inv(q)))

function convert_platform_target_scene_coordinates(Np,Nst,Nt,p_xyz,t_xyz_3xN,targets_loc,s_xyz_3xN,s_loc_3xN,avg_peg, params)
    @unpack display_geometry_coord, ts_coord_sys = params

    Nsc=size(s_loc_3xN,2)
    p_loc=zeros(3,Np,Nst)
    t_loc=zeros(3,Nt)
    s_loc=zeros(Nsc)
    if display_geometry_coord=="LLH"
        for i=1:Np
            p_xyz_i=p_xyz[:,i,:]
            p_xyz_i=reshape(p_xyz_i,3,Nst)
            p_loc[:,i,:]=xyz_to_geo(p_xyz_i)
        end
        t_loc=xyz_to_geo(t_xyz_3xN)
        s_loc=xyz_to_geo(s_xyz_3xN)
    elseif display_geometry_coord=="SCH"
        for i=1:Np
            p_xyz_i=p_xyz[:,i,:]
            p_xyz_i=reshape(p_xyz_i,3,Nst)
            p_loc[:,i,:]=xyz_to_sch(p_xyz_i,avg_peg)
        end
        if ts_coord_sys=="SCH"
            t_loc=targets_loc
            s_loc=s_loc_3xN
        else # ts_coord_sys == LLH or XYZ
            t_loc=xyz_to_sch(t_xyz_3xN,avg_peg)
            s_loc=xyz_to_sch(s_xyz_3xN,avg_peg)
        end
    elseif display_geometry_coord=="XYZ"
        p_loc=p_xyz
        t_loc=t_xyz_3xN
        s_loc=s_xyz_3xN
    end
    return p_loc,t_loc,s_loc
end

"""
Calculates average peg point and average platform heights over all platforms and all slow-time (pulse) locations
"""
function avg_peg_h(p_xyz)
    # Average Platform Heading
      Np=size(p_xyz)[2] # number of platforms
      Nst=size(p_xyz)[3] # number of slow-time samples (pulses processed)
      p_headings=zeros(1,Np)
      p_geo=zeros(3,Np,Nst)
      for i=1:Np
        p_xyz_i=p_xyz[:,i,:]
        p_xyz_i=reshape(p_xyz_i,3,Nst)
        p_geo[:,i,:]=xyz_to_geo(Float64.(p_xyz_i))
        p_lat=reshape(p_geo[1,i,:],1,Nst)
        p_lon=reshape(p_geo[2,i,:],1,Nst)
        p_headings[i]=mean(compute_heading(p_lat,p_lon))
      end
      p_avg_heading=mean(p_headings)
      p_avg_geo=mean(mean(p_geo,dims=2),dims=3) # average LLH of platforms over platforms and slow-time locations
      p_h_avg=p_avg_geo[3]
      avg_peg=PegPoint(p_avg_geo[1],p_avg_geo[2],p_avg_heading)
      return avg_peg,p_h_avg
end
function avg_peg_h(p_xyz,p_vel)
    # Average Platform Heading
      Np=size(p_xyz)[2] # number of platforms
      Nst=size(p_xyz)[3] # number of slow-time samples (pulses processed)
      p_headings=zeros(1,Np)
      p_geo=zeros(3,Np,Nst)
      for i=1:Np
        p_xyz_i=p_xyz[:,i,:]
        p_xyz_i=reshape(p_xyz_i,3,Nst)
        p_geo[:,i,:]=xyz_to_geo(Float64.(p_xyz_i))
        p_vel_i=reshape(p_vel[:,i,:],3,Nst)
        p_lat=reshape(p_geo[1,i,:],1,Nst)
        p_lon=reshape(p_geo[2,i,:],1,Nst)
        p_headings[i]=mean(compute_heading(p_lat,p_lon,p_vel_i))
      end
      p_avg_heading=mean(p_headings)
      p_avg_geo=mean(mean(p_geo,dims=2),dims=3) # average LLH of platforms over platforms and slow-time locations
      p_h_avg=p_avg_geo[3]
      avg_peg=PegPoint(p_avg_geo[1],p_avg_geo[2],p_avg_heading)
      # Median heading
      #p_vel_i=p_vel[:,Int(round(Np/2)),Int(round(Nst/2))]
      #p_median_heading=compute_heading(p_geo[1,Int(round(Np/2)),Int(round(Nst/2))],p_geo[2,Int(round(Np/2)),Int(round(Nst/2))],p_vel_i)
      #p_h_avg=p_geo[3,Int(round(Np/2)),Int(round(Nst/2))]
      #avg_peg=PegPoint(p_geo[1,Int(round(Np/2)),Int(round(Nst/2))],p_geo[2,Int(round(Np/2)),Int(round(Nst/2))],p_median_heading[1])
      return avg_peg,p_h_avg
end

"""
Converts Lat/Lon/Height to ECEF XYZ position

# Arguments
 - `geo::3x1 Float Array`, Lat/Long/Height (deg,deg,m)
 - `earth_radius::Float64`, optional, earth equatorial radius (in meters)
 - `earth_eccentricity::Float64`, optional, earth eccentricity

# Output
 - `xyz::3x1 Float Array`, ECEF xyz positions (in meters)
"""
function geo_to_xyz(geo,earth_radius::Float64=6.378137e6,earth_eccentricity::Float64=0.08181919084267456)
    xyz=zeros(size(geo))
    θ=geo[1,:] # latitude [deg]
    ϕ=geo[2,:] # longitude [deg]
    h=geo[3,:] # height [m]
    re=earth_radius./(float(1).-earth_eccentricity^2*sind.(θ).^2).^0.5
    xyz[1,:]=(re+h).*cosd.(θ).*cosd.(ϕ)
    xyz[2,:]=(re+h).*cosd.(θ).*sind.(ϕ)
    xyz[3,:]=(re.*(float(1).-earth_eccentricity^2)+h).*sind.(θ)
    return xyz
end

function find_min_max_range(t_xyz_grid,p_xyz)
    ranges=zeros(1,size(t_xyz_grid,2)*size(p_xyz,2)*size(p_xyz,3))
    m=1;
    for i=1:size(t_xyz_grid,2)
        for j=1:size(p_xyz,2)
            for k=1:size(p_xyz,3)
                #ranges[m]=distance(t_xyz_grid[:,i],p_xyz[:,j,k])
                @views ranges[m] = distance(t_xyz_grid[:, i], p_xyz[:, j, k])
                m=m+1;
            end
        end
    end
    min_range=minimum(ranges)
    max_range=maximum(ranges)
    return min_range,max_range
end

function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
end

"""
Converts ECEF XYZ position to Lat/Log/Height

# Arguments
 - `xyz::3xN Float Array`, ECEF xyz positions (in meters)
 - `earth_radius::Float64`, optional, earth equatorial radius (in meters)
 - `earth_eccentricity::Float64`, optional, earth eccentricity

# Output
 - `geo::3xN Float Array`, Lat/Long/Height (deg,deg,m)
"""
function xyz_to_geo(xyz::Array{Float64,},earth_radius::Float64=6.378137e6,earth_eccentricity::Float64=0.08181919084267456)
    ecc_sq = earth_eccentricity^2; #eccentricity squared
    geo=zeros(size(xyz))
    x=xyz[1,:]
    y=xyz[2,:]
    z=xyz[3,:]

    #latitude computation
    semi_minor_axis=earth_radius*sqrt(1-earth_eccentricity^2)
    p=sqrt.(x.^2+y.^2)
    α=atan.(z, p*sqrt(1-ecc_sq))
    θ=atan.(z+ecc_sq/(1-ecc_sq)*semi_minor_axis*sin.(α).^3, p-ecc_sq*earth_radius*cos.(α).^3)

    #longitude
    ϕ=atan.(y,x)

    #height computation
    re=earth_radius./sqrt.(1.0.-ecc_sq.*sin.(θ).^2)
    h=p./cos.(θ)-re

    geo[1,:] = θ*180/π
    geo[2,:] = ϕ*180/π
    geo[3,:] = h
    return geo
end

"""
Converts SCH coordinates to ECEF XYZ

# Arguments
 - `SCH::3xN Float Array`, SCH coordinates in meters
 - `peg::PegPointType`, see PegPoint for more details

# Output
 - `XYZ:3xN Float Array`, ECEF XYZ coordinates
"""
function sch_to_xyz(sch::Array{Float64,1},peg::PegPoint)

    #conversion from S,C,H to Stheta, Ctheta, H coordinates
    Stheta=sch[1]/peg.Ra
    Clamda=sch[2]/peg.Ra

    #convert [Stheta, Clamda, h] vector to [X',Y',Z'] vector
    XYZPrime=[(peg.Ra+sch[3])*cos(Clamda)*cos(Stheta), (peg.Ra+sch[3])*cos(Clamda)*sin(Stheta),(peg.Ra+sch[3])*sin(Clamda)]

    #compute and return the xyz value
    return peg.Mxyzprime_xyz*XYZPrime+peg.O
end
function sch_to_xyz(sch::Array{Float64,2},peg::PegPoint)
    @assert size(sch,1)==3 "SCH vector needs to be 3xN"

    #set up xyz output
    xyz=zeros(size(sch))

    #compute the xyz value per SCH triplet
    for ipt=1:size(sch)[2]
        xyz[:,ipt]=sch_to_xyz(sch[:,ipt],peg)
    end
    return xyz
end
"""
Converts ECEF XYZ coordinates to SCH

# Arguments
 - `XYZ::3xN Float Array`, XYZ (ECEF) coordinates in meters
 - `peg::PegPointType`, see PegPoint for more details

# Output
 - `SCH:3xN Float Array`, SCH coordinates
"""
function xyz_to_sch(xyz::Array{Float64,1},peg)

    #compute x'y'z' coordinates from peg and xyz
    xyzp = (peg.Mxyzprime_xyz)'*(xyz-peg.O)

    #norm of the x'y'z' vector (this is equal to peg.Ra + h)
    r = norm(xyzp);

    #compute s and c coordinates
    s = peg.Ra*atan(xyzp[2], xyzp[1]);
    c = peg.Ra*asin(xyzp[3]/r);
    h = r - peg.Ra;
    return [s,c,h]
end
function xyz_to_sch(xyz::Array{Float64,2},peg)
    @assert size(xyz,1)==3 "SCH vector needs to be 3xN"

    #set up xyz output
    sch=zeros(size(xyz))

    #compute the sch value per xyz triplet
    for ipt=1:size(xyz,2)
        sch[:,ipt]=xyz_to_sch(xyz[:,ipt],peg)
    end
    return sch
end

"""
Compute range from ray-ellipse intersection
 - Usage: ρ = get_rho(position, look vector, Earth Radius, Earth Eccentricity)

# Arguments
- `P::3xN Float Array`: position vector
- `lv::3xN Float Array`: look vector
- `Ra::Float`: [optional] Ellipsoid Equatorial Radius
- `e::Float`: [optional] Ellipsoid eccentricity (usually Earth Eccentricity)

# Return
- `ρ::Float`: range to surface of ellipsoid from P in the direction of lv
"""
function get_rho(P::Array{Float64,1}, lv::Array{Float64,1}, Ra::Float64=6.378137e6, e::Float64=0.08181919084267456)
    @assert size(P,1)==3 "POS vector needs to be 3 x 1"
    @assert size(lv,1)==3 "Look vector needs to be 3 x 1"
    @assert size(P)==size(lv) "POS and Look Vector need to have same size"

    #normalize look vector
    lv = lv./norm(lv);

    #estimate rho based on ray-ellipse intersection
    ra = (((lv[1]^2) + (lv[2]^2)) / (Ra^2)) + ((lv[3]^2) / ((Ra^2) * (1 - (e^2))))
    rb = 2 * ((((lv[1] * P[1]) + (lv[2] * P[2])) / (Ra^2)) + ((lv[3] * P[3]) / ((Ra^2) * (1 - (e^2)))))
    rc = (((P[1]^2) + (P[2]^2)) / (Ra^2)) + ((P[3]^2) / ((Ra^2) * (1 - (e^2)))) - 1
    if (((rb^2) - (4 * ra * rc))) < 0
        error("Cannot compute rho, input correct platform co-ordinates or look vector")
    else
        ρ = (- rb - sqrt((rb^2) - (4 * ra * rc))) / (2 * ra)
    end
    return ρ
end
function get_rho(pos::Array{Float64,2}, lv::Array{Float64,2}, Ra::Float64=6.378137e6, e::Float64=0.08181919084267456)
    @assert size(pos,1)==3 "POS vector needs to be 3 x N"
    @assert size(lv,1)==3 "Look vector needs to be 3 x N"
    @assert size(pos)==size(lv) "POS and Look Vector need to have same size"

    #initialize variables
    ρ = zeros(size(pos,2))

    #compute t,c,n triplet for each position and velocity
    for itp=1:size(pos)[2]
        ρ[itp] = get_rho(pos[:,itp], lv[:,itp], Ra, e)
    end
    return ρ
end

"""
 Compute TCN frame based on SC position and velocity
 - Usage: t,c,n = get_tcn(pos, vel)

# Arguments
- `pos::3xN Float Array`: position vector (usually satellite position in XYZ)
- `vel::3xN Float Array`: velocity vector (usually satellite velocity in XYZ)

# Return
- `t-hat::3xN Flot Array`: definition of the t_hat (tangential velocity) vector in base frame (usually XYZ)
- `c-hat::3xN Float Array`: definition of the c_hat (cross-track) vector in base frame (usually XYZ)
- `n-hat::3xN Float Array`: definition of the n_hat (nadir) vector in base frame (usually XYZ)

 """
function get_tcn(pos::Array{Float64,1}, vel::Array{Float64,1})
    @assert length(pos)==3 "POS needs to be 3 x 1"
    @assert length(vel)==3 "VEL needs to be 3 x 1"
    @assert size(pos)==size(vel) "POS and VEL need to have same size"

    #create geocetric TCN frame for reference satellite described
    #in coordinates of input 'position' and 'velocity' assumed to be ECEF
    llh = xyz_to_geo(pos);
    e,n,u = enu_from_geo(llh[1], llh[2])

    nhat = -dropdims(u,dims=2); # unit nadir vector
    chat = cross(nhat, vel/norm(vel)); #cross-track vector
    that = cross(chat, nhat); #track vector
    return that, chat, nhat
end
function get_tcn(pos::Array{Float64,2}, vel::Array{Float64,2})
    @assert size(pos,1)==3 "POS needs to be 3 x N"
    @assert size(vel,1)==3 "VEL needs to be 3 x N"
    @assert size(pos)==size(vel) "POS and VEL need to have same size"

    #initialize variables
    that = zeros(size(pos))
    chat = zeros(size(pos))
    nhat = zeros(size(pos))

    #compute t,c,n triplet for each position and velocity
    for itp=1:size(pos)[2]
        that[:,itp], chat[:,itp], nhat[:,itp] = get_tcn(pos[:,itp], vel[:,itp])
    end
    return that, chat, nhat
end


"""
Compute orientation quaternion based on TCN basis
- Usage: quat = tcn_quat(pos, vel)

# Arguments
- `pos::3xN Float Array`: position vector (usually satellite position in XYZ)
- `vel::3xN Float Array`: velocity vector (usually satellite velocity in XYZ)

# Return
- `quat::

"""
function tcn_quat(pos::Array{Float64,1}, vel::Array{Float64,1})
    #get TCN basis from position and velocity
    t,c,n = Geometry.get_tcn(pos, vel)
    return q_tcn = [dcm_to_quat(DCM([t';c';n']))]
end
function tcn_quat(pos::Array{Float64,2}, vel::Array{Float64,2})
    @assert size(pos,1)==3 "POS needs to be 3 x N"
    @assert size(vel,1)==3 "VEL needs to be 3 x N"
    @assert size(pos)==size(vel) "POS and VEL need to have same size"

    #initialize variables
    quat = Array{Quaternion{Float64},1}(undef, size(pos)[2])
    #quat = zeros(4,size(pos)[2])

    #compute tcn quaternion for each position and velocity
    for itp=1:size(pos)[2]
        quat[itp] = tcn_quat(pos[:,itp], vel[:,itp])[1]
    end
    return quat
end



"""
 Compute look vector based on TCN frame
 - Usage: look_vector = tcn_lvec(t, c, n, θ, ϕ)

# Arguments
- `t::3x1 Float Array`: t-axis of the TCN frame in parent frame (typically ECEF)
- `c::3x1 Float Array`: c-axis of the TCN frame in parent frame (typically ECEF)
- `n::3x1 Float Array`: n-axis of the TCN frame in parent frame (typically ECEF)
- `θ::Float`: elevation angle (rad)
- `ϕ::Float`: azimuth angle (rad)

# Return
- `look_vec::3x1 Flot Array`: unit look vector in TCN parent frame (typically ECEF)
 """
function tcn_lvec(t, c, n, θ_el, ϕ_az)
    return  t*sin(θ_el)*sin(ϕ_az) + c*sin(θ_el)*cos(ϕ_az) + n*cos(θ_el)
end

"""
Compute heading based on lat/long or lat/long/velocity. If only lat/long are provided
the heading is computed based on difference in successive geodetic points (needs at
least two geodetic points). Thre resulting heading is of length N-1.
In case velocity vector (ECEF) is also given heading is computed for every point.
- Usage:
    + heading = compute_heading(lat, lon)
    + heading = compute_heading(lat, lon, velocity)

# Arguments
    - `lat::N-element Array`: geodetic latitude (deg)
    - `lon::N-element Array`: geodetic longitude (deg)
    - `vel::3xN-element Array`: (optional) velocity in ECEF (m/s)

# Return
    - `heading::N (or N-1) element Float Array`: heading relative to north (0 is due north, 90 is due east) (deg)
"""
function compute_heading(lat, lon)
    @assert length(lat) == length(lon) "lat and lon inputs should be same size"
    @assert length(lat) > 1 "need at least two geodetic points to compute heading"
    #compute heading based on geographic coordinates
    Δlon = diff(lon[:]); #difference in longitudes
    lat_start = lat[1:end-1]; #start latitutde of
    lat_end   = lat[2:end];
    X = cosd.(lat_end).*sind.(Δlon);
    Y = cosd.(lat_start).*sind.(lat_end) - sind.(lat_start).*cosd.(lat_end).*cosd.(Δlon);
    heading = atan.(X,Y)*180/π;
    heading[findall(x->x<0.0, heading)] = heading[findall(x->x<0.0, heading)] .+360;
    heading[findall(x->x>359.9, heading)] .= 0;
    return heading
end
function compute_heading(lat, lon, vel)
    @assert length(lat) == length(lon) == size(vel,2) "lat, lon, vel inputs should be same size"
    @assert size(vel,1)==3 "velocity should be a 3xN vector"
    @assert size(lat,1)==1 "lat/lon should be 1xN vectors"

    #compute ENU basis for each lat, lon
    e,n,u = Geometry.enu_from_geo(lat, lon)

    #space for headding
    heading = zeros(size(lat));

    for ii=1:length(lat)
        v_enu = cat(e[:,ii],n[:,ii],u[:,ii], dims=2)'*vel[:,ii];
        v_enu = v_enu./norm(v_enu);
        heading[ii] = atan(v_enu[1], v_enu[2])*180/π;
        if heading[ii] < 0 && heading[ii] < -0.5
            heading[ii] = heading[ii] + 360.0;
        elseif heading[ii] < 0 && heading[ii] >= -0.5
            heading[ii] = 0.0;
        end
    end

    return heading
end


"""
Compute ENU basis vectors based on lat, lon
# Usage
    - e,n,u = enu_from_geo(lat, lon)

# Arguments
    - `lat::N-element Array`: geodetic latitude (deg)
    - `lon::N-element Array`: geodetic longitude (deg)

# Return
- `E::N-element Array`: Easting axis in ECEF (m)
- `N::N-element Array`: Northing axis in ECEF (m)
- `U::N-element Array`: Up axis in ECEF (m)
"""
function enu_from_geo(lat, lon)
    @assert length(lat) == length(lon) "lat and lon inputs should be same size"

    #e, n, u basis vectors
    e = zeros(3,length(lat));
    n = zeros(3,length(lat));
    u = zeros(3,length(lat));

    for ii=1:length(lat)
        e[:,ii] = cat(-sind.(lon[ii]),                  cosd.(lon[ii]),                 0, dims =1)
        n[:,ii] = cat(-cosd.(lon[ii]).*sind.(lat[ii]), -sind.(lon[ii]).*sind.(lat[ii]), cosd.(lat[ii]), dims =1)
        u[:,ii] = cat( cosd.(lon[ii]).*cosd.(lat[ii]),  sind.(lon[ii]).*cosd.(lat[ii]), sind.(lat[ii]), dims =1)
    end
    return e,n,u
end

"""
Compute ENU basis vectors from LLH with an origin based on new location
# Usage
    - e,n,u = geo_to_enu_new_org(lat, lon, orgllh)
# Arguments
    - `lat::N-element Array`: geodetic latitude (deg)
    - `lon::N-element Array`: geodetic longitude (deg)
    - `height::N-element Array`: height (m)
    - `orgllh::3 x 1 Array`: new origin of ENU in LLH (m)
# Return
- `E::N-element Array`: Easting axis in ECEF (m) w.r.t org_ece
- `N::N-element Array`: Northing axis in ECEF (m) w.r.t org_ece
- `U::N-element Array`: Up axis in ECEF (m) w.r.t org_ece
"""
function llh_to_enu_new_org(lat, lon, height, orgllh,earth_radius::Float64=6.378137e6,earth_eccentricity::Float64=0.08181919084267456)
   
    ESQ = earth_eccentricity^2
    orgece = geo_to_xyz(orgllh);
    numEntries = length(lat)
    enu = zeros(3,numEntries)
    for i = 1:numEntries
        ecef = zeros(3,1) 
        
        SP = sind(lat[i])
        CP = cosd(lat[i])
        SL = sind(lon[i])
        CL = cosd(lon[i])
        GSQ = 1.0 - (ESQ*SP*SP);
        EN = earth_radius / sqrt(GSQ);
        Z = (EN + height[i]) * CP;
        ecef[1,1] = Z * CL;
        ecef[2,1] = Z * SL;
        EN = EN - (ESQ * EN);
        ecef[3,1] = (EN + height[i]) * SP;
        
        difece = ecef - orgece;   # difference between coordinate origins
        
        #Rotate the difference vector into ENU coordinates
        
        sla = sind(orgllh[1,1]); cla = cosd(orgllh[1,1]);
        slo = sind(orgllh[2,1]); clo = cosd(orgllh[2,1]);
        
        enu[:,i] = [  -slo      clo      0 ; 
            -sla*clo  -sla*slo  cla;
                cla*clo   cla*slo  sla] * difece;       
    end#for

    # E = enu[1,:]; N = enu[2,:]; U = enu[3,:];
    return enu
end#function


"""
Compute Spherical Coordinates from ENU base system
# Usage
    - theta,phi,r = geo_to_enu_new_org(e,n,u)
# Arguments
    - `east::N-element Array`: Easting axis in ECEF (m)
    - `north::N-element Array`: Northing axis in ECEF (m)
    - `up::N-element Array`: Up axis in ECEF (m)
# Return
- `sph::3xN-element Array`: spherical coordinates with [1,:] = range, [2,:] = elevation, [3,:] = azimuth
- `elev::N-element Array`: elev angle (deg)
- `phi::N-element Array`: azimuth angle (deg)
- `range::N-element Array`: range from origin to spherical point
"""
function enu_to_sph(east, north, up)
    numEntries = length(east);
    sph = zeros(3,numEntries)
    for i = 1:numEntries
        
        az,elev,range = cart2sph(east[i],north[i],up[i])
            
        sph[1,i] = range#radius
        sph[2,i] = elev#elevation
        sph[3,i] = az#azimuth

        # println("range: $range")
        # println("elev: $elev")
        # println("az: $az")
        @assert  ((elev >= -90) & (elev <= 90) & (az >= -180) & (az <= 180) & (range >= 0)) "Invalid Range"
    end
   
    return sph
end#function

"""
Compute Spherical Coordinates from cartesian
# Usage
    - theta,phi,r = cart2sph(x,y,z)
# Arguments
    - `x::Float32`: x coord (m)
    - `y::Float32`: y coord (m)
    - `z::Float32`: z coord (m)
# Return
- `elev::N-element Array`: elev angle (deg)
- `az::N-element Array`: azimuth angle (deg)
- `r::N-element Array`: range (m)
"""
function cart2sph(x,y,z)
    hypotxy = sqrt(x^2+y^2)
    r = sqrt(hypotxy^2+z^2)
    elev = atand(z,hypotxy)
    az = atand(y,x)
    return az,elev,r
end#function
"""
Compute ENU basis vectors from LLH with an origin based on new location
# Usage
    - e,n,u = geo_to_enu_new_org(lat, lon, orgllh)

# Arguments
    - `lat::N-element Array`: geodetic latitude (deg)
    - `lon::N-element Array`: geodetic longitude (deg)
    - `height::N-element Array`: height (m)
    - `orgllh::3 x 1 Array`: new origin of ENU in LLH (m)

# Return
- `E::N-element Array`: Easting axis in ECEF (m) w.r.t org_ece
- `N::N-element Array`: Northing axis in ECEF (m) w.r.t org_ece
- `U::N-element Array`: Up axis in ECEF (m) w.r.t org_ece
"""
function llh_to_enu_new_org(lat, lon, height, orgllh,earth_radius::Float64=6.378137e6,earth_eccentricity::Float64=0.08181919084267456)
   
    ESQ = earth_eccentricity^2
    orgece = geo_to_xyz(orgllh);
    numEntries = length(lat)
    enu = zeros(3,numEntries)
    for i = 1:numEntries
        ecef = zeros(3,1) 
        
        SP = sind(lat[i])
        CP = cosd(lat[i])
        SL = sind(lon[i])
        CL = cosd(lon[i])
        GSQ = 1.0 - (ESQ*SP*SP);
        EN = earth_radius / sqrt(GSQ);
        Z = (EN + height[i]) * CP;
        ecef[1,1] = Z * CL;
        ecef[2,1] = Z * SL;
        EN = EN - (ESQ * EN);
        ecef[3,1] = (EN + height[i]) * SP;
        
        difece = ecef - orgece;   # difference between coordinate origins
        
        #Rotate the difference vector into ENU coordinates
        
        sla = sind(orgllh[1,1]); cla = cosd(orgllh[1,1]);
        slo = sind(orgllh[2,1]); clo = cosd(orgllh[2,1]);
        
        enu[:,i] = [  -slo      clo      0 ; 
            -sla*clo  -sla*slo  cla;
                cla*clo   cla*slo  sla] * difece;       
    end#for

    # E = enu[1,:]; N = enu[2,:]; U = enu[3,:];
    return enu
end#function


"""
Compute Spherical Coordinates from ENU base system
# Usage
    - theta,phi,r = geo_to_enu_new_org(e,n,u)

# Arguments
    - `east::N-element Array`: Easting axis in ECEF (m)
    - `north::N-element Array`: Northing axis in ECEF (m)
    - `up::N-element Array`: Up axis in ECEF (m)

# Return
- `sph::3xN-element Array`: spherical coordinates with [1,:] = range, [2,:] = elevation, [3,:] = azimuth
- `elev::N-element Array`: elev angle (deg)
- `phi::N-element Array`: azimuth angle (deg)
- `range::N-element Array`: range from origin to spherical point
"""
function enu_to_sph(east, north, up)
    numEntries = length(east);
    sph = zeros(3,numEntries)
    for i = 1:numEntries
        
        az,elev,range = cart2sph(east[i],north[i],up[i])
            
        sph[1,i] = range#radius
        sph[2,i] = elev#elevation
        sph[3,i] = az#azimuth

        # println("range: $range")
        # println("elev: $elev")
        # println("az: $az")
        @assert  ((elev >= -90) & (elev <= 90) & (az >= -180) & (az <= 180) & (range >= 0)) "Invalid Range"
    end
   
    return sph
end#function

"""
Compute Spherical Coordinates from cartesian
# Usage
    - theta,phi,r = cart2sph(x,y,z)

# Arguments
    - `x::Float32`: x coord (m)
    - `y::Float32`: y coord (m)
    - `z::Float32`: z coord (m)

# Return
- `elev::N-element Array`: elev angle (deg)
- `az::N-element Array`: azimuth angle (deg)
- `r::N-element Array`: range (m)
"""
function cart2sph(x,y,z)
    hypotxy = sqrt(x^2+y^2)
    r = sqrt(hypotxy^2+z^2)
    elev = atand(z,hypotxy)
    az = atand(y,x)
    return az,elev,r
end#function

end #end module
