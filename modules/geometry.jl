module Geometry
using ReferenceFrameRotations
using LinearAlgebra

" Create Quaternion based on rotation angle and axis"
quat(rot_angle, rot_ax) = Quaternion(cosd(rot_angle/2.0), rot_ax*sind(rot_angle/2.0))
" Rotate frame given a rotation quaternion "
rotate_frame(v,q) = convert(Array{Float64,1}, vect(inv(q)*v*q))
" Rotate vector, given a rotation quaternion "
rotate_vec(v,q) = convert(Array{Float64,1}, vect(q*v*inv(q)))

"""
Converts Lat/Log/Height to ECEF XYZ position

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
                ranges[m]=distance(t_xyz_grid[:,i],p_xyz[:,j,k])
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
Creates a Peg point based on peg coordinates

# Arguments
 - `peg::3xN Float Array`, peg coordinates [peglat, pegLon, pegHeading] in (deg,deg,deg)
 - `earth_radius::Float64`, optional, earth equatorial radius (in meters)
 - `earth_eccentricity::Float64`, optional, earth eccentricity

# Output
 - `peg::PegPoint`, in addition to peg coordinates also contains other parameters necessary for peg calculations
"""
mutable struct PegPoint
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




abstract type AbstractFrame end
abstract type AbstractPlanetFrame <: AbstractFrame end
abstract type AbstractSpaceCraftFrame <: AbstractPlanetFrame end
abstract type AbstractAntennaFrame <: AbstractSpaceCraftFrame end

mutable struct Planet_frame <:AbstractFrame
    org
    qtrn::Quaternion
    a
    e_sqr
end

mutable struct  SpaceCraft_frame <:AbstractSpaceCraftFrame
    org #origin relative to planet frame
    qtrn::Quaternion# with respect to the planet
    antenna::AbstractAntennaFrame # antenna frame defined relative to sc frame
end

mutable struct Antenna_frame <:AbstractAntennaFrame
    org #origin relative to spacecraft frame
    qtrn::Quaternion # with respect to the SpaceCraft
end

function rotate_spacecraft(sc_frame, rot_angle, rot_ax )
    sc_quat = quat(rot_angle, rot_ax)
    sc_frame.qtrn =  sc_quat * sc_frame.qtrn
    sc_frame.antenna.qtrn =  sc_frame.qtrn * sc_frame.antenna.qtrn
    return sc_frame
end

function rotate_spacecraft(sc_frame, quat::Quaternion )
    sc_frame.qtrn =  quat * sc_frame.qtrn
    sc_frame.antenna.qtrn =  sc_frame.qtrn * sc_frame.antenna.qtrn
    return sc_frame
end

function rotate_antenna(ant_frame, rot_angle, rot_ax )
    ant_quat = quat(rot_angle, rot_ax)
    ant_frame.qtrn =  ant_quat * ant_frame.qtrn
    return ant_frame
end

function rotate_antenna(ant_frame, quat::Quaternion )
    ant_frame.qtrn =  quat * ant_frame.qtrn
    return ant_frame
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
    nhat = -pos/norm(pos); # unit nadir vector
    vrad = dot(vel, -nhat)*(-nhat); #radial velocity
    vtan = vel  - vrad; #tangential velocity
    that = vtan/norm(vtan); #track vector
    chat = cross(nhat, that); #cross-track vector
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







end
