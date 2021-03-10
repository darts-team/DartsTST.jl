module Geometry
using ReferenceFrameRotations
using LinearAlgebra

quat(rot_angle, rot_ax) = Quaternion(cosd(rot_angle/2.0), rot_ax*sind(rot_angle/2.0))
rotate_frame(v,q) = vect(inv(q)*v*q)
rotate_vec(v,q) = vect(q*v*inv(q))

function geo_to_xyz(geo,earth_radius,earth_eccentricity) # geo (θϕh) is a 3xN array for N  points
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

function xyz_to_geo(xyz,earth_radius,earth_eccentricity)
    geo=zeros(3)
    x=xyz[1]
    y=xyz[2]
    z=xyz[3]
    b=earth_radius*(1-earth_eccentricity^2)^0.5
    if x>=0
        ϕ=atand(y/x)
    elseif x<0
        ϕ=sign(y)*atand(y/x)+180
    end
    p=(x^2+y^2)^0.5
    alpha=atand((z/p)*(1/(1-e^2))^0.5)
    θ=atand((z+(earth_eccentricity^2/(1-earth_eccentricity^2))*b*sind(alpha)^3)/(p-earth_eccentricity^2*earth_radius*cosd(alpha)^3))
    re=earth_radius/(1-earth_eccentricity^2*sind(θ)^2)^0.5
    h=p/cosd(θ)-re
    geo=[θ, ϕ, h]
    return geo
end

function peg_calculations(peg,earth_radius,earth_eccentricity)
    e2  = earth_eccentricity^2 #eccentricity squared
    #break out peg parameters
    pegθ  = peg[1]*π/180
    pegϕ  = peg[2]*π/180
    peghed  = peg[3]*π/180
    repeg = earth_radius/sqrt(1-e2*sin(pegθ)^2)
    rnpeg = earth_radius*(1-e2)/sqrt((1-e2*sin(pegθ)^2)^3)
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
    Mxyzprime_xyz=Menu_xyz*Mxyzprime_enu
    return Mxyzprime_xyz,O,ra
end

function sch_to_xyz_2(sch,Mxyzprime_xyz,O,ra)
    xyz = zeros(3)
    #break out SCH vectors
    s   = sch[1]
    c   = sch[2]
    h   = sch[3]
    #conversion from S,C,H to Stheta, Ctheta, H coordinates
    Stheta=s/ra
    Clamda=c/ra
    #convert [Stheta, Clamda, h] vector to [X',Y',Z'] vector
    XYZPrime=[(ra+h)*cos(Clamda)*cos(Stheta), (ra+h)*cos(Clamda)*sin(Stheta),(ra+h)*sin(Clamda)]
    #compute the xyz value
    xyz=Mxyzprime_xyz*XYZPrime+O;
    xyz=[xyz[1],xyz[2],xyz[3]]
    return xyz
end

function sch_to_xyz(sch,peg,earth_radius,earth_eccentricity) # works with multiple points (array inputs) #TODO works only for grid
    xyz=zeros(size(sch))
    XYZPrime=zeros(size(sch))
    e2  = earth_eccentricity^2 #eccentricity squared
    #break out SCH vectors
    s   = sch[1,:]
    c   = sch[2,:]
    h   = sch[3,:]
    #break out peg parameters
    pegθ  = peg[1]*π/180
    pegϕ  = peg[2]*π/180
    peghed  = peg[3]*π/180
    repeg = earth_radius/sqrt(1-e2*sin(pegθ)^2)
    rnpeg = earth_radius*(1-e2)/sqrt((1-e2*sin(pegθ)^2)^3)
    ra = repeg*rnpeg/(repeg*cos(peghed)^2+rnpeg*sin(peghed)^2)

    #conversion from S,C,H to Stheta, Ctheta, H coordinates
    Stheta=s/ra
    Clamda=c/ra

    #convert [Stheta, Clamda, h] vector to [X',Y',Z'] vector
    XYZPrime[1,:]=(ra.+h).*cos.(Clamda).*cos.(Stheta)
    XYZPrime[2,:]=(ra.+h).*cos.(Clamda).*sin.(Stheta)
    XYZPrime[3,:]=(ra.+h).*sin.(Clamda)

    #ENU to XYZ transformation matrix
    Menu_xyz = [-sin(pegϕ) -sin(pegθ)*cos(pegϕ) cos(pegθ)*cos(pegϕ);
                 cos(pegϕ) -sin(pegθ)*sin(pegϕ) cos(pegθ)*sin(pegϕ);
                  0            cos(pegθ)             sin(pegθ)]

    # X'Y'Z' to ENU transformation matrix
    Mxyzprime_enu = [0 sin(peghed) -cos(peghed);
                     0 cos(peghed) sin(peghed);
                     1    0           0]
    #Up vector in XYZ
    Uxyz = Menu_xyz*[0 0 1]';

    #vector from center of ellipsoid to pegpoint
    P = [repeg*cos(pegθ)*cos(pegϕ), repeg*cos(pegθ)*sin(pegϕ), repeg*(1-e2)*sin(pegθ)]

    #translation vector
    O = P-ra*Uxyz

    #compute the xyz value
    for i=1:size(sch)[2]
        xyz[:,i]=Menu_xyz*Mxyzprime_enu*XYZPrime[:,i]+O;
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
- `P::3x1 Float Array`: position vector
- `lv::3x1 Float Array`: look vector
- `Ra::Float`: Ellipsoid Radius at P (usually Earth Radius)
- `e::Float`: Ellipsoid eccentricity (usually Earth Eccentricity)

# Return
- `ρ::Float`: range to surface of ellipsoid from P in the direction of lv
"""
function get_rho(P, lv, Ra, e)
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

"""
 Compute TCN frame based on SC position and velocity
 - Usage: t,c,n = get_tcn(pos, vel)

# Arguments
- `pos::3x1 Float Array`: position vector (usually satellite position in XYZ)
- `vel::3x1 Float Array`: velocity vector (usually satellite velocity in XYZ)

# Return
- `t-hat::3x1 Flot Array`: definition of the t_hat (tangential velocity) vector in base frame (usually XYZ)
- `c-hat::3x1 Float Array`: definition of the c_hat (cross-track) vector in base frame (usually XYZ)
- `n-hat::3x1 Float Array`: definition of the n_hat (nadir) vector in base frame (usually XYZ)

 """
function get_tcn(pos, vel)
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


end
