module Geometry
using ReferenceFrameRotations

quat(rot_angle, rot_ax) = Quaternion(cosd(rot_angle/2.0), rot_ax*sind(rot_angle/2.0))
rotate_frame(v,q) = vect(inv(q)*v*q)
rotate_vec(v,q) = vect(q*v*inv(q))

function geo_to_xyz(geo,a,e) # geo (θϕh) is a 3xN array for N  points
    xyz=zeros(size(geo))
    θ=geo[1,:] # latitude
    ϕ=geo[2,:] # longitude
    h=geo[3,:] # height
    re=a./(float(1).-e^2*sind.(θ).^2).^0.5
    xyz[1,:]=(re+h).*cosd.(θ).*cosd.(ϕ)
    xyz[2,:]=(re+h).*cosd.(θ).*sind.(ϕ)
    xyz[3,:]=(re.*(float(1).-e^2)+h).*sind.(θ)
    return xyz
end

function xyz_to_geo(xyz,a,e)
    geo=zeros(3)
    x=xyz[1]
    y=xyz[2]
    z=xyz[3]
    b=a*(1-e^2)^0.5
    if x>=0
        ϕ=atand(y/x)
    elseif x<0
        ϕ=sign(y)*atand(y/x)+180
    end
    p=(x^2+y^2)^0.5
    alpha=atand((z/p)*(1/(1-e^2))^0.5)
    θ=atand((z+(e^2/(1-e^2))*b*sind(alpha)^3)/(p-e^2*a*cosd(alpha)^3))
    re=a/(1-e^2*sind(θ)^2)^0.5
    h=p/cosd(θ)-re
    geo=[θ, ϕ, h]
    return geo
end

function peg_calculations(peg,a,e)
    e2  = e^2 #eccentricity squared
    #break out peg parameters
    pegθ  = peg[1]*π/180
    pegϕ  = peg[2]*π/180
    peghed  = peg[3]*π/180
    repeg = a/sqrt(1-e2*sin(pegθ)^2)
    rnpeg = a*(1-e2)/sqrt((1-e2*sin(pegθ)^2)^3)
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

function sch_to_xyz(sch,peg,a,e) # works with multiple points (array inputs) #TODO works only for grid
    xyz=zeros(size(sch))
    XYZPrime=zeros(size(sch))
    e2  = e^2 #eccentricity squared
    #break out SCH vectors
    s   = sch[1,:]
    c   = sch[2,:]
    h   = sch[3,:]
    #break out peg parameters
    pegθ  = peg[1]*π/180
    pegϕ  = peg[2]*π/180
    peghed  = peg[3]*π/180
    repeg = a/sqrt(1-e2*sin(pegθ)^2)
    rnpeg = a*(1-e2)/sqrt((1-e2*sin(pegθ)^2)^3)
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




end
