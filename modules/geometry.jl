module Geometry

function geo_to_xyz(geo,a,e)
    xyz=zeros(3)
    lat=geo[1]
    lon=geo[2]
    h=geo[3]
    re=a/(1-e^2*sind(lat)^2)^0.5
    xyz[1]=(re+h)*cosd(lat)*cosd(lon)
    xyz[2]=(re+h)*cosd(lat)*sind(lon)
    xyz[3]=(re*(1-e^2)+h)*sind(lat)
    return xyz
end

function xyz_to_geo(xyz,a,e)
    geo=zeros(3)
    x=xyz[1]
    y=xyz[2]
    z=xyz[3]
    b=a*(1-e^2)^0.5
    if x>=0
        lon=atand(y/x)
    elseif x<0
        lon=sign(y)*atand(y/x)+180
    end
    p=(x^2+y^2)^0.5
    alpha=atand((z/p)*(1/(1-e^2))^0.5)
    lat=atand((z+(e^2/(1-e^2))*b*sind(alpha)^3)/(p-e^2*a*cosd(alpha)^3))
    re=a/(1-e^2*sind(lat)^2)^0.5
    h=p/cosd(lat)-re
    geo=[lat, lon, h]
    return geo
end

function sch_to_xyz(sch,peg,a,e)
    xyz=zeros(3)
    s=sch[1]
    c=sch[2]
    h=sch[3]
    lat0=peg[1]
    lon0=peg[2]
    heading=peg[3]
    re=a/(1-e^2*sind(lat0)^2)^0.5
    rn=a*(1-e^2)/(1-e^2*sind(lat0)^2)^1.5
    ra=re*rn/(re*cosd(heading)^2+rn*sind(heading)^2)
    xp=(ra+h)*cosd(c/ra)*cosd(s/ra)
    yp=(ra+h)*cosd(c/ra)*sind(s/ra)
    zp=(ra+h)*sind(c/ra)
    xyzp=[xp,yp,zp]
    M_ENU_to_xyz=[-sind(lon0) -sind(lat0)*cosd(lon0) cosd(lat0)*cosd(lon0);cosd(lon0) -sind(lat0)*sind(lon0) cosd(lat0)*sind(lon0);0 cosd(lat0) sind(lat0)]
    M_xyzp_to_ENU=[0 sind(heading) -cosd(heading);0 cosd(heading) sind(heading);1 0 0]
    P=[re*cosd(lat0)*cosd(lon0),re*cosd(lat0)*sind(lon0),re*(1-e^2)*sind(lat0)]
    ENU=M_xyzp_to_ENU*xyzp
    E=ENU[1]
    N=ENU[2]
    U=ENU[3]
    O=P.-ra*U
    xyz=M_ENU_to_xyz*ENU+O
    return xyz
end

end
