module Scene

function form3Dgrid_for(t_θ,t_ϕ,t_h) # for loop method
  t_geo_grid=zeros(3,length(t_θ)*length(t_ϕ)*length(t_h))
  m=0
  for i=1:length(t_θ)
    for j=1:length(t_ϕ)
      for k=1:length(t_h)
        m=m+1
        t_geo_grid[:,m]=[t_θ[i],t_ϕ[j],t_h[k]]
      end
    end
  end
  return t_geo_grid
end

function form3Dgrid_array(t_θ,t_ϕ,t_h) # array method
  t_θ1=Array{Float64}(undef,1,length(t_θ))
  t_ϕ1=Array{Float64}(undef,1,length(t_ϕ))
  t_h1=Array{Float64}(undef,1,length(t_h))
  t_θ1[:]=t_θ
  t_ϕ1[:]=t_ϕ
  t_h1[:]=t_h
  t_θ_all=repeat(t_θ1,inner=[1,1],outer=[1,length(t_ϕ)*length(t_h)])
  t_ϕ_all=repeat(t_ϕ1,inner=[1,length(t_θ)],outer=[1,length(t_h)])
  t_h_all=repeat(t_h1,inner=[1,length(t_θ)*length(t_ϕ)],outer=[1,1])
  t_geo_grid=[t_θ_all;t_ϕ_all;t_h_all]
  return t_geo_grid
end

function lookangle_to_range(a,el,p_h) # look angle to slant/ground range,  p_h is platform height, spherical planet
    inc=asind(sind(el)*(a+p_h)/a); # incidence angle
    alpha=inc-el; # deg planet-central angle
    rg=a*alpha*pi/180; # ground range
    rs=(((a+h_p)*sind(alpha)).^2+((a+h_p).*cosd(alpha)-a).^2).^0.5; # slant range
    return rs, rg
end

function slantrange_to_lookangle(a,rs,p_h) # slant range to look angle and ground range,  p_h is platform height, spherical planet
  el=acos((rs.^2+(a+p_h).^2-a^2)./(2*rs.*(a+p_h))); # rad look angle
  inc=asin((a+p_h)./a.*sin(el)); # rad incidence angle
  alpha=inc-el; # rad planet-central angle
  rg=a*alpha; # ground range
  return rg,el*180/pi
end

function groundrange_to_lookangle(a,rg,p_h) # ground range to look angle and slant range,  p_h is platform height, spherical planet
  alpha=rg/a # rad planet-central angle
  rs=(((a+h_p)*sind(alpha)).^2+((a+h_p).*cosd(alpha)-a).^2).^0.5; # slant range
  el=acos((rs.^2+(a+p_h).^2-a^2)./(2*rs.*(a+p_h))); # rad look angle
  return rs,el*180/pi
end

function azelh_to_xyz(azelh,p_geo,a,e)  # azelh is a 3xN array for N  points, p_h is platform height

    az=azelh[1,:] #  azimuth angles
    el=azelh[2,:] # elevation angles
    h=azelh[3,:] # heights

    p_h=p_geo[3] # platform height

    rs,rg=lookangle_to_range(a,el,p_h) # slant range between target and platform, assumes spherical planet

    vL_sch=rs*[sin(el)*sin(az),sin(el)*cos(az),-cos(el)] # look vectors in geo (for each target at each look angle
    peg=[p_geo[1],p_geo[2],0] # peg point in geo
    vL_xyz=sch_to_xyz(vL_sch,peg,a,e) # look vectors in xyz (for each target at each look angle)

    vP=geo_to_xyz(p_geo,a,e) # platform vector in xyz

    vT=vL+vP # target vectors in xyz

    # rotate look vector around platform vector by azimuth degrees

    return vT
end

end
