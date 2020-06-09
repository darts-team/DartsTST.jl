module Scene

include("geometry.jl")

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

function lookangle_to_range(a,θ_l,p_h) # look angle to slant/ground range,  p_h is platform height, spherical planet
    inc=asind.(sind.(θ_l)*(a+p_h)/a) # deg incidence angle
    α=inc-θ_l # deg planet-central angle
    rg=a*α*pi/180 # ground range
    rs=(((a+p_h)*sind.(α)).^2+((a+p_h).*cosd.(α).-a).^2).^0.5 # slant range
    return rs, rg
end

function slantrange_to_lookangle(a,rs,p_h) # slant range to look angle and ground range,  p_h is platform height, spherical planet
  θ_l=acos.((rs.^2+(a+p_h).^2-a^2)./(2*rs.*(a+p_h))) # rad look angle
  inc=asin.((a+p_h)./a.*sin.(θ_l)) # rad incidence angle
  α=inc-θ_l # rad planet-central anglesin
  rg=a*α # ground range
  return rg,θ_l*180/pi
end

function groundrange_to_lookangle(a,rg,p_h) # ground range to look angle and slant range,  p_h is platform height, spherical planet
  α=rg/a # rad planet-central angle
  rs=(((a+p_h)*sind.(α)).^2+((a+p_h).*cosd.(α)-a).^2).^0.5 # slant range
  θ_l=acos.((rs.^2+(a+p_h).^2-a^2)./(2*rs.*(a+p_h))) # rad look angle
  return rs,θ_l*180/pi
end

function lookh_to_xyz(lookh,p_geo,peg,a,e)  # look is a 3xN array for N  points, p_h is platform height
    vL2_sch=zeros(size(lookh))
    peg_xyz_grid=zeros(size(lookh))
    vL_xyz=zeros(size(lookh))
    vT=zeros(size(lookh))
    ϕ_l=lookh[1,:] # deg  azimuth angles
    θ_l=lookh[2,:] # deg elevation angles
    h=lookh[3,:] # heights

    p_h=p_geo[3] # platform height
    p_xyz=Geometry.geo_to_xyz(p_geo,a,e) # platform position in xyz

    rs,rg=lookangle_to_range(a,θ_l,p_h) # slant range between target and platform, assumes spherical planet

    vL2_sch[1,:]=sind.(θ_l).*sind.(ϕ_l)
    vL2_sch[2,:]=sind.(θ_l).*cosd.(ϕ_l)
    vL2_sch[3,:]=-cosd.(θ_l) # look vectors in geo (for each target at each look angle
    vL2_xyz=Geometry.sch_to_xyz(vL2_sch,peg,a,e) # look vectors in xyz (for each target at each look angle)
    peg_geo=[peg[1],peg[2],0]
    peg_xyz=Geometry.geo_to_xyz(peg_geo,a,e)
    peg_xyz_grid=repeat(peg_xyz,1,size(lookh)[2])
    vL_xyz=rs.*(vL2_xyz.-peg_xyz_grid)
    vT=vL_xyz.+p_xyz # target vectors in xyz (for azimuth angle of 0 deg)
    #vT=[vT[1,1],vT[2,1],vT[3,1]]

    # rotate look vector around platform vector by azimuth degrees
  #  q = Geometry.quat(ϕ_l, vP) #create a quaternion to rotate a vector by ϕ_l degrees about zaxis
  #  vT_rotated = Geometry.rotate_vec(vT, q) #rotate a vector aligned with the x-axis, by q1

    return vT#_rotated
end

end
