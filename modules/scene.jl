module Scene
#TODO add function definitions, comments, define input types, remove unused functions

include("geometry.jl")
using LinearAlgebra

"generates 3D (volumetric) grid from three 1D (linear) arrays using nested for-loops\n
Inputs\n
  array1,array2,array3: three 1D arrays of size N1,N2,N3 representing the principal axes of the 3D volume\n
Outputs\n
  vol_grid: coordinates of 3D (volumetric) grid of points as 3xN matrix where N=N1xN2xN3"
function form3Dgrid_for(array1,array2,array3)
  vol_grid=zeros(3,length(array1)*length(array2)*length(array3))
  m=0
  for i=1:length(array1)
    for j=1:length(array2)
      for k=1:length(array3)
        m=m+1
        vol_grid[:,m]=[array1[i],array2[j],array3[k]]
      end
    end
  end
  return vol_grid
end

"generates 3D (volumetric) grid from three 1D (linear) arrays using array processing\n
Inputs\n
  array1,array2,array3: three 1D arrays of size N1,N2,N3 representing the principal axes of the 3D volume\n
Outputs\n
  vol_grid: coordinates of 3D (volumetric) grid of points as 3xN matrix where N=N1xN2xN3"
function form3Dgrid_array(array1,array2,array3)
  array11=Array{Float64}(undef,1,length(array1))
  array22=Array{Float64}(undef,1,length(array2))
  array33=Array{Float64}(undef,1,length(array3))
  array11[:].=array1
  array22[:].=array2
  array33[:].=array3
  array1_all=repeat(array11,inner=[1,1],outer=[1,length(array2)*length(array3)])
  array2_all=repeat(array22,inner=[1,length(array1)],outer=[1,length(array3)])
  array3_all=repeat(array33,inner=[1,length(array1)*length(array2)],outer=[1,1])
  vol_grid=[array1_all;array2_all;array3_all]
  return vol_grid
end

"look angle to slant/ground range conversion (spherical planet assumed)\n
Inputs\n
  p_h: platform height\n
  ra: radius of spherical planet\n
  t_h: target heights vector\n
Outputs\n
  rs: slant ranges to the targets\n
  rg: ground ranges to the targets"
function lookangle_to_range(ra,θ_l,p_h,t_h)
    ra=ra.+t_h
    p_h=p_h.-t_h
    inc=asind.(sind.(θ_l).*(ra+p_h)./ra) # deg incidence angle
    α=inc-θ_l # deg planet-central angle
    rg=ra.*α*pi/180 # ground range
    rs=(((ra+p_h).*sind.(α)).^2+((ra.+p_h).*cosd.(α).-ra).^2).^0.5 # slant range
    return rs, rg
end

"slant range to look angle and ground range conversion (spherical planet assumed)\n
Inputs\n
  p_h: platform height\n
  ra: radius of spherical planet\n
  rs: slant ranges to the targets\n
Outputs\n
  θ_l: look angles to the targets\n
  rg: ground ranges to the targets"
function slantrange_to_lookangle(ra,rs,p_h) # target height is assumed 0 TODO add target height
  θ_l=acos.((rs.^2+(ra+p_h).^2-ra^2)./(2*rs.*(ra+p_h))) # rad look angle
  inc=asin.((ra+p_h)./ra.*sin.(θ_l)) # rad incidence angle
  α=inc-θ_l # rad planet-central anglesin
  rg=ra*α # ground range
  return rg,θ_l*180/pi
end

"ground range to look angle and slant range conversion (spherical planet assumed)\n
Inputs\n
  p_h: platform height\n
  ra: radius of spherical planet\n
  rg: ground ranges to the targets\n
Outputs\n
  θ_l: look angles to the targets\n
  rs: slant ranges to the targets"
function groundrange_to_lookangle(ra,rg,p_h) # # target height is assumed 0 TODO add target height
  α=rg/ra # rad planet-central angle
  rs=(((ra+p_h)*sind.(α)).^2+((ra+p_h).*cosd.(α)-ra).^2).^0.5 # slant range
  θ_l=acos.((rs.^2+(ra+p_h).^2-ra^2)./(2*rs.*(ra+p_h))) # rad look angle
  return rs,θ_l*180/pi
end

"look angles and target heights to target xyz conversion (spherical planet assumed)\n
Inputs\n
  earth_radius: radius of earth \n
  earth_eccentricity: eccentricity of earth\n
  peg: peg point coordinates for sch to xyz conversion \n
  p_θϕh: platform location in geographic coordinates (lat/lon/height)\n
  lookh: 3xN array for N point targets (1st dimension: azimuth look angles, 2nd dimension: elevation look angles, 3rd dimension: target heights) \n
Outputs\n
  vT_xyz: target positions in xyz"
function lookh_to_xyz(lookh,p_θϕh,peg,earth_radius,earth_eccentricity)
    vL2_sch=zeros(size(lookh))
    peg_xyz_grid=zeros(size(lookh))
    vL_xyz=zeros(size(lookh))
    vT=zeros(size(lookh))
    ϕ_l=lookh[1,:] # deg  azimuth angles
    θ_l=lookh[2,:] # deg elevation angles
    t_h=lookh[3,:] # target heights
    p_h=p_θϕh[3] # platform height
    p_xyz=Geometry.geo_to_xyz(p_θϕh,earth_radius,earth_eccentricity) # platform position in xyz
    rs,rg=lookangle_to_range(earth_radius,θ_l,p_h,t_h) # slant range between target and platform, assumes spherical planet TODO t_h can be a smaller vector and then repeat
    vL2_sch[1,:]=sind.(θ_l).*sind.(ϕ_l)
    vL2_sch[2,:]=sind.(θ_l).*cosd.(ϕ_l) #TODO add d for left/right looking
    vL2_sch[3,:]=-cosd.(θ_l) # look vectors in geo (for each target at each look angle
    vL2_xyz=Geometry.sch_to_xyz(vL2_sch,peg,earth_radius,earth_eccentricity) # look vectors in xyz (for each target at each look angle)
    peg_geo=[peg[1],peg[2],0]
    peg_xyz=Geometry.geo_to_xyz(peg_geo,earth_radius,earth_eccentricity)
    peg_xyz_grid=repeat(peg_xyz,1,size(lookh)[2])
    rs_grid=repeat(rs',3,1)
    vL_xyz=rs_grid.*(vL2_xyz.-peg_xyz_grid)
    vT_xyz=vL_xyz.+p_xyz # target vectors in xyz (for azimuth angle of 0 deg)
    return vT_xyz
end

"\n
Inputs\n
  earth_radius: radius of earth \n
  earth_eccentricity: eccentricity of earth\n
  peg: peg point coordinates for sch to xyz conversion \n
  p_θϕh: platform location in geographic coordinates (lat/lon/height)\n
  t_c: cross-track coordinates of target points \n
  t_h: height coordinates of target points\n
  rot_P: amount of rotations (as an array) in degrees of look vector about platform position vector\n
Outputs\n
  t_xyz_grid_rot: target positions in xyz as a 3xN array"
function chP_to_xyz_grid(t_c,t_h,rot_P,p_θϕh,peg,earth_radius,earth_eccentricity)
  t_s=0 # before rotation around P, targets are at the same along-track position with platform
  p_xyz=Geometry.geo_to_xyz(p_θϕh,earth_radius,earth_eccentricity) # platform position in xyz
  p_xyz=p_xyz/norm(p_xyz)
  t_sch_grid=Scene.form3Dgrid_for(t_s,t_c,t_h) # using 3 nested for loops
  #t_geo_grid=Scene.form3Dgrid_array(t_s,t_c,t_h) # using array processing
  t_xyz_grid=Geometry.sch_to_xyz(t_sch_grid,peg,earth_radius,earth_eccentricity)
  t_xyz_grid_rot=zeros(3,size(t_c)[1]*size(t_h)[1]*size(rot_P)[1])
  # rotate look vector around platform vector by azimuth degrees
  for i=1:size(rot_P)[1]
    q = Geometry.quat(rot_P[i], p_xyz) #create a quaternion to rotate a vector by rot_P degrees about platform position vector
    for j=1:size(t_xyz_grid)[2]
      n=(i-1)*size(t_xyz_grid)[2]+j
      t_xyz_grid_rot[:,n] = Geometry.rotate_vec(t_xyz_grid[:,j], q) #rotate target position vectors around platform position vector
    end
  end
  return t_xyz_grid_rot
end

"Converts 2D scene array of size 3xN to 3D scene array of size Ns1xNs2xNs3 which is useful for displaying tomograms\n
Inputs\n
  image_3xN: 2D scene array of size 3xN (N=Ns1xNs2xNs3) \n
  Ns1,Ns2,Ns3: lengths of scene vectors in each principle axis \n
Outputs\n
  image_3D: 3D scene array of size Ns1xNs2xNs3"
function convert_image_3xN_to_3D(image_3xN,Ns1,Ns2,Ns3)
  image_3D=zeros(Ns1,Ns2,Ns3)
  for i=1:Ns1
    for j=1:Ns2
      for k=1:Ns3
        indx=(i-1)*Ns2*Ns3+(j-1)*Ns3+k
        image_3D[i,j,k]=image_3xN[indx] # square for power?
      end
    end
  end
  return image_3D
end

end
