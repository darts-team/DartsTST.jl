module Scene
#TODO add function definitions, comments, define input types, remove unused functions

#include("geometry.jl")
using ..Geometry
using LinearAlgebra

"""
Generates 3D (volumetric) grid from three 1D (linear) arrays using nested for-loops
## Arguments
  - array1,1D arrays of size N1 representing first principal axes of the 3D volume
  - array2,1D arrays of size N2 representing second principal axes of the 3D volume
  - array3:1D arrays of size N3 representing third principal axes of the 3D volume
## Outputs
  vol_grid: coordinates of 3D (volumetric) grid of points as 3xN matrix where N=N1xN2xN3
"""
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

"""
Generates 3D (volumetric) grid from three 1D (linear) arrays using array processing
## Arguments
  - array1,1D arrays of size N1 representing first principal axes of the 3D volume
  - array2,1D arrays of size N2 representing second principal axes of the 3D volume
  - array3:1D arrays of size N3 representing third principal axes of the 3D volume
## Outputs
  vol_grid: coordinates of 3D (volumetric) grid of points as 3xN matrix where N=N1xN2xN3
"""
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

"""
Convert look angle to slant/ground range (spherical planet assumed)
## Inputs
  - θ_l: look angles (deg)
  - p_h: platform height (m)
  - ra: radius of spherical planet (m)
  - t_h: target heights vector (m)
## Outputs
  - rs: slant ranges to the targets
  - rg: ground ranges to the targets
"""
function lookangle_to_range(θ_l,p_h,t_h, ra)
    ra=ra.+t_h
    p_h=p_h.-t_h
    inc=asind.(sind.(θ_l).*(ra+p_h)./ra) # deg incidence angle
    α=inc-θ_l # deg planet-central angle
    rg=ra.*α*pi/180 # ground range
    rs=(((ra+p_h).*sind.(α)).^2+((ra.+p_h).*cosd.(α).-ra).^2).^0.5 # slant range
    return rs, rg
end

"""
Slant range to look angle and ground range conversion (spherical planet assumed)
## Inputs
  - p_h: platform height
  - ra: radius of spherical planet
  - rs: slant ranges to the targets
## Outputs
  - θ_l: look angles to the targets
  - rg: ground ranges to the targets
"""
function slantrange_to_lookangle(ra,rs,p_h) # target height is assumed 0 TODO add target height
  θ_l=acos.((rs.^2+(ra+p_h).^2-ra^2)./(2*rs.*(ra+p_h))) # rad look angle
  inc=asin.((ra+p_h)./ra.*sin.(θ_l)) # rad incidence angle
  α=inc-θ_l # rad planet-central anglesin
  rg=ra*α # ground range
  return rg,θ_l*180/pi
end

"""
Ground range to look angle and slant range conversion (spherical planet assumed)
## Arguments
  p_h: platform height
  ra: radius of spherical planet
  rg: ground ranges to the targets
## Outputs
  - θ_l: look angles to the targets
  - rs: slant ranges to the targets
"""
function groundrange_to_lookangle(ra,rg,p_h) # # target height is assumed 0 TODO add target height
  α=rg/ra # rad planet-central angle
  rs=(((ra+p_h)*sind.(α)).^2+((ra+p_h).*cosd.(α)-ra).^2).^0.5 # slant range
  θ_l=acos.((rs.^2+(ra+p_h).^2-ra^2)./(2*rs.*(ra+p_h))) # rad look angle
  return rs,θ_l*180/pi
end

"""
Convert Azimuth/Elevation look angles and target heights to target xyz (spherical approximation)
## Arguments
  - lookh: 3xN array for N point targets (1st dimension: azimuth look angles [deg], 2nd dimension: elevation look angles [deg], 3rd dimension: target heights [m])
  - platform_geo: platform location in geographic coordinates (lat/lon/height)
  - peg: peg point struct (see Geometry.PegPoint for details)
  - earth_radius: radius of earth (optional)
  - earth_eccentricity: eccentricity of earth (optional)
## Outputs
  - T_xyz: vector of target positions in xyz
"""
function lookh_to_xyz(lookh,platform_geo,peg::Geometry.PegPoint,earth_radius::Float64=6.378137e6,earth_eccentricity::Float64=0.08181919084267456)
    vL2_sch=zeros(size(lookh))
    peg_xyz_grid=zeros(size(lookh))
    vL_xyz=zeros(size(lookh))
    vT=zeros(size(lookh))
    ϕ_l=lookh[1,:] # deg  azimuth angles
    θ_l=lookh[2,:] # deg elevation angles
    t_h=lookh[3,:] # target heights
    p_h=platform_geo[3] # platform height
    p_xyz=Geometry.geo_to_xyz(platform_geo,earth_radius,earth_eccentricity) # platform position in xyz
    rs,rg=lookangle_to_range(θ_l,p_h,t_h, earth_radius) # slant range between target and platform, assumes spherical planet TODO t_h can be a smaller vector and then repeat
    vL2_sch[1,:]=sind.(θ_l).*sind.(ϕ_l)
    vL2_sch[2,:]=sind.(θ_l).*cosd.(ϕ_l) #TODO add d for left/right looking
    vL2_sch[3,:]=-cosd.(θ_l) # look vectors in geo (for each target at each look angle
    #peg_point = Geometry.PegPoint(peg[1], peg[2], peg[3],earth_radius,earth_eccentricity) # create peg point struct
    vL2_xyz=Geometry.sch_to_xyz(vL2_sch,peg) # look vectors in xyz (for each target at each look angle)
    peg_geo=[peg.pegLat,peg.pegLon,0]
    peg_xyz=Geometry.geo_to_xyz(peg_geo,earth_radius,earth_eccentricity)
    peg_xyz_grid=repeat(peg_xyz,1,size(lookh,2))
    rs_grid=repeat(rs',3,1)
    vL_xyz=rs_grid.*(vL2_xyz.-peg_xyz_grid)
    vT_xyz=vL_xyz.+p_xyz # target vectors in xyz (for azimuth angle of 0 deg)
    return vT_xyz
end

"""
Converts target positions defined in chP (C,H- of SCH) and azimuth angle vector around platform position vector to ECEF xyz
## Arguments
  - t_c: cross-track coordinates of target points
  - t_h: height coordinates of target points
  - rot_P: amount of rotations (as an array) in degrees of look vector about platform position vector
  - platform_geo: platform location in geographic coordinates (lat/lon/height)
  - peg: peg point coordinates for sch to xyz conversion (see Geometry.PegPoint)
  - earth_radius: radius of earth (optional)
  - earth_eccentricity: eccentricity of earth (optional)
## Outputs
  - XYZ-grid: target positions in ECEF-xyz as a 3xN array
"""
function chP_to_xyz_grid(t_c,t_h,rot_P,platform_geo,peg::Geometry.PegPoint,earth_radius::Float64=6.378137e6,earth_eccentricity::Float64=0.08181919084267456)
  t_s=0 # before rotation around P, targets are at the same along-track position with platform
  p_xyz=Geometry.geo_to_xyz(platform_geo,earth_radius,earth_eccentricity) # platform position in xyz
  p_xyz=p_xyz/norm(p_xyz)
  t_sch_grid=Scene.form3Dgrid_for(t_s,t_c,t_h) # using 3 nested for loops
  #t_geo_grid=Scene.form3Dgrid_array(t_s,t_c,t_h) # using array processing
  t_xyz_grid=Geometry.sch_to_xyz(t_sch_grid,peg)
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

"""
Converts 1D scene array of size 1xN to 3D scene array of size Ns1xNs2xNs3 which is useful for displaying tomograms
## Arguments
  - image_1xN: 1D scene array of size 1xN (N=Ns1xNs2xNs3)
  - Ns1,Ns2,Ns3: lengths of scene vectors in each principle axis
## Outputs
  - image_3D: 3D scene array of size Ns1xNs2xNs3
"""
function convert_image_1xN_to_3D(image_1xN,Ns1,Ns2,Ns3)
  image_3D=zeros(Ns1,Ns2,Ns3)
  for i=1:Ns1
    for j=1:Ns2
      for k=1:Ns3
        indx=(i-1)*Ns2*Ns3+(j-1)*Ns3+k
        image_3D[i,j,k]=image_1xN[indx] # square for power?
      end
    end
  end
  return image_3D
end

""" Compute look vector based based on TCN basis vectors,
 # Arguments
  - t-hat: 3x1 definition of the t_hat (tangential velocity) vector in base frame (usually XYZ)
  - c-hat: 3x1 definition of the c_hat (cross-track) vector in base frame (usually XYZ)
  - n-hat: 3x1 definition of the n_hat (nadir) vector in base frame (usually XYZ)
  - θ-rng: scalar look angle (as measured from the nadir vector) [radians]
  - ϕ-az:  azimuth angle (angle projected on the T-C plane measured from the C-vector) [radians]
# Return
  + 3x1 normalized look vector in base frame (usually XYZ)
- to get TCN frame  see modules/geometry.jl->get_tcn()
"""
lookvec_fromTCN(t_hat,c_hat,n_hat,θ_rng,ϕ_az) =  sin(θ_rng) * sin(ϕ_az) * t_hat + sin(θ_rng)* cos(ϕ_az) * c_hat + cos(θ_rng) * n_hat

end
