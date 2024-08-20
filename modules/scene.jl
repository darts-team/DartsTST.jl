module Scene
#TODO add function definitions, comments, define input types, remove unused functions

#external packages
using LinearAlgebra
using Optim
using Parameters
using StaticArrays
using Interpolations
using Plots

#local packages
include("geometry.jl")
using .Geometry

c               = 299792458
earth_radius    = 6378.137e3 # Earth semi-major axis at equator
earth_eccentricity = sqrt(0.00669437999015) # Earth eccentricity

mutable struct target_str
  #loc # target location
  #ref # target reflectivity
end

function construct_targets_str(params)
  @unpack target_pos_mode, t_loc_1, t_loc_2, t_loc_3, t_ref = params

  if target_pos_mode=="layered-grid"
    t_loc_3xN = Scene.form3Dgrid_for(t_loc_1, t_loc_2, t_loc_3) # using 3 nested for loops

    # create and flatten 3D grid of layered reflectivities from arbitrary t_ref profile
    itp = LinearInterpolation( range(1, stop=length(t_loc_3), length=length(t_ref)), t_ref )
    t_ref_interpolatedProfile = itp( range(1, stop=length(t_loc_3), length=length(t_loc_3)) ) # t_ref interpolated on the t_loc_3 axes
    t_ref_3d  = repeat( t_ref_interpolatedProfile, length(t_loc_1), length(t_loc_2), 1) # repeat the reflectivity in "horizontal S-C layers"
    t_ref_1xN = Scene.convert_3D_to_1xN(t_ref_3d)
    Nt = size(t_loc_3xN, 2) # number of targets
    #@info "Number of targets and interpolated profile" Nt, t_ref_interpolatedProfile

  elseif target_pos_mode=="layered-grid-GEDIL2"
    t_loc_3xN = Scene.form3Dgrid_for(t_loc_1, t_loc_2, t_loc_3) # using 3 nested for loops

    # create and flatten 3D grid of layered reflectivities from arbitrary t_ref profile
    t_ref_3d  = repeat( t_ref, length(t_loc_1), length(t_loc_2), 1) # repeat the reflectivity in "horizontal S-C layers"
    t_ref_1xN = Scene.convert_3D_to_1xN(t_ref_3d)
    Nt = size(t_loc_3xN, 2) # number of targets
    #@info "Number of targets and  profile" Nt, t_ref

  elseif target_pos_mode=="shaped-grid" # target positions are defined as a volumetric grid (useful for distributed target)
    @warn "Target position mode shaped-grid not implemented yet"

  elseif target_pos_mode=="grid" # target positions are defined as a volumetric grid (useful for distributed target)
    @warn "Grid target_pos_mode is deprecated, use layered-grid or shaped-grid instead"
    t_loc_3xN = Scene.form3Dgrid_for(t_loc_1, t_loc_2, t_loc_3) # using 3 nested for loops
    #t_loc_3xN=Scene.form3Dgrid_array(trg_prm.loc_1,trg_prm.loc_2,trg_prm.loc_3) # using array processing
    t_ref_1xN = Scene.convert_3D_to_1xN(t_ref)

  elseif target_pos_mode=="CR" # target positions are defined as 3xN (useful for a few discrete targets)
    t_loc_3xN = vcat(t_loc_1, t_loc_2, t_loc_3)
    t_ref_1xN = t_ref
  end


Nt = size(t_loc_3xN, 2) # number of targets
  #targets=Array{target_str}(undef,Nt)
  #for i=1:Nt;
  #  targets[i]=target_str(t_loc_3xN[:,i],t_ref_1xN[i])
  #end

  return t_loc_3xN, t_ref_1xN, Nt
end

"""
Generate Input Target Scene in 3D (scene limited by input scene arrays)
"""
function generate_input_scene_3D(targets_ref, Nt, params)
    @unpack s_loc_1, s_loc_2, s_loc_3,
            t_loc_1, t_loc_2, t_loc_3, target_pos_mode = params

    inputscene_3D=zeros(length(s_loc_1),length(s_loc_2),length(s_loc_3))
    # TODO check if target is inside the scene. gives error if target is outside the scene.
    ind_1=Int.(ones(Nt,1));ind_2=Int.(ones(Nt,1));ind_3=Int.(ones(Nt,1))
    if length(s_loc_1)>1;ind_1=round.(Int64,(t_loc_1.-s_loc_1[1])/(s_loc_1[2]-s_loc_1[1]).+1);else;end
    if length(s_loc_2)>1;ind_2=round.(Int64,(t_loc_2.-s_loc_2[1])/(s_loc_2[2]-s_loc_2[1]).+1);else;;end
    if length(s_loc_3)>1;ind_3=round.(Int64,(t_loc_3.-s_loc_3[1])/(s_loc_3[2]-s_loc_3[1]).+1);else;end
    if target_pos_mode=="grid" # TODO check if there are targets outside the scene
        ind_3xN=Int64.(Scene.form3Dgrid_for(ind_1,ind_2,ind_3))
        for i=1:Nt
            inputscene_3D[ind_3xN[1,i],ind_3xN[2,i],ind_3xN[3,i]]=abs.(targets_ref[i])
        end
    elseif target_pos_mode=="CR"
        for i=1:Nt
            inputscene_3D[ind_1[i],ind_2[i],ind_3[i]]=abs.(targets_ref[i])
        end
    end
    return inputscene_3D
end

"""
Convert target and scene coordinates to XYZ
"""
# calculate avg heading from platform positions
function convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, orbit_vel, params)
  @unpack ts_coord_sys, look_angle, left_right_look = params
  if ts_coord_sys=="LLH" # convert LLH to XYZ
      t_xyz_3xN=Geometry.geo_to_xyz(targets_loc,earth_radius,earth_eccentricity)
      s_xyz_3xN=Geometry.geo_to_xyz(s_loc_3xN,earth_radius,earth_eccentricity)
      avg_peg=[]
  elseif ts_coord_sys=="SCH" # convert SCH to XYZ
      avg_peg,p_h_avg=Geometry.avg_peg_h(orbit_pos,orbit_vel)
      slant_range,ground_range=Scene.lookangle_to_range(look_angle,p_h_avg,0,avg_peg.Ra) # slant_range (equal to ref_range?)
      targets_loc_sch=targets_loc
      if left_right_look == "left";C_dir=1;elseif left_right_look == "right";C_dir=-1;end
      targets_loc_sch[2,:]=targets_loc_sch[2,:].+C_dir*ground_range # note: this also changes targets_loc!
      t_xyz_3xN=Geometry.sch_to_xyz(targets_loc_sch,avg_peg)
      scene_loc_sch=s_loc_3xN
      scene_loc_sch[2,:]=scene_loc_sch[2,:].+C_dir*ground_range # note: this also changes s_loc_3xN
      s_xyz_3xN=Geometry.sch_to_xyz(scene_loc_sch,avg_peg)
  elseif ts_coord_sys=="XYZ" # no conversion needed
      t_xyz_3xN=targets_loc
      s_xyz_3xN=s_loc_3xN
      avg_peg=[]
  end
  return t_xyz_3xN,s_xyz_3xN,avg_peg
end

# calculate avg heading from platform velocities
function convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params)
  @unpack ts_coord_sys, look_angle, left_right_look = params

  if ts_coord_sys=="LLH" # convert LLH to XYZ
      t_xyz_3xN=Geometry.geo_to_xyz(targets_loc,earth_radius,earth_eccentricity)
      s_xyz_3xN=Geometry.geo_to_xyz(s_loc_3xN,earth_radius,earth_eccentricity)
      avg_peg=[]
  elseif ts_coord_sys=="SCH" # convert SCH to XYZ
      avg_peg,p_h_avg=Geometry.avg_peg_h(orbit_pos)
      slant_range,ground_range=Scene.lookangle_to_range(look_angle,p_h_avg,0,avg_peg.Ra) # slant_range (equal to ref_range?)
      targets_loc_sch=targets_loc
      if left_right_look == "left";C_dir=1;elseif left_right_look == "right";C_dir=-1;end
      targets_loc_sch[2,:]=targets_loc_sch[2,:].+C_dir*ground_range # note: this also changes targets_loc!
      t_xyz_3xN=Geometry.sch_to_xyz(targets_loc_sch,avg_peg)
      scene_loc_sch=s_loc_3xN
      scene_loc_sch[2,:]=scene_loc_sch[2,:].+C_dir*ground_range # note: this also changes s_loc_3xN
      s_xyz_3xN=Geometry.sch_to_xyz(scene_loc_sch,avg_peg)
  elseif ts_coord_sys=="XYZ" # no conversion needed
      t_xyz_3xN=targets_loc
      s_xyz_3xN=s_loc_3xN
      avg_peg=[]
  end
  return t_xyz_3xN,s_xyz_3xN,avg_peg
end

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
function slantrange_to_lookangle(ra,rs,p_h,t_h) # t_h: target height
  ra=ra.+t_h
  p_h=p_h.-t_h
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
  rs=(((ra+p_h)*sin.(α)).^2+((ra+p_h).*cos.(α)-ra).^2).^0.5 # slant range
  θ_l=acos.((rs.^2+(ra+p_h).^2-ra^2)./(2*rs.*(ra+p_h))) # rad look angle
  return rs,θ_l*180/pi
end


"""
Convert look angle to incidence angle
## Inputs
- θ_l: look angles (deg)
- p_h: platform height (m)
- ra: radius of spherical planet (m)
- t_h: target heights vector (m)
## Outputs
- inc: incidence angles (deg)
"""
function lookangle_to_incangle(θ_l,p_h,t_h, ra)
  ra=ra.+t_h
  p_h=p_h.-t_h
  inc=asind.(sind.(θ_l).*(ra+p_h)./ra) 
  return inc
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
          image_3D[i,j,k]=image_1xN[indx]
        end
      end
    end
    return image_3D
  end

  """
  Converts 3D array of size Ns1xNs2xNs3 to 1D array of size 1xN which is used to convert target amplitudes defined in 3D to 1xN for generate_raw_data
    ## Arguments
    - array_3D: 3D array of size Ns1xNs2xNs3
    ## Outputs
    - array_1xN: 1D array of size 1xN (N=Ns1xNs2xNs3)
    """
  function convert_3D_to_1xN(array_3D)
    N1,N2,N3=size(array_3D)
    N=N1*N2*N3
    array_1xN=zeros(Float64,1,N) ##?????
    for i=1:N1
      for j=1:N2
        for k=1:N3
          indx=(i-1)*N2*N3+(j-1)*N3+k
          array_1xN[indx]=array_3D[i,j,k]
        end
      end
    end
    return array_1xN
  end

  #TODO add function definition
  function take_1D_cuts(image_3D, params)
    @unpack s_loc_1, s_loc_2, s_loc_3, t_loc_1, t_loc_2, t_loc_3, res_dB, PSF_image_point, PSF_cuts, PSF_direction, PSF_image_point, display_1D_cuts = params
    target_location = [t_loc_1 t_loc_2 t_loc_3[1]]

    if PSF_cuts == 1
      scene_axis11=s_loc_1;scene_axis22=s_loc_2;scene_axis33=s_loc_3
      image_1D_1, image_1D_2, image_1D_3 = obtain_1D_slices(image_3D, target_location, s_loc_1, s_loc_2, s_loc_3, PSF_image_point)
      # Calculate Scene Resolutions and Plot 1D cuts
      if length(image_1D_1)>1
        scene_res1=s_loc_1[2]-s_loc_1[1] # scene resolution along the 1st axis
        if display_1D_cuts
          display(plot(s_loc_1,20*log10.(image_1D_1/maximum(image_1D_1)),xaxis=("scene axis 1 in scene units"),ylabel=("amplitude (dB)"),size=(1600,900),leg=false)) # plot the cut along axis 1
        end
      else;scene_res1=NaN;end
      if length(image_1D_2)>1
        scene_res2=s_loc_2[2]-s_loc_2[1] # scene resolution along the 2nd axis
        if display_1D_cuts
           display(plot(s_loc_2,20*log10.(image_1D_2/maximum(image_1D_2)),xaxis=("scene axis 2 in scene units"),ylabel=("amplitude (dB)"),size=(1600,900),leg=false)) # plot the cut along axis 2
        end
      else;scene_res2=NaN;end
      if length(image_1D_3)>1
        scene_res3=s_loc_3[2]-s_loc_3[1] # scene resolution along the 3rd axis
        if display_1D_cuts
          display(plot(s_loc_3,20*log10.(image_1D_3/maximum(image_1D_3)),xaxis=("scene axis 3 in scene units"),ylabel=("amplitude (dB)"),size=(1600,900),leg=false)) # plot the cut along axis 3
        end
      else;scene_res3=NaN;end
      scene_res=[scene_res1 scene_res2 scene_res3]
    elseif PSF_cuts == 2 # tilted cut is taken from the scene center
      image_1D_1, scene_axis11, scene_axis22, scene_axis33 = obtain_1D_slice_tilted(image_3D, s_loc_1, s_loc_2, s_loc_3, PSF_direction)
      image_1D_2=NaN
      image_1D_3=NaN
      # Plot tilted 1D cut
      scene_res=((scene_axis11[2]-scene_axis11[1])^2+(scene_axis22[2]-scene_axis22[1])^2+(scene_axis33[2]-scene_axis33[1])^2)^0.5 # scene resolution along the PSF direction
      scene_axis=(0:scene_res:(length(image_1D_1)-1)*scene_res).-(length(image_1D_1)-1)*scene_res/2
      if display_1D_cuts
        display(plot(scene_axis,20*log10.(abs.(image_1D_1)/maximum(abs.(image_1D_1))),xaxis=("scene axis along specified direction"),ylim=(-50,0),ylabel=("amplitude (dB)"),size=(900,900),leg=false)) # plot the tilted cut
      end
    end

    return scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res
  end

  #TODO add function definition
  function obtain_1D_slice_tilted(image_3D,scene_axis1,scene_axis2,scene_axis3,PSF_direction)
      Ns_1=length(scene_axis1);Ns_2=length(scene_axis2);Ns_3=length(scene_axis3)
      slice_index1=Int(ceil(Ns_1/2));slice_index2=Int(ceil(Ns_2/2));slice_index3=Int(ceil(Ns_3/2))
      xc=scene_axis1[slice_index1];yc=scene_axis2[slice_index2];zc=scene_axis3[slice_index3]
      x2=PSF_direction[1];y2=PSF_direction[2];z2=PSF_direction[3]
      # 3D line equation: (x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1) and x1=0;y1=0;z1=0 is the origin at center of scene
      scene_axis10=scene_axis1.-xc
      scene_axis20=scene_axis2.-yc
      scene_axis30=scene_axis3.-zc
      if x2!=0
          scene_axis1_int=scene_axis10
          scene_axis2_int=scene_axis1_int*(y2/x2)
          scene_axis3_int=scene_axis1_int*(z2/x2)
      elseif y2!=0
          scene_axis2_int=scene_axis20
          scene_axis1_int=scene_axis2_int*(x2/y2)
          scene_axis3_int=scene_axis2_int*(z2/y2)
      elseif z2!=0
          scene_axis3_int=scene_axis30
          scene_axis1_int=scene_axis3_int*(x2/z2)
          scene_axis2_int=scene_axis3_int*(y2/z2)
      else;println("The PSF direction relative to center/origin cannot be [0,0,0]!");end
      if Ns_1>1;scene_axis1_ind=round.(Int64,(scene_axis1_int.-scene_axis10[1])/(scene_axis10[2]-scene_axis10[1]).+1);else;image_3D=dropdims(image_3D,dims=1);end
      if Ns_2>1;scene_axis2_ind=round.(Int64,(scene_axis2_int.-scene_axis20[1])/(scene_axis20[2]-scene_axis20[1]).+1);else;image_3D=dropdims(image_3D,dims=2);end
      if Ns_3>1;scene_axis3_ind=round.(Int64,(scene_axis3_int.-scene_axis30[1])/(scene_axis30[2]-scene_axis30[1]).+1);else;image_3D=dropdims(image_3D,dims=3);end
      itp=interpolate(image_3D,BSpline(Cubic(Free(OnGrid()))))
      image_1D=zeros(Float64,length(scene_axis1_int))
      for i=1:length(scene_axis1_int)
          if Ns_1>1 && Ns_2>1 && Ns_3>1
              if scene_axis1_ind[i]>0 && scene_axis1_ind[i]<=Ns_1 && scene_axis2_ind[i]>0 && scene_axis2_ind[i]<=Ns_2 && scene_axis3_ind[i]>0 && scene_axis3_ind[i]<=Ns_3
                  image_1D[i]=itp(scene_axis1_ind[i],scene_axis2_ind[i],scene_axis3_ind[i])
              else;image_1D[i]=NaN;end
          elseif Ns_1==1
              if scene_axis2_ind[i]>0 && scene_axis2_ind[i]<=Ns_2 && scene_axis3_ind[i]>0 && scene_axis3_ind[i]<=Ns_3
                  image_1D[i]=itp(scene_axis2_ind[i],scene_axis3_ind[i])
              else;image_1D[i]=NaN;end
          elseif Ns_2==1
              if scene_axis1_ind[i]>0 && scene_axis1_ind[i]<=Ns_1 && scene_axis3_ind[i]>0 && scene_axis3_ind[i]<=Ns_3
                  image_1D[i]=itp(scene_axis1_ind[i],scene_axis3_ind[i])
              else;image_1D[i]=NaN;end
          elseif Ns_3==1
              if scene_axis2_ind[i]>0 && scene_axis2_ind[i]<=Ns_2 && scene_axis1_ind[i]>0 && scene_axis1_ind[i]<=Ns_1
                  image_1D[i]=itp(scene_axis1_ind[i],scene_axis2_ind[i])
              else;image_1D[i]=NaN;end
          end
      end
      image_1D=deleteat!(vec(image_1D),findall(isnan,vec(image_1D)))
      scene_axis11=scene_axis1_int.+xc
      scene_axis22=scene_axis2_int.+yc
      scene_axis33=scene_axis3_int.+zc
      #NaN_ind=findall(image_1D==NaN)
      return image_1D,scene_axis11,scene_axis22,scene_axis33
  end

  #TODO add function definition
  function obtain_1D_slices(image_3D,target_location,scene_axis1,scene_axis2,scene_axis3,PSF_peak_target)
      if PSF_peak_target==1 # 1D slices from PSF peak
          max_ind=findall(image_3D .==maximum(image_3D))
          slice_index1=max_ind[1][1]
          slice_index2=max_ind[1][2]
          slice_index3=max_ind[1][3]
      elseif PSF_peak_target==2 # PSF slices from actual target location
          slice_index1=findall(target_location[1] .==scene_axis1)
          slice_index2=findall(target_location[2] .==scene_axis2)
          slice_index3=findall(target_location[3] .==scene_axis3)
      elseif PSF_peak_target==3 # PSF slices from 3D scene center
          Ns_1=length(scene_axis1)
          Ns_2=length(scene_axis2)
          Ns_3=length(scene_axis3)
          slice_index1=Int(ceil(Ns_1/2))
          slice_index2=Int(ceil(Ns_2/2))
          slice_index3=Int(ceil(Ns_3/2))
      end
      if PSF_peak_target==2 && (isempty(slice_index1) || isempty(slice_index2) || isempty(slice_index3))
          image_1D_1=NaN;image_1D_2=NaN;image_1D_3=NaN
      else
          gr() # or plotly()
          if length(scene_axis1)>1
              image_slice=image_3D[:,slice_index2,slice_index3]
              image_1D_1=zeros(Float64,length(image_slice))
              image_1D_1[:]=image_slice
          else;image_1D_1=NaN;end #println("PSF metrics along 1st dimension cannot be calculated since image has no 1st dimension.");end
          if length(scene_axis2)>1
              image_slice=image_3D[slice_index1,:,slice_index3]
              image_1D_2=zeros(Float64,length(image_slice))
              image_1D_2[:]=image_slice
          else;image_1D_2=NaN;end #println("PSF metrics along 2nd dimension cannot be calculated since image has no 2nd dimension.");end
          if length(scene_axis3)>1
              image_slice=image_3D[slice_index1,slice_index2,:]
              image_1D_3=zeros(Float64,length(image_slice))
              image_1D_3[:]=image_slice
          else;image_1D_3=NaN;end #println("PSF metrics along 3rd dimension cannot be calculated since image has no 3rd dimension.");end
      end
      return image_1D_1,image_1D_2,image_1D_3
  end

  """
  Compute target points along an iso-range line on the ellipsoid
  # Arguments
  - `pos::3x1 Float Array`: Cartesian Position of Spacecraft (typically ECEF)
  - `vel::3x1 Float Array`: Cartesian Velocity of Spacecraft (typically ECEF)
  - `ϕ::1xN Float Array`: azimuth angles (rad)
  - `ρ_0::Float64`: slant range (m)
  - `θ_0::Float64`: loook angle associated with ρ_0 and ϕ_0 (rad)
  - `earth_radius::Float64`: radius of earth (optional)
  - `earth_eccentricity::Float64`: eccentricity of earth (optional)

  # Return
  - `Targets::3xN Flot Array`: Targets along iso-range line in SC(pos/vel) parent frame (typically ECEF)
  """
  function isorange_ellipsoid(pos::Array{Float64,1}, vel::Array{Float64,1},
    ϕ_in::Array{Float64,1},ρ_0::Float64, θ_0::Float64,
    ellipsoid_radius::Float64=6.378137e6,ellipsoid_eccentricity::Float64=0.08181919084267456)

    @assert size(pos,1)==3 "POS needs to be 3 x 1"
    @assert size(vel,1)==3 "VEL needs to be 3 x 1"
    @assert ρ_0 >0 "Slant Range needs to be positive"

    #get TCN frame from position and velocity
    that,chat,nhat = Geometry.get_tcn(pos, vel)

    #set aside space for output variable
    Targ_vec = zeros(3,length(ϕ_in))

    #iterate over various ϕ values, est theta for that that suits
    for i_ϕ = 1:length(ϕ_in)

      #define function that is difference between slant-range at (θ, ϕ_in) and slant_range at (θ_in, 0)
      mylhatf(θ) = (Geometry.get_rho(pos, Geometry.tcn_lvec(that, chat, nhat, θ, ϕ_in[i_ϕ]), ellipsoid_radius, ellipsoid_eccentricity) - ρ_0)^2

      #find minimum of that function by iterating over look angle
      res = optimize(mylhatf, θ_0 - π/180, θ_0 + π/180, Brent(), rel_tol = 1e-10);
      θ_est = Optim.minimizer(res); #estimated look angle

      #estimate slant range for this look vector
      ρ_est = Geometry.get_rho(pos, Geometry.tcn_lvec(that, chat, nhat, θ_est, ϕ_in[i_ϕ]), ellipsoid_radius, ellipsoid_eccentricity)

      #compute Target coordinates (using T = P + l)
      Targ_vec[:,i_ϕ] = pos + Geometry.tcn_lvec(that, chat, nhat, θ_est, ϕ_in[i_ϕ])*ρ_est;
    end

    return Targ_vec
  end
  """
  Compute target points along an iso-range line on the ellipsoid (at various heights)
  # Arguments
  - `pos::3x1 Float Array`: Cartesian Position of Spacecraft (typically ECEF)
  - `vel::3x1 Float Array`: Cartesian Velocity of Spacecraft (typically ECEF)
  - `ϕ::1xNϕ Float Array`: azimuth angles (rad)
  - `hgt::1xNh Float Array`: heights above ellipsoid (m)
  - `ρ_0::Float64`: slant range (m)
  - `θ_0::Float64`: loook angle associated with ρ_0 and ϕ_0 (rad)
  - `earth_radius::Float64`: radius of earth (optional)
  - `earth_eccentricity::Float64`: eccentricity of earth (optional)

  # Return
  - `Targets::[3 x Nϕ x Nh] Float Array`: Targets along iso-range line in SC(pos/vel) parent frame (typically ECEF)
  """
  function isorange_ellipsoid(pos::Array{Float64,1}, vel::Array{Float64,1},
    ϕ_in::Array{Float64,1},hgt_in::Array{Float64,1}, ρ_0::Float64, θ_0::Float64,
    ellipsoid_radius::Float64=6.378137e6,ellipsoid_eccentricity::Float64=0.08181919084267456)

    @assert size(pos,1)==3 "POS needs to be 3 x 1"
    @assert size(vel,1)==3 "VEL needs to be 3 x 1"
    @assert ρ_0 >0 "Slant Range needs to be positive"

    #get TCN frame from position and velocity
    that,chat,nhat = Geometry.get_tcn(pos, vel)

    #look vector at zero azimuth, zero height
    lhat_0 = Geometry.tcn_lvec(that, chat, nhat, θ_0, 0.0); #unit look vector
    lvec_0 = lhat_0*ρ_0; #full look vector

    #set aside space for output variable
    Targ_vec = zeros(3,length(ϕ_in), length(hgt_in))

    #iterate over various ϕ values, est theta for that that suits
    for i_hgt = 1:length(hgt_in)

      #rotate look_vec about t_hat to get target at ith hgt
      qr = Geometry.quat(-hgt_in[i_hgt]/ρ_0*180/π, that)
      lvec_hgt = Geometry.rotate_vec(lvec_0, qr)

      #compute target vector at ith height
      Targ_hgt = pos + lvec_hgt;

      #Geodetic coordinates at zero azimuth
      targ_geo_hgt = Geometry.xyz_to_geo(Targ_hgt, ellipsoid_radius,ellipsoid_eccentricity);

      #get targets for various phi
      # add height rotation to θ_0 and Target geodetic height to ellipsoid
      # keep slant range (ρ_0) constant
      Targ_vec[:,:,i_hgt] = isorange_ellipsoid(pos, vel, ϕ_in,
          ρ_0, θ_0+hgt_in[i_hgt]/ρ_0, ellipsoid_radius+targ_geo_hgt[3], ellipsoid_eccentricity);
    end

    return Targ_vec
  end
  """
  Compute target points along an iso-range line on the ellipsoid (at various heights and look angles)
  # Arguments
  - `pos::3x1 Float Array`: Cartesian Position of Spacecraft (typically ECEF)
  - `vel::3x1 Float Array`: Cartesian Velocity of Spacecraft (typically ECEF)
  - `ϕ::1xNϕ Float Array`: azimuth angles (rad)
  - `hgt::1xNh Float Array`: heights above ellipsoid (m)
  - `θ::1xNθ Float64 Array`: loook angles to evaluate (rad)
  - `earth_radius::Float64`: radius of earth (optional)
  - `earth_eccentricity::Float64`: eccentricity of earth (optional)

  # Return
  - `Targets::[3 x Nϕ x Nh x Nθ] Float Array`: Targets along iso-range line in SC(pos/vel) parent frame (typically ECEF)
  """
  function isorange_ellipsoid(pos::Array{Float64,1}, vel::Array{Float64,1},
    ϕ_in::Array{Float64,1},hgt_in::Array{Float64,1}, θ_in::Array{Float64,1},
    ellipsoid_radius::Float64=6.378137e6,ellipsoid_eccentricity::Float64=0.08181919084267456)

    @assert size(pos,1)==3 "POS needs to be 3 x 1"
    @assert size(vel,1)==3 "VEL needs to be 3 x 1"

    #get TCN frame from position and velocity
    that,chat,nhat = Geometry.get_tcn(pos, vel)

    #set aside space for output variable
    Targ_vec = zeros(3,length(ϕ_in), length(hgt_in), length(θ_in))

    #iterate over various ϕ values, est theta for that that suits
    for i_θ = 1:length(θ_in)

      #look vector at zero azimuth, zero height and the ith look angle
      lhat_0 = Geometry.tcn_lvec(that, chat, nhat, θ_in[i_θ], 0.0); #unit look vector
      # slant range to ellipsoid for ith look angle
      ρ_0 = Geometry.get_rho(pos, lhat_0, ellipsoid_radius,ellipsoid_eccentricity);

      # get targets for various ϕ and heights
      Targ_vec[:,:,:,i_θ] = isorange_ellipsoid(pos, vel, ϕ_in, hgt_in, ρ_0, θ_in[i_θ],
          ellipsoid_radius, ellipsoid_eccentricity);
    end

    return Targ_vec
  end



end #end module
