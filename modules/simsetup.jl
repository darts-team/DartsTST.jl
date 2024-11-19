module SimSetup

#external packages
using ReferenceFrameRotations
using LinearAlgebra
using Parameters
using Statistics
using Distributions

#local packages
using ..Geometry
using ..Antenna
using ..Scene
using ..DEM


"""
Antenna structure. Contains pattern, rotation and origin wrt to SC
# Usage example
    - ant = SimSetup.sc_ant(pattern, rot = Geometry.quat(-90, [0,0,1]), org=[0.,0.,0.])
# Arguments
    - `pattern::Antenna.AntGrid or Antenna.AntCut`, pattern from netCDF file
    - `kw(rot)::Quaternion` rotation quaternion describing rotation from SC (IJK) to Pattern frame
    - `kw(org)::Array` origin of antenna frame in SC frame
"""
mutable struct sc_ant
    pattern
    rot # rotation of antenna wrt the SpaceCraft frame
    org #origin relative to spacecraft frame
    function sc_ant(pattern; rot=Geometry.quat(-90, [0,0,1]), org=[0.,0.,0.])
        new(pattern, rot, org)
    end
end
"""
SpaceCraft structure. Contains sc position, velocity, orientation and antenna structure
# Usage example
    - ant = SimSetup.spacecraft(pos, vel, ant = pattern, rot = Geometry.quat(-90, [0,0,1]))
# Arguments
    - `position::3xN Array`, orbit position (ECEF, m)
    - `velocity::3xN Array`, orbit velocity (ECEF, m/s)
    - `kw(ant)::sc_ant`, spacecraft antenna structure
    - `kw(rot)::Quaternion` rotation quaternion describing rotation from ECEF to SC (nominally TCN)
"""
mutable struct spacecraft
    pos #position in ECEF (m,m,m)
    vel #velocity relative in ECEF (m/s,m/s, m/s)
    ant # antenna frame defined relative to sc frame
    rot # with respect to the planet
    look_angle # look angle in degrees (+ve )
    side # which side is the spacecraft looking
    function spacecraft(pos, vel; ant = [], rot=Geometry.tcn_quat(pos, vel), look_angle = 35, side = "left")
        @assert size(rot,1) == size(pos,2) "SC rotation quaternion must be N-element vector for a position vector of size 3xN"
        @assert size(pos) == size(vel) "SC position and velocity must have same size"
        @assert side == "left" || side == "l" || side == "right" || side == "r" "side should either be left, l, right, r"

        #set up the antenna with correct orientation
        lsgn(sd) = (sd == "l" || sd == "left") ? 1 : -1
        q_ant_look = Geometry.quat(lsgn(side)*look_angle, [0,1,0]) #rotate about antenna y-hat (the long side) to the correct look angle
        ant = rotate_antenna(ant, q_ant_look)

        #create structure
        new(pos, vel, ant, rot, look_angle, side)
    end
end

"""
Rotate Spacecraft about its native axes. Nominally SC is aligned with TCN.
Rotation about T is equivalent to roll, rotation about C is pitch, rotation about N is yaw.
 # Usage
    - sc = rotate_spacecraft(sc, Geometry.quat(10, [0,1,0]))
    - rotate_spacecraft!(sc, Geometry.quat(10, [0,1,0]))
    will rotate sc by 10 degrees about it's y-axis (or C-axis)
    equivalent to 10 degree pitch

# Arguments
    - `sc::spacecraft` SimSetup.spacecraft structure
    - `quat::Quaternion` option1: quaternion describing rotation
    - `rpy::Array{Float64,2}` option 2: RPY vector [r,p,y] x Npts
# Output
    - `sc::spacecraft` SimSetup.spacecraft structure with modified rotation
"""
function rotate_spacecraft(sc::spacecraft, quat::Quaternion)
    sc_r = deepcopy(sc);
    sc_r.rot =  sc.rot.*[quat]
    return sc_r
end
function rotate_spacecraft(sc::spacecraft, quat::Array{Quaternion{Float64},1})
    @assert length(quat) == length(sc.rot)
    sc_r = deepcopy(sc);
    sc_r.rot =  sc.rot.*quat
    return sc_r
end
function rotate_spacecraft(sc::spacecraft, rpy::Array{Float64,})
    @assert size(rpy,1) == 3 "RPY must be a 3xN vector"
    @assert size(rpy,2) == size(sc.pos,2)

    #construct quat
    quat = Array{Quaternion{Float64},1}(undef, size(sc.pos,2))
    for itp = 1:size(sc.pos,2)
        quat[itp] = Geometry.quat(rpy[1,itp], [1, 0, 0]);
        quat[itp] = quat[itp]*Geometry.quat(rpy[2,itp], [0, 1, 0]);
        quat[itp] = quat[itp]*Geometry.quat(rpy[3,itp], [0, 0, 1]);
    end
    return rotate_spacecraft(sc, quat)
end
function rotate_spacecraft!(sc::spacecraft, quat::Quaternion )
    sc.rot =  sc.rot.*[quat]
end
function rotate_spacecraft!(sc::spacecraft, quat::Array{Quaternion{Float64},1})
    @assert length(quat) == length(sc.rot)
    sc.rot =  sc.rot.*quat
end
function rotate_spacecraft!(sc::spacecraft, rpy::Array{Float64,})
    @assert size(rpy,1) == 3 "RPY must be a 3xN vector"
    @assert size(rpy,2) == size(sc.pos,2)

    #construct quat
    quat = Array{Quaternion{Float64},1}(undef, size(sc.pos,2))
    for itp = 1:size(sc.pos,2)
        quat[itp] = Geometry.quat(rpy[1,itp], [1, 0, 0]);
        quat[itp] = quat[itp]*Geometry.quat(rpy[2,itp], [0, 1, 0]);
        quat[itp] = quat[itp]*Geometry.quat(rpy[3,itp], [0, 0, 1]);
    end
    rotate_spacecraft!(sc, quat)
end

"""
Rotate Antenna about its native axes. Antenna Y-axis [0,1,0] is aligned
with the long-side of the antenna, while X-axis [1,0,0] is aligned with
the short-side of the antenna. Antenna Z-axis is the boresight.
 # Usage
    - sc_ant = rotate_antenna(sc_ant, Geometry.quat(90, [1,0,0]))
    - rotate_antenna!(sc_ant, Geometry.quat(90, [1,0,0]))
    will rotate sc_ant by 90 degrees about it's x-axis (short side)

# Arguments
    - `ant::sc_ant` SimSetup.sc_ant structure
    - `quat::Quaternion` quaternion describing rotation

# Output
    - `ant::sc_ant` SimSetup.sc_ant structure with modified rotation
"""
function rotate_antenna(ant, quat::Quaternion)
    ant_r = deepcopy(ant);
    ant_r.rot =  ant.rot*quat
    return ant_r
end
function rotate_antenna!(ant, quat::Quaternion )
    ant.rot =  ant.rot*quat
end

"""
Interpolate antenna pattern given spacecraft structure and a set of target points
 # Usage
    - cp,xp = interpolate_pattern(sc, targets)
    cp (co-pol) and xp (cross-pol) antenna patterns

# Arguments
    - `sc::spacecraft` SimSetup.spacecraft structure
    - `targ::3xN float array` target points in ECEF xyz
    - `ind::Int` index of the sc position, quaternion to use
    - `freq::Float64`, [Optional] antenna freqency to use, in GHz

# Output
    - `cp` co-pol antenn pattern of size 1xN
    - `xp` cross-pol antenna pattern of size 1xN
"""
function interpolate_pattern(sc, targ::Array{Float64,}, ind::Int = 1, freq::Float64=1.25)
    @assert size(targ,1) == 3 "targ_xyz must be a 3xN vector"
    @assert ind <= size(sc.pos,2) "index must be less than the size of the spacecraft position/orientation array"

    #create look vector from targets and spacecraft
    lvec = targ .- sc.pos;
    lhat = lvec./norm.(eachcol(lvec))';

    #compute the quaternion to rotate from ECEF to antenna frame
    q_ecef_ant = sc.rot[ind] * sc.ant.rot


    #project look vector onto the antenna frame
    lhat_ant = zeros(3,size(lhat,2))
    for ii=1:size(lhat,2)
        lhat_ant[:,ii] = Geometry.rotate_frame(lhat[:,ii], q_ecef_ant);
    end

    #inteprolate pattern
    cp,xp = Antenna.interpolate_pattern(sc.ant.pattern, lhat_ant, freq);
    return cp,xp
end


function applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, t_xyz_3xN, params)
    @unpack ts_coord_sys, t_loc_1, t_loc_2, t_loc_3, user_defined_orbit, look_angle, antennaFile, target_pos_mode = params
 
    if user_defined_orbit > 0
       @warn "Antenna pattern not tested with custom orbits. Skipping."
       return
    end
 
    avg_p_xyz=reshape(mean(mean(p_xyz,dims=2),dims=3),3)
    avg_p_vel=reshape(mean(mean(orbit_vel,dims=2),dims=3),3)
    if ts_coord_sys=="SCH"
        look_ang=look_angle
    elseif ts_coord_sys=="XYZ" || ts_coord_sys=="LLH"
        platform_heights=zeros(Np);slant_ranges=zeros(Np)
        avg_t_xyz=mean(t_xyz_3xN,dims=2) # average target location in XYZ
        avg_rs=Geometry.distance(avg_t_xyz,avg_p_xyz) # average slant range
        for i=1:Np
            p_xyz_i=p_xyz[:,i,:] # p_xyz: 3 x Np x Nst
            p_xyz_i=reshape(p_xyz_i,3,Nst) # p_xyz: 3 x Nst
            p_LLH=Geometry.xyz_to_geo(p_xyz_i)
            platform_heights[i]=mean(p_LLH[3,:]) # average platform heights over slow-time for each platform
        end
        avg_p_h=mean(platform_heights) # average platform height over platforms and slow-time
        if avg_rs<avg_p_h;avg_rs=avg_p_h;end
        avg_rg,look_ang=Scene.slantrange_to_lookangle(earth_radius,avg_rs,avg_p_h,0) # assuming target height is 0 (negligible effect), look_ang: average look angle
    end
    vgrid = Antenna.AntGrid(antennaFile) # read in vpol grid, takes time to load
    ant = sc_ant(vgrid); #create antenna structure, additional arguments are rotation and origin
    sc = spacecraft(avg_p_xyz, Float64.(avg_p_vel), ant = ant, look_angle = look_ang, side = "right"); ##create spacecraft structure; ant, look_angle, side are optional
    co_pol,cross_pol = interpolate_pattern(sc, t_xyz_3xN);#inteprolate pattern (cp:co-pol, xp: cross-pol), outputs are 1xNt complex vectors
 
    targets_ref=targets_ref.*transpose(co_pol).^2/maximum(abs.(co_pol))^2 #TODO separate TX and RX, include range effect?
 
    if target_pos_mode=="grid" # plotting projected pattern only for grid type target distribution
        projected_pattern_3D=Scene.convert_image_1xN_to_3D(abs.(co_pol),length(t_loc_1),length(t_loc_2),length(t_loc_3))#take magnitude and reshape to 3D
        include("modules/plotting.jl");coords_txt=Plotting.coordinates(ts_coord_sys)
        Nt_1=length(t_loc_1)
        Nt_2=length(t_loc_2)
        Nt_3=length(t_loc_3)
        gr()
        if Nt_2>1 && Nt_1>1;display(heatmap(t_loc_2,t_loc_1, 20*log10.(projected_pattern_3D[:,:,1]),ylabel=coords_txt[1],xlabel=coords_txt[2],title = "Antenna Pattern Projected on Targets (V-copol)", fill=true,size=(1600,900)));end #, clim=(-80,40),aspect_ratio=:equal
        if Nt_3>1 && Nt_2>1;display(heatmap(t_loc_3,t_loc_2, 20*log10.(projected_pattern_3D[1,:,:]),ylabel=coords_txt[2],xlabel=coords_txt[3],title = "Antenna Pattern Projected on Targets (V-copol)", fill=false,size=(1600,900)));end #, clim=(-80,40),aspect_ratio=:equal
        if Nt_3>1 && Nt_1>1;display(heatmap(t_loc_3,t_loc_1, 20*log10.(projected_pattern_3D[:,1,:]),ylabel=coords_txt[1],xlabel=coords_txt[3],title = "Antenna Pattern Projected on Targets (V-copol)", fill=false,size=(1600,900)));end #, clim=(-80,40),aspect_ratio=:equal
    end
 end


 function segment_simulation_grid(trg_ref_lat, trg_ref_lon, lat_extent, lon_extent, NB_lat, NB_lon)

    B_lat_extent = lat_extent / NB_lat
    B_lon_extent = lon_extent / NB_lon

    Add_B_lat_extent = 0 #B_lat_extent*0.1 
    Add_B_lon_extent = 0 #B_lon_extent*0.1 #10% additional region added 


    trg_ref_lat_list = zeros(NB_lat * NB_lon)
    trg_ref_lon_list = zeros(NB_lat * NB_lon)

    k = 1
    for i=1:NB_lat
        for j=1:NB_lon
            trg_ref_lat_list[k] = (trg_ref_lat - Add_B_lat_extent) + ((i-1) * (B_lat_extent  + Add_B_lat_extent))
            trg_ref_lon_list[k] = (trg_ref_lon - Add_B_lon_extent) + ((j-1) * (B_lon_extent  + Add_B_lon_extent))
            k=k+1
        end
    end

    return trg_ref_lat_list, trg_ref_lon_list, B_lat_extent+(2*Add_B_lat_extent), B_lon_extent+(2*Add_B_lon_extent)

end


function uniform_random(a, b, n)
    rand_nums = [(b - a) * rand() + a for _ in 1:n]
    return rand_nums
end

function define_target_pixels(t_loc_1_range, t_loc_2_range, t_loc_3_range, DEM_region, num_targ_vol::Int=1, target_mode::Int=1, var_1::Float64=0.0, var_2::Float64=0.0, var_3::Float64=0.0)
    t_loc_1                 =  zeros(length(t_loc_1_range) * length(t_loc_2_range) * length(t_loc_3_range) * num_targ_vol)
    t_loc_2                 =  zeros(length(t_loc_1_range) * length(t_loc_2_range) * length(t_loc_3_range) * num_targ_vol)
    t_loc_3                 =  zeros(length(t_loc_1_range) * length(t_loc_2_range) * length(t_loc_3_range) * num_targ_vol)
    
    m=1
    for i=1:length(t_loc_1_range)
        for j= 1:length(t_loc_2_range)
            for k=1:length(t_loc_3_range)
                for l=1:num_targ_vol
                    if target_mode == 1
                        t_loc_1[m] = t_loc_1_range[i]
                        t_loc_2[m] = t_loc_2_range[j]
                        t_loc_3[m] = t_loc_3_range[k] + DEM_region[i,j]
                        m=m+1
                    elseif target_mode == 2
                        t_loc_1[m] = t_loc_1_range[i] + (rand(Uniform(-1,1)) .* var_1)
                        t_loc_2[m] = t_loc_2_range[j] + (rand(Uniform(-1,1)) .* var_2)
                        t_loc_3[m] = (t_loc_3_range[k] + DEM_region[i,j]) + (rand(Uniform(-1,1)) .* var_3)
                        m=m+1
                    elseif target_mode == 3
                        if l==1
                            t_loc_1[m] = t_loc_1_range[i]
                            t_loc_2[m] = t_loc_2_range[j]
                            t_loc_3[m] = t_loc_3_range[k] + DEM_region[i,j]
                            m=m+1
                        else
                            t_loc_1[m] = t_loc_1_range[i] + (rand(Uniform(-1,1)) .* var_1)
                            t_loc_2[m] = t_loc_2_range[j] + (rand(Uniform(-1,1)) .* var_2)
                            t_loc_3[m] = (t_loc_3_range[k] + DEM_region[i,j])  + (rand(Uniform(-1,1)) .* var_3)
                            m=m+1
                        end
                    end
                end
            end
        end
    end

    t_loc_3xN               = vcat(t_loc_1', t_loc_2', t_loc_3')
    return t_loc_3xN
end


function define_target_pixels(t_loc_1_range, t_loc_2_range, t_loc_3_range, DEM_region,  slope_lat_region, slope_lon_region, num_targ_vol::Int=1, target_mode::Int=1, var_1::Float64=0.0, var_2::Float64=0.0, var_3::Float64=0.0)
    t_loc_1                 =  zeros(length(t_loc_1_range) * length(t_loc_2_range) * length(t_loc_3_range) * num_targ_vol)
    t_loc_2                 =  zeros(length(t_loc_1_range) * length(t_loc_2_range) * length(t_loc_3_range) * num_targ_vol)
    t_loc_3                 =  zeros(length(t_loc_1_range) * length(t_loc_2_range) * length(t_loc_3_range) * num_targ_vol)
    
    m=1
    for i=1:length(t_loc_1_range)
        for j= 1:length(t_loc_2_range)
            for k=1:length(t_loc_3_range)
                for l=1:num_targ_vol
                    if target_mode == 1
                        t_loc_1[m] = t_loc_1_range[i]
                        t_loc_2[m] = t_loc_2_range[j]
                        t_loc_3[m] = t_loc_3_range[k] + DEM_region[i,j]
                        m=m+1
                    elseif target_mode == 2
                        r_var_1 = (rand(Uniform(-1,1)) .* var_1)
                        r_var_2 = (rand(Uniform(-1,1)) .* var_2)
                        t_loc_1[m] = t_loc_1_range[i] + r_var_1
                        t_loc_2[m] = t_loc_2_range[j] + r_var_2
                        t_loc_3[m] = (t_loc_3_range[k] + DEM_region[i,j]) + ( (DEM.deg_to_m_lat(r_var_1) * tand(slope_lat_region[i,j])) + (DEM.deg_to_m_lon(r_var_2, t_loc_1_range[i]) * tand(slope_lon_region[i,j])) ) + (rand(Uniform(-1,1)) .* var_3)
                        m=m+1
                    elseif target_mode == 3
                        if l==1
                            t_loc_1[m] = t_loc_1_range[i]
                            t_loc_2[m] = t_loc_2_range[j]
                            t_loc_3[m] = t_loc_3_range[k] + DEM_region[i,j]
                            m=m+1
                        else
                            t_loc_1[m] = t_loc_1_range[i] + (rand(Uniform(-1,1)) .* var_1)
                            t_loc_2[m] = t_loc_2_range[j] + (rand(Uniform(-1,1)) .* var_2)
                            t_loc_3[m] = (t_loc_3_range[k] + DEM_region[i,j])  + (rand(Uniform(-1,1)) .* var_3)
                            m=m+1
                        end
                    end
                end
            end
        end
    end

    t_loc_3xN               = vcat(t_loc_1', t_loc_2', t_loc_3')
    return t_loc_3xN
end


function define_scene_pixels(s_loc_1_range, s_loc_2_range, s_loc_3_range, DEM_region)
    s_loc_1                 =  zeros(length(s_loc_1_range) * length(s_loc_2_range) * length(s_loc_3_range) )
    s_loc_2                 =  zeros(length(s_loc_1_range) * length(s_loc_2_range) * length(s_loc_3_range) )
    s_loc_3                 =  zeros(length(s_loc_1_range) * length(s_loc_2_range) * length(s_loc_3_range) )
    
    m=1
    for i=1:length(s_loc_1_range)
        for j= 1:length(s_loc_2_range)
            for k=1:length(s_loc_3_range)
                    s_loc_1[m] = s_loc_1_range[i]
                    s_loc_2[m] = s_loc_2_range[j]
                    s_loc_3[m] = s_loc_3_range[k] + DEM_region[i,j]
                    m=m+1               
            end
        end
    end

    s_loc_3xN               = vcat(s_loc_1', s_loc_2', s_loc_3')
    return s_loc_3xN
end

 
function get_platform_ref_from_trg_ref(trg_ref, p_height, p_heading, p_la, p_look_dir, earth_radius=6378.137e3) # LLH coords

    platform_ref_point      = [0.0;0.0;0.0]
    ground_range            = Scene.lookangle_to_range(p_la, p_height, trg_ref[3], earth_radius)[2]  
    displacement_N          = ground_range .* cosd(p_heading + 90) 
    displacement_E          = ground_range .* sind(p_heading + 90)
    if p_look_dir == "right"
        platform_ref_point[1]             = trg_ref[1] - ((displacement_N/earth_radius)*180/pi)
        platform_ref_point[2]             = trg_ref[2] - ((displacement_E/(earth_radius*cosd(platform_ref_point[1])))*180/pi) 
    elseif p_look_dir == "left"
        platform_ref_point[1]             = trg_ref[1] + ((displacement_N/earth_radius)*180/pi)
        platform_ref_point[2]             = trg_ref[2] + ((displacement_E/(earth_radius*cosd(platform_ref_point[1])))*180/pi) 
    end 
    platform_ref_point[3]   = p_height

    return platform_ref_point

end

function get_trg_ref_from_platform_ref(p_ref, trg_height, p_heading, p_la, p_look_dir, earth_radius=6378.137e3) # LLH coords

    trg_ref_point           = [0.0;0.0;0.0]
    ground_range            = Scene.lookangle_to_range(p_la,p_ref[3],trg_height, earth_radius)[2]  
    displacement_N          = ground_range .* cosd(p_heading + 90) 
    displacement_E          = ground_range .* sind(p_heading + 90)
    if p_look_dir == "right"
        trg_ref_point[1]             = p_ref[1] + ((displacement_N/earth_radius)*180/pi)
        trg_ref_point[2]             = p_ref[2] + ((displacement_E/(earth_radius*cosd(p_ref[1])))*180/pi) 
    elseif p_look_dir == "left"
        trg_ref_point[1]             = p_ref[1] - ((displacement_N/earth_radius)*180/pi)
        trg_ref_point[2]             = p_ref[2] - ((displacement_E/(earth_radius*cosd(p_ref[1])))*180/pi) 
    end 
    trg_ref_point[3]        = trg_height

    return trg_ref_point

end


end #end module
