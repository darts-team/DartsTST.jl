module SimSetup

#external packages
using ReferenceFrameRotations
using LinearAlgebra

#local packages
using ..Geometry
using ..Antenna


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


end #end module
