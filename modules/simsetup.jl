module SimSetup

#external packages
using ReferenceFrameRotations
using LinearAlgebra

#local packages
include("geometry.jl")
using .Geometry
include("antenna.jl")
using .Antenna

#abstract type AbstractFrame end
#abstract type AbstractPlanetFrame <: AbstractFrame end
#abstract type AbstractSpaceCraftFrame <: AbstractPlanetFrame end
#abstract type AbstractAntennaFrame <: AbstractSpaceCraftFrame end

#mutable struct Planet_frame <:AbstractFrame
#    org
#    qtrn::Quaternion
#    a
#    e_sqr
#end
mutable struct sc_ant
    pattern
    orientation # with respect to the SpaceCraft
    org #origin relative to spacecraft frame
    function sc_ant(pattern, orientation=Geometry.quat(-90, [0,0,1]), org=[0.,0.,0.])
        new(pattern, orientation, org)
    end
end

mutable struct spacecraft
    org #origin in ECEF (m,m,m)
    vel #velocity relative in ECEF (m/s,m/s, m/s)
    orientation # with respect to the planet
    antenna::sc_ant # antenna frame defined relative to sc frame
    #function SpaceCraft(org::Float)
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


function interpolate_antenna(sc::spacecraft, targ_xyz::Array{Float64,}, freq::Float64=1.25)
    println("SpaceCraft ")


end


end #end module
