include("modules/geometry.jl")

using Plots
#gr()
pyplot()
#plotly()
#plotlyjs()
## planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)
# TODO calculate ra
## coordinate transformations
geo=[35.38987,-111.8116,9748.89523]
#geo=[0,0,0]
xyz=Geometry.geo_to_xyz(geo,a,e)
geo2=Geometry.xyz_to_geo(xyz,a,e)
println("XYZ: ", xyz)
println("GEO: ", geo)
println("GEO2: ", geo2)

sch=[-19766.4,23.145535442,9748.895229822]
println(sch)
peg=[35.2117072245,-111.8112805579,179.8535529463]
xyz2=Geometry.sch_to_xyz(sch,peg,a,e)
println("SCH to XYZ: ", xyz2)

Mxyzprime_xyz,O,ra=Geometry.peg_calculations(peg,a,e) # TODO use structure
xyz3=Geometry.sch_to_xyz_2(sch,Mxyzprime_xyz,O,ra)
println("SCH to XYZ: ", xyz3)
## rotations
# rotating vector with quaternion
q = Geometry.quat(45, [0,1,0]) #create a quaternion to rotate a vector by 45 degrees about yaxis [0,1,0]
rotated_vec = Geometry.rotate_vec([1,0,0], q) #rotate a vector aligned with the x-axis, by q
println("\nRotated Vector: ", rotated_vec)

# combine rotations q3*q2*q1 means first q1 then q2 then q3
q1 = Geometry.quat(90, [1,0,0]) #create a quaternion to rotate a vector by 90 degrees about xaxis
q2 = Geometry.quat(90, [0,1,0]) #create a quaternion to rotate a vector by 90 degrees about yaxis
q3 = Geometry.quat(90, [0,0,1]) #create a quaternion to rotate a vector by 90 degrees about zaxis
rotated_vec = Geometry.rotate_vec([1,0,0], q1) #rotate a vector aligned with the x-axis, by q1
println("\nRotated Vector: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q2) #rotate a vector aligned with the x-axis, by q2
println("Rotated Vector: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q3) #rotate a vector aligned with the x-axis, by q3
println("Rotated Vector: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q2*q1) #rotate a vector aligned with the x-axis, by q2*q1
println("Rotated Vector (q2*q1): ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q1*q2) #rotate a vector aligned with the x-axis, by q1*q2
println("Rotated Vector (q1*q2): ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q1*q2*q3) #rotate a vector aligned with the x-axis, by q1*q2*q3
println("Combined Rotated Vector, q1*q2*q3: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q1*q3*q2) #rotate a vector aligned with the x-axis, by q1*q3*q2
println("Combined Rotated Vector, q1*q3*q2: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q2*q1*q3) #rotate a vector aligned with the x-axis, by q2*q1*q3
println("Combined Rotated Vector, q2*q1*q3: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q2*q3*q1) #rotate a vector aligned with the x-axis, by q2*q3*q1
println("Combined Rotated Vector, q2*q3*q1: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q3*q2*q1) #rotate a vector aligned with the x-axis, by q3*q2*q1
println("Combined Rotated Vector, q3*q2*q1: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q3*q1*q2) #rotate a vector aligned with the x-axis, by q3*q1*q2
println("Combined Rotated Vector, q3*q1*q2: ", rotated_vec)

# rotate a frame, then describe vector in the rotated frame
qy = Geometry.quat(90, [0,1,0]) #create a quaternion to rotate a frame by 90 degrees about yaxis [0,1,0]
qz= Geometry.quat(90, [0,0,1]) #create a quaternion to rotate a frame by 90deg about z-axis
projected_vec = Geometry.rotate_frame([1,0,0], qy) #project a vector aligned with the x-axis, to a rotated frame described by q
println("\nProjected Vector: ", projected_vec)
projected_vec = Geometry.rotate_frame([1,0,0], qz) #project a vector aligned with the x-axis, to a rotated frame described by q
println("Projected Vector: ", projected_vec)
# EXtrinsic vs. intrinsic rotations:
# intrinsic rotation: second rotation is about rotated frame
# extrinsic rotation: second rotation is about the original frmae
# q1*q2 can be thought of in two ways:
# (A) rotate q1 and then rotate q2 about rotated frame (intrinsic)
# (B) rotate q2 first and then rotate q1 about original frame (extrinsic)
projected_vec = Geometry.rotate_frame([1,0,0], qy*qz) #project a vector aligned with the x-axis, to a rotated frame described by q
println("\nIntrinsic projected vector: ", projected_vec)
projected_vec = Geometry.rotate_frame([1,0,0], qz*qy) #project a vector aligned with the x-axis, to a rotated frame described by q
println("Extrinsic projected vector: ", projected_vec)