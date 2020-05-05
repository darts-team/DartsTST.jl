include("modules/geometry.jl")


a=6378.137e3
e=sqrt(0.00669437999015)
geo=[35.38987,-111.8116,9748.89523]
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



# rotating vector with quaternion
q = Geometry.quat(45, [0,1,0]) #create a quaternion to rotate a vector by 45 degrees about yaxis [0,1,0]
rotated_vec = Geometry.rotate_vec([1,0,0], q) #rotate a vector aligned with the x-axis, by q
println("\nRotated Vector: ", rotated_vec)

# rotate a frame, then describe vector in the rotated frame
q = Geometry.quat(45, [0,1,0]) #create a quaternion to rotate a frame by 45 degrees about yaxis [0,1,0]
projected_vec = Geometry.rotate_frame([1,0,0], q) #project a vector aligned with the x-axis, to a rotated frame described by q
println("\nProjected Vector: ", projected_vec)


# combine rotations
q1 = Geometry.quat(35, [0,1,0]) #create a quaternion to rotate a vector by 10 degrees about yaxis [0,1,0]
q2 = Geometry.quat(8, [0,0,1]) #create a quaternion to rotate a vector by 35 degrees about yaxis [0,1,0]
q3 = Geometry.quat(2, [0,1,0]) #create a quaternion to rotate a vector by 35 degrees about yaxis [0,1,0]
rotated_vec = Geometry.rotate_vec([1,0,0], q1) #rotate a vector aligned with the x-axis, by q1
println("\nRotated Vector: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q2*q1) #rotate a vector aligned with the x-axis, by q2
println("\nRotated Vector (q2*q1): ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q1*q2) #rotate a vector aligned with the x-axis, by q2
println("\nRotated Vector (q1*q2): ", rotated_vec)

rotated_vec = Geometry.rotate_vec(rotated_vec, q3) #rotate a vector aligned with the x-axis, by q2
println("\nCombined Rotated Vector: ", rotated_vec)
rotated_vec = Geometry.rotate_vec([1,0,0], q1*q2*q3) #rotate a vector aligned with the x-axis, by q1 and q2
println("Combined Rotated Vector, qmult: ", rotated_vec)



# EXtrinsic vs. intrinsic rotations
qy = Geometry.quat(90, [0,1,0]) #create a quaternion to rotate a frame by 45 degrees about yaxis [0,1,0]
qz= Geometry.quat(90, [0,0,1]) #create a quaternion to rotate a frame by 90deg about z-axis

#intrinsic rotation: second rotation is about rotated frame
projected_vec = Geometry.rotate_frame([1,0,0], qy*qz) #project a vector aligned with the x-axis, to a rotated frame described by q
println("\nIntrinsic projected vector: ", projected_vec)

#extrinsic rotation: second rotation is about the original frmae
projected_vec = Geometry.rotate_frame([1,0,0], qz*qy) #project a vector aligned with the x-axis, to a rotated frame described by q
println("\nExtrinsic projected vector: ", projected_vec)
