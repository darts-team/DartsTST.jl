include("modules/geometry.jl")
a=6378.137e3
e=0.00669437999015^0.5
geo=[35.38987,-111.8116,9748.89523]
xyz=Geometry.geo_to_xyz(geo,a,e)
geo2=Geometry.xyz_to_geo(xyz,a,e)
println(xyz)
println(geo)
println(geo2)

sch=[-19766.4,23.145535442,9748.895229822]
println(sch)
peg=[35.2117072245,-111.8112805579,179.8535529463]
xyz2=Geometry.sch_to_xyz(sch,peg,a,e)
println(xyz2)
