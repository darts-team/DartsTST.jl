using Test
include("geometry.jl")

@testset "Geometry Module" begin
    @test [0,0,-1] == Geometry.rotate_frame([1,0,0], Geometry.quat(-90, [0,1,0]))
    @test [6.378137e6,0,0] == Geometry.geo_to_xyz([0.,0.,0.],6.378137e6  , 0.)
    @test [0,0,6.378137e6] == Geometry.geo_to_xyz([90.,0.,0.],6.378137e6  , 0.)
    @test [0,6.378137e6,0] == Geometry.geo_to_xyz([0.,90.,0.],6.378137e6  , 0.)
    @test [0.,1.,0.] == Geometry.geo_to_xyz( Geometry.xyz_to_geo([0.,1.,0.]))
    @test [15.,13.,12.] ≈ Geometry.xyz_to_geo(Geometry.geo_to_xyz([15.,13.,12.]))

end
