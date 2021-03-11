using Test
include("geometry.jl")
include("../modules/sync.jl")
include("../inputs/input_parameters_sync.jl")
# Testing some of the deterministic functions.
pos = cat([1000 2000 1000; 3000 4000 5000; 0 0 0], [5000 6000 7000; 7000 8000 9000; 1000 2000 3000], [-1000 -2000 -1000; -3000 -4000 -5000; 0 0 0], dims=3)

# define parameters in order to run sync module
parameters.mode=1
parameters.fc=1e9
parameters.fp=1
parameters.sync_pri=1
parameters.sync_processing_time=.005
parameters.sync_signal_len=1024
parameters.sync_fc=1e9
parameters.sync_fs=25e6
parameters.sync_fbw=25e6
parameters.sync_fmin=0.01
parameters.f_osc=10e6
parameters.sync_clk_fs=1e3
parameters.master=1
osc_coeffs = repeat([1 1 1 1 1],3)
parameters.osc_coeffs=osc_coeffs
parameters.sigma_freq_offsets=zeros(3)
parameters.no_sync_flag=false
slow_time=(0:1/parameters.fp:2)

@testset "simtool" begin
    @testset "Geometry Module" begin
        @test [0,0,-1] == Geometry.rotate_frame([1,0,0], Geometry.quat(-90, [0,1,0]))
        @test [6.378137e6,0,0] == Geometry.geo_to_xyz([0,0,0],6.378137e6  , 0.)
        @test [0,0,6.378137e6] == Geometry.geo_to_xyz([90,0,0],6.378137e6  , 0.)
        @test [0,6.378137e6,0] == Geometry.geo_to_xyz([0,90,0],6.378137e6  , 0.)
        @test [0,1,0] == Geometry.geo_to_xyz( Geometry.xyz_to_geo([0,1,0], 1 , 0) , 1 , 0)
        @test [0,cosd(30),sind(30)] == Geometry.geo_to_xyz(Geometry.xyz_to_geo([0,cosd(30),sind(30)], 1 , 0) , 1 , 0)
    end
    @testset "Sync Module" begin
        @test ([0.9084405771506712, 0.8548258968972741, 13.848179529735841, 39.02668776561919, 3.2794229862786333], [-2.5, -1.5, -0.5, 0.5, 1.5]) == Sync.osc_psd_twosided(5,5,[1 1 1 1 1])
        @test (3,3) == size(Sync.get_sync_phase(slow_time,pos, parameters))
        @test (3,3) == size(Sync.getSensorCRLB_network(pos, 1024,1e9,25e6, 25e6)[1])
        @test [4 5 1 2 3]' == Sync.shift([1 2 3 4 5]',2)
        @test [0.0, 0.1, 0.2, 0.3, 0.4] == Sync.osc_timeseries_from_psd_twosided([1 2 3 4 5], 10)[2]
        @test 5 == length(Sync.osc_timeseries_from_psd_twosided([1 2 3 4 5], 10)[1])
        @test 5 == length(Sync.sync_effects_on_PSD([1 2 3 4 5], [-2 -1 0 1 2], 0.005, 1, 1, 1, 1))
    end
end