module TomoWorkflow

include("../modules/generate_raw_data.jl")
include("../modules/process_raw_data.jl")
include("../modules/geometry.jl")
include("../modules/scene.jl")
include("../modules/range_spread_function.jl") # as RSF
include("../modules/orbits.jl")
include("../modules/sync.jl")
include("../modules/error_sources.jl")
include("../modules/performance_metrics.jl")
include("../modules/antenna.jl")
include("../modules/simsetup.jl")
include("../modules/user_parameters.jl")
include("../modules/plotting.jl")
using NCDatasets
using Statistics
using Parameters
using Dates
using StaticArrays
using .UserParameters

"""
Uses the input parameters to generate complex tomogram data

# Arguments
- `params::user_parameters struct`: structure of input parameters

# Output
- `image_3D::Array with dimensions of given scene`: Complex 3D tomographic image
"""
function generate_tomo(params)

    # Compute orbits time, position, and velocity
    orbit_time, orbit_pos, orbit_vel = Orbits.computeTimePosVel(params)
    # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
    p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)

    # Create target/scene location 
    targets_loc, targets_ref, Nt = Scene.construct_targets_str(params) # Nt: number of targets, targets: structure array containing target locations and reflectivities
    s_loc_3xN  = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops
    t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions

    # Read number of platforms (todo: move into a struct)
    Np  = size(orbit_pos)[2] # number of platforms

    # Apply antenna pattern
    if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
        Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
    end

    # Generate range spread function (matched filter output)
    min_range, max_range = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
    Trx = 2*(max_range-min_range)/c + 5*params.pulse_length # s duration of RX window
    Srx, MF, ft, t_rx = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
    # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

    # Generate TomoSAR raw data
    ref_range = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
    rawdata = Generate_Raw_Data.main_RSF_slowtime(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
    if params.enable_thermal_noise # adding random noise based on SNR after range (fast-time) processing
        rawdata = Error_Sources.random_noise(rawdata, params)
    end 

    # Add phase error
    sync_osc_coeffs = repeat(params.sync_a_coeff_dB, Np)
    if params.enable_sync_phase_error
        rawdata = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, s_xyz_3xN, sync_osc_coeffs, params)
    end

    # Process raw data to generate image
    if params.processing_steps === :bp3d # 1-step processing
        image_3D = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
    elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
        SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
        image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
    end

return image_3D

end#function

end#module