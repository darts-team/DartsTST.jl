module UserParameters

using Parameters
using StaticArrays

export earth_radius
export earth_eccentricity
export c

@consts begin
    c = 299792458 # speed of light (m/s)
    earth_radius = 6378.137e3 # Earth semi-major axis at equator
    earth_eccentricity = sqrt(0.00669437999015) # Earth eccentricity
end

# Input parameters structure with default values
# To run simulations with different values of input parameters,
# create a predefined struct in predefined-input-parameters.jl or
# or add parameter=value pairs directly in the workflow when creating the inputParameter object.
@with_kw struct inputParameters

    mode::Int  = 2 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
    tx_el::Int = 1 # which element transmits for SIMO (max value N)
    processing_steps = :bp3d  # :bp3d --> 1-step, :bp2d3d --> 2-step for SAR and tomographic processing

    # radar parameters
    fc::Float64  = 1.25e9 # center frequency (Hz) L-band; fc=3.2e9 # center frequency (Hz) S-band; fc=6e9 # center frequency (Hz) C-band
    fp::Float64  = 10 # pulse repetition frequency (Hz)
    SNR::Float64 = 50 # SNR for single platform and single pulse before fast-time processing dB (for additive random noise only) TODO calculate based on sigma-zero (which depends on target type, wavelength, look angle, polarization) and NESZ (which depends on radar specs and processing)
    SAR_duration::Float64   = 3 # synthetic aperture duration (s)
    SAR_start_time::Float64 = 0 # SAR imaging start time (s)

    # platform locations in xyz (including slow-time locations)
    user_defined_orbit::Int = 2 # 0: use orbits file; 1: user defined orbits in SCH; 2: user defined orbits in TCN
    orbit_filename::String = "orbit_output_062021.nc" # position in km, time in sec; "orbitOutput_082020.nc" --> TODO: convert to :file, :sch, :tcn

    # User defined orbits, set either SCH or TCN, see user_defined_orbit
    p_t0_LLH::Array{Float64,1} = [0;0;750e3] # initial lat/lon (deg) and altitude (m) of reference platform (altitude is assumed constant over slow-time if SCH option)
    Torbit::Float64    = 10*60 # orbital duration (s) (should be larger than 2 x (SAR_start_time+SAR_duration) )
    dt_orbits::Float64 = 0.5 # orbit time resolution (s)
    p_heading::Float64 = 0 # heading (deg), all platforms assumed to have the same heading, 0 deg is north
    pos_n   = [-7.5 -5 -2 0 3.7 5.5 6.5]*1e3 # SCH option: relative position of each platform along n (m), 0 is the reference location
    pos_TCN = [0 -6 0; 0 -5 0; 0 -2 0; 0 0 0; 0 3.5 0; 0 5 0]*1e3   # TCN option: Np x 3 matrix; each row is the TCN coordinate of each platform relative to reference

    # target locations and reflectvities
    target_pos_mode::String="CR" #  targets are defined as three 1D arrays forming either a volumetric grid ("layered-grid" or "shaped-grid") or a 3xN array ("CR" for corner reflectors)
    ts_coord_sys::String="SCH" # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene
    display_geometry_coord::String="SCH" # platform/target/scene geometry (scatter plot) coordinate system: "LLH", "SCH", "XYZ"
    look_angle=30 # in cross-track direction, required only if SCH coordinates, using same look angle for targets and scene (deg)
    t_loc_1 = [0.] # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2 = [0.] # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3 = [0.] # m heights if LLH or SCH, Z if XYZ
    t_ref   = [1.] # reflectivities: a list of CRs in CR mode; an arbitrary vertical profile that will be interpolated on t_loc_3 axis in *grid modes


    # image/scene pixel coordinates
    s_loc_1 = 0 # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2 = -40:1:40 # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3 = -40:1:40 # m  heights if LLH or SCH, Z if XYZ

    # range spread function (RSF) parameters
    pulse_length::Float64 = 10e-6 # s pulse length
    Δt::Float64 = 1e-9 # s fast-time resolution (ADC sampling rate effect is excluded for now)
    bandwidth::Float64 = 40e6 # bandwidth (Hz)

    # derived parameters (some are needed further below)
    λ = c/fc # wavelength (m)
    h = p_t0_LLH[3] # default altitude
    inc_angle = asind((earth_radius+h)./earth_radius.*sind(look_angle))
    Ns_1 = length(s_loc_1)
    Ns_2 = length(s_loc_2)
    Ns_3 = length(s_loc_3)

    # performance metrics
    res_dB::Float64 = 4 # dB two-sided resolution relative power level (set to 0 for peak-to-null Rayleigh resolution), positive value needed
    PSF_image_point::Int = 3 # 1: peak location, 2: target location, 3: center of 3D scene
    PSF_cuts::Int = 2 # 1: principal axes (SCH, LLH, XYZ based on ts_coord_sys), 2: a single cut along PSF_direction_xyz in scene coordinates relative to center of scene
    PSF_direction = [0 1 tand(inc_angle)] # # direction (in ts_coord_sys) relative to scene center to take 1D PSF cut along a line which goes through center of scene (used only if PSF_cuts=2), direction along non-existing scene dimension is ignored; default cut is along n. For cut along r, use [0 1 -1/tand(inc_angle)]

    # antenna settings
    antennaFile::String = "inputs/darts_ant_03192021.nc"

    # synchronization parameters
    sync_pri::Float64 = 1.0 # (s) repetition interval of sync
    sync_processing_time::Float64 = 0.001 # processing time between stage 1 and stage 2 sync
    sync_signal_len::Int64 = 1024 # waveform length
    sync_fc = fc # waveform center frequency - set to same as radar freq
    sync_fs::Float64 = 25e6 # sync receiver sampling rate
    sync_fbw = sync_fs # LFM bandwidth
    sync_fmin::Float64 = 1.0 # minimum frequency > 0 in Hz to window PSD
    sync_clk_fs::Float64 = 1e3; # sample rate of clock error process
    sync_master::Int = 1; # selection of master transmitter for sync (assumes a simplified communication achitecture- all talking with one master platform)

    sync_osc_type::String = "USO"
    sync_a_coeff_dB = [-95 -90 -200 -130 -155] # [USO: Krieger]
    sync_f_osc::Float64 = 10e6 # local oscillator frequency --> depends on oscillator type

    # positioning parameters

    # simulation options
    enable_thermal_noise::Bool    = false # whether to enable or disable random additive noise (e.g. thermal noise)
    display_geometry::Bool        = false # whether to display geometry plots
    display_1D_cuts::Bool         = false # whether to 1D cuts from Scene module
    display_custom_orbit::Bool    = false # whether to show custom orbit on Earth sphere (for a duration of Torbit)
    display_RSF_rawdata::Bool     = false # whether to display RSF and rawdata plots
    display_tomograms::Int        = 1 # how to display tomograms, 0: do not display, 1: display only 3 slices at the reference point, 2: display all slices in each dimension, 3: display as 3D scatter plot
    include_antenna::Bool         = false # whether to include projected antenna pattern
    display_input_scene::Bool     = false # display input scene (targets) and delta between input/output scenes (3 slices at the center of scene) with same scene size as output tomogram scene
    no_sync_flag::Bool            = false # if flag == true, no sync is used. flag == false results in normal sync process estimation
    enable_sync_phase_error::Bool = false # if flag == true, oscillator phase errors considered. If false, ideal oscillators used

    # logging level


end


"""
    validateInputParams(params)

Check consistency of user parameters in inputParameters object `params`.
Return `true` if `params` contains valid parameters.

    # Examples
```julia-repl
julia> paramsIsValid = validateInputParams(params)
true
```
"""
function validateInputParams(params)

    if params.target_pos_mode == "CR"
        @assert length(params.t_loc_1) == length(params.t_loc_2) == length(params.t_loc_3) "Size of target location arrays must be equal for target_pos_mode=CR"
    end

    @assert params.target_pos_mode in ["CR", "layered-grid", "shaped-grid", "grid"]  "Target position mode is not valid"

    # Add more @assert's here

    return true
end

end