c=299792458 # speed of light (m/s)
# planetary shape constants
earth_radius=6378.137e3 # semi-major axis at equator
earth_eccentricity=sqrt(0.00669437999015)
# MIMO parameters
mode=2 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)
# radar parameters
fc=1.25e9 # center frequency (Hz) L-band
# fc=3e9 # center frequency (Hz) S-band
# fc=6e9 # center frequency (Hz) C-band
fp=10 # pulse repetition frequency (Hz)
SNR=50 # SNR for single platform and single pulse before fast-time processing dB (for additive random noise only) TODO calculate based on sigma-zero (which depends on target type, wavelength, look angle, polarization) and NESZ (which depends on radar specs and processing)
# platform locations in xyz taken from orbits (including slow-time)
orbit_filename="orbitOutput_082020.nc" # position in km, time in sec
SAR_duration=2 # synthetic aperture duration (s)
SAR_start_time=0 # SAR imaging start time (s)
# target/scene locations and target reflectvities
target_pos_mode="CR" #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
ts_coord_sys="XYZ" # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene
display_geometry_coord="LLH" # platform/target/scene geometry (scatter plot) coordinate system: "LLH", "SCH", "XYZ"
# length(t_loc_1)==length(t_loc_2)==length(t_loc_3) should hold
if ts_coord_sys=="XYZ"
    t_loc_1=[earth_radius] # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=[0] # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=[0] # m  heights if LLH or SCH, Z if XYZ
    t_ref=  [1] # reflectivities
    s_loc_1=earth_radius # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2=-350:2:350 # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3=-54:0.2:54 # m  heights if LLH or SCH, Z if XYZ
    PSF_direction=[0 13 2] # direction (in ts_coord_sys) relative to scene center to take 1D PSF cut along a line which goes through center of scene (used only if PSF_cuts=2), direction along non-existing scene dimension is ignored
elseif ts_coord_sys=="SCH"
    t_loc_1=[0] # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=[0] # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=[0] # m  heights if LLH or SCH, Z if XYZ
    t_ref=  [1] # reflectivities
    s_loc_1=-20:0.5:20 # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2=-610:5:610 # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3=0 # m  heights if LLH or SCH, Z if XYZ
    look_angle=0 # in cross-track direction, required only if SCH coordinates, using same look angle for targets and scene (deg)
    PSF_direction=[20 610 0] # direction (in ts_coord_sys) relative to scene center to take 1D PSF cut along a line which goes through center of scene (used only if PSF_cuts=2), direction along non-existing scene dimension is ignored
elseif ts_coord_sys=="LLH"
    t_loc_1=[0] # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=[0] # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=[0] # m  heights if LLH or SCH, Z if XYZ
    t_ref=  [1] # reflectivities
    s_loc_1=-0.0004:0.000005:0.0004 # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2=-0.0029:0.000025:0.0029 # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3=0 # m  heights if LLH or SCH, Z if XYZ
    PSF_direction=[1 7 0] # direction (in ts_coord_sys) relative to scene center to take 1D PSF cut along a line which goes through center of scene (used only if PSF_cuts=2), direction along non-existing scene dimension is ignored
end
# range spread function (RSF) parameters
pulse_length=10e-6 # s pulse length
Δt=1e-8 # s fast-time resolution (ADC sampling rate effect is excluded for now)
bandwidth=40e6 # bandwidth (Hz)
# performance metrics
res_dB=3 # dB two-sided resolution relative power level (set to 0 for peak-to-null Rayleigh resolution), positive value needed
PSF_image_point=1 # 1: peak location, 2: target location, 3: center of 3D scene
PSF_cuts=2 # 1: principal axes (SCH, LLH, XYZ based on ts_coord_sys), 2: a single cut along PSF_direction_xyz shown in xyz TODO (works only if scene is defined in XYZ, SCH/LLH to be added)
# simulation options
enable_thermal_noise=false # whether to enable or disable random additive noise (e.g. thermal noise)
enable_fast_time=true # whether to enable or disable fast-time axis, 0:disable, 1: enable
display_geometry=false # whether to display geometry plots
display_RSF_rawdata=false # whether to display RSF and rawdata plots
display_tomograms=1 # how to display tomograms, 0: do not display, 1: display only 3 slices at the reference point, 2: display all slices in each dimension, 3: display as 3D scatter plot
include_antenna=false # whether to include projected antenna pattern
display_input_scene=false # display input scene (targets) and delta between input/output scenes (3 slices at the center of scene) with same scene size as output tomogram scene