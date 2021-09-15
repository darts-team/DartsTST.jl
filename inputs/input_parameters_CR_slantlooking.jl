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
#orbit_filename="orbitOutput_082020.nc" # position in km, time in sec
orbit_filename="orbit_output_062021.nc" # position in km, time in sec
SAR_duration=0 # synthetic aperture duration (s)
SAR_start_time=0 # SAR imaging start time (s)
# target locations and reflectvities
target_pos_mode="CR" #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
ts_coord_sys="SCH" # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene
display_geometry_coord="SCH" # platform/target/scene geometry (scatter plot) coordinate system: "LLH", "SCH", "XYZ"
look_angle=30 # in cross-track direction, required only if SCH coordinates, using same look angle for targets and scene (deg)
# length(t_loc_1)==length(t_loc_2)==length(t_loc_3) should hold
t_loc_1=[0] # deg latitude if LLH, along-track if SCH, X if XYZ
t_loc_2=[0] # deg longitude if LLH, cross-track if SCH, Y if XYZ
t_loc_3=[100] # m  heights if LLH or SCH, Z if XYZ
t_ref=  [1] # reflectivities
# image/scene pixel coordinates
s_loc_1=0 # deg latitude if LLH, along-track if SCH, X if XYZ
s_loc_2=-900:2:900 # deg longitude if LLH, cross-track if SCH, Y if XYZ
s_loc_3=-500:2:700 # m  heights if LLH or SCH, Z if XYZ
# range spread function (RSF) parameters
pulse_length=10e-6 # s pulse length
Δt=1e-8 # s fast-time resolution (ADC sampling rate effect is excluded for now)
bandwidth=40e6 # bandwidth (Hz)
# performance metrics
res_dB=3 # dB two-sided resolution relative power level (set to 0 for peak-to-null Rayleigh resolution), positive value needed
PSF_image_point=1 # 1: peak location, 2: target location, 3: center of 3D scene
# simulation options
enable_thermal_noise=false # whether to enable or disable random additive noise (e.g. thermal noise)
enable_fast_time=true # whether to enable or disable fast-time axis, 0:disable, 1: enable
display_geometry=false # whether to display geometry plots
display_RSF_rawdata=false # whether to display RSF and rawdata plots
display_tomograms=1 # how to display tomograms, 0: do not display, 1: display only 3 slices at the reference point, 2: display all slices in each dimension, 3: display as 3D scatter plot
include_antenna=false # whether to include projected antenna pattern
display_input_scene=false # display input scene (targets) and delta between input/output scenes (3 slices at the center of scene) with same scene size as output tomogram scene
