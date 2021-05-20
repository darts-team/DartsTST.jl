c=299792458 # speed of light (m/s)
# planetary shape constants
earth_radius=6378.137e3 # semi-major axis at equator
earth_eccentricity=sqrt(0.00669437999015)
# MIMO parameters
mode=1 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)
# radar parameters
fc=1.25e9 # center frequency (Hz)
fp=3 # pulse repetition frequency (Hz)
SNR=50 # SNR for single platform and single pulse before fast-time processing dB (for additive random noise only) TODO calculate based on sigma-zero (which depends on target type, wavelength, look angle, polarization) and NESZ (which depends on radar specs and processing)
# platform locations in xyz taken from orbits (including slow-time)
orbit_filename="orbitOutput_082020.nc" # position in km, time in sec
SAR_duration=2 # synthetic aperture duration (s)
SAR_start_time=0 # SAR imaging start time (s)
# target locations (volumetric grid) defined in geo (θϕh)
target_pos_mode="CR" #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
ts_coord_sys="LLH" # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene
display_geometry_coord="XYZ" # platform/target/scene geometry (scatter plot) coordinate system: "LLH", "SCH", "XYZ"
t_loc_1=7 # deg latitude
t_loc_2=0 # deg longitude
t_loc_3=0 # m  heights
t_ref=1 # reflectivities
# image/scene pixel coordinates
s_loc_1=7-0.0003:0.000025:7+0.0003 # deg latitude
s_loc_2=-0.0012:0.0001:0.0012 # deg longitude
s_loc_3=-45:5:45 # m  heights
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
display_tomograms=1 # how to display tomograms, 0: do not display, 1: display only 3 slices at the scene center, 2: display all slices in each dimension, 3: display as 3D scatter plot
