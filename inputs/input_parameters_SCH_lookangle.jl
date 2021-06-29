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
# target locations and reflectvities
target_pos_mode="CR" #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
ts_coord_sys="SCH" # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene
display_geometry_coord="SCH" # platform/target/scene geometry (scatter plot) coordinate system: "LLH", "SCH", "XYZ"
if ts_coord_sys=="SCH" # if SCH, target and scene locations are defined relative to the point where look angle vector intersects the surface
    look_angle=30 # in cross-track direction, required only if SCH coordinates, using same look angle for targets and scene (deg)
    #p_avg_heading=0.1 # average heading of platforms, due North is 0, due East is 90 (deg), required only if SCH coordinates TODO we should get this from orbits!
end
if target_pos_mode=="grid" # target positions are defined as a volumetric grid (useful for distributed target)
    t_loc_1=-10:10:10 # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=-20:20:20 # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=0:30:30 # m  heights if LLH or SCH, Z if XYZ
    t_ref=rand(Float64,length(t_loc_1),length(t_loc_2),length(t_loc_3)) # uniform random reflectivities between 0 and 1, a 3D input array (e.g. 3D image) can be used instead
elseif target_pos_mode=="CR" # ("CR" for corner reflector) target positions are defined as 3xN array (useful for a few discrete targets)
    # length(t_loc_1)==length(t_loc_2)==length(t_loc_3) should hold
    t_loc_1=[0] # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=[0] # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=[40] # m  heights if LLH or SCH, Z if XYZ
    t_ref=  [1] # reflectivities
end
# image/scene pixel coordinates
s_loc_1=-40:2:40 # deg latitude if LLH, along-track if SCH, X if XYZ
s_loc_2=-60:2:60 # deg longitude if LLH, cross-track if SCH, Y if XYZ
s_loc_3=  0:2:80 # m  heights if LLH or SCH, Z if XYZ
# s_loc_1=-10:.25:10 # deg latitude if LLH, along-track if SCH, X if XYZ
# s_loc_2=-40:.5:40 # deg longitude if LLH, cross-track if SCH, Y if XYZ
# s_loc_3=  (-15:.25:15) .+ 40 # m  heights if LLH or SCH, Z if XYZ

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
display_geometry=true # whether to display geometry plots
display_RSF_rawdata=false # whether to display RSF and rawdata plots
display_tomograms=1 # how to display tomograms, 0: do not display, 1: display only 3 slices at the reference point, 2: display all slices in each dimension, 3: display as 3D scatter plot
