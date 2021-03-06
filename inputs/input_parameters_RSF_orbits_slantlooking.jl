c=299792458 # speed of light (m/s)
# planetary shape constants
earth_radius=6378.137e3 # semi-major axis at equator
earth_eccentricity=sqrt(0.00669437999015)
# MIMO parameters
mode=1 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)
# radar parameters
fc=1e9 # center frequency (Hz)
fp=3 # pulse repetition frequency (Hz)
SNR=50 # SNR for single platform and single pulse before fast-time processing dB (for additive random noise only) TODO calculate based on sigma-zero (which depends on target type, wavelength, look angle, polarization) and NESZ (which depends on radar specs and processing)
# platform locations in xyz taken from orbits (including slow-time)
orbit_filename="orbitOutput_082020.nc" # position in km, time in sec
SAR_duration=3 # synthetic aperture duration (s)
SAR_start_time=0 # SAR imaging start time (s)
# target locations (volumetric grid) defined in geo (θϕh)
t_θ=7 # deg latitude
t_ϕ=0 # deg longitude
t_h=0 # m  heights
# image/scene pixel coordinates
s_θ=7-0.0002:0.000004:7+0.0002 # deg latitude
s_ϕ=-0.0008:0.00004:0.0008 # deg longitude
s_h=-30:1:30 # m  heights
# range spread function (RSF) parameters
enable_fast_time=true # whether to enable or disable fast-time axis, 0:disable, 1: enable
enable_thermal_noise=false # whether to enable or disable random additive noise (e.g. thermal noise)
Trx=300e-6 # s duration of RX window (may need to be increased if aperture or scene is large) TODO (adjust based on max/min range)
pulse_length=10e-6 # s pulse length
Δt=1e-8 # s fast-time resolution (ADC sampling rate effect is excluded for now)
bandwidth=10e6 # bandwidth (Hz)
# performance metrics
res_dB=3 # dB two-sided resolution relative power level (set to 0 for peak-to-null Rayleigh resolution), positive value needed
