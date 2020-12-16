c=299792458 # speed of light (m/s)
# planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)
# radar parameters
mode=1 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)
fc=1e9 # center frequency (Hz)
fp=1 # pulse repetition frequency (Hz)
# platform locations in xyz taken from orbits (including slow-time)
orbit_filename="orbitOutput_082020.nc" # position in km, time in sec
SAR_duration=10 # synthetic aperture duration (s)
SAR_start_time=0 # SAR imaging start time (s)
# target locations (volumetric grid) defined in geo (θϕh)
t_θ=0 # deg latitude
t_ϕ=0 # deg longitude
t_h=0 # m  heights
# image/scene pixel coordinates
s_θ=-0.00001:0.0000002:0.00001 # deg latitude
s_ϕ=-0.00001:0.0000002:0.00001 # deg longitude
s_h=0 # m  heights
# range spread function (RSF) parameters
Trx=300e-6 # s duration of RX window (may need to be increased if aperture or scene is large) TODO (adjust based on max/min range)
τ=10e-6 # s pulse length
Δt=1e-8 # s fast-time resolution (ADC sampling rate effect is excluded for now)
B=10e6 # bandwidth (Hz)
