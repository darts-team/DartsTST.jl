c=299792458 # speed of light (m/s)
# planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)
# radar parameters
mode=2 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)
fc=1e9 # center frequency (Hz)
# platform locations (volumetric grid) defined in geo (θϕh)
p_θ=-0.04:0.01:0.04 # deg latitude
p_ϕ=-0.08:0.01:0.08 # deg longitude
p_h=500e3#495e3:1e2:505e3 # m  heights
# target locations (volumetric grid) defined in geo (θϕh)
t_θ=0 # deg latitude
t_ϕ=0#-0.001:0.0005:0.001 # deg longitude
t_h=0 # m  heights
# image/scene pixel coordinates
s_θ=0-0.0004:0.00001:0+0.0004 # deg latitude
s_ϕ=-0.0004:0.00001:0.0004 # deg longitude
s_h=0 # m  heights
