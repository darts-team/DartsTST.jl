#include("modules/geometry.jl")
#include("modules/scene.jl")
include("modules/raw_data.jl")
using Plots
pyplot()
pygui(true)

c=299792458 # speed of light (m/s)
#planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)
# radar parameters
mode=1 #1: SAR, 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)

fc=1e9 # center frequency (Hz)

# platform locations (volumetric grid) defined in geo (θϕh)
p_θ=0 # deg latitude
p_ϕ=-0.1:0.02:0.1 # deg longitude
p_h=500e3 # m  heights

# target locations (volumetric grid) defined in geo (θϕh)
t_θ=-0.005:0.0001:0.005 # deg latitude
t_ϕ=-0.005:0.0001:0.005 # deg longitude
t_h=0:10:1000 # m  heights

# generate raw data
rawdata=Raw_Data.main(t_θ,t_ϕ,t_h,p_θ,p_ϕ,p_h,mode,tx_el,fc,a,e)

# plot raw data
plot(abs.(rawdata),leg=false,title="raw data magnitude")
plot(angle.(rawdata),leg=false,title="raw data phase")
# save raw data
