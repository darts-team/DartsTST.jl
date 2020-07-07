include("modules/raw_data.jl")
using Plots
#pyplot()

c=299792458 # speed of light (m/s)
#planetary shape constants
a=6378.137e3
e=sqrt(0.00669437999015)
# radar parameters
mode=3 #1: SAR (ping-pong), 2:SIMO, 3:MIMO
tx_el=1 # which element transmits for SIMO (max value N)

fc=1e9 # center frequency (Hz)

# platform locations (volumetric grid) defined in geo (θϕh)
p_θ=0 # deg latitude
p_ϕ=-0.1:0.01:0.1 # deg longitude
p_h=500e3 # m  heights

# target locations (volumetric grid) defined in geo (θϕh)
t_θ=-0.005:0.0005:0.005 # deg latitude
t_ϕ=-0.005:0.0005:0.005 # deg longitude
t_h=0:50:1000 # m  heights

# generate raw data
rawdata=Raw_Data.main(t_θ,t_ϕ,t_h,p_θ,p_ϕ,p_h,mode,tx_el,fc,a,e)

# plot raw data
if mode==3
    display(heatmap(abs.(rawdata),c=cgrad([:black,:white]),xlabel="RX platform",ylabel="TX platform",title="raw data magnitude"))
    heatmap(angle.(rawdata)*180/pi,c=cgrad([:black,:white]),xlabel="RX platform",ylabel="TX platform",title="raw data phase (deg)")
else
    display(plot(abs.(rawdata),leg=false,title="raw data magnitude"))
    plot(angle.(rawdata)*180/pi,leg=false,title="raw data phase (deg)")
end
# save raw data
#savefig("yyyyyy.xxx")
