include("../modules/tomographic_ISLR.jl")
using Plots

mode=1 # tomographic mode, 1: SAR, 2: SIMO, 3: MIMO
#plat_distr=-10000:2500:10000 # locations of platforms (m) (uniform spacing)
plat_distr=[-10000;-5000;-2000;-1000;3000;5500;7000;7700;10000] # locations of platforms (m) (non-uniform distribution)
# plat_distr is relative platform locations projected perpendicular to look angle direction relative to the reference platform location
look_angle=30 # perpendicular to along-track and with respect to nadir (deg)
local_slope=5 # local slope of terrain (deg)
altitude=700e3 # of the reference platform (m)
max_veg_H=20 # maximum vegetation height (m)
freq=1.2e9 # radar center frequency (Hz)
scene_res=0.2 # pixel resolution of the scene (m)

PSF,ISLR,xsc=tomographic_ISLR.main(mode,plat_distr,look_angle,local_slope,altitude,max_veg_H,freq,scene_res)

if mode==1;mode_text="SAR";elseif mode==2;mode_text="SIMO";elseif mode==3;mode_text="MIMO";end
display(plot(xsc,PSF,leg=false,title="Tomographic PSF ("*mode_text*") (ISLR: "*string(round(ISLR,digits=2))*" dB)",xaxis=("scene (m)"),ylabel=("amplitude (linear)"),size=(1600,900)))
display(plot(xsc,20*log10.(PSF),leg=false,title="Tomographic PSF ("*mode_text*") (ISLR: "*string(round(ISLR,digits=2))*" dB)",xaxis=("scene (m)"),ylabel=("amplitude (dB)"),size=(1600,900)))
display("ISLR: "*string(round(ISLR,digits=2))*" dB")
