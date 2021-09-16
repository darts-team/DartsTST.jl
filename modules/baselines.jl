using NCDatasets
using Plots
include("orbits.jl")

#orbit_filename="orbitOutput_082020.nc" # position in km, time in sec
orbit_filename="orbit_output_062021.nc" # position in km, time in sec
look_angle=50 # deg
max_time=20000 # seconds

time12=orbit_dataset["time"][1:2]
dt=time12[2]-time12[1]
max_time_index=Int(round(max_time/dt))
time_index=1:1:max_time_index

orbit_dataset=Dataset("inputs/"*orbit_filename) # Read orbits data in NetCDF format
time=orbit_dataset["time"][time_index]
p=orbit_dataset["position"][:,:,time_index]
v=orbit_dataset["velocity"][:,:,time_index]
baselines=Orbits.get_perp_baselines(p,v,look_angle)
maxbaselines=reshape(maximum(maximum(baselines,dims=1),dims=2),max_time_index)
maxbaseline,maxb_time_ind=findmax(maxbaselines)
println("Maximum Baseline Along Orbit: ",round.(maxbaseline,digits=1)," km")
println("Maximum Baseline Time Index: ",maxb_time_ind)
println("Maximum Baseline Time: ",Int((maxb_time_ind-1)*dt)," s")
plotly();plot(time/60,maxbaselines,leg=false,size=(1000,600),xlabel="Orbit Time (min)",ylabel="Perpendicular Baseline (km)",title="Orbit: "*orbit_filename*", look angle: "*string(look_angle)*" deg")
