using NCDatasets
using Plots
using Dates
include("orbits.jl")

#orbit_filename="orbitOutput_082020.nc" # position in km, time in sec
orbit_filename="orbit_output_062021.nc" # position in km, time in sec
orbit_dataset=Dataset("inputs/"*orbit_filename) # Read orbits data in NetCDF format

look_angle=30 # deg
max_time=20000 # seconds

time12=orbit_dataset["time"][1:2]
dt=time12[2]-time12[1]
max_time_index=Int(round(max_time/dt))
time_index=1:1:max_time_index

orbit_time=orbit_dataset["time"][time_index]
orbit_pos_ECI=orbit_dataset["position"][:,:,time_index]*1e3
orbit_vel_ECI=orbit_dataset["velocity"][:,:,time_index]*1e3

try #does file have dcm already?
    global dcm=orbit_dataset["dcm"];
catch #if not generate from Orbits
    dv = orbit_dataset.attrib["epoch"];
    local epoch = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
    global dcm = Orbits.eci_dcm(orbit_time, epoch);
end
p,v=Orbits.ecef_orbitpos(orbit_pos_ECI,orbit_vel_ECI,dcm) # ECI to ECEF

baselines=Orbits.get_perp_baselines(p,v,look_angle)
maxbaselines=reshape(maximum(maximum(baselines,dims=1),dims=2),max_time_index)
maxbaseline,maxb_time_ind=findmax(maxbaselines)
println("Maximum Baseline Along Orbit: ",round.(maxbaseline/1e3,digits=3)," km")
println("Maximum Baseline Time Index: ",maxb_time_ind)
println("Maximum Baseline Time: ",Int((maxb_time_ind-1)*dt)," s")
plotly();plot(orbit_time/60,maxbaselines/1e3,leg=false,size=(1000,600),xlabel="Orbit Time (min)",ylabel="Perpendicular Baseline (km)",title="Orbit: "*orbit_filename*", look angle: "*string(look_angle)*" deg")
