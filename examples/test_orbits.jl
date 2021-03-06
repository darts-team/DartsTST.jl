# Test file to highlight the functions in the orbits.jl module
# Note: run this script where the orbit file exists

include("modules/orbits.jl")

using NCDatasets
using Plots
using Dates

# Read in NetCDF orbit file
ds = Dataset("orbitOutput_082020.nc");

# read in time and position data
read_samps = 1:15000;
orbit_time = ds["time"][read_samps];
orbit_pos  = ds["position"][:,:,read_samps];

# plot the ECI orbit
plot(orbit_time, orbit_pos[:,1,:]', xaxis=("time (sec)"), ylabel=("ECEF position (km)"))

# read in the DCM and convert to ECEF position and plot
dcm = ds["dcm"];
orbit_pos_ecef = Orbits.ecef_orbitpos(orbit_pos, dcm);
# plot the ECEF orbit computed using dcm from NetCDF file
scatter!(orbit_time[1:100:end], orbit_pos_ecef[:,1,1:100:end]', marker = (:+))

# read in epoch, get DCM and convert to ECEF and plot
dv = ds.attrib["epoch"];
epoch = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
dcm = Orbits.eci_dcm(orbit_time, epoch);
orbit_pos_ecef2 = Orbits.ecef_orbitpos(orbit_pos, dcm);
# plot the ECEF orbit computed using DCM based on epoch
scatter!(orbit_time[1:100:end], orbit_pos_ecef2[:,1,1:100:end]', marker = (:circle))


# interpolate orbit to new timebase and plot it
fs = 1e2;
tv = collect(0:1/fs:50).+12.0;
orbit_pos_interp = Orbits.interp_orbit(orbit_time, orbit_pos, tv);
plot(tv, orbit_pos_interp[:,1,:]', xaxis=("time (sec)"), ylabel=("ECEF position (km)"))
ind = findall(x->x .>= tv[1] && x .<= tv[end], orbit_time);
scatter!(orbit_time[ind], orbit_pos[:,1,ind]', marker = (:+))

# get perpendicular baselines for the interpolated orbit and plot them
orbit_vel = ds["velocity"][:,:,read_samps];
bperp = Orbits.get_perp_baselines(orbit_pos, orbit_vel, 30);
plot(orbit_time, reshape(bperp, size(bperp,1)*size(bperp,2), length(orbit_time))', xaxis=("time (sec)"), ylabel=("Perp baseline (km)"))
