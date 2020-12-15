module Sync
include("input_parameters_sync.jl") # clock/sync system parameters
include("input_parameters_1") # radar system parameters
include("sync_helper_functions.jl") # this script houses the functions to describe parts of sync module

#-start-function--------------------------------------------------------------------------------------------
"takes the position vectors and time vector to produce a phase error matrix output"
function main(time_vector, pos )

# - steps -
# 1) determine number of platforms from input position vector
# 2) estimate SNR of inter-platform connections through Friis transmission formula, calculate CRLB?
# 3) determine the PSD of the phase error from system, clock, and oscillator_quality
# 4) generate a set (# of platforms) of random phase sequences using PSDs. Make sure the time vector is sufficiently long (longer than entire time vector input)
# 5) for each platform, interpolate from the internal time base (generated with the phase error sequences) to the global time base (the time_vector input)
# note: we may need to keep the output as a time vector with the time base. Then within the raw data generator, interpolate to the path delay times
# 6) create output matrix comprised of the phase errors for each platform and epoch

# Questions to consider and things to add:
# 1) time since sync: do we need to include a transmit schedule and apply different delay times to each transmitter (in multi tx modes)
# 2) do new phase error curves have to be created every time the sync is completed (sync_pri)?

using Interpolations
using LinearAlgebra
using FFTW

#TODO input time vector from orbits module, rename to time_vector

# find total elapsed time over the course of the orbits
t_elapse = time_vector[end] - time_vector[1]

fs = round(clk_args_N/t_elapse) # upper bound of frequency in the PSD generation. Limited (i.e. selected) by the sync_pri
clk_args_N = ceil(fs*t_elapse)  # Number of sample points in the PSD. Needs to meet the length of min {1/sync_pri,t_elapse}
m = f/clk_args_f_osc;           # frequency up-conversion factor (scale factor from LO to RF)

# verify size of position input
szp = size(pos)
@assert szp[1]==3 "POS needs to be 3 x [Np x Nt] or 3 x Nt"
if ndims(pos) == 3
    nplat = szp[2]
    ntimes = szp[3]
elseif ndims(pos) == 2
    nplat = 1
    ntimes = szp[2]
end

if nplat ==1
phase_err = zeros(size(pos)) # this assumes no sync is used. Does not include any oscillator error effects
return phase_err_s

else # multi platform case

# find CRLB based on SNR of sync connection
sig_crlb = getSensorCRLB(pos,Ns,sdradar_args_N,sdradar_args_fc,sdradar_args_fs,sdradar_args_fbw,master)

# TODO: Does the sig_crlb need to be updated at every sync time (with changing SNR) - yes, we should change the PSDs every SRI


# create matrix of phase error values
phase_err = Array{ComplexF64}(szp[2],szp[3]) # Nplatforms x Ntimes

for i = 1:nplat #for each platform, generate oscillator phase errors at each time point

a_coeff_db = osc_coeffs[i,:] # grab the clock coefficients for each platform
[Sphi,f_psd] = osc_psd_twosided(fs,clk_args_N,a_coeff_db)

#for j = 1 : ntimes # for looping over segments of times with sync_pri period
#Sphi_sync = sync_effects_on_PSD(Sphi,f_psd,sync_radar_offset,sig_crlb[i,j],sync_pri,sdradar_args_fs)

Sphi_sync = sync_effects_on_PSD(Sphi,f_psd,sync_radar_offset,sig_crlb[i],sync_pri,sdradar_args_fs)


#------- GPS only model----------------------------
#TODO GPS only model for clock error corrections? How would this look? (ask Sam)
#---------------------------------------------------

# generate time series of phase error
[r,t] = osc_timeseries_from_psd_twosided(Sphi_sync,sdradar_args_fs)

# interpolate phase error to the measurement times
itp = CubicSplineInterpolation(t, r)    # create interpolant
phase = itp(time_vector)                # interpolates to the time vector given by orbits module
phase_err[i,:] = phase.*m               # scale by frequency up-conversion factor

end # loop over time
end # loop over platforms
end # main function

end #module
