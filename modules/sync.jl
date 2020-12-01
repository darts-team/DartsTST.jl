module Sync
include("input_parameters_sync.jl")

#-start-function--------------------------------------------------------------------------------------------
"takes the position vectors and time vector to produce a phase error matrix output"
function get_phase_error(time_vector (slow time), pos )

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

# find total elapsed time
t_elapse = time_vector[end] - time_vector[1]

fs = round(clk_args_N/t_elapse)

#eci_pos = 3 x N_plat x N_time

szp = size(pos)
@assert szp[1]==3 "POS needs to be 3 x [Np x Nt] or 3 x Nt"
if ndims(pos) == 3
    nplat = szp[2]
elseif ndims(pos) == 2
    nplat = 1
end

sig_crlb = getSensorCRLB(pos,N,,fc,fs)
#find CRLB based on SNR of sync connection
getSensorCRLB(pos,Ns,clk_args_N,sdradar_args_fc,sdradar_args_fs,sdradar_args_fbw,master)

# define PSD array
put matrix here [N x num freq points(clk_args_N)]


for i = 1:nplat #for each platform, generate oscillator phase errors at each time point

a_coeff_db = osc_coeffs[i,:]
[Sphi_2S,f] = osc_psd_twosided(fs,clk_args_N,a_coeff_db)




#-------- sync effects on PSD------------------------
# PSD of clk phase after finite offset sync
Sphi_tilde = (2*Sphi.*(1 - cos(2*pi*sync_radar_offset*f_psd)))

# PSD of AWGN with variance determined by CRLB of sync algorithm
S_n = (((sig_crlb[i,j])*sdradar_args_fs*(2*pi)*2./sqrt(2)).^2).*ones(size(Sphi))

# PSD of clock phase error after synchronization
Sphi_sync = (Sphi_tilde+S_n).*rect(f_psd,-1/sync_pri,1/sync_pri,1)

#------- GPS only model----------------------------
%TODO GPS only model for clock error corrections?
#---------------------------------------------------

end #for

m = f/clk_args_f_osc; # frequency up-conversion factor (scale factor from LO to RF)


end # end get_sync_phase
#--end--function-------------------------------------------------------------------------------------------


#-start-function--------------------------------------------------------------------------------------------
function [Sphi_2S,f] = osc_psd_twosided(fs,N,a_coeff_db)
#   Generate PSD of clock phase error()
#INPUTS
#     fs = 2000; # max PSD frequency [sample rate of clock phase error process]
#     N = 600000; # number of PSD sample points
#     a_coeff_db = [-28 -40 -200 -130 -155]; # coefficients of the noise characteristic asymptotes. (optional)
#     fmin = .1; # minimum frequency .> 0 in Hz to window PSD [optional, default = .1 Hz].
#
# OUTPUTS
#     Sphi_2S: # two sided oscillator phase error PSD
#     f:    # frequency vector
#     Sphi_2S: # one sided oscillator phase error PSD [optional]

# for 3 seconds: use fs = 200000; N = 600000
# for 30 seconds: use fs = 20000; N = 600000
# for 300 seconds: use fs = 2000; N = 600000
# for 3000 seconds: use fs = 2000; N = 6000000
# above comments from Sam - use fs = N / time ?

f = ((-N/2):1:(N/2-1))*fs/N

a_coeff = 10.^(.1*a_coeff_db)

fmin = .1
#if (nargin>3)
#    fmin = varargin[2]
#end

# oscillator phase error PSD
Sphi_2S = a_coeff*[f.^(-4) f.^(-3) f.^(-2) f.^(-1) f.^0]'; # two sided
# Sphi_2S[f==0] = numel(Sphi_2S)
# Sphi_2S[f==0] = 1
Sphi_2S[f==0] = 0
# try flattening out lower frequencies
# Sphi_2S[abs(f)<fmin] = 1

fmin1 = min(abs(f[f>0]))
if (fmin1<fmin)
    Sphi_2S[abs(f)<fmin] = a_coeff*[fmin.^(-4) fmin.^(-3) fmin.^(-2) fmin.^(-1) fmin.^0]'
#     Sphi_2S[f==0] = 0
else()
#     Sphi_2S[f==0] = (a_coeff*[fmin1.^(-4) fmin1.^(-3) fmin1.^(-2) fmin1.^(-1) fmin1.^0]')*(1)
    # conserve power around dc
    Sphi_2S[abs(f)<fmin] = (a_coeff*[fmin.^(-4) fmin.^(-3) fmin.^(-2) fmin.^(-1) fmin.^0]')*(fmin/fmin1)
#     Sphi_2S[f==0] = (a_coeff*[fmin1.^(-4) fmin1.^(-3) fmin1.^(-2) fmin1.^(-1) fmin1.^0]')*(fmin1/fmin)
end
# Sphi_2S[abs(f)<fmin] = a_coeff*[fmin.^(-4) fmin.^(-3) fmin.^(-2) fmin.^(-1) fmin.^0]'
# Sphi_2S[f==0] = 0

Sphi_1S = zeros(size(Sphi_2S))
Sphi_1S[f>=0] = 2*Sphi_2S[f>=0]; # one sided

end
end #function
#--end--function-------------------------------------------------------------------------------------------


#-start-function--------------------------------------------------------------------------------------------
function osc_timeseries_from_psd_twosided(Sphi,fs)
#   Generate time series from two sided PSD of clock phase error.
#   Will first calculate the one sided PSD which = 0 for f<0
#INPUTS
#     Sphi: # 2 sided oscillator phase error PSD
#     fs = 2000; # max PSD frequency [sample rate of clock phase error process]
# OUTPUTS
#     r:    # random time sequence realization
#     t = (0:(N-1))*1/fs; # time vector

N = numel(Sphi)
f = ((-N/2):1:(N/2-1))*fs/N
Sphi_1S = zeros(size(Sphi))
Sphi_1S[f>=0] = 2*Sphi[f>=0]; # one sided with f=0 included
# Sphi_1S[f>0] = 2*Sphi[f>0]; # one sided

# A_Sphi = sqrt(Sphi_1S)
A_Sphi = sqrt(Sphi_1S.*N)
# A_Sphi = sqrt(Sphi_1S).*N
# A_Sphi = sqrt(Sphi_1S/N)
# A_Sphi = sqrt(Sphi_1S)./N

# random phase vector with amplitude of PSD
phi_u = rand(1,N)*2*pi
# amp_u = rand(1,N); # random amplitude
amp_u = 1; # constant amplitude

# sequence with random phase & random amplitude of PSD
z = A_Sphi.*exp(1i*phi_u).*amp_u
# z[f==0]=abs(A_Sphi[f==0])

z2 = ifftshift(z)

# random time sequence realization
r = fftshift(ifft(z2))

# time vector
t = (0:(N-1))*1/fs

end
return [r,t]
end
#--end--function-------------------------------------------------------------------------------------------


#-start-function-------------------------------------------------------------------------------------------
function getSensorCRLB(pos,Ns,N,fc,fs,fbw,master)
# pos:              gives locations of all platforms in network
# Ns:               number of platforms
# M:                number of sample times in time_vec
# N:                number of samples in waveform (Period (s) = N/fs)
# fs:               sampling frequency
# fbw:              sync waveform bandwidth
# master:           select which Tx is the "master" for the sync


k = 1.380649e-23; # Joule/Kelvin -  Boltzmann constant
T = 25 + 273.15; # Kelvin -  let approximate temperature be 25C
Np = k*T*BW; # (W) noise power in receiver - NOTE: this doesn't consider the noise figure in receiver electronics
c = 299792458; # m/s SOL

# initialize SNR and CRLB matrices
snr_dB = Array{Float64,2}(undef, M, Ns)
crlb_var_syncclk = Array{Float64,2}(undef, M, Ns)

for i = 1 : M # loop over number of sample times

for j = 1 : Ns
  # need Pt Gt Gr
  # find R, distance between the platform and master
  master_pos =
  R = sqrt()

  # use Friis formula to calculate received power, SNR
  Pr = POL * Pt*Gt*Gr*c^2/(4*pi*fc*R)^2;
  snr_linear = Pr/Np;
  snr_dB[i,j] = 10*log10(snr_linear); # save SNR for debug

  crlb_var_syncclk[i,j] = 1./(sqrt(2*Ns)*(2*pi*fbw)^2*((10.^(.1*snr_db[i,j]))*(N/fs)*fs)) # clk error after sync for Ns sensors
end #for

return sig_crlb
#--end--function-------------------------------------------------------------------------------------------

#-start-function-------------------------------------------------------------------------------------------
function rect(t,t1,t2,height) # does Julia have a built in rect function?
# Unit energy rect pulse either from -T/2 to T/2 with height 1/sqrt(T)
# | from t1 to t2 with height 1/sqrt(t2-t1)
# example: r=rect(t,tmin,tmax,1) # rect from tmin to tmax with height 1

t1 = max(t[1],t1)
t2 = min(t[end],t2)
T = abs(t2-t1)

delta = (t[end]-t[1])/(numel(t))
ind1 = floor((t1-t[1])/delta)+1
ind2 = floor((t2-t[1])/delta)
r = zeros(1,numel(t))
r[ind1:ind2] = height*ones(1,ind2-ind1+1)
return r
end
#--end--function-------------------------------------------------------------------------------------------

end #module
