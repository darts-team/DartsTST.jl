
#-start-function--------------------------------------------------------------------------------------------
function osc_psd_twosided(fs,N,a_coeff_db)
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

f = collect( ( (-N/2):1:(N/2-1) ) ).*fs./N

temp = 0.1 .* a_coeff_db
a_coeff = [ 10^temp[1], 10^temp[2], 10^temp[3], 10^temp[4], 10^temp[5] ]

fmin = .1 # minimum PSD frequency. Zeros out PSD below this value


Sphi_2S = a_coeff.*[f.^(-4); f.^(-3); f.^(-2); f.^(-1); f.^0]'

#fmin1 = min(abs(f[f>0])) # original matlab code

#idx = findall(f -> f > 0,f)
#fmin1 = f[idx[1]] # takes first index of frequency vector that is greater than 0

#if fmin1<fmin
#    idx = findall(f -> abs(f) < fmin,f)
#      Sphi_2S[idx] =  a_coeff[1].*fmin^(-4).+a_coeff[2]*fmin^(-3).+a_coeff[3]*fmin^(-2).+a_coeff[4]*fmin^(-1).+a_coeff[5]*fmin^0
#else
#    # conserve power around dc
#    idx = findall(f -> abs(f) < fmin,f)
#    Sphi_2S[idx] = ones(length(idx)) .* sum([a_coeff[1]*fmin.^(-4) a_coeff[2]*fmin.^(-3) a_coeff[3]*fmin.^(-2) a_coeff[4]*fmin.^(-1) a_coeff[5]*fmin.^0 ]).*(fmin/fmin1)
#end

return Sphi_2S, f

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

N = length(Sphi)
f = collect( ( (-N/2) : 1 : (N/2-1) ) ) .* fs ./N

Sphi_1S = zeros(N)

findall(idx->idx>=0, Sphi)
idx = findall(Sphi -> Sphi >= 0,Sphi)
#Sphi_1S[f>=0] = 2 .* Sphi[f>=0]; # one sided with f=0 included

Sphi_1S[idx] = 2* Sphi[idx] # one sided with f=0 included

A_Sphi = sqrt(Sphi_1S.*N)


# random phase vector with amplitude of PSD
phi_u = rand(N) .*2 .* pi
# amp_u = rand(1,N); # random amplitude
amp_u = 1; # constant amplitude

# sequence with random phase & random amplitude of PSD
z = A_Sphi.*exp(1i*phi_u).*amp_u
# z[f==0]=abs(A_Sphi[f==0])

z2 = ifftshift(z)

# random time sequence realization
r = real(fftshift(ifft(z2))) # take real component of time series

# time vector
t = collect(0:(N-1)) * 1/fs

return r,t
end#function
#--end--function-------------------------------------------------------------------------------------------


#-start-function-------------------------------------------------------------------------------------------
"generates a matrix of estimated CRLB values based on relative positions of platforms and Friis transmission"
function getSensorCRLB(pos,N,fc,fs,fbw,master)
# pos:              gives locations of all platforms in network (3 x N_plat x N_time)
# nplat:            number of platforms
# ntimes:           number of sample times in time_vec
# N:                number of samples in waveform (Period (s) = N/fs)
# fs:               sampling frequency
# fbw:              sync waveform bandwidth
# master:           select which Tx is the "master" for the sync

if ndims(pos) == 3
    nplat = szp[2]
    ntimes = szp[3]
elseif ndims(pos) == 2
    nplat = 1
    ntimes = szp[2]
end

k = 1.380649e-23; # Joule/Kelvin -  Boltzmann constant
T = 25 + 273.15; # Kelvin -  let approximate temperature be 25C
Np = k*T*fbw; # (W) noise power in receiver - NOTE: this doesn't consider the noise figure in receiver electronics
c = 299792458; # m/s SOL

# initialize SNR and CRLB matrices
snr_dB = Array{Float64,2}(undef, nplat, ntimes)
sig_crlb = Array{Float64,2}(undef, nplat, ntimes)

for i = 1 : nplat # number of platforms
if i == master
  sig_crlb[i,:] = zeros(ntimes) # ignore oscillator effects if only single platform
  snr_dB[i,:] = zeros(ntimes) # setting SNR to zero just to satisfy output
else

for j = 1 : ntimes # loop over number of sample times
master_pos = pos[:,master,j]

# findall R, distance between the platform and master
receiver_pos = pos[:,i,:]

R = sqrt((receiver_pos[1]-master_pos[1])^2+(receiver_pos[2]-master_pos[2])^2+(receiver_pos[3]-master_pos[3])^2)


  # need Pt, Gt, Gr, POL. These could be calculated in the orbits module using antenna pattern/platform rotations/
  Pt = 1e-3     # let transmit power be 1 mW (needs to be evaluated in trade study, more likely a consequence of orbit design)
  Gt = 1        # Isotropic rx/tx antennas for sync
  Gr = 1        # Isotropic rx/tx antennas for sync
  POL = 1       # Polarization mismatch factor (let antennas be aligned here) - 0 is total misalignment, 1 is total alignment


  # use Friis formula to calculate received power, SNR
  Pr = POL * Pt*Gt*Gr*c^2/(4*pi*fc*R)^2;
  snr_linear = Pr/Np;
  snr_dB[i,j] = 10*log10(snr_linear); # save SNR for debug/trade study

  sig_crlb[i,j] = 1/(sqrt(2*nplat)*(2*pi*fbw)^2*(snr_linear)*(N/fs)*fs) # clk error after sync for nplat sensors
end # if
end #for i
end #for j
return [sig_crlb, snr_dB] # return tuple of CRLB and SNR values
end#function
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
r = zeros(1,length(t))
r[ind1:ind2] = height*ones(1,ind2-ind1+1)
return r
end
#--end--function-------------------------------------------------------------------------------------------

##-start-function-------------------------------------------------------------------------------------------
function sync_effects_on_PSD(Sphi,f_psd,sync_radar_offset,sig_crlb,sync_pri,sdradar_args_fs)
# Function describes PSD of clk phase after finite offset sync
# Sphi is input PSD

Sphi_tilde = (2*Sphi.*(1 - cos(2*pi*sync_radar_offset.*f_psd) ) )

# PSD of AWGN with variance determined by CRLB of sync algorithm
S_n = ( ( (sig_crlb) * sdradar_args_fs*(2*pi)*2/sqrt(2)).^2) .*ones(length(Sphi))

# PSD of clock phase error after synchronization
Sphi_sync = (Sphi_tilde+S_n).*rect(f_psd,-1/sync_pri,1/sync_pri,1) # low pass filter the PSD by the Sync process

return Sphi_sync

end#function
#--end--function-------------------------------------------------------------------------------------------
