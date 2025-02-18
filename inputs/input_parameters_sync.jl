using NCDatasets

use_orbits_flag = true # true if using an orbit file to inform number of platforms
disable_freq_offset = true  # true = no linear phase ramp (ideal osc frequency), false = linear phase ramp error
if !use_orbits_flag
    setNumPlatforms = 3 # manually select number of Rx platforms

sync_pri = 2 # (s) repetition interval of sync

sync_processing_time = 0.001 # processing time between stage 1 and stage 2 sync
sync_signal_len = 1024 # waveform length
sync_fc = 1.25e9 # waveform center frequency
sync_fs = 25e6; # sync receiver sampling rate
sync_fbw = sync_fs # LFM bandwidth

# osc_type = "USO" # putting a oscillator type variable here to auto-name save files
# osc_type = "USRP"
# osc_type = "Wenzel5MHz"
# osc_type = "Wenzel100MHz"
osc_type = "MicroSemi" #Microsemi GPS-3500




#defines oscillator quality. Either leave as single row to use across all platforms, or define values for each platform as a new row
#define oscillator quality and frequency
if osc_type == "USO"
    a_coeff_dB = [-95 -90 -200 -130 -155] # [USO: Krieger]
    f_osc = 10e6 # local oscillator frequency
elseif osc_type == "USRP"
    a_coeff_dB = [-66 -62 -80 -110 -153] # [USRP E312]
    f_osc = 10e6 # local oscillator frequency
elseif osc_type == "Wenzel5MHz"
    a_coeff_dB = [-1000 -128 -1000 -150 -178] # [Wenzel 5MHz oscillator] - NOTE: fractional dB values were rounded up for Wenzel oscillators (to keep as Int64 values)
    f_osc = 5e6 # local oscillator frequency
elseif osc_type == "Wenzel100MHz"
    a_coeff_dB = [-1000 -73 -1000 -104 -181] # [Wenzel 100MHz oscillator]
    f_osc = 100e6 # local oscillator frequency
elseif osc_type == "MicroSemi"
    a_coeff_dB = [-120 -114 -999 -134 -166 ] # [Microsemi GPS-3500 oscillator]
    f_osc = 10e6 # local oscillator frequency #TODO(right center freq?)
end

# here we assume all platforms are the same quality. however, we can redefine the osc_coeffs to have different values. (likely scenario in SIMO mode with master transmitter)
## TODO need to figure out a way to get the number of platforms without hardcoding this orbit file name - use OrbitFileName variable?
# ------this is ugly code that needs to be replaced------
# orbit_dataset=Dataset("inputs/orbitOutput_082020.nc") # Read orbits data in NetCDF format
if use_orbits_flag
    orbit_dataset=Dataset("inputs/"*orbit_filename) 
    orbit_pos=orbit_dataset["position"][:,:,orbit_time_index] # read in position data
    nplat=size(orbit_pos)[2] # number of platforms
else
    nplat = setNumPlatforms
end
##-----------
if size(a_coeff_dB)[1] == 1
    osc_coeffs = repeat(a_coeff_dB,nplat)
end #if

if disable_freq_offset == true # option to remove linear phase drift due to osc frequency offset
    sigma_freq_offsets = zeros(nplat)
else
    sigma_freq_offsets = 1.5e-3 # Hz - std. dev. of the frequency offset of the oscillator. This is the linear phase ramp value
    sigma_freq_offsets = sigma_freq_offsets .* ones(nplat) # convert to matrix form, one value for each oscillator
end



sync_fmin = 1.0 # minimum frequency > 0 in Hz to window PSD -- Float64
# f_osc = 10e6 # local oscillator frequency
sync_clk_fs = 1e3; # sample rate of clock error process
master = 1; # selection of master transmitter for sync (assumes a simplified communication achitecture- all talking with one master platform)


no_sync_flag = false; # if flag == true, no sync is used. flag == false results in normal sync process estimation
## make a struct of important input parameters
#list key parameters in here, they will get passed to most(?) modules
mutable struct keyParameters
    # radar parameters
    mode #1: SAR (ping-pong), 2:SIMO, 3:MIMO
    tx_el # which element transmits for SIMO (max value N)
    fc # center frequency (Hz)
    fp # pulse repetition frequency (Hz)
    # TODO: add others... set defaults

    # clock parameters
    sync_pri
    sync_processing_time
    sync_signal_len
    sync_fc
    sync_fs
    sync_fbw
    sync_fmin
    f_osc
    sync_clk_fs
    master
    osc_coeffs
    sigma_freq_offsets

    no_sync_flag
end

## define instance of parameter structure
parameters = keyParameters(mode,
tx_el,
fc,
fp,
sync_pri,
sync_processing_time,
sync_signal_len,
sync_fc,
sync_fs,
sync_fbw,
sync_fmin,
f_osc,
sync_clk_fs,
master,
osc_coeffs,
sigma_freq_offsets,
no_sync_flag)