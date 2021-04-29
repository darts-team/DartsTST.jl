using NCDatasets

use_orbits_flag = true # true if using an orbit file to inform number of platforms
disable_freq_offset = false # false = frequency mismatch + phase ramp. true = no phase ramp

if !use_orbits_flag
    setNumPlatforms = 3 # manually select number of Rx platforms
end 
sync_pri = .1 # (s) repetition interval of sync

sync_processing_time = 0.001 # processing time between stage 1 and stage 2 sync
sync_signal_len = 1024 # waveform length
sync_fc = 1e9 # waveform center frequency
sync_fs = 25e6; # sync receiver sampling rate
sync_fbw = sync_fs # LFM bandwidth

# osc_type = "USO" # putting a oscillator type variable here to auto-name save files
osc_type = "USRP"

#defines oscillator quality. Either leave as single row to use across all platforms, or define values for each platform as a new row

if osc_type == "USO"
    a_coeff_dB = [-95 -90 -200 -130 -155] # [USO: Krieger]
elseif osc_type == "USRP"
    a_coeff_dB = [-28 -40 -200 -130 -155] # [USRP E312]
end

# here we assume all platforms are the same quality. however, we can redefine the osc_coeffs to have different values. (likely scenario in SIMO mode with master transmitter)
## TODO need to figure out a way to get the number of platforms without hardcoding this orbit file name - use OrbitFileName variable?
# ------this is ugly code that needs to be replaced------
# orbit_dataset=Dataset("inputs/orbitOutput_082020.nc") # Read orbits data in NetCDF format
if use_orbits_flag
    orbit_dataset=Dataset(orbit_filename) # Read orbits data in NetCDF format
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
    # value above comes from: (1.5e-3 * 2*pi * 1sec) * (180/pi rad/deg) ~= .54 deg/s linear phase drift on LO (roughly matching Krieger paper)
    # at RF: .54 deg/s * (1GHz/10MHz upconversion) = 54 deg/s linear phase drift
    sigma_freq_offsets = sigma_freq_offsets .* ones(nplat) # convert to matrix form, one value for each oscillator
end



sync_fmin = 0.01 # minimum frequency > 0 in Hz to window PSD
f_osc = 10e6 # local oscillator frequency
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
