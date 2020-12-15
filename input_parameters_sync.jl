# add new clock based parameters
# need to have selection of oscillator_quality for each receiver (for non-uniform constellation)
include("input_parameters_1.jl") # radar system parameters

sync_pri = 1; # repition interval of sync

sync_processing_time = .0005; # processing time between stage 1 and stage 2 sync

sync_radar_offset = 0.001 # time between sync epoch and radar pulse
sdradar_args_N = 1024 # waveform length
sdradar_args_fc = 1e9 # waveform center frequency
sdradar_args_fs = 50e6; # sync receiver sampling rate
sdradar_args_fbw = sdradar_args_fs # LFM bandwidth

a_coeff_db = [-28 -40 -200 -130 -155] # defines oscillator quality

# here we assume all platforms are the same quality. however, we can redefine the osc_coeffs to have different values. (likely scenario in SIMO mode with master transmitter)
nplat = length(p_θ)
osc_coeffs = repeat(a_coeff_db,nplat)

clk_args_fmin = 0.01 # minimum frequency > 0 in Hz to window PSD
#clk_args_fs = 1e3; # 
clk_args_f_osc = 10e6 # local oscillator frequency
clk_args_N = 60000; # number of sample points in PSD
master = 1; # selection of master transmitter for sync