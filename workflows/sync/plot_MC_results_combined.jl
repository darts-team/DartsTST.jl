# begin
using JLD2, Plots, SharedArrays, StatsPlots, PyCall, Statistics, CurveFit, MAT

##
  # filename1 =  "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_1_Measured_sync_pri_sweep_along_n_expandedScene.jld2"
  filename1 =  "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_1_Measured_sync_pri_sweep_along_n.jld2"
  @load filename1 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals1 = ideal_peak_dB .- (20 .* log10.(peaks))
  resVals1 = resolutions .- ideal_res
  PSLRvals1 = PSLRs .- ideal_PSLR
  ISLRvals1 = ISLRs .- ideal_ISLR
  # locVals1 = loc_errors

  filename2 =  "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_pri_sweep_along_n.jld2"
  @load filename2 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals2 = ideal_peak_dB .- (20 .* log10.(peaks))
  resVals2 = resolutions .- ideal_res
  PSLRvals2 = PSLRs .- ideal_PSLR
  ISLRvals2 = ISLRs .- ideal_ISLR
  filename3 =  "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_3_Measured_sync_pri_sweep_along_n.jld2"
  @load filename3 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals3 = ideal_peak_dB .- (20 .* log10.(peaks))
  resVals3 = resolutions .- ideal_res
  PSLRvals3 = PSLRs .- ideal_PSLR
  ISLRvals3 = ISLRs .- ideal_ISLR
  # gr()


  coeff_plot = sync_PRIs.*ones(size(peaks,2))

  file = matopen("syncModule_MonteCarlo_RFSOC_GPSDO_PRI_sweep.mat","w")
  write(file, "peakvals1", collect(peakvals1))
  write(file, "peakvals2", collect(peakvals2))
  write(file, "peakvals3", collect(peakvals3))
  write(file, "resVals1", collect(resVals1))
  write(file, "resVals2", collect(resVals2))
  write(file, "resVals3", collect(resVals3))
  write(file, "PSLRvals1", collect(PSLRvals1))
  write(file, "PSLRvals2", collect(PSLRvals2))
  write(file, "PSLRvals3", collect(PSLRvals3))
  write(file, "ISLRvals1", collect(ISLRvals1))
  write(file, "ISLRvals2", collect(ISLRvals2))
  write(file, "ISLRvals3", collect(ISLRvals3))
  write(file, "sync_PRIs", collect(sync_PRIs))
  close(file)


##
    # filename1 =  "sync data/Lband/GPSDO Measured/syncModule_MonteCarlo_mode_1_Measured_sync_pri_sweep_along_n.jld2"
    # @load filename1 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
    # ideal_peak_dB = 20 .* log10.(ideal_peak)
    # peakvals1 = ideal_peak_dB .- (20 .* log10.(peaks))
    # resVals1 = resolutions .- ideal_res
    # PSLRvals1 = PSLRs .- ideal_PSLR
    # ISLRvals1 = ISLRs .- ideal_ISLR
    # # locVals1 = loc_errors

    # filename2 =  "sync data/Lband/GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_pri_sweep_along_n.jld2"
    # @load filename2 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
    # ideal_peak_dB = 20 .* log10.(ideal_peak)
    # peakvals2 = ideal_peak_dB .- (20 .* log10.(peaks))
    # resVals2 = resolutions .- ideal_res
    # PSLRvals2 = PSLRs .- ideal_PSLR
    # ISLRvals2 = ISLRs .- ideal_ISLR
    # filename3 =  "sync data/Lband/GPSDO Measured/syncModule_MonteCarlo_mode_3_Measured_sync_pri_sweep_along_n.jld2"
    # @load filename3 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
    # ideal_peak_dB = 20 .* log10.(ideal_peak)
    # peakvals3 = ideal_peak_dB .- (20 .* log10.(peaks))
    # resVals3 = resolutions .- ideal_res
    # PSLRvals3 = PSLRs .- ideal_PSLR
    # ISLRvals3 = ISLRs .- ideal_ISLR
    # # gr()
    #
    #
    # coeff_plot = sync_PRIs.*ones(size(peaks,2))
    #
    # file = matopen("syncModule_MonteCarlo_GPSDO_PRI_sweep.mat","w")
    # write(file, "peakvals1", collect(peakvals1))
    # write(file, "peakvals2", collect(peakvals2))
    # write(file, "peakvals3", collect(peakvals3))
    # write(file, "resVals1", collect(resVals1))
    # write(file, "resVals2", collect(resVals2))
    # write(file, "resVals3", collect(resVals3))
    # write(file, "PSLRvals1", collect(PSLRvals1))
    # write(file, "PSLRvals2", collect(PSLRvals2))
    # write(file, "PSLRvals3", collect(PSLRvals3))
    # write(file, "ISLRvals1", collect(ISLRvals1))
    # write(file, "ISLRvals2", collect(ISLRvals2))
    # write(file, "ISLRvals3", collect(ISLRvals3))
    # write(file, "sync_PRIs", collect(sync_PRIs))
    # close(file)
    #
    #
##
    filename1 =  "sync data/Lband/syncModule_MonteCarlo_mode_1_RoseL_sync_pri_sweep_along_n.jld2"
    @load filename1 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
    ideal_peak_dB = 20 .* log10.(ideal_peak)
    peakvals1 = ideal_peak_dB .- (20 .* log10.(peaks))
    resVals1 = resolutions .- ideal_res
    PSLRvals1 = PSLRs .- ideal_PSLR
    ISLRvals1 = ISLRs .- ideal_ISLR
    # locVals1 = loc_errors

    filename2 =  "sync data/Lband/syncModule_MonteCarlo_mode_2_RoseL_sync_pri_sweep_along_n.jld2"
    @load filename2 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
    ideal_peak_dB = 20 .* log10.(ideal_peak)
    peakvals2 = ideal_peak_dB .- (20 .* log10.(peaks))
    resVals2 = resolutions .- ideal_res
    PSLRvals2 = PSLRs .- ideal_PSLR
    ISLRvals2 = ISLRs .- ideal_ISLR
    filename3 =  "sync data/Lband/syncModule_MonteCarlo_mode_3_RoseL_sync_pri_sweep_along_n.jld2"
    @load filename3 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
    ideal_peak_dB = 20 .* log10.(ideal_peak)
    peakvals3 = ideal_peak_dB .- (20 .* log10.(peaks))
    resVals3 = resolutions .- ideal_res
    PSLRvals3 = PSLRs .- ideal_PSLR
    ISLRvals3 = ISLRs .- ideal_ISLR
    # gr()


    coeff_plot = sync_PRIs.*ones(size(peaks,2))

    file = matopen("syncModule_MonteCarlo_RoseL_PRI_sweep.mat","w")
    write(file, "peakvals1", collect(peakvals1))
    write(file, "peakvals2", collect(peakvals2))
    write(file, "peakvals3", collect(peakvals3))
    write(file, "resVals1", collect(resVals1))
    write(file, "resVals2", collect(resVals2))
    write(file, "resVals3", collect(resVals3))
    write(file, "PSLRvals1", collect(PSLRvals1))
    write(file, "PSLRvals2", collect(PSLRvals2))
    write(file, "PSLRvals3", collect(PSLRvals3))
    write(file, "ISLRvals1", collect(ISLRvals1))
    write(file, "ISLRvals2", collect(ISLRvals2))
    write(file, "ISLRvals3", collect(ISLRvals3))
    write(file, "sync_PRIs", collect(sync_PRIs))
    close(file)


## COMPARE OSCILLATORS FOR SIMO MODE
filename1 =  "sync data/Lband/GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_pri_sweep_along_n.jld2"
@load filename1 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
ideal_peak_dB = 20 .* log10.(ideal_peak)
peakvals_gpsdo = ideal_peak_dB .- (20 .* log10.(peaks))
resVals_gpsdo = resolutions .- ideal_res
PSLRvals_gpsdo = PSLRs .- ideal_PSLR
ISLRvals_gpsdo = ISLRs .- ideal_ISLR

filename2 =  "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_pri_sweep_along_n.jld2"
@load filename2 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
ideal_peak_dB = 20 .* log10.(ideal_peak)
peakvals_rfsoc = ideal_peak_dB .- (20 .* log10.(peaks))
resVals_rfsoc = resolutions .- ideal_res
PSLRvals_rfsoc = PSLRs .- ideal_PSLR
ISLRvals_rfsoc = ISLRs .- ideal_ISLR

filename3 =  "sync data/Lband/syncModule_MonteCarlo_mode_2_RoseL_sync_pri_sweep_along_n.jld2"
@load filename3 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
ideal_peak_dB = 20 .* log10.(ideal_peak)
peakvals_rosel = ideal_peak_dB .- (20 .* log10.(peaks))
resVals_rosel = resolutions .- ideal_res
PSLRvals_rosel = PSLRs .- ideal_PSLR
ISLRvals_rosel = ISLRs .- ideal_ISLR

filename4 =  "sync data/Lband/syncModule_MonteCarlo_mode_2_USO_sync_pri_sweep_along_n.jld2"
@load filename4 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
ideal_peak_dB = 20 .* log10.(ideal_peak)
peakvals_uso = ideal_peak_dB .- (20 .* log10.(peaks))
resVals_uso = resolutions .- ideal_res
PSLRvals_uso = PSLRs .- ideal_PSLR
ISLRvals_uso = ISLRs .- ideal_ISLR

filename5 =  "sync data/Lband/syncModule_MonteCarlo_mode_2_USRP_sync_pri_sweep_along_n.jld2"
@load filename5 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
ideal_peak_dB = 20 .* log10.(ideal_peak)
peakvals_usrp = ideal_peak_dB .- (20 .* log10.(peaks))
resVals_usrp = resolutions .- ideal_res
PSLRvals_usrp = PSLRs .- ideal_PSLR
ISLRvals_usrp = ISLRs .- ideal_ISLR

file = matopen("syncModule_MonteCarlo_oscillator_compare_PRI_sweep.mat","w")
write(file, "peakvals_gpsdo", collect(peakvals_gpsdo))
write(file, "resVals_gpsdo", collect(resVals_gpsdo))
write(file, "PSLRvals_gpsdo", collect(PSLRvals_gpsdo))
write(file, "ISLRvals_gpsdo", collect(ISLRvals_gpsdo))

write(file, "peakvals_rfsoc", collect(peakvals_rfsoc))
write(file, "resVals_rfsoc", collect(resVals_rfsoc))
write(file, "PSLRvals_rfsoc", collect(PSLRvals_rfsoc))
write(file, "ISLRvals_rfsoc", collect(ISLRvals_rfsoc))

write(file, "peakvals_rosel", collect(peakvals_rosel))
write(file, "resVals_rosel", collect(resVals_rosel))
write(file, "PSLRvals_rosel", collect(PSLRvals_rosel))
write(file, "ISLRvals_rosel", collect(ISLRvals_rosel))

write(file, "peakvals_uso", collect(peakvals_uso))
write(file, "resVals_uso", collect(resVals_uso))
write(file, "PSLRvals_uso", collect(PSLRvals_uso))
write(file, "ISLRvals_uso", collect(ISLRvals_uso))

write(file, "peakvals_usrp", collect(peakvals_usrp))
write(file, "resVals_usrp", collect(resVals_usrp))
write(file, "PSLRvals_usrp", collect(PSLRvals_usrp))
write(file, "ISLRvals_usrp", collect(ISLRvals_usrp))
write(file, "sync_PRIs", collect(sync_PRIs))

close(file)
