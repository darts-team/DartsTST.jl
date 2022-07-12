# # plot results of Sync PRI Monte Carlo Study
# using JLD2, Plots, SharedArrays, StatsPlots, PyCall, Statistics
# 
# # define a couple functions to filter NaNs and make plots
# function filterNaNs2D(vals)
#   #this function filters out the NaNs, but returns an array of vectors. Unfortunately, not helpful for plotting boxplots.
#   np = pyimport("numpy")
#   mask = np.isnan(vals) # finds NaNs
#   mask = [!i for i in mask] # find non-nan indices
# 
#   valsF = Array{Vector{Float64},1}(undef, size(mask,1))
#   for i = 1 : size(mask,1)
#     valsF[i] = vals[i, mask[i,:]]
#   end#for
# 
#   return valsF
# end#function
# 
# function makeBoxPlot(DSSs,vals,titleString::String,saveFilename::String)
#   len = size(vals,1)
# 
#   boxplot(delays[1].*ones(len),vals[1], title = titleString,
#   xlabel = "Delay Since Last Sync (s)", legend = false)
# 
#   for i = 2 : len - 1
#     boxplot!(delays[i].*ones(len),vals[i], title = titleString,
#     xlabel = "Delay Since Last Sync (s)", legend = false)
#   end#for
# 
#   display(boxplot!(delays[len].*ones(len),vals[len], title = titleString,
#   xlabel = "Delay Since Last Sync (s)", legend = false))
# 
#   savefig(saveFilename)
# end#function
# 
# 
# # run the rest of the code, calling an input .jld2 file with the monte carlo data stored
# begin 
# 
# 
#   # filename = "sync data/Lband/syncModule_MonteCarlo_mode_1_Measured_sync_DSS_sweep.jld2"
#   filename = "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_DSS_sweep.jld2"
# 
#           #TODO find a way to make log x axis: try "xscale=:log10, xticks=(1:10, 1:10))" or something to set x axis ticks as needed
# 
#   @load filename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors delays
# 
#   gr()
# 
#   #matrices are either (numDSS x numTrials) or (3 x numDSS x numTrials)
#   num2plot = 14 #allows us to cut off later delay values if too big
#   #for 3 axis values
#   # peaks = peaks[1:num2plot,:]; resolutions = resolutions[:,1:num2plot,:]; PSLRs = PSLRs[:,1:num2plot,:];
#   # ISLRs = ISLRs[:,1:num2plot,:]; loc_errors = loc_errors[:,1:num2plot,:]; delays = delays[1:num2plot]';
#   # for 1 axis values
#   peaks = peaks[1:num2plot,:]; resolutions = resolutions[1:num2plot,:]; PSLRs = PSLRs[1:num2plot,:];
#   ISLRs = ISLRs[1:num2plot]; loc_errors = loc_errors[1:num2plot,:]; delays = delays[1:num2plot]';
# 
# 
#   DSS_plot = delays  .* ones(size(peaks,2))
# 
# ## trying box plots # ---- Peak loss ------
# 
#   ideal_peak_dB = 20 .* log10.(ideal_peak)
#   peakvals = ideal_peak_dB .- (20 .* log10.(peaks))
#   display(boxplot(DSS_plot, peakvals', title = "Peak Power Loss (dB)",
#     xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loss.png")
#   savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loss.png")
# 
#   # ---- PSLR ----
#   vals1 = filterNaNs2D(PSLRs[1,:,:] .- ideal_PSLR[1])
#   vals2 = filterNaNs2D(PSLRs[2,:,:] .- ideal_PSLR[2])
#   vals3 = filterNaNs2D(PSLRs[3,:,:] .- ideal_PSLR[3])
# 
# 
#   makeBoxPlot(delays,vals1,"PSLR Change Along-Track (dB)","RFSoC_meas_SIMO_DSS_along_n_sweep_pslr_change_AT.png")
#   makeBoxPlot(delays,vals2,"PSLR Change Cross-Track (dB)","RFSoC_meas_SIMO_DSS_along_n_sweep_pslr_change_CT.png")
#   makeBoxPlot(delays,vals3,"PSLR Change Height (dB)","RFSoC_meas_SIMO_DSS_along_n_sweep_pslr_change_hgt.png")
# 
#   # makeStdevMeanPlot(delays, vals1,"PSLR Change Cross-Track (dB)", true)
#   # display(boxplot(delays,vals1',title = "PSLR Change Cross-Track (dB)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_pslr_change_CT.png")
#   # display(boxplot(delays,vals2',title = "PSLR Change Along-Track (dB)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_pslr_change_AT.png")
#   # display(boxplot(delays,vals3',title = "PSLR Change Height (dB)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_pslr_change_hgt.png")
# 
# 
#   # ---- ISLR ------
#   vals1 = filterNaNs2D(ISLRs[1,:,:] .- ideal_ISLR[1])
#   vals2 = filterNaNs2D(ISLRs[2,:,:] .- ideal_ISLR[2])
#   vals3 = filterNaNs2D(ISLRs[3,:,:] .- ideal_ISLR[3])
#   makeBoxPlot(delays,vals1,"ISLR Change Along-Track (dB)","RFSoC_meas_SIMO_DSS_along_n_sweep_islr_change_AT.png")
#   makeBoxPlot(delays,vals2,"ISLR Change Cross-Track (dB)","RFSoC_meas_SIMO_DSS_along_n_sweep_islr_change_CT.png")  
#   makeBoxPlot(delays,vals3,"ISLR Change Height (dB)","RFSoC_meas_SIMO_DSS_along_n_sweep_islr_change_hgt.png")
# 
#   # display(boxplot(delays,vals1',title = "ISLR Change Cross-Track (dB)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_islr_change_CT.png")
#   # display(boxplot(delays,vals2',title = "ISLR Change Along-Track (dB)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_islr_change_AT.png")
#   # display(boxplot(delays,vals3',title = "ISLR Change Height (dB)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_islr_change_hgt.png")
# 
#   # ---- Resolution ------
#   vals1 = filterNaNs2D(resolutions[1,:,:] .- ideal_res[1])
#   vals2 = filterNaNs2D(resolutions[2,:,:] .- ideal_res[2])
#   vals3 = filterNaNs2D(resolutions[3,:,:] .- ideal_res[3])
# 
#   makeBoxPlot(delays,vals1,"Resolution Change Along-Track (m)","RFSoC_meas_SIMO_DSS_along_n_sweep_resolution_change_AT.png")
#   makeBoxPlot(delays,vals2,"Resolution Change Cross-Track (m)","RFSoC_meas_SIMO_DSS_along_n_sweep_resolution_change_CT.png")
#   makeBoxPlot(delays,vals3,"Resolution Change Height (m)","RFSoC_meas_SIMO_DSS_along_n_sweep_resolution_change_hgt.png")
# 
# 
#   # display(boxplot(delays,vals1',title = "Resolution Change Cross-Track",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_resolution_change_CT.png")
#   # display(boxplot(delays,vals2',title = "Resolution Change Along-Track (m)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_resolution_change_AT.png")
#   # display(boxplot(delays,vals3',title = "Resolution Change Height (m)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_resolution_change_hgt.png")
# 
#   # ---- Peak location error ------
#   vals1 = filterNaNs2D(loc_errors[1,:,:])
#   vals2 = filterNaNs2D(loc_errors[2,:,:])
#   vals3 = filterNaNs2D(loc_errors[3,:,:])
#   makeBoxPlot(delays,vals1,"Peak Location Error Along-Track (m)","RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_AT.png")
#   makeBoxPlot(delays,vals2,"Peak Location Error Cross-Track (m)","RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_CT.png")
#   makeBoxPlot(delays,vals3,"Peak Location Error Height (m)","RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_hgt.png")
# 
#   # display(boxplot(DSS_plot, vals1', title = "Peak Location Error Cross-Track (m)",
#   #  xlabel = "Delay Since Last Sync (s)", legend = false))
#   #  savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_lat.png")
#   # display(boxplot(DSS_plot, vals2', title = "Peak Location Error Along-Track (m)",
#   #   xlabel = "Delay Since Last Sync (s)", legend = false))
#   #   savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_lon.png")
#   # display(boxplot(DSS_plot, vals3', title = "Peak Location Error Height (m)",
#   #    xlabel = "Delay Since Last Sync (s)", legend = false))
#   # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_hgt.png")
# 
# end#begin
# 
# 
# 
# # i=6
# # mean(peakvals[i,:])
# # mean(PSLRs[1,i,:] .- ideal_PSLR[1])
# # mean(PSLRs[2,i,:] .- ideal_PSLR[2])
# # mean(PSLRs[3,i,:] .- ideal_PSLR[3])
# # mean(ISLRs[1,i,:] .- ideal_ISLR[1])
# # mean(ISLRs[2,i,:] .- ideal_ISLR[2])
# # mean(ISLRs[3,i,:] .- ideal_ISLR[3])
# # mean(resolutions[1,i,:] .- ideal_res[1])
# # mean(resolutions[2,i,:] .- ideal_res[2])
# # mean(resolutions[3,i,:] .- ideal_res[3])
# # mean(loc_errors[1,i,:])
# # mean(loc_errors[2,i,:])
# # mean(loc_errors[3,i,:])
# 
# 

# plot results of Sync PRI Monte Carlo Study - along n
using JLD2, Plots, SharedArrays, StatsPlots, PyCall, Statistics, MAT

# define a couple functions to filter NaNs and make plots
function filterNaNs2D(vals)
  #this function filters out the NaNs, but returns an array of vectors. Unfortunately, not helpful for plotting boxplots.
  np = pyimport("numpy")
  mask = np.isnan(vals) # finds NaNs
  mask = [!i for i in mask] # find non-nan indices

  valsF = Array{Vector{Float64},1}(undef, size(mask,1))
  for i = 1 : size(mask,1)
    valsF[i] = vals[i, mask[i,:]]
  end#for

  return valsF
end#function

function makeBoxPlot(DSSs,vals,titleString::String,saveFilename::String)
  len = size(vals,1)

  boxplot(delays[1].*ones(len),vals[1], title = titleString,
  xlabel = "Delay Since Last Sync (s)", legend = false)

  for i = 2 : len - 1
    boxplot!(delays[i].*ones(len),vals[i], title = titleString, 
    xlabel = "Delay Since Last Sync (s)", legend = false)
  end#for

  display(boxplot!(delays[len].*ones(len),vals[len], title = titleString,
  xlabel = "Delay Since Last Sync (s)", legend = false))

  savefig(saveFilename)
end#function


# run the rest of the code, calling an input .jld2 file with the monte carlo data stored
begin 

  # filename = "sync data/Lband/syncModule_MonteCarlo_mode_1_Measured_sync_DSS_sweep.jld2"
  filename = "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_DSS_sweep_along_n.jld2"

          #TODO find a way to make log x axis: try "xscale=:log10, xticks=(1:10, 1:10))" or something to set x axis ticks as needed

  @load filename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors delays

  gr()

  #matrices are either (numDSS x numTrials) or (3 x numDSS x numTrials)
  num2plot = 14 #allows us to cut off later delay values if too big
  peaks = peaks[1:num2plot,:]; resolutions = resolutions[1:num2plot,:]; PSLRs = PSLRs[1:num2plot,:];
  ISLRs = ISLRs[1:num2plot]; loc_errors = loc_errors[1:num2plot,:]; delays = delays[1:num2plot]';


  DSS_plot = delays  .* ones(size(peaks,2))

## trying box plots # ---- Peak loss ------

  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals = ideal_peak_dB .- (20 .* log10.(peaks))
  display(boxplot(DSS_plot, peakvals', title = "Peak Power Loss (dB)",
    xlabel = "Delay Since Last Sync (s)", legend = false))
  # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loss.png")
  savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loss.png")

  # ---- PSLR ----
  vals1 = filterNaNs2D(PSLRs .- ideal_PSLR)
  makeBoxPlot(delays,vals1,"PSLR Change (dB)","RFSoC_meas_SIMO_DSS_along_n_sweep_pslr_change_N.png")
  
  # plot(delays, mean(peakvals,dims=2)')

  # ---- ISLR ------
  vals1 = filterNaNs2D(ISLRs .- ideal_ISLR)
  makeBoxPlot(delays,vals1,"ISLR Change Along-Track (dB)","RFSoC_meas_SIMO_DSS_along_n_sweep_islr_change_N.png")

  # ---- Resolution ------
  vals1 = filterNaNs2D(resolutions .- ideal_res)
  makeBoxPlot(delays,vals1,"Resolution Change (m)","RFSoC_meas_SIMO_DSS_along_n_sweep_resolution_change_N.png")
  
  # ---- Peak location error ------
  # vals1 = filterNaNs2D(loc_errors[:,:])
  # makeBoxPlot(delays,vals1,"Peak Location Error (m)","RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_N.png")
  

  # display(boxplot(DSS_plot, vals1', title = "Peak Location Error Cross-Track (m)",
  #  xlabel = "Delay Since Last Sync (s)", legend = false))
  #  savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_lat.png")
  # display(boxplot(DSS_plot, vals2', title = "Peak Location Error Along-Track (m)",
  #   xlabel = "Delay Since Last Sync (s)", legend = false))
  #   savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_lon.png")
  # display(boxplot(DSS_plot, vals3', title = "Peak Location Error Height (m)",
  #    xlabel = "Delay Since Last Sync (s)", legend = false))
  # savefig("RFSoC_meas_SIMO_DSS_along_n_sweep_peak_loc_hgt.png")


  filename1 = "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_1_Measured_sync_DSS_sweep_along_n.jld2"
  filename2 = "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_2_Measured_sync_DSS_sweep_along_n.jld2"
  filename3 = "sync data/Lband/RFSOC+GPSDO Measured/syncModule_MonteCarlo_mode_3_Measured_sync_DSS_sweep_along_n.jld2"
  
  @load filename1 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors delays
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals1 = ideal_peak_dB .- (20 .* log10.(peaks))
  PSLRs1 = PSLRs .- ideal_PSLR; ISLRs1 = ISLRs .- ideal_ISLR; res1 = resolutions .- ideal_res; loc_errors1 = loc_errors; delays1 = delays
  
  @load filename2 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors delays
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals2 = ideal_peak_dB .- (20 .* log10.(peaks))
  PSLRs2 = PSLRs .- ideal_PSLR; ISLRs2 = ISLRs .- ideal_ISLR; res2 = resolutions .- ideal_res; loc_errors2 = loc_errors; delays2 = delays

  @load filename3 peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors delays
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals3 = ideal_peak_dB .- (20 .* log10.(peaks))
  PSLRs3 = PSLRs .- ideal_PSLR; ISLRs3 = ISLRs .- ideal_ISLR; res3 = resolutions .- ideal_res; loc_errors3 = loc_errors; delays3 = delays


  file = matopen("syncModule_MonteCarlo_DSS_along_n_all_Modes_LBand.mat","w")
  write(file, "peakvals1", collect(peakvals1))
  write(file, "peakvals2", collect(peakvals2))
  write(file, "peakvals3", collect(peakvals3))
  write(file, "PSLRs1", collect(PSLRs1))
  write(file, "PSLRs2", collect(PSLRs2))
  write(file, "PSLRs3", collect(PSLRs3))
  write(file, "ISLRs1", collect(ISLRs1))
  write(file, "ISLRs2", collect(ISLRs2))
  write(file, "ISLRs3", collect(ISLRs3))
  write(file, "res1", collect(res1))
  write(file, "res2", collect(res2))
  write(file, "res3", collect(res3))
  write(file, "loc_errors1", collect(loc_errors1))
  write(file, "loc_errors2", collect(loc_errors2))
  write(file, "loc_errors3", collect(loc_errors3))
  write(file, "delays1", collect(delays1))
  write(file, "delays2", collect(delays2))
  write(file, "delays3", collect(delays3))
  close(file)
  
end#begin
