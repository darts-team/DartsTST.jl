# plot results of Sync PRI Monte Carlo Study
using JLD2, Plots, SharedArrays, StatsPlots, PyCall, Statistics

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

function makeBoxPlot(SRIs,vals,titleString::String,saveFilename::String)
  len = size(vals,1)
  
  boxplot(sync_PRIs[1].*ones(len),vals[1], title = titleString,
  xlabel = "Sync Repetition Interval (s)", legend = false)
  
  for i = 2 : len - 1
    boxplot!(sync_PRIs[i].*ones(len),vals[i], title = titleString,
    xlabel = "Sync Repetition Interval (s)", legend = false)
  end#for
  
  display(boxplot!(sync_PRIs[len].*ones(len),vals[len], title = titleString,
  xlabel = "Sync Repetition Interval (s)", legend = false))
  
  savefig(saveFilename)
end#function


# run the rest of the code, calling an input .jld2 file with the monte carlo data stored
begin 
  
  # filename = "sync data/syncModule_MonteCarlo_mode_2_USO_sync_pri_sweep_wFreq.jld2"# this uses sync_data/*, a folder which isn't used on the master branch
  # filename = "sync data/syncModule_MonteCarlo_mode_2_USO_NOsync_pri_sweep_noFreq.jld2"
  filename = "sync data/syncModule_MonteCarlo_mode_2_USO_sync_pri_sweep_noFreq.jld2"
  
  @load filename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs
  
  gr()

  #matrices are either (numSRI x numTrials) or (3 x numSRI x numTrials)

  SRI_plot = sync_PRIs.*ones(size(peaks,2))

## trying box plots
  # ---- Peak loss ------
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals = ideal_peak_dB .- (20 .* log10.(peaks))
  display(boxplot(SRI_plot, peakvals', title = "Peak Power Loss (dB)",
    xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_peak_loss.png")
  savefig("USO_SIMO_SRI_sweep_peak_loss.png")

  # ---- PSLR ----
  vals1 = filterNaNs2D(PSLRs[1,:,:] .- ideal_PSLR[1])
  vals2 = filterNaNs2D(PSLRs[2,:,:] .- ideal_PSLR[2])
  vals3 = filterNaNs2D(PSLRs[3,:,:] .- ideal_PSLR[3])

  
  makeBoxPlot(sync_PRIs,vals1,"PSLR Change Along-Track (dB)","USO_SIMO_SRI_sweep_with_phase_ramp_pslr_change_AT.png")
  makeBoxPlot(sync_PRIs,vals2,"PSLR Change Cross-Track (dB)","USO_SIMO_SRI_sweep_with_phase_ramp_pslr_change_CT.png")
  makeBoxPlot(sync_PRIs,vals3,"PSLR Change Height (dB)","USO_SIMO_SRI_sweep_with_phase_ramp_pslr_change_hgt.png")
    
  # makeStdevMeanPlot(sync_PRIs, vals1,"PSLR Change Cross-Track (dB)", true)
  # display(boxplot(sync_PRIs,vals1',title = "PSLR Change Cross-Track (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_pslr_change_CT.png")
  # display(boxplot(sync_PRIs,vals2',title = "PSLR Change Along-Track (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_pslr_change_AT.png")
  # display(boxplot(sync_PRIs,vals3',title = "PSLR Change Height (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_pslr_change_hgt.png")


  # ---- ISLR ------
  vals1 = filterNaNs2D(ISLRs[1,:,:] .- ideal_ISLR[1])
  vals2 = filterNaNs2D(ISLRs[2,:,:] .- ideal_ISLR[2])
  vals3 = filterNaNs2D(ISLRs[3,:,:] .- ideal_ISLR[3])
  makeBoxPlot(sync_PRIs,vals1,"ISLR Change Along-Track (dB)","USO_SIMO_SRI_sweep_with_phase_ramp_islr_change_AT.png")
  makeBoxPlot(sync_PRIs,vals2,"ISLR Change Cross-Track (dB)","USO_SIMO_SRI_sweep_with_phase_ramp_islr_change_CT.png")  
  makeBoxPlot(sync_PRIs,vals3,"ISLR Change Height (dB)","USO_SIMO_SRI_sweep_with_phase_ramp_islr_change_hgt.png")
  
  # display(boxplot(sync_PRIs,vals1',title = "ISLR Change Cross-Track (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_islr_change_CT.png")
  # display(boxplot(sync_PRIs,vals2',title = "ISLR Change Along-Track (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_islr_change_AT.png")
  # display(boxplot(sync_PRIs,vals3',title = "ISLR Change Height (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_islr_change_hgt.png")

  # ---- Resolution ------
  vals1 = filterNaNs2D(resolutions[1,:,:] .- ideal_res[1])
  vals2 = filterNaNs2D(resolutions[2,:,:] .- ideal_res[2])
  vals3 = filterNaNs2D(resolutions[3,:,:] .- ideal_res[3])
  
  makeBoxPlot(sync_PRIs,vals1,"Resolution Change Along-Track (m)","USO_SIMO_SRI_sweep_with_phase_ramp_resolution_change_AT.png")
  makeBoxPlot(sync_PRIs,vals2,"Resolution Change Cross-Track (m)","USO_SIMO_SRI_sweep_with_phase_ramp_resolution_change_CT.png")
  makeBoxPlot(sync_PRIs,vals3,"Resolution Change Height (m)","USO_SIMO_SRI_sweep_with_phase_ramp_resolution_change_hgt.png")
  
  
  # display(boxplot(sync_PRIs,vals1',title = "Resolution Change Cross-Track",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_resolution_change_CT.png")
  # display(boxplot(sync_PRIs,vals2',title = "Resolution Change Along-Track (m)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_resolution_change_AT.png")
  # display(boxplot(sync_PRIs,vals3',title = "Resolution Change Height (m)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_resolution_change_hgt.png")

  # ---- Peak location error ------
  vals1 = filterNaNs2D(loc_errors[1,:,:])
  vals2 = filterNaNs2D(loc_errors[2,:,:])
  vals3 = filterNaNs2D(loc_errors[3,:,:])
  makeBoxPlot(sync_PRIs,vals1,"Peak Location Error Along-Track (m)","USO_SIMO_SRI_sweep_with_phase_ramp_peak_loc_AT.png")
  makeBoxPlot(sync_PRIs,vals2,"Peak Location Error Cross-Track (m)","USO_SIMO_SRI_sweep_with_phase_ramp_peak_loc_CT.png")
  makeBoxPlot(sync_PRIs,vals3,"Peak Location Error Height (m)","USO_SIMO_SRI_sweep_with_phase_ramp_peak_loc_hgt.png")
  
  # display(boxplot(SRI_plot, vals1', title = "Peak Location Error Cross-Track (m)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  #  savefig("USO_SIMO_SRI_sweep_with_phase_ramp_peak_loc_lat.png")
  # display(boxplot(SRI_plot, vals2', title = "Peak Location Error Along-Track (m)",
  #   xlabel = "Sync Repetition Interval (s)", legend = false))
  #   savefig("USO_SIMO_SRI_sweep_with_phase_ramp_peak_loc_lon.png")
  # display(boxplot(SRI_plot, vals3', title = "Peak Location Error Height (m)",
  #    xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("USO_SIMO_SRI_sweep_with_phase_ramp_peak_loc_hgt.png")

end#begin
