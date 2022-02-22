# plot results of Sync PRI Monte Carlo Study
using JLD2, Plots, SharedArrays, StatsPlots, PyCall, Statistics
np = pyimport("numpy")
# define a couple functions to filter NaNs and make plots
function filterNaNs2D(vals)
  #this function filters out the NaNs, but returns an array of vectors. Unfortunately, not helpful for plotting boxplots.

  mask = np.isnan(vals) # finds NaNs
  mask = [!i for i in mask] # find non-nan indices

  valsF = Array{Vector{Float64},1}(undef, size(mask,1))
  for i = 1 : size(mask,1)
    valsF[i] = vals[i, mask[i,:]]
  end#for

  return valsF
end#function

function makeBoxPlot(coeff1s,vals,titleString::String,saveFilename::String)
  len = size(vals,1)

  boxplot(osc_coeff_sweep[1].*ones(len),vals[1], title = titleString,
  xlabel = "PSD Coefficient Value (dB)", legend = false)

  for i = 2 : len - 1
    boxplot!(osc_coeff_sweep[i].*ones(len),vals[i], title = titleString,
    xlabel = "PSD Coefficient Value (dB)", legend = false)
  end#for

  display(boxplot!(osc_coeff_sweep[len].*ones(len),vals[len], title = titleString,
  xlabel = "PSD Coefficient Value (dB)", legend = false))
  plot!(ylim=(0, 20))
  savefig(saveFilename)
end#function


# run the rest of the code, calling an input .jld2 file with the monte carlo data stored
begin
  filename = "sync data/syncModule_MonteCarlo_mode_2_coeff_number1_sync_osc_sweep_wSync.jld2"
  @load filename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors osc_coeff_sweep
  gr()

  #matrices are either (numcoeff1 x numTrials) or (3 x numcoeff1 x numTrials)

  coeff_plot = osc_coeff_sweep.*ones(size(peaks,2))

## trying box plots
  # ---- Peak loss ------
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals = ideal_peak_dB .- (20 .* log10.(peaks))
  display(boxplot(coeff_plot, peakvals', title = "Peak Power Loss Oscillator Coefficient #1", ylabel = "Power Loss (dB)",
    xlabel = "PSD Coefficient Value (dB)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_peak_loss.png")
  savefig("SIMO_coeff1_sweep_peak_loss.png")

  # ---- PSLR ----
  vals1 = filterNaNs2D(PSLRs[1,:,:] .- ideal_PSLR[1])
  vals2 = filterNaNs2D(PSLRs[2,:,:] .- ideal_PSLR[2])
  vals3 = filterNaNs2D(PSLRs[3,:,:] .- ideal_PSLR[3])


  makeBoxPlot(osc_coeff_sweep,vals1,"PSLR Change Along-Track (dB)","SIMO_coeff1_sweep_no_phase_ramp_pslr_change_AT.png")
  makeBoxPlot(osc_coeff_sweep,vals2,"PSLR Change Cross-Track (dB)","SIMO_coeff1_sweep_no_phase_ramp_pslr_change_CT.png")
  makeBoxPlot(osc_coeff_sweep,vals3,"PSLR Change Height (dB)","SIMO_coeff1_sweep_no_phase_ramp_pslr_change_hgt.png")

  # makeStdevMeanPlot(osc_coeff_sweep, vals1,"PSLR Change Cross-Track (dB)", true)
  # display(boxplot(osc_coeff_sweep,vals1',title = "PSLR Change Cross-Track (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_pslr_change_CT.png")
  # display(boxplot(osc_coeff_sweep,vals2',title = "PSLR Change Along-Track (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_pslr_change_AT.png")
  # display(boxplot(osc_coeff_sweep,vals3',title = "PSLR Change Height (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_pslr_change_hgt.png")


  # ---- ISLR ------
  vals1 = filterNaNs2D(ISLRs[1,:,:] .- ideal_ISLR[1])
  vals2 = filterNaNs2D(ISLRs[2,:,:] .- ideal_ISLR[2])
  vals3 = filterNaNs2D(ISLRs[3,:,:] .- ideal_ISLR[3])
  makeBoxPlot(osc_coeff_sweep,vals1,"ISLR Change Along-Track (dB)","SIMO_coeff1_sweep_no_phase_ramp_islr_change_AT.png")
  makeBoxPlot(osc_coeff_sweep,vals2,"ISLR Change Cross-Track (dB)","SIMO_coeff1_sweep_no_phase_ramp_islr_change_CT.png")
  makeBoxPlot(osc_coeff_sweep,vals3,"ISLR Change Height (dB)","SIMO_coeff1_sweep_no_phase_ramp_islr_change_hgt.png")

  # display(boxplot(osc_coeff_sweep,vals1',title = "ISLR Change Cross-Track (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_islr_change_CT.png")
  # display(boxplot(osc_coeff_sweep,vals2',title = "ISLR Change Along-Track (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_islr_change_AT.png")
  # display(boxplot(osc_coeff_sweep,vals3',title = "ISLR Change Height (dB)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_islr_change_hgt.png")

  # ---- Resolution ------
  vals1 = filterNaNs2D(resolutions[1,:,:] .- ideal_res[1])
  vals2 = filterNaNs2D(resolutions[2,:,:] .- ideal_res[2])
  vals3 = filterNaNs2D(resolutions[3,:,:] .- ideal_res[3])

  makeBoxPlot(osc_coeff_sweep,vals1,"Resolution Change Along-Track (m)","SIMO_coeff1_sweep_no_phase_ramp_resolution_change_AT.png")
  makeBoxPlot(osc_coeff_sweep,vals2,"Resolution Change Cross-Track (m)","SIMO_coeff1_sweep_no_phase_ramp_resolution_change_CT.png")
  makeBoxPlot(osc_coeff_sweep,vals3,"Resolution Change Height (m)","SIMO_coeff1_sweep_no_phase_ramp_resolution_change_hgt.png")


  # display(boxplot(osc_coeff_sweep,vals1',title = "Resolution Change Cross-Track",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_resolution_change_CT.png")
  # display(boxplot(osc_coeff_sweep,vals2',title = "Resolution Change Along-Track (m)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_resolution_change_AT.png")
  # display(boxplot(osc_coeff_sweep,vals3',title = "Resolution Change Height (m)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_resolution_change_hgt.png")

  # ---- Peak location error ------
  vals1 = filterNaNs2D(loc_errors[1,:,:])
  vals2 = filterNaNs2D(loc_errors[2,:,:])
  vals3 = filterNaNs2D(loc_errors[3,:,:])
  makeBoxPlot(osc_coeff_sweep,vals1,"Peak Location Error Along-Track (m)","SIMO_coeff1_sweep_no_phase_ramp_peak_loc_AT.png")
  makeBoxPlot(osc_coeff_sweep,vals2,"Peak Location Error Cross-Track (m)","SIMO_coeff1_sweep_no_phase_ramp_peak_loc_CT.png")
  makeBoxPlot(osc_coeff_sweep,vals3,"Peak Location Error Height (m)","SIMO_coeff1_sweep_no_phase_ramp_peak_loc_hgt.png")

  # display(boxplot(coeff_plot, vals1', title = "Peak Location Error Cross-Track (m)",
  #  xlabel = "Sync Repetition Interval (s)", legend = false))
  #  savefig("SIMO_coeff1_sweep_no_phase_ramp_peak_loc_lat.png")
  # display(boxplot(coeff_plot, vals2', title = "Peak Location Error Along-Track (m)",
  #   xlabel = "Sync Repetition Interval (s)", legend = false))
  #   savefig("SIMO_coeff1_sweep_no_phase_ramp_peak_loc_lon.png")
  # display(boxplot(coeff_plot, vals3', title = "Peak Location Error Height (m)",
  #    xlabel = "Sync Repetition Interval (s)", legend = false))
  # savefig("SIMO_coeff1_sweep_no_phase_ramp_peak_loc_hgt.png")

end#begin
