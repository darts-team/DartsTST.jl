# plot results of Sync PRI Monte Carlo Study
using JLD2, Plots, SharedArrays
# old code with fewer vars saved
# filename = "sync data/syncModule2_MonteCarlo_mode_3_USO_sync_pri_sweep_noFreq.jld2"
# @load filename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors 
# sync_PRIs = [.1 1 2 5 10]


# newer code with updated vars saved
# filename = "syncModule2_MonteCarlo_mode_3_USO_sync_pri_sweep_wFreq.jld2"
filename = "syncModule_MonteCarlo_mode_3_USRP_sync_pri_sweep_wFreq.jld2"

@load filename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_θ s_ϕ s_h

gr()

#matrices are either (numSRI x numTrials) or (3 x numSRI x numTrials)

SRI_plot = sync_PRIs.*ones(size(peaks,2))

# # ---- Peak loss ------
# ideal_peak_dB = 10 .* log10.(ideal_peak)
# peakvals = ideal_peak_dB .- (10 .* log10.(peaks))
# display(scatter(SRI_plot, peakvals', title = "Peak Power Loss (dB)",
#   xlabel = "SRI (s)", legend = false))
# savefig("USO_SRI_sweep_no_phase_ramp_peak_loss.png")
# 
# # ---- ISLR ------
# vals1 = ideal_ISLR[1] .- ISLRs[1,:,:]
# vals2 = ideal_ISLR[2] .- ISLRs[2,:,:]
# vals3 = ideal_ISLR[3] .- ISLRs[3,:,:]
# p1 = scatter(SRI_plot, vals1', title = "ISLR Change Lat (dB)",
#  xlabel = "SRI (s)", legend = false)
# p2 = scatter(SRI_plot, vals2', title = "ISLR Change Lon (dB)",
#   xlabel = "SRI (s)", legend = false)
# p3 = scatter(SRI_plot, vals3', title = "ISLR Change Height (dB)",
#    xlabel = "SRI (s)", legend = false)
# display(plot(p1, p2, p3, layout = (3,1), legend = false))
# savefig("USO_SRI_sweep_no_phase_ramp_islr_change.png")
# 
# # ---- PSLR ------
# vals1 = ideal_PSLR[1] .- PSLRs[1,:,:]
# vals2 = ideal_PSLR[2] .- PSLRs[2,:,:]
# vals3 = ideal_PSLR[3] .- PSLRs[3,:,:]
# p1 = scatter(SRI_plot, vals1', title = "PSLR Change Lat (dB)",
#  xlabel = "SRI (s)", legend = false)
# p2 = scatter(SRI_plot, vals2', title = "PSLR Change Lon (dB)",
#   xlabel = "SRI (s)", legend = false)
# p3 = scatter(SRI_plot, vals3', title = "PSLR Change Height (dB)",
#    xlabel = "SRI (s)", legend = false)
# display(plot(p1, p2, p3, layout = (3,1), legend = false))
# savefig("USO_SRI_sweep_no_phase_ramp_pslr_change.png")
# 
# # ---- Resolution ------
# vals1 = (ideal_res[1] .- resolutions[1,:,:]) .*110e3
# vals2 = (ideal_res[2] .- resolutions[2,:,:]) .*110e3
# vals3 = (ideal_res[3] .- resolutions[3,:,:])
# p1 = scatter(SRI_plot, vals1', title = "Resolution Change Lat (deg*110km -> m)",
#  xlabel = "SRI (s)", legend = false)
# p2 = scatter(SRI_plot, vals2', title = "Resolution Change Lon (deg*110km -> m)",
#   xlabel = "SRI (s)", legend = false)
# p3 = scatter(SRI_plot, vals3', title = "Resolution Change Height (m)",
#    xlabel = "SRI (s)", legend = false)
# display(plot(p1, p2, p3, layout = (3,1), legend = false, ylims = (0,1)))
# savefig("USO_SRI_sweep_no_phase_ramp_resolution_change.png")
# 
# # ---- Peak location error ------
# vals1 = loc_errors[1,:,:] .*110e3
# vals2 = loc_errors[2,:,:] .*110e3
# vals3 = loc_errors[3,:,:] 
# p1 = scatter(SRI_plot, vals1', title = "Peak Location Error Lat (deg*110km -> m)",
#  xlabel = "SRI (s)", legend = false)
# p2 = scatter(SRI_plot, vals2', title = "Peak Location Error Lon (deg*110km -> m)",
#   xlabel = "SRI (s)", legend = false)
# p3 = scatter(SRI_plot, vals3', title = "Peak Location Error Height (m)",
#    xlabel = "SRI (s)", legend = false)
# display(plot(p1, p2, p3, layout = (3,1), legend = false, ylims = (0,1)))
# savefig("USO_SRI_sweep_no_phase_ramp_peak_loc.png")

## trying box plots
using StatsPlots

# ---- Peak loss ------
ideal_peak_dB = 10 .* log10.(ideal_peak)
peakvals = ideal_peak_dB .- (10 .* log10.(peaks))
display(boxplot(SRI_plot, peakvals', title = "Peak Power Loss (dB)",
  xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_peak_loss.png")

# ---- PSLR ----
vals1 = PSLRs[1,:,:] .- ideal_PSLR[1]
vals2 = PSLRs[2,:,:] .- ideal_PSLR[2]
vals3 = PSLRs[3,:,:] .- ideal_PSLR[3]
display(boxplot(sync_PRIs,vals1',title = "PSLR Change Cross-Track (dB)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_pslr_change_lat.png")
display(boxplot(sync_PRIs,vals2',title = "PSLR Change Along-Track (dB)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_pslr_change_lon.png")
display(boxplot(sync_PRIs,vals3',title = "PSLR Change Height (dB)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_pslr_change_hgt.png")


# ---- ISLR ------
vals1 = ISLRs[1,:,:] .- ideal_ISLR[1]
vals2 = ISLRs[2,:,:] .- ideal_ISLR[2]
vals3 = ISLRs[3,:,:] .- ideal_ISLR[3]
display(boxplot(sync_PRIs,vals1',title = "ISLR Change Cross-Track (dB)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_islr_change_lat.png")
display(boxplot(sync_PRIs,vals2',title = "ISLR Change Along-Track (dB)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_islr_change_lon.png")
display(boxplot(sync_PRIs,vals3',title = "ISLR Change Height (dB)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_islr_change_hgt.png")

# ---- Resolution ------
vals1 = (resolutions[1,:,:] .- ideal_res[1]) .*110e3
vals2 = (resolutions[2,:,:] .- ideal_res[2]) .*110e3
vals3 = (resolutions[3,:,:] .- ideal_res[3])
display(boxplot(sync_PRIs,vals1',title = "Resolution Change Cross-Track",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_resolution_change_lat.png")
display(boxplot(sync_PRIs,vals2',title = "Resolution Change Along-Track (m)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_resolution_change_lon.png")
display(boxplot(sync_PRIs,vals3',title = "Resolution Change Height (m)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_resolution_change_hgt.png")

# ---- Peak location error ------
vals1 = loc_errors[1,:,:] .*110e3
vals2 = loc_errors[2,:,:] .*110e3
vals3 = loc_errors[3,:,:] 
display(boxplot(SRI_plot, vals1', title = "Peak Location Error Cross-Track (m)",
 xlabel = "Sync Repetition Interval (s)", legend = false))
 savefig("USRP_SRI_sweep_with_phase_ramp_peak_loc_lat.png")
display(boxplot(SRI_plot, vals2', title = "Peak Location Error Along-Track (m)",
  xlabel = "Sync Repetition Interval (s)", legend = false))
  savefig("USRP_SRI_sweep_with_phase_ramp_peak_loc_lon.png")
display(boxplot(SRI_plot, vals3', title = "Peak Location Error Height (m)",
   xlabel = "Sync Repetition Interval (s)", legend = false))
savefig("USRP_SRI_sweep_with_phase_ramp_peak_loc_hgt.png")