# plot results of Sync PRI Monte Carlo Study
using JLD2, Plots, SharedArrays, StatsPlots, PyCall, Statistics, CurveFit, MAT
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

function makeBoxPlot(coeff5s,vals,titleString::String,saveFilename::String)
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

begin
  # filename1 = "sync data/Lband/syncModule_MonteCarlo_mode_2_coeff_number1_sync_osc_sweep_wSync.jld2"
  # @load filename1 peaks ideal_peak osc_coeff_sweep #resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR loc_errors 
  # ideal_peak_dB = 20 .* log10.(ideal_peak)
  # peakvals1 = ideal_peak_dB .- (20 .* log10.(peaks))
  # filename2 = "sync data/Lband/syncModule_MonteCarlo_mode_2_coeff_number2_sync_osc_sweep_wSync.jld2"
  # @load filename2 peaks ideal_peak#resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR loc_errors osc_coeff_sweep
  # peakvals2 = ideal_peak_dB .- (20 .* log10.(peaks))
  # filename3 = "sync data/Lband/syncModule_MonteCarlo_mode_2_coeff_number3_sync_osc_sweep_wSync.jld2"
  # @load filename3 peaks ideal_peak#resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR loc_errors osc_coeff_sweep
  # peakvals3 = ideal_peak_dB .- (20 .* log10.(peaks))
  # filename4 = "sync data/Lband/syncModule_MonteCarlo_mode_2_coeff_number4_sync_osc_sweep_wSync.jld2"
  # @load filename4 peaks ideal_peak #resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR loc_errors osc_coeff_sweep
  # peakvals4 = ideal_peak_dB .- (20 .* log10.(peaks))
  # filename5 = "sync data/Lband/syncModule_MonteCarlo_mode_2_coeff_number5_sync_osc_sweep_wSync.jld2"
  # @load filename5 peaks ideal_peak #resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR loc_errors osc_coeff_sweep
  # peakvals5 = ideal_peak_dB .- (20 .* log10.(peaks))
  # gr()
  # 
  # 
  # coeff_plot = osc_coeff_sweep.*ones(size(peaks,2))
  # 
  # boxplot(coeff_plot, peakvals1', title = "Peak Power Loss for Oscillator Power Law Coefficients", ylabel = "Power Loss (dB)", marker_color="red",
  #   xlabel = "PSD Coefficient Value (dB)", legend = false)
  # boxplot(coeff_plot, peakvals2', title = "Peak Power Loss for Oscillator Power Law Coefficients", ylabel = "Power Loss (dB)", marker_color="green",
  #   xlabel = "PSD Coefficient Value (dB)", legend = false)
  # boxplot(coeff_plot, peakvals3', title = "Peak Power Loss for Oscillator Power Law Coefficients", ylabel = "Power Loss (dB)", marker_color="blue",
  #   xlabel = "PSD Coefficient Value (dB)", legend = false)
  # boxplot(coeff_plot, peakvals4', title = "Peak Power Loss for Oscillator Power Law Coefficients", ylabel = "Power Loss (dB)", marker_color="orange",
  #   xlabel = "PSD Coefficient Value (dB)", legend = false)
  # boxplot(coeff_plot, peakvals5', title = "Peak Power Loss for Oscillator Power Law Coefficients", ylabel = "Power Loss (dB)", marker_color="black",
  #   xlabel = "PSD Coefficient Value (dB)", legend = false)
  # 
  # meanvals1 = mean(peakvals1,dims=2)
  # # a,b = exp_fit(osc_coeff_sweep, meanvals1) # this one fails, used matlab
  # a = 4402; b = 0.186
  # x_temp = -100:.1:-30 
  # y_temp = a.*x_temp .^b
  # plot(osc_coeff_sweep', meanvals1)
  # plot!(x_temp,y_temp)
  # 
  # meanvals2 = mean(peakvals2,dims=2)
  # meanvals3 = mean(peakvals3,dims=2)
  # meanvals4 = mean(peakvals4,dims=2)
  # meanvals5 = mean(peakvals5,dims=2)
  # 
  # plot(osc_coeff_sweep',meanvals1')
      
      
      
      
  
  # # @load filename peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors osc_coeff_sweep
  # file = matopen("syncModule_MonteCarlo_mode_2_sync_osc_sweep_.mat","w")
  # write(file, "peakvals1", collect(peakvals1))
  # write(file, "peakvals2", collect(peakvals2))
  # write(file, "peakvals3", collect(peakvals3))
  # write(file, "peakvals4", collect(peakvals4))
  # write(file, "peakvals5", collect(peakvals5))
  # write(file, "osc_coeff_sweep", collect(osc_coeff_sweep))
  # close(file)
      
      
      
      
  filename1 = "sync data/Lband/syncModule_MonteCarlo_mode_1_coeff_number1_sync_osc_sweep_noSync.jld2"
  @load filename1 peaks ideal_peak osc_coeff_sweep #resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR loc_errors 
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals1 = ideal_peak_dB .- (20 .* log10.(peaks))
  filename2 = "sync data/Lband/syncModule_MonteCarlo_mode_2_coeff_number1_sync_osc_sweep_wSync.jld2"
  @load filename2 peaks ideal_peak#resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR loc_errors osc_coeff_sweep
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals2 = ideal_peak_dB .- (20 .* log10.(peaks))
  filename3 = "sync data/Lband/syncModule_MonteCarlo_mode_3_coeff_number1_sync_osc_sweep_wSync.jld2"
  @load filename3 peaks ideal_peak#resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR loc_errors osc_coeff_sweep
  ideal_peak_dB = 20 .* log10.(ideal_peak)
  peakvals3 = ideal_peak_dB .- (20 .* log10.(peaks))
  gr()


  coeff_plot = osc_coeff_sweep.*ones(size(peaks,2))
  
  file = matopen("syncModule_MonteCarlo_coeff1_all_Modes_sync_osc_sweep_.mat","w")
  write(file, "peakvals1", collect(peakvals1))
  write(file, "peakvals2", collect(peakvals2))
  write(file, "peakvals3", collect(peakvals3))
  write(file, "osc_coeff_sweep", collect(osc_coeff_sweep))
  close(file)
      
end#begin