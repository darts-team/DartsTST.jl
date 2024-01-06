# writing a wrapper function to take the .jld2 files on epstein and resave them into .mat files so that we can read them on other machines
using JLD2, MAT, SharedArrays

list_of_osc_names = ["syncModule_MonteCarlo_mode_2_MeasuredRFSoCwGPSDO_sync_pri_sweep_along_n","syncModule_MonteCarlo_mode_2_RFSoc_sync_pri_sweep_along_n",
"syncModule_MonteCarlo_mode_2_USO_sync_pri_sweep_along_n","syncModule_MonteCarlo_mode_2_RoseL_sync_pri_sweep_along_n",
 "syncModule_MonteCarlo_mode_2_MeasuredGPSDO_sync_pri_sweep_along_n"]
 
 for i = 1 : length(list_of_osc_names)
     fname = list_of_osc_names[i]*".jld2"
    
    @load fname peaks resolutions PSLRs ISLRs ideal_res ideal_PSLR ideal_ISLR ideal_peak loc_errors sync_PRIs s_loc_1 s_loc_2 s_loc_3
    
    
    file = matopen(list_of_osc_names[i]*".mat","w")
    write(file, "peaks", collect(peaks))
    write(file, "resolutions", collect(resolutions))
    write(file, "PSLRs", collect(PSLRs))
    write(file, "ISLRs", collect(ISLRs))
    write(file, "ideal_res", collect(ideal_res))
    write(file, "ideal_PSLR", collect(ideal_PSLR))
    write(file, "ideal_peak", collect(ideal_peak))
    write(file, "ideal_ISLR", collect(ideal_ISLR))
    write(file, "loc_errors", collect(loc_errors))
    write(file, "sync_PRIs", collect(sync_PRIs))
    write(file, "s_loc_1", collect(s_loc_1))
    write(file, "s_loc_2", collect(s_loc_2))
    write(file, "s_loc_3", collect(s_loc_3))
    close(file)
     
 end