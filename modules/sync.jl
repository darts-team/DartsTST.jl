module Sync

using Interpolations
using LinearAlgebra
using FFTW
using Plots
using Distributed
using StaticArrays
using XLSX
using StatsFuns
# using DataInterpolations
#-start-function--------------------------------------------------------------------------------------------

"""
takes the position vectors and slow time vector to produce a phase error matrix output

# Arguments
- `time_vector::N times x 1 Array`: slow time vector
- `pos::3xN Array`: positions of platforms at slow time points (ECEF)
- `parameters::Parameters`: structure of simulation configuration parameters
- `sync_PSD::Nplatform x Npulses x sync_clk_fs Array`: OPTIONAL INPUT, precalculated synchronization PSDs

"""

function get_sync_phase(time_vector::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, pos::Array{Float64,3}, osc_coeffs::Matrix{Int64}, parameters)
    # - Inputs -
    # time_vector   : a vector of time values at which the orbit positions for each platform are sampled
    # pos           : a matrix (3 x [Np x Nt] or 3 x Nt) of the orbit positions for each platform
    # parameters    : structure with list of radar and clock parameters

    # - sync algorithm steps -
    # 1) determine number of platforms from input position vector
    # 2) estimate SNR of inter-platform connections through Friis transmission formula, calculate CRLB
    # 3) determine the PSD of the phase error from system, clock, and oscillator_quality
    #    Note: PSDs will change with every sync (based on inter-platform SNR)
    # 4) generate a set (# of platforms) of random phase sequences using PSDs.
    # 5) create phase error vectors over each sync PRI. Collect and put together into one large vector
    # 6) for each platform, interpolate from the internal time base (generated with the phase error sequences) to the global time base (the time_vector input)
    # 7) create output matrix comprised of the phase errors for each platform and epoch

    # unpack variables from input structure required in sync module
    sync_pri                = parameters.sync_pri
    sync_processing_time    = parameters.sync_processing_time
    sync_signal_len         = parameters.sync_signal_len
    sync_fc                 = parameters.sync_fc
    sync_fs                 = parameters.sync_fs
    sync_fbw                = parameters.sync_fbw
    sync_fmin               = parameters.sync_fmin
    f_osc                   = parameters.sync_f_osc
    sync_clk_fs             = parameters.sync_clk_fs
    master                  = parameters.sync_master
    #osc_coeffs              = parameters.sync_osc_coeffs
    sync_processing_time    = parameters.sync_processing_time
    mode                    = parameters.mode
    no_sync_flag            = parameters.no_sync_flag
    fc                      = parameters.fc
    phase_offset_flag       = parameters.phase_offset_flag
    use_measured_psd_flag   = parameters.use_measured_psd_flag
    delay_since_sync        = parameters.delay_since_sync

    if use_measured_psd_flag
        osc_psd_meas_filename = parameters.osc_psd_meas_filename
        f_psd_meas, osc_psd_meas = read_PSD_excel_data(osc_psd_meas_filename)
    end

    # find total elapsed time over the course of the orbits
    t_elapse = time_vector[end] - time_vector[1]
    dt = time_vector[end]-time_vector[end-1] # assumes equal spacing
    time_vec =  collect(time_vector) # change to an array from a range. This is the output time vector
    push!(time_vec,time_vec[end]+dt) # adds an additional sample to the end for interpolation

    if no_sync_flag
        clk_args_N = ceil(4 * sync_clk_fs * t_elapse )
    else
        sync_prf = 1/sync_pri
        clk_args_N = ceil(4 * sync_clk_fs * sync_pri) # Number of sample points in the PSD. Needs to have the frequency resolution for the PSDs to caputure SRI
    end #if
    if delay_since_sync > 0 #enforces that minimum frequency of PSD covers the total time period of the aperture+delay_since_sync
        clk_args_N = ceil(sync_clk_fs * (t_elapse + delay_since_sync) )
    end #if
    if clk_args_N < 40e3 # enfore a minimum number of points. Needed for PSD accuracy
        clk_args_N = 40e3
    end#if
    up_convert =  sqrt(2) * fc / f_osc        # frequency up-conversion factor (scale factor from LO to RF)
    up_convert_psd = 1 #2*(fc / f_osc)^2 # frequency up-conversion factor (scale factor from LO to RF) for use on PSD #TODO reverted back to upscaling in phase domain instead of PSD because PSD domain doesn't account for frequency offset error upscaling

    # verify size of position input
    szp = size(pos)
    @assert szp[1] == 3 "POS needs to be 3 x [Np x Nt] or 3 x Nt"
    if ndims(pos) == 3
        nplat = szp[2]
        ntimes = szp[3]
    elseif ndims(pos) == 2
        nplat = 1
        ntimes = szp[2]
    end

## Determine TDMA Schedule TODO: change to match Ilgin's algorithm
# there's a sync TDMA schedule that we'll lump together into one event with length (sync_processing_time*nplat)
# Then, calculate TDMA schedule for radar pulses. The phase error state of each platform is calculated at each of these time points.
# Following, the total phase error (output of sync module) is dependent on the radar mode (SIMO/MIMO, etc.)
# Ex. total phase error = phase error of Tx platform at Tx time + phase error of Rx.

# Create sync schedule (in slow time domain)
t_sync = collect(time_vector[1]:sync_pri:time_vector[end]+sync_pri) # create sync PRI time base (the schedule of syncs)

# create the radar TDMA schedule
if mode == 1 # SAR
    slot_len    = dt/nplat # split the PRI evenly between each tx //TODO: replace with Ilgin's schedule
    tdma_radar  = collect(time_vector[1]:slot_len:time_vector[end]+dt) # creates schedule with radar time slots
    tx_map      = round.(mod.(tdma_radar./slot_len,nplat)).+1 # map of which platform is transmitting at each TDMA time point
    # phase_state keeps track of phase error for each individual platform at each transmit time
    phase_state = Array{Float64}(undef, nplat, length(tdma_radar)) # Nplatforms x N TDMA times

elseif mode == 2 # SIMO
    # only the one transmitter. So the schedule is given by the input time vector
    tdma_radar = time_vec[1:end-1]
    tx_map = master.*ones(size(tdma_radar)) #only single transmitter
    # phase_state keeps track of phase error for each individual platform at each transmit time
    phase_state = Array{Float64}(undef, nplat, ntimes) # Nplatforms x Ntimes

elseif mode == 3 # MIMO, # code for this mode could get more complex if we decide to have only N transmitters out of M platforms
    slot_len    = dt/nplat # split the PRI evenly between each tx
    tdma_radar  = collect(time_vector[1]:slot_len:time_vector[end]+dt) # creates schedule with radar time slots
    tx_map      = round.(mod.(tdma_radar./slot_len,nplat)).+1 # map of which platform is transmitting at each time
    phase_state = Array{Float64}(undef, nplat, length(tdma_radar)) # Nplatforms x N TDMA times (collects phase error at each platform at each TDMA slot)

else
    error("Undefined radar mode")
end

## find the phase error

    # find CRLB based on SNR of sync connection at each time point. Will interpolate to SRI time vector later (This assumes that SNR doesn't vary rapidly as compared to the orbit sampling)
    (sig_crlb, snr_dB) = getSensorCRLB_network(pos, sync_signal_len, sync_fc, sync_fs, sync_fbw)
    #(sig_crlb, snr_dB) = getSensorCRLB_network(pos, sync_signal_len, sync_fc, sync_fs, sync_fbw,parameters)
    #display(plot(time_vector,snr_dB',xlabel="Slow time (s)",ylabel="SNR (dB)",title="Link to Master Platform SNR"))
    #display(plot(time_vector,crlb,xlabel="time (s)",title="CRLB"))

    # create matrix of phase error values
    if mode == 1 || mode == 2
         phase_err   = Array{Float64}(undef, szp[2], szp[3]) # Nplatforms x Ntimes
    elseif mode == 3
         phase_err   = Array{Float64}(undef, szp[2], szp[2], szp[3]) # Nplatforms x Nplatforms x Ntimes
    end#if

     # save the PSD outputs so they can be reused. Saves time when multiple trials are used and is more important at longer SRI
     #size is a mode dependent
     if mode == 1 || mode == 2
         pulse_idx = findall(tx_map->tx_map == 1,tx_map) #fixed to first transmitters number of pulses. will be greater than any others by default, letting array be either equal to or greater than needed size
     else
         pulse_idx = collect(1:length(tx_map))
     end#if
     platform_pulse_times = tdma_radar[pulse_idx]
     sync_PSDs   = Array{Float64}(undef, nplat, length(platform_pulse_times), convert(Int64,clk_args_N)) # Nplatforms x Ntimes x (length of PSD = sync_clk_fs)

    for i = 1:nplat #for each platform, generate oscillator phase errors at each time point

        a_coeff_db = osc_coeffs[i,:]   # grab the clock coefficients for each platform
        crlbs      = sig_crlb[i,:]      # grab the CRLB previously calculated
        #interpolate CRLB values to sync PRI times
        push!(crlbs,crlbs[end])
        itp_crlb = LinearInterpolation(time_vec, crlbs, extrapolation_bc = Line()) # create interpolant with CRLB values sampled at orbit sample times. Repeats last value to avoid extrapolation errors


        if use_measured_psd_flag #use measured PSD values
            (Sphi, f_psd) = osc_psd_measured(sync_clk_fs,clk_args_N,f_psd_meas, osc_psd_meas)

        else # use PSD from coefficient values
            (Sphi, f_psd) = osc_psd_twosided(sync_clk_fs, clk_args_N, a_coeff_db, sync_fmin)    # this gives the basic PSD of the oscillator (no sync involved)
        end

        phase_err_internal = zeros(0)       # initialize an empty array to collect the phase error values at internal time base
        internal_time_vec = zeros(0)        # initialize an empty array to keep track of internal time base

        if mode == 2 # SIMO
            pulse_idx = collect(1:length(tx_map))
        else
            pulse_idx = findall(tx_map->tx_map == i,tx_map) # gives all radar pulse times for given transmitter
        end#if SIMO mode
        platform_pulse_times = tdma_radar[pulse_idx] # this gives the times for each pulse transmitted by given platform

        ## this approach loops over the sync interval times, generates a phase error PSD as each time point
        # 1) get CRLB value, 2) find post-sync PSD...?
        for j = 1 : length(platform_pulse_times) # loops over each pulse time

            if no_sync_flag # no sync case
                Sphi_uc = up_convert_psd .* Sphi
                (r, t)    = osc_timeseries_from_psd_twosided(Sphi_uc, sync_clk_fs, phase_offset_flag,delay_since_sync)
                #store PSD
                sync_PSDs[i,j,:] = Sphi
            else
                pulse_time = platform_pulse_times[j] # time of current pulse
                sync_idx = findlast(t_sync-> t_sync <= pulse_time, t_sync) # find most recent sync time index
                sync_time = t_sync[sync_idx] # most recent sync time
                # sync_radar_offset = pulse_time - sync_time # time difference between the pulse and most recent sync
                sync_radar_offset = pulse_time - sync_time + delay_since_sync # 3/29/22 adding arbitrary delay since sync delay here

                if sync_radar_offset < (sync_processing_time*nplat) # enforce minimum sync processing time for offset
                    sync_radar_offset = (sync_processing_time*nplat)
                end

                # calculate post-sync PSD from crlb, sync_radar_offset
                crlb = itp_crlb(sync_time) # this tells us what the crlb was at the sync
                Sphi_sync = sync_effects_on_PSD(Sphi, f_psd, sync_radar_offset, crlb, sync_prf, sync_fs, sync_clk_fs, dt, f_osc)
                # generate time series of phase error
                Sphi_sync_uc = up_convert_psd .* Sphi
                (r, t)    = osc_timeseries_from_psd_twosided(Sphi_sync_uc, sync_clk_fs,phase_offset_flag, delay_since_sync)

                #store PSD
                sync_PSDs[i,j,:] = Sphi_sync
            end # if no_sync_flag



            # t and r are longer than PRI, cut to PRI
            idx_pri = findfirst(t -> t >= dt,t)
            r = r[1:idx_pri]
            t = t[1:idx_pri]

            # append phase errors to internal matrices. Will interpolate to all radar pulse times after loop

            if length(internal_time_vec)>0
                append!(internal_time_vec, internal_time_vec[end].+ t) # adding previous ending time to continue with time base (not reset from 0)
                append!(phase_err_internal,r .+ phase_err_internal[end]) # makes phase error continuous
            else
                append!(internal_time_vec, t .+ time_vector[1]) # adding from starting time from input time vector
                append!(phase_err_internal,r) # makes phase error continuous
            end#if start of internal time vec

        end # for platform pulse times

        #= add in linear phase ramp term (uses frequency offset of oscillator)
        sigma_freq_offset = sigma_freq_offsets[i]
        freq_offset = randn(1)*sigma_freq_offset # zero mean random
        phase_ramp = 2*pi* freq_offset .* internal_time_vec
        phase_err_internal = phase_err_internal + phase_ramp
        =#

        # reset phase error to 0 at each sync time. Do this by finding the start/stop indices at each of the sync times
        if no_sync_flag
            #no reset to 0 error
        else
            for nsync = 1 : length(t_sync)-1
                sync_time       = t_sync[nsync]
                next_sync_time  = t_sync[nsync+1]
                start_idx       = findfirst(internal_time_vec -> internal_time_vec >= sync_time, internal_time_vec)
                stop_idx        = findfirst(internal_time_vec -> internal_time_vec >= next_sync_time, internal_time_vec)
                # grab phase values during this interval
                if start_idx != nothing
                    if stop_idx == nothing
                        stop_idx        = length(internal_time_vec)
                        phase_vals      = phase_err_internal[start_idx : stop_idx]
                        phase_err_internal[start_idx:stop_idx] = phase_vals .- phase_vals[1]
                    else
                        phase_vals      = phase_err_internal[start_idx : stop_idx - 1]
                        phase_err_internal[start_idx:stop_idx-1] = phase_vals .- phase_vals[1]
                    end #if stop_idx
                else
                    # println("Start idx = nothing, nsync = ", nsync)
                end#if start_idx
            end # forloop
        end# no sync_flag

        # Now, create interpolant and interpolate to the TDMA schedule
        itp_tdma = LinearInterpolation(internal_time_vec, phase_err_internal,extrapolation_bc = Flat()) # create interpolant, hold same value at end if out of bounds
        phase_state[i,:] = itp_tdma(tdma_radar);
        # if i==1
        #     # display(plot(tdma_radar,phase_state[i,:].*180/pi,xlabel="time (s)",ylabel="Oscillator phase (deg)",title="Phase Error at TDMA Sampling"))
        #     display(plot(tdma_radar,phase_state[i,:].*180/pi.*up_convert,xlabel="time (s)",ylabel="RF phase (deg)",title="Phase Error at TDMA Sampling"))
        # else
        #     # display(plot!(tdma_radar,phase_state[i,:].*180/pi,xlabel="time (s)",ylabel="Oscillator phase (deg)",title="Phase Error at TDMA Sampling"))
        #     display(plot!(tdma_radar,phase_state[i,:].*180/pi.*up_convert,xlabel="time (s)",ylabel="RF phase (deg)",title="Phase Error at TDMA Sampling"))
        # end

    end # loop over platforms

    if mode == 1 # ping-pong
            for m = 1 : nplat
                pulse_idx = findall(tx_map->tx_map == m,tx_map) # gives all radar pulse times for given transmitter
                phase_err[m,:] = up_convert.* phase_state[m,pulse_idx[1:ntimes]]
            end#for nplat

    elseif mode == 2 #SIMO
        phase_err = up_convert.* phase_state

    elseif mode == 3 # MIMO
        for m = 1 : nplat
            pulse_idx = findall(tx_map->tx_map == m,tx_map) # gives all radar pulse times for given transmitter
            for n = 1 : nplat
                phase_err[m,n,:] = up_convert.* phase_state[n,pulse_idx[1:ntimes]] # gives phase error state for rx platformat tx times
            end#for nplat
        end#for nplat
    end#if

    return phase_err, sync_PSDs

end #get_sync_phase
# ------------------------------------------------------------------------------------------------------
# overload function with case where Sync PSDs are given
function get_sync_phase(time_vector::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, pos::Array{Float64,3}, parameters, sync_PSDs::Array{Float64,3})


    # unpack variables from input structure required in sync module
    sync_pri                = parameters.sync_pri
    sync_processing_time    = parameters.sync_processing_time
    sync_signal_len         = parameters.sync_signal_len
    sync_fc                 = parameters.sync_fc
    sync_fs                 = parameters.sync_fs
    sync_fbw                = parameters.sync_fbw
    sync_fmin               = parameters.sync_fmin
    f_osc                   = parameters.sync_f_osc
    sync_clk_fs             = parameters.sync_clk_fs
    master                  = parameters.sync_master
    # osc_coeffs              = parameters.osc_coeffs
    sync_processing_time    = parameters.sync_processing_time
    mode                    = parameters.mode
    no_sync_flag            = parameters.no_sync_flag
    fc                      = parameters.fc
    phase_offset_flag       = parameters.phase_offset_flag
    # sigma_freq_offsets      = parameters.sigma_freq_offsets
    delay_since_sync        = parameters.delay_since_sync


    # find total elapsed time over the course of the orbits
    t_elapse = time_vector[end] - time_vector[1]
    dt = time_vector[end]-time_vector[end-1] # assumes equal spacing
    time_vec =  collect(time_vector) # change to an array from a range. This is the output time vector
    push!(time_vec,time_vec[end]+dt) # adds an additional sample to the end for interpolation

    if no_sync_flag
        clk_args_N = ceil(4 * sync_clk_fs * t_elapse )
    else
        sync_prf = 1/sync_pri
        clk_args_N = ceil(4 * sync_clk_fs * sync_pri) # Number of sample points in the PSD. Needs to have the frequency resolution for the PSDs to caputure SRI

    end #if
    if delay_since_sync > 0 #enforces that minimum frequency of PSD covers the total time period of the aperture+delay_since_sync
        clk_args_N = ceil(sync_clk_fs * (t_elapse + delay_since_sync) )
    end #if
        
    if clk_args_N < 80e3 # enfore a minimum number of points. Needed for PSD accuracy
        clk_args_N = 80e3
    end#if
    up_convert =  fc / f_osc        # frequency up-conversion factor (scale factor from LO to RF)

    # verify size of position input
    szp = size(pos)
    @assert szp[1] == 3 "POS needs to be 3 x [Np x Nt] or 3 x Nt"
    if ndims(pos) == 3
        nplat = szp[2]
        ntimes = szp[3]
    elseif ndims(pos) == 2
        nplat = 1
        ntimes = szp[2]
    end

## Determine TDMA Schedule TODO: change to match Ilgin's algorithm
# there's a sync TDMA schedule that we'll lump together into one event with length (sync_processing_time*nplat)
# Then, calculate TDMA schedule for radar pulses. The phase error state of each platform is calculated at each of these time points.
# Following, the total phase error (output of sync module) is dependent on the radar mode (SIMO/MIMO, etc.)
# Ex. total phase error = phase error of Tx platform at Tx time + phase error of Rx.

# Create sync schedule (in slow time domain)
t_sync = collect(time_vector[1]:sync_pri:time_vector[end]+sync_pri) # create sync PRI time base (the schedule of syncs)

# create the radar TDMA schedule
if mode == 1 # ping-pong
    slot_len    = dt/nplat # split the PRI evenly between each tx //TODO: replace with Ilgin's schedule
    tdma_radar  = collect(time_vector[1]:slot_len:time_vector[end]+dt) # creates schedule with radar time slots
    tx_map      = round.(mod.(tdma_radar./slot_len,nplat)).+1 # map of which platform is transmitting at each TDMA time point
    # phase_state keeps track of phase error for each individual platform at each transmit time
    phase_state = Array{Float64}(undef, nplat, length(tdma_radar)) # Nplatforms x N TDMA times

elseif mode == 2 #SIMO
    # only the one transmitter. So the schedule is given by the input time vector
    tdma_radar = time_vec[1:end-1]
    tx_map = master.*ones(size(tdma_radar)) #only single transmitter
    # phase_state keeps track of phase error for each individual platform at each transmit time
    phase_state = Array{Float64}(undef, nplat, ntimes) # Nplatforms x Ntimes

elseif mode == 3 # MIMO, # code for this mode could get more complex if we decide to have only N transmitters out of M platforms
    slot_len    = dt/nplat # split the PRI evenly between each tx //TODO: replace with Ilgin's schedule
    tdma_radar  = collect(time_vector[1]:slot_len:time_vector[end]+dt) # creates schedule with radar time slots
    tx_map      = round.(mod.(tdma_radar./slot_len,nplat)).+1 # map of which platform is transmitting at each time
    phase_state = Array{Float64}(undef, nplat, length(tdma_radar)) # Nplatforms x N TDMA times (collects phase error at each platform at each TDMA slot)

else
    error("Undefined radar mode")
end

## find the phase error

    # find CRLB based on SNR of sync connection at each time point. Will interpolate to SRI time vector later (This assumes that SNR doesn't vary rapidly as compared to the orbit sampling)
    (sig_crlb, snr_dB) = getSensorCRLB_network(pos, sync_signal_len, sync_fc, sync_fs, sync_fbw)
    #(sig_crlb, snr_dB) = getSensorCRLB_network(pos, sync_signal_len, sync_fc, sync_fs, sync_fbw,parameters)
    #display(plot(time_vector,snr_dB',xlabel="Slow time (s)",ylabel="SNR (dB)",title="Link to Master Platform SNR"))
    #display(plot(time_vector,crlb,xlabel="time (s)",title="CRLB"))

    # create matrix of phase error values
    if mode == 1 || mode == 2
         phase_err   = Array{Float64}(undef, szp[2], szp[3]) # Nplatforms x Ntimes
    elseif mode == 3
         phase_err   = Array{Float64}(undef, szp[2], szp[2], szp[3]) # Nplatforms x Nplatforms x Ntimes
    end#if

     # save the PSD outputs so they can be reused. Saves time when multiple trials are used and is more important at longer SRI
     #size is a mode dependent
     if mode == 1 || mode == 2
         pulse_idx = findall(tx_map->tx_map == 1,tx_map) #fixed to first transmitters number of pulses. will be greater than any others by default, letting array be either equal to or greater than needed size
     else
         pulse_idx = collect(1:length(tx_map))
     end#if
     platform_pulse_times = tdma_radar[pulse_idx]

    for i = 1:nplat #for each platform, generate oscillator phase errors at each time point

        phase_err_internal = zeros(0)       # initialize an empty array to collect the phase error values at internal time base
        internal_time_vec = zeros(0)        # initialize an empty array to keep track of internal time base

        if mode == 2 # SIMO
            pulse_idx = collect(1:length(tx_map))
        else
            pulse_idx = findall(tx_map->tx_map == i,tx_map) # gives all radar pulse times for given transmitter
        end#if SIMO mode
        platform_pulse_times = tdma_radar[pulse_idx] # this gives the times for each pulse transmitted by given platform

        ## this approach loops over the sync interval times, generates a phase error PSD as each time point
        # 1) get CRLB value, 2) find post-sync PSD
        for j = 1 : length(platform_pulse_times) # loops over each pulse time

            Sphi_sync = sync_PSDs[i,j,:] # load precalculated PSD
            # generate time series of phase error
            Sphi_sync_uc = up_convert_psd .* Sphi
            (r, t)    = osc_timeseries_from_psd_twosided(Sphi_sync_uc, sync_clk_fs,phase_offset_flag, delay_since_sync)

            # t and r are longer than PRI, cut to PRI
            idx_pri = findfirst(t -> t >= dt,t)
            r = r[1:idx_pri]
            t = t[1:idx_pri]

            # append phase errors to internal matrices. Will interpolate to all radar pulse times after loop

            if length(internal_time_vec)>0
                append!(internal_time_vec, internal_time_vec[end].+ t) # adding previous ending time to continue with time base (not reset from 0)
                append!(phase_err_internal,r .+ phase_err_internal[end]) # makes phase error continuous
            else
                append!(internal_time_vec, t .+ time_vector[1]) # adding from starting time from input time vector
                append!(phase_err_internal,r) # makes phase error continuous
            end#if start of internal time vec

        end # for platform pulse times

        # add in linear phase ramp term (uses frequency offset of oscillator)
        # sigma_freq_offset = sigma_freq_offsets[i]
        # freq_offset = randn(1)*sigma_freq_offset # zero mean random
        # phase_ramp = 2*pi* freq_offset .* internal_time_vec
        # phase_err_internal = phase_err_internal + phase_ramp

        # reset phase error to 0 at each sync time. Do this by finding the start/stop indices at each of the sync times
        if no_sync_flag

        else
            for nsync = 1 : length(t_sync)-1
                sync_time       = t_sync[nsync]
                next_sync_time  = t_sync[nsync+1]
                start_idx       = findfirst(internal_time_vec -> internal_time_vec >= sync_time, internal_time_vec)
                stop_idx        = findfirst(internal_time_vec -> internal_time_vec >= next_sync_time, internal_time_vec)
                # grab phase values during this interval
                if start_idx != nothing
                    if stop_idx == nothing
                        stop_idx        = length(internal_time_vec)
                        phase_vals      = phase_err_internal[start_idx : stop_idx]
                        phase_err_internal[start_idx:stop_idx] = phase_vals .- phase_vals[1]
                    else
                        phase_vals      = phase_err_internal[start_idx : stop_idx - 1]
                        phase_err_internal[start_idx:stop_idx-1] = phase_vals .- phase_vals[1]
                    end #if stop_idx
                else

                end#if start_idx
            end # forloop
        end# no sync_flag

        # Now, create interpolant and interpolate to the TDMA schedule
        itp_tdma = LinearInterpolation(internal_time_vec, phase_err_internal,extrapolation_bc = Flat()) # create interpolant, hold same value at end if out of bounds
        phase_state[i,:] = itp_tdma(tdma_radar);
         # display(plot(tdma_radar,phase_state[i,:].*180/pi,xlabel="time (s)",ylabel="phase (deg)",title="Phase Error at TDMA Sampling"))

    end # loop over platforms

    if mode == 1 # SAR
            for m = 1 : nplat
                pulse_idx = findall(tx_map->tx_map == m,tx_map) # gives all radar pulse times for given transmitter
                phase_err[m,:] = up_convert.* phase_state[m,pulse_idx[1:ntimes]]
            end#for nplat

    elseif mode == 2 #SIMO
        phase_err = up_convert.* phase_state

    elseif mode == 3 # MIMO
        for m = 1 : nplat
            pulse_idx = findall(tx_map->tx_map == m,tx_map) # gives all radar pulse times for given transmitter
            for n = 1 : nplat
                phase_err[m,n,:] = up_convert.* phase_state[n,pulse_idx[1:ntimes]] # gives phase error state for rx platformat tx times
            end#for nplat
        end#for nplat
    end#if

    return phase_err, sync_PSDs

end #get_sync_phase
# --------end of phase error function---------------------------------------------------------------------------------------


#-------------Helper functions-----------------
"""
calculate the two-sided PSD of the clock phase error

# Arguments
- `fs::Integer`: max PSD frequency [sample rate of clock phase error process]
- `N::Integer`: number of PSD sample points
- `a_coeff_db::1 x 5 Array`: coefficients of the noise characteristic asymptotes

"""
#-start-function--------------------------------------------------------------------------------------------
function osc_psd_twosided(fs::Float64,N::Float64,a_coeff_db::Array{Int64,1}, sync_fmin::Float64)
#   Generate PSD of clock phase error
#INPUTS
#     fs = 2000; # max PSD frequency [sample rate of clock phase error process]
#     N = 600000; # number of PSD sample points
#     a_coeff_db = [-28 -40 -200 -130 -155]; # coefficients of the noise characteristic asymptotes.
#     fmin = .1; # minimum frequency .> 0 in Hz to window PSD [optional, default = .01 Hz].
#
# OUTPUTS
#     Sphi_2S: # two sided oscillator phase error PSD
#     f:    # frequency vector

    f = collect( ( (-N/2):1:(N/2-1) ) ).*fs./N

    temp = 0.1 .* a_coeff_db
    a_coeff = [ 10^temp[1], 10^temp[2], 10^temp[3], 10^temp[4], 10^temp[5] ]

    fmin = sync_fmin # minimum PSD frequency. Reduces PSD below this value


    Sphi_2S = a_coeff' * [f.^(-4), f.^(-3), f.^(-2), f.^(-1), f.^0]

    # find zero frequency value, set power to 0
    idx0 = findall(f -> f == 0,f)
    Sphi_2S[idx0] .= 0

    # Previous code 7/15/21 This code set frequency values below fmin to the PSD amplitude at fmin.
    # idx = findall(f -> f > 0,f)
    # fmin1 = f[idx[1]] # takes first index of frequency vector that is greater than 0, grabs frequency value
    #
    # if fmin1<fmin
    #     idx = findall(f -> abs(f) < fmin,f)
    #       #Sphi_2S[idx] .=  a_coeff[1].*fmin^(-4).+a_coeff[2]*fmin^(-3).+a_coeff[3]*fmin^(-2).+a_coeff[4]*fmin^(-1).+a_coeff[5]*fmin^0
    #     temp = a_coeff' * [fmin.^(-4), fmin.^(-3), fmin.^(-2), fmin.^(-1), fmin.^(0)]
    #     Sphi_2S[idx] .= temp
    # else
    #     # conserve power around dc
    #     idx = findall(f -> abs(f) < fmin,f)
    #     Sphi_2S[idx] .= sum([a_coeff[1]*fmin.^(-4) a_coeff[2]*fmin.^(-3) a_coeff[3]*fmin.^(-2) a_coeff[4]*fmin.^(-1) a_coeff[5]*fmin.^0 ]).*(fmin/fmin1)
    # end#if

    # code added 7/15/21, This creates a "high-pass filter" with fmin as the cutoff frequency
    # create a roll-off from fmin to 0 frequency
    temp = a_coeff' * [fmin.^(-4), fmin.^(-3), fmin.^(-2), fmin.^(-1), fmin.^(0)] # this is the PSD value at fmin
    idx_fmin = findlast(f -> abs(f) < fmin,f)
    if !isnothing(idx_fmin)
        idx = idx0[1]:idx_fmin
        # Sphi_2S[idx] .= range(0,stop=temp,length=length(idx)) # linear function from 0 freq to fmin
        
        # linear roll-off filter changed to logistic function #10/18/22
        # logistic function: f(x) = L/ [1 + e^-k(x-x0)]. L = final amplitude, k = growth rate, x0 = 0.5L crossing
        k = 10 # selected to have a less steep slope, but not linear (k=1 shows linear curve)
        # f[idx] are the x values, need to shift by fmin/2
        temp_logistic = temp .* logistic.( k.*(f[idx].-(fmin/2) ) )
        Sphi_2S[idx] = temp_logistic
    end

    # take magnitude to remove any negative PSD values.
    Sphi_2S = abs.(Sphi_2S)
    return Sphi_2S, f

end #function
#--end--function-------------------------------------------------------------------------------------------
"""
Use measured one-sided PSD of the clock phase error and convert to a 2-sided PSD for later use in sync module. Interpolates to required frequency vector

# Arguments
- `fs::Integer`: max PSD frequency [sample rate of clock phase error process]
- `N::Integer`: number of PSD sample points
- `f_psd_meas::1 x N_psd Array`: frequency vector of measured PSD
- `osc_psd_meas::1 x N_psd Array`: measured PSD amplitudes (1-sided)

"""
#-start-function--------------------------------------------------------------------------------------------
function osc_psd_measured(fs::Float64, N::Float64, f_psd_meas::Vector{Float64}, osc_psd_meas::Vector{Float64})
    #mirror osc_psd to make two-sided
    osc_psd = collect(hcat(-reverse(osc_psd_meas)', osc_psd_meas')./2) #halve magnitude to go from 1-sided to 2-sided
    f_psd = collect(hcat(-reverse(f_psd_meas)', f_psd_meas'))

    # create frequency vector to interpolate to
    f = collect( ( (-N/2):1:(N/2-1) ) ).*fs./N

    freq_itp = LinearInterpolation(f_psd[1,:], osc_psd[1,:], extrapolation_bc=Line());

    Sphi_2S = 10 .^ (freq_itp.(f)./10)

    # find zero frequency value, set power to 0
    idx0 = findall(f -> f == 0,f)
    Sphi_2S[idx0] .= 0

    # take magnitude to remove any negative PSD values.
    Sphi_2S = abs.(Sphi_2S)
    return Sphi_2S, f

end #function
#--end--function-------------------------------------------------------------------------------------------

#-start-function--------------------------------------------------------------------------------------------
"""
generates realizations of phase error given a two-sided oscillator PSD

# Arguments
- `Sphi::Nx1 Array`: 2 sided oscillator phase error PSD
- `fs::Integer`: max PSD frequency [sample rate of clock phase error process]
- `a_coeff_db:: 1x5 Array`: coefficients of the noise characteristic asymptotes
- `no_sync_flag:: Bool`: Flag if sync is not used: Phase error resets to 0 at sync event if flag = 0
- `delay_since_sync:: Bool`:

"""
function osc_timeseries_from_psd_twosided(Sphi::Array{Float64,1},fs::Float64,phase_offset_flag::Bool,delay_since_sync::Float64)
#   Generate time series from two sided PSD of clock phase error.
#   Will first calculate the one sided PSD which = 0 for f<0
#INPUTS
#     Sphi: # 2 sided oscillator phase error PSD
#     fs = 2000; # max PSD frequency [sample rate of clock phase error process]
# OUTPUTS
#     r:    # random time sequence realization
#     t = (0:(N-1))*1/fs; # time vector

N = length(Sphi)
f = collect( ( (-N/2) : 1 : (N/2-1) ) ) .* fs ./N

Sphi_1S = zeros(N)
idx = findall(f -> f >= 0,f)
Sphi_1S[idx] = 2 .* Sphi[idx] # one sided with f=0 included
A_Sphi = sqrt.(abs.(Sphi_1S.*N))

# random phase vector with amplitude of PSD
phi_u = rand(N) .*2 .* pi
amp_u = 1; # constant amplitude

# sequence with random phase & random amplitude of PSD
z = A_Sphi.*exp.(1im*phi_u).*amp_u

z2 = ifftshift(z)

# random time sequence realization
r = real(fftshift(ifft(z2))) # take real component of time series

# time vector
t = collect(0:(N-1)) .* 1/fs

#introducing the time delay since last synchronization
if delay_since_sync > 0 # generate extended phase value for time before aperture start. This will be "t=0" to "t=delay_since_sync". Aperture starts at "delay_since_sync"
    seq_len = N/fs # length of the random phase instance that is generated from PSD
    numSeqNeeded = ceil(delay_since_sync / seq_len) + 1
    time_vec = collect(0:(N*numSeqNeeded)) .* 1/fs
    phase_hist = Array{Float64}(undef,length(time_vec))
    for i = 1 : numSeqNeeded
        phi_u = rand(N) .*2 .* pi
        z = A_Sphi.*exp.(1im*phi_u).*amp_u
        z2 = ifftshift(z)
        r_int = real(fftshift(ifft(z2)))# phase values internal to loop

        if i>1 # make phase continuous
            r_int = r_int .+ phase_hist[floor(Int,(i-1)*N)]
        else
            r_int = r_int.-r_int[1]  # makes sync event at t=0 in this frame
        end
        phase_hist[ floor(Int,(i-1)*N+1) : floor(Int,N*i) ] = r_int



    end

    # plot(time_vec, phase_hist.*180 ./pi,xaxis=("time (sec)"),ylabel=("Oscillator Phase (deg) @10MHz"))
    # savefig("phase_history_wDelay")
    # now we will interpolate the phases to the pulse times of the aperture
    t_int = t .+ delay_since_sync # shift the aperture time vector "t" to the time vector that now includes the delay_since_sync
    itp_phase_hist = LinearInterpolation(time_vec, phase_hist) # create interpolant
    r = itp_phase_hist(t_int);


end



if phase_offset_flag && delay_since_sync == 0 # only considered if delay_since_sync is "0" meaning it is not being used. If delay_since_sync < 0, then phase values will be set to 0 at that sync time
    r = r.-r[1] # this zeros out the error at t=0. Sync "resets" the clock drift to 0. Will be added in sequence to keep phase errors continuous
end



return r,t
end#function
#--end--function-------------------------------------------------------------------------------------------


#-start-function-------------------------------------------------------------------------------------------
"""
generates a matrix of estimated CRLB values at each position based on relative positions of platforms and Friis transmission formula with respect to the master platform

# Arguments
- `pos::3 x N_plat x N_time Array`: vector of platform positions
- `N::Integer`: number of samples in waveform (Period (s) = N/fs)
- `fc::Integer`: center frequency of transmitted sync waveform
- `fs::Integer`: sampling frequency of sync receiver
- `fbw::Integer`: sync waveform bandwidth
- `master::Integer`: master platform number

"""
function getSensorCRLB(pos::Array{Float64,3},N::Float64,fc::Float64,fs::Float64,fbw::Float64,master::Int64)
#function getSensorCRLB(pos,N,fc,fs,fbw,master, Parameters)
# pos:              gives locations of all platforms in network (3 x N_plat x N_time)
# nplat:            number of platforms
# ntimes:           number of sample times in time_vec
# N:                number of samples in waveform (Period (s) = N/fs)
# fs:               sampling frequency
# fbw:              sync waveform bandwidth
# master:           select which Tx is the "master" for the sync

# TODO Pass in variable structure with sync power budget information

    szp = size(pos)
    if ndims(pos) == 3
        nplat = szp[2]
        ntimes = szp[3]
    elseif ndims(pos) == 2
        nplat = 1
        ntimes = szp[2]
    end

    k = 1.380649e-23; # Joule/Kelvin -  Boltzmann constant
    T = 25 + 273.15; # Kelvin -  let approximate temperature be 25C
    Np = k*T*fbw; # (W) noise power in receiver - NOTE: this doesn't consider the noise figure in receiver electronics
    c = 299792458; # m/s SOL

    # initialize SNR and CRLB matrices
    snr_dB = Array{Float64,2}(undef, nplat, ntimes)
    sig_crlb = Array{Float64,2}(undef, nplat, ntimes)

    for i = 1 : nplat # number of platforms
      if i == master
        sig_crlb[i,:] = zeros(ntimes) # ignore oscillator effects if only single platform
        snr_dB[i,:] = zeros(ntimes) # setting SNR to zero just to satisfy output
      else

      for j = 1 : ntimes # loop over number of sample times
      master_pos = pos[:,master,j]

      # findall R, distance between the platform and master
      receiver_pos = pos[:,i,:]

      R = sqrt((receiver_pos[1]-master_pos[1])^2+(receiver_pos[2]-master_pos[2])^2+(receiver_pos[3]-master_pos[3])^2)


      # need Pt, Gt, Gr, POL. These could be calculated in the orbits module using antenna pattern/platform rotations/
      Pt = 1e-3     # let transmit power be 1 mW (needs to be evaluated in trade study, more likely a consequence of orbit design)
      Gt = 1        # Isotropic rx/tx antennas for sync
      Gr = 1        # Isotropic rx/tx antennas for sync
      POL = 1       # Polarization mismatch factor (let antennas be aligned here) - 0 is total misalignment, 1 is total alignment


      # use Friis formula to calculate received power, SNR
      Pr = POL * Pt*Gt*Gr*c^2/(4*pi*fc*R)^2;
      snr_linear = Pr/Np;
      snr_dB[i,j] = 10*log10(snr_linear); # save SNR for debug/trade study

      sig_crlb[i,j] = 1/(sqrt(2*nplat)*(2*pi*fbw)^2*(snr_linear)*(N/fs)*fs) # clk error after sync for nplat sensors
      end # for ntimes
      end #if
    end #for nplat

    return (sig_crlb, snr_dB) # return tuple of CRLB and SNR values
end#function
#--end--function-------------------------------------------------------------------------------------------

#-start-function-------------------------------------------------------------------------------------------
"""
rectangle function for filter

# Arguments
- `t::Float`: vector of platform positions
- `t1::Float`: start time of rect
- `t2::Float`: end time of rect
- `height::Float`: amplitude of rect

"""
function rect(t,t1,t2,height) # does Julia have a built in rect function?
# Unit energy rect pulse either from -T/2 to T/2 with height 1/sqrt(T)
# | from t1 to t2 with height 1/sqrt(t2-t1)
# example: r=rect(t,tmin,tmax,1) # rect from tmin to tmax with height 1

    t1 = max(t[1],t1)
    t2 = min(t[end],t2)
    T = abs(t2-t1)

    delta = (t[end]-t[1])/(length(t))
    ind1 = floor(Int, (t1-t[1])/delta)+1
    ind2 = floor(Int, (t2-t[1])/delta)
    r = zeros(length(t))
    r[ind1:ind2] .= height

    return r
end
#--end--function-------------------------------------------------------------------------------------------

##-start-function-------------------------------------------------------------------------------------------
"""
calculates the post-synchronization phase error PSD

# Arguments
- `Sphi::Nsamplesx1 Array`: unsync'd phase error PSD
- `f_psd::Nsamplesx1 Array`: frequency vector for PSD
- `sync_radar_offset::Float`: time elapsed since last sync
- `sig_crlb::Float`: Cramer-Rao lower bound of Gaussian noise in sync
- `sync_prf::Float`: repetition frequency of sync process (1/SRI)
- `sync_fs::Float`: sync receiver sampling rate
- `sync_clk_fs::Float`: sample rate of clock error process
- `PRI::Float`: radar pulse repetition interval


"""
function sync_effects_on_PSD(Sphi::Array{Float64,1},f_psd::Array{Float64,1},sync_radar_offset::Float64,sig_crlb::Float64,sync_prf::Float64,sync_fs::Float64, sync_clk_fs::Float64, PRI::Float64, f_osc::Float64)
# Function describes PSD of clk phase after finite offset sync
# Sphi is input PSD

    sync_pri = 1 / sync_prf

    S_n = (2*(sig_crlb*sync_fs.*2*pi).^2).*ones(length(Sphi)).*(sync_pri*1.5).*(sync_fs./f_osc)^2; # scale by (sync_fs*sync_pri) to account for lower bandwidth/downsampling to conserve PSD power

    #adding code 9/25/21 this is the least squares linear predictor
    n_hist = 3 # number of previous samples used in prediction
    cb = sync_pri*(n_hist^2*(n_hist^2-1))/12
    # avoid divide by 0. Will still reduce to original expression for Sphi_tilde
    if (n_hist==1)
        cb = 1
    end

    dtp = sync_radar_offset
    sincpwr = 3 # exponent that sinc functions get raised to

    S_sync1 = zeros(length(f_psd))
    for i=0:(n_hist-1)
        S_sync1 = S_sync1 .+ n_hist.* (i.-(n_hist-1)./2).*(cos.(2 .*pi.*f_psd.*((n_hist-1-i).*sync_pri.+dtp))-cos.(2 .*pi.*f_psd.*((n_hist-1-i).*sync_pri)))
    end
    S_sync2 = zeros(length(f_psd))
    for i=0:(n_hist-1)
        for j=0:(n_hist-1)
            S_sync2 = S_sync2 .+ (n_hist.^2).*(i.*j.+(((n_hist-1).^2)/4)).*exp.(1i.*2 .*pi*f_psd.*(i-j)*sync_pri) .- ((n_hist.^2).*(n_hist-1)).*i.*cos.(2 .*pi.*f_psd.*(i-j).*sync_pri)
        end
    end
    S_n_ls = S_n .* (-1*(dtp/cb).*n_hist.*(n_hist-1) .+ ((dtp/cb).^2).*(n_hist.^3 .* (n_hist.^2 .-1 ))./12)

    S_phi_sync_tilde = Sphi.*((2 .- 2 .*cos.(2 .*pi.*dtp.*f_psd) .-1 .*(2 .* dtp./cb) .*S_sync1 .+1 .*((dtp./cb).^2).*S_sync2))

    Sphi_tilde_alias = downsampled_spectrum(S_phi_sync_tilde,round(sync_clk_fs/(2/sync_pri)),sync_prf).*rect(f_psd,-1/sync_pri, 1/sync_pri+(f_psd[2]-f_psd[1]) , 1)
    Sphi_sync = (Sphi_tilde_alias.+S_n.+S_n_ls) .*abs.(sinc.(f_psd.*PRI)) .^(sincpwr+.25) ./ abs.(sinc.(f_psd.*sync_pri)).^.25


    # # this is the finite offset filter stage
    # Sphi_tilde = (2 .* Sphi .* (1 .- cos.(2*pi*sync_radar_offset.*f_psd) ) )
    #
    #
    # #here is the aliasing error stage from long PRI
    # # Sphi_tilde_alias = downsampled_spectrum(Sphi_tilde,sync_clk_fs,sync_prf)
    # Sphi_tilde_alias = downsampled_spectrum(Sphi_tilde, sync_clk_fs, sync_pri) .*rect(f_psd,-1/sync_pri, 1/sync_pri+(f_psd[2]-f_psd[1]) , 1) # updated from Sam's code 6/21/21
    #
    #
    # # PSD of AWGN with variance determined by CRLB of sync algorithm
    # # updated noise floor power calculation
    # S_n = (2*(sig_crlb*sync_fs .*2*pi ).^2 ) .* ones(length(Sphi)) .* (sync_clk_fs * sync_pri*1.5) # scale by (clk_args.fs*sync_pri) to account for lower bandwidth/downsampling to conserve PSD power
    #
    #
    # # #debugging negative Sphi
    # # idx_test = findall(S_n -> S_n < 0,S_n)
    # # if !isempty(idx_test)
    # #     println("negative psd values")
    # #     println(S_n[idx_test])
    # #
    # #     println(f_psd[idx_test])
    # # end
    #
    # # PSD of clock phase error after synchronization assumes SRI and PRI are different, and SRI is a multiple of PRI
    # Sphi_sync = abs.( (Sphi_tilde_alias .+ S_n) .* ( abs.(sinc.(f_psd .* PRI)).^4.25 ) ./ abs.(sinc.(f_psd.* sync_pri)).^0.25 )

    return Sphi_sync

end#function
#--end--function-------------------------------------------------------------------------------------------


#-start-function-------------------------------------------------------------------------------------------

"""
generates a matrix of estimated CRLB values at each position based on relative positions of platforms and Friis transmission using a mean-network based value

# Arguments
- `pos::3 x N_plat x N_time Array`: vector of platform positions
- `N::Integer`: number of samples in waveform (Period (s) = N/fs)
- `fc::Integer`: center frequency of transmitted sync waveform
- `fs::Integer`: sampling frequency of sync receiver
- `fbw::Integer`: sync waveform bandwidth

"""
function getSensorCRLB_network(pos::Array{Float64,3},N::Int64,fc::Float64,fs::Float64,fbw::Float64)
#function getSensorCRLB_network(pos,N,fc,fs,fbw,Parameters) #TODO add in support for parameter structure
# pos:              gives locations of all platforms in network (3 x N_plat x N_time)
# nplat:            number of platforms
# ntimes:           number of sample times in time_vec
# N:                number of samples in waveform (Period (s) = N/fs)
# fs:               sampling frequency
# fbw:              sync waveform bandwidth
# unpack extra Parameters variables
# T       = Parameters.RF_temp
# c       = Parameters.SOL
# Pt      = Parameters.sync_Pt # sync link transmit power
# Gt      = Parameters.sync_Gt
# Gr      = Parameters.sync_Gr
# POL     = Parameters.sync_pol_loss # sync polarization mismatch factor



# TODO Pass in variable structure with sync power budget information

    szp = size(pos)
    if ndims(pos) == 3
        nplat = szp[2]
        ntimes = szp[3]
    elseif ndims(pos) == 2
        nplat = 1
        ntimes = szp[2]
    end

    k = 1.380649e-23; # Joule/Kelvin -  Boltzmann constant
    T = 25 + 273.15; # Kelvin -  let approximate temperature be 25C
    Np = k*T*fbw; # (W) noise power in receiver - NOTE: this doesn't consider the noise figure in receiver electronics
    c = 299792458; # m/s SOL

    # initialize SNR and CRLB matrices
    snr_dB = Array{Float64,3}(undef, nplat, nplat, ntimes)
    sig_crlb = Array{Float64,2}(undef, nplat, ntimes)

    for i = 1 : nplat # number of platforms
        for j = 1 : ntimes # loop over number of sample times
            for l = 1 : nplat
                if i!=l
                  tx_pos = pos[:,i,j] # position of tx platform
                  rx_plat = pos[:,l,j] # position of a rx platform

                  # find R, distance between the platforms
                  R = sqrt((rx_plat[1]-tx_pos[1])^2+(rx_plat[2]-tx_pos[2])^2+(rx_plat[3]-tx_pos[3])^2)

                  # need Pt, Gt, Gr, POL. These could be calculated in the orbits module using antenna pattern/platform rotations/
                  Pt = 1e-6     # let transmit power be 1 uW (needs to be evaluated in trade study, more likely a consequence of orbit design)
                  Gt = 1        # Isotropic rx/tx antennas for sync
                  Gr = 1        # Isotropic rx/tx antennas for sync
                  POL = 1       # Polarization mismatch factor (let antennas be aligned here) - 0 is total misalignment, 1 is total alignment


                  # use Friis formula to calculate received power, SNR
                  Pr = POL * Pt*Gt*Gr*c^2/(4*pi*fc*R)^2;
                  snr_linear = Pr/Np;
              else
                  snr_linear = 0; # doesn't contribute to mean SNR (divide by nplat-1 later)
              end#if
              snr_dB[i,l,j] = 10*log10(snr_linear); # save SNR for debug/trade study
            end #nplat

            snr_avg = sum(10 .^(snr_dB[i,:,j]./10) )/ (nplat-1) # find mean SNR to all platforms

            sigma_2TOF = sqrt( 3/ ( (2*pi*fbw)^2 * snr_avg * (N/fs) * fs ) ) # this may become a useful parameter for the positioning module/task

            sig_crlb[i,j] = sqrt( ( (nplat-1) / nplat ) * sigma_2TOF^2) # latest equation from Sam 06/17/2021
        end # for ntimes
    end #for nplat

    return (sig_crlb, snr_dB) # return tuple of CRLB and SNR values
end#function
#--end--function-------------------------------------------------------------------------------------------

#-start-function-------------------------------------------------------------------------------------------
"""
calculates the effects of downsampling on the phase error PSD

# Arguments
- `Sphi::Nx1 Array`: 2 sided oscillator phase error PSD
- `sync_clk_fs::Float`: sample rate of clock error process
- `sync_prf::Float`: repetition frequency of sync process (1/SRI)

"""
function downsampled_spectrum(Sphi::Array{Float64,1},sync_clk_fs::Float64,sync_prf::Float64)
    # M: downsampling factor

    #TODO: set maximum M value in order to reduce computational time? High SRI values get very excessive..
    M = round(sync_clk_fs/(2.0*sync_prf))
    
    #testing: set maximum vaue for M at 1000. The extra spectral folding doesn't change much in the end result but makes it way slower
    if M > 1000
        M = 1000
    end

    N = length(Sphi)
    sMfft=zeros(N) #initialize

    for i = 1 : M
        sMfft = sMfft + shift(Sphi,(i-1)*N/M)
    end


    #testing distributed approach
    ## set-up distributed  TODO: add in processing variables into input parameters? something like setting the number of cores/threads etc.
    # nproc = nprocs()
    # if nproc ==1
    #     addprocs(1) # testing 2 processes
    # end
    # sMfft   = @distributed (+) for i in 1:M
    #      shift(Sphi,(i-1)*N/M)
    # end

    return sMfft

end#function
#--end--function-------------------------------------------------------------------------------------------

#-start-function-------------------------------------------------------------------------------------------
"""
this function produces an all-real output of a shift of input vector x by n samples

# Arguments
- `x::Nx1 Array`: input vector to be shifted
- `n::Integer`: number of samples to shift by

"""
function shift(x,n)
# this function produces an all-real output of a shift of input vector x by n samples

    if round(n)==n    # this will almost always be the case since we're passing in a round() value
        z = circshift(x,n)
    else
        N = length(x)
        t = collect(-1/2 : 1/N : 1/2 - 1/N)
        maxx = maximum(abs.(x))
        xfft = fftshift(fft(x))
        xfft = xfft.*exp.(-1im*2*pi*n.*t)
        z = ifft(ifftshift(xfft))
        z = real(z)
        scale = sqrt( sum(abs.(x).^2) / sum(z.^2) )
        z = scale.*z
    end
     return z
end#function
#--end--function-------------------------------------------------------------------------------------------

"""
This function reads measured PSD data in from an excel file, returns frequency vector and PSD amplitude

# Arguments
- `filename::String`: filename of PSD data
"""
function read_PSD_excel_data(filename::String)
    xf           = XLSX.readxlsx(filename)
    sh           = xf["Sheet1"]
    
    rn=0; # this loop finds number of rows in sheet
    for r in XLSX.eachrow(sh)
        rn = XLSX.row_number(r)
    end
    
    f_psd_meas = Array{Float64,1}(undef,rn)
    osc_psd_meas = Array{Float64,1}(undef,rn)
    
    for r in XLSX.eachrow(sh)
        rn = XLSX.row_number(r)
        f_psd_meas[rn] = r[1]
        osc_psd_meas[rn] = r[2]
    end

    f_psd_meas   = convert(Array{Float64},f_psd_meas)
    osc_psd_meas = convert(Array{Float64},osc_psd_meas)

    return f_psd_meas, osc_psd_meas
end #function
end #module


for r in XLSX.eachrow(sh)
    rn = XLSX.row_number(r)
    println(rn)
end