module Error_Sources
using ..Sync
using Parameters
using JLD2, Plots


function random_noise_raw(rawdata, params) # based on SNR input (obsolete)
    @unpack SNR,enable_fast_time,mode = params
    
    Np_RX=size(rawdata)[2] # number of RX platforms
    if mode==3;Np_TX=size(rawdata)[3];end # number of TX platforms
    Nst=size(rawdata)[1] # number of slow-time samples
    add_noise_amp=10^(-SNR/20); # relative additive noise amplitude (set to 0 for no additive random noise)
    if enable_fast_time
        Nft=size(rawdata)[end] # number of fast-time samples
        add_rnd_noise=(add_noise_amp/2^0.5)*(randn(Nst,Np_RX,Nft)+im*randn(Nst,Np_RX,Nft)); # additive complex random noise, unique value for each platform,  fast-time sample, and slow-time sample
    else
        add_rnd_noise=(add_noise_amp/2^0.5)*(randn(Nst,Np_RX)+im*randn(Nst,Np_RX)); # additive complex random noise, unique value for each platform and slow-time sample
    end
    for s=1:Nst # slow-time (pulses)
        for i=1:Np_RX # RX platform
            if enable_fast_time
                if mode==1 || mode==2 # SAR (ping-pong) or SIMO
                    rawdata[s,i,:]=rawdata[s,i,:]+add_rnd_noise[s,i,:]
                elseif mode==3 # MIMO
                    for k=1:Np_TX # TX platform for MIMO
                        rawdata[s,i,k,:]=rawdata[s,i,k,:]+add_rnd_noise[s,i,:]
                    end
                end
            else
                if mode==1 || mode==2 # SAR (ping-pong) or SIMO
                    rawdata[s,i]=rawdata[s,i]+add_rnd_noise[s,i]
                elseif mode==3 # MIMO
                    for k=1:Np_TX # TX platform for MIMO
                        rawdata[s,i,k]=rawdata[s,i,k]+add_rnd_noise[s,i]
                    end
                end
            end
        end
    end
    return rawdata
end

function random_noise_raw_new(rawdata, params) # obsolete
    @unpack SNR,mode = params
    
    Noise = 30

    Np_RX=size(rawdata)[2] # number of RX platforms
    if mode==3;Np_TX=size(rawdata)[3];end # number of TX platforms
    Nst=size(rawdata)[1] # number of slow-time samples
    add_noise_amp=10^(Noise/10); # relative additive noise amplitude (set to 0 for no additive random noise)
    Nft=size(rawdata)[end] # number of fast-time samples
    add_rnd_noise=(add_noise_amp/2^0.5)*(randn(Nst,Np_RX,Nft)+im*randn(Nst,Np_RX,Nft)); # additive complex random noise, unique value for each platform,  fast-time sample, and slow-time sample
    rawdata2 = zeros(ComplexF64, size(rawdata))

    for s=1:Nst # slow-time (pulses)
        for i=1:Np_RX # RX platform
            if mode==1 || mode==2 # SAR (ping-pong) or SIMO
                rawdata2[s,i,:]=rawdata[s,i,:]+add_rnd_noise[s,i,:] 
            elseif mode==3 # MIMO
                for k=1:Np_TX # TX platform for MIMO
                    rawdata2[s,i,k,:]=rawdata[s,i,k,:]+add_rnd_noise[s,i,:]
                end
            end
        end
    end
    return rawdata2, add_rnd_noise
end

function random_noise_image(image_3D, params, ref_pix_ind, Lsa, Nst, Gtx_ref, Grx_ref, sig0_ref, R_ref, θ_ref, α_ref)
    @unpack pulse_length,bandwidth,T_noise,pow_tx,Rx_bw_factor,L_sys,fp,λ,c,kBolt = params
    # Note: image_3D must be a complex image! TODO remove abs() from process_raw_data functions
    # new inputs
    # ref_pix_ind: image indices for the reference pixel
    # Lsa: synthetic aperture length
    # Nst: number of slow-time pulses in Lsa
    # Gtx_ref & Grx_ref: Tx/Rx antenna gains for the reference scene pixel
    # sig0_ref: sigma-0 for the reference scene pixel
    # R_ref: slant range for the reference scene pixel
    # θ_ref & α_ref: incidence angle and slope for the reference scene pixel
    
    # Calculate final SNR after processing for the reference scene pixel
    NESZ_ref = (4^4 * pi^3 * kBolt * T_noise * bandwidth * Rx_bw_factor * L_sys * Lsa * R_ref^3 * sin((θ_ref-α_ref)*pi/180)) / (pow_tx * Nst * pulse_length * λ^3 * Gtx_ref * Grx_ref * c)
    SNR_ref = sig0_ref / NESZ_ref # linear scale

    # Calculate noise level for the reference scene pixel
    ref_pix_sig =  image_3D[ref_pix_ind]
    Npix=size(image_3D) # number of pixels
    noise_lvl = (ref_pix_sig/2^0.5/SNR_ref^0.5)*(randn(Npix)+im*randn(Npix))

    # Add noise to all image pixels
    noisy_image = image_3D + noise_lvl
    return noisy_image
end

function ADC_sampling(rawdata,params)
    @unpack fs_ADC,Δt = Parameters
    ds_rate = 1 / (fs_ADC * Δt) # downsampling rate (for now requires fs_ADC to be an integer multiple of Δt)
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata_ds = rawdata[:,:,1:ds_rate:end]
        Nft_ds = size(rawdata_ds,3)
    elseif mode==3 # MIMO
        rawdata_ds = rawdata[:,:,:,1:ds_rate:end]
        Nft_ds = size(rawdata_ds,4)
    end
    return rawdata_ds, Nft_ds
end

"""
Calculates oscillator phase error and sychronization effects. Returns the raw data with phase errors.
# Arguments
- `rawdata::Num slow time points x Num platforms Array x Num fast time points`: complex valued raw data. Num slow time points x Num platforms Array if enable_fast_time = false
- `orbit_pos_interp::3xNplatform x N time Array`: positions of platforms at pulse times
- `t_xyz_3xN::3xN Array`: array of target locations
- `enable_fast_time::Boolean`: flag if fast time is used
- `parameters::Parameters`: structure of simulation configuration parameters
- `sync_PSD::Nplatform x Npulses x sync_clk_fs Array`: OPTIONAL INPUT, precalculated synchronization PSDs
"""
function synchronization_errors!(rawdata,slow_time,orbit_pos_interp, t_xyz_3xN, osc_coeffs, params)
    @unpack mode, sync_master = params
    
    Np_RX=size(rawdata)[2] # number of RX platforms
    Nst=size(rawdata)[1] # number of slow-time samples
    if mode==3
        Np_TX=size(rawdata)[3]# number of TX platforms
    end 
    
    ## Synchronization Effects
    (phase_err, sync_PSDs) = Sync.get_sync_phase(slow_time,orbit_pos_interp, t_xyz_3xN, osc_coeffs, params) 
    # note: phase_err is (Nplat x N slow-time) for modes 1 & 2, but (Nplat x Nplat x N slow-time) for MIMO
    # for MIMO, first axis is the transmitting platform number, 2nd is receive platform, 3rd is slow-time number

    ## combine with raw data

        Nft=size(rawdata)[end] # number of fast-time sampless
        if mode == 1 #ping-pong
            for s = 1 : Nst # slow-time (pulses)
                for i = 1 : Np_RX
                    rawdata[s,i,:] = rawdata[s,i,:] .* exp(im*(phase_err[i,s,1] - phase_err[i,s,2])) # because tx and rx are the same platform
                end#N platforms
            end#slow time
        elseif mode == 2 # SIMO
            for s = 1 : Nst # slow-time (pulses)
                tx_phase = phase_err[sync_master, s , 1] # transmitter phase error state 
                for i = 1 : Np_RX
                    rawdata[s,i,:] = rawdata[s,i,:].*exp(im*(phase_err[i,s,2] - tx_phase) )
                end#N platforms
            end#slow time
        elseif mode == 3 #MIMO
            for s = 1 : Nst # slow-time (pulses)
                for i = 1 : Np_TX # Tx platform for MIMO
                    for k = 1 : Np_RX # Rx platform for MIMO
                        rawdata[s,k,i,:] = rawdata[s,k,i,:].*exp(im*(phase_err[i,i,s,1] - phase_err[i,k,s,2]) ) #tx phase at tx time + rx phase at rx time
                    end#N Rx platforms
                end#N Tx platforms
            end#slow time
        end#mode
    
    return rawdata
end#sync_error function

#this is the same function but is overwritten with the synchronization PSD input
function synchronization_errors(rawdata,slow_time,orbit_pos_interp,enable_fast_time,parameters,sync_PSDs)
    mode   = parameters.mode
    master = parameters.master
    
    Np_RX=size(rawdata)[2] # number of RX platforms
    Nst=size(rawdata)[1] # number of slow-time samples
    if mode==3
        Np_TX=size(rawdata)[3]# number of TX platforms
    end 
    
    ## Synchronization Effects
    phase_err = Sync.get_sync_phase(slow_time,orbit_pos_interp,parameters,sync_PSDs) 
    # note: phase_err is (Nplat x N slow-time) for modes 1 & 2, but (Nplat x Nplat x N slow-time) for MIMO
    # for MIMO, first axis is the transmitting platform number, 2nd is receive platform, 3rd is slow-time number
    
    ## combine with raw data
    if enable_fast_time
        
        Nft=size(rawdata)[end] # number of fast-time samples
        if mode == 1 #ping-pong
            for s = 1 : Nst # slow-time (pulses)
                for i = 1 : Np_RX
                    rawdata[s,i,:] = rawdata[s,i,:] .* exp(im*2*phase_err[i,s])' # 2x because tx and rx are the same platform
                end#N platforms
            end#slow time
        elseif mode == 2 # SIMO
            for s = 1 : Nst # slow-time (pulses)
                tx_phase = phase_err[master,s] # transmitter phase error state 
                for i = 1 : Np_RX
                    rawdata[s,i,:] = rawdata[s,i,:].*exp(im*(phase_err[i,s] + tx_phase) )
                end#N platforms
            end#slow time
        elseif mode == 3 #MIMO
            for s = 1 : Nst # slow-time (pulses)
                for i = 1 : Np_TX # Tx platform for MIMO
                    for k = 1 : Np_RX # Rx platform for MIMO
                        rawdata[s,k,i,:] = rawdata[s,k,i,:].*exp(im*(phase_err[i,i,s] + phase_err[i,k,s]) ) #tx phase + rx phase at tx time
                    end#N Rx platforms
                end#N Tx platforms
            end#slow time
        end#mode
        
    else # no fast time
        if mode == 1 #ping-pong
            for s = 1 : Nst # slow-time (pulses)
                for i = 1 : Np_RX
                    rawdata[s,i] = rawdata[s,i]*exp(im*2*phase_err[i,s]) # 2x because tx and rx are the same platform
                end#N platforms
            end#slow time
        elseif mode == 2 # SIMO
            for s = 1 : Nst # slow-time (pulses)
                tx_phase = phase_err[master,s] # transmitter phase error state 
                for i = 1 : Np_RX
                    rawdata[s,i] = rawdata[s,i]*exp(im*(phase_err[i,s] + tx_phase) )
                end#N platforms
            end#slow time
        elseif mode == 3 #MIMO
            for s = 1 : Nst # slow-time (pulses)
                for i = 1 : Np_TX # Tx platform for MIMO
                    for k = 1 : Np_RX # Rx platform for MIMO
                        rawdata[s,k,i] = rawdata[s,k,i]*exp(im*(phase_err[i,i,s] + phase_err[i,k,s]) ) #tx phase + rx phase at tx time
                    end#N Rx platforms
                end#N Tx platforms
            end#slow time
        end#mode    
end#if fast time

    return rawdata
end#sync_error function #2

function position_errors
end

end#module