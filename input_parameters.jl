module input_parameters


struct Parameter
    
    ##### planet parameters
    radius_a      # equitorial radius
    eccentricity  # eccentricity
    SOL           # (m/s) - speed of light
    ##### radar parameters
    mode        #1: SAR (ping-pong), 2:SIMO, 3:MIMO
    tx_el       # which element transmits for SIMO (max value N)
    fc          # center frequency (Hz)
    fp          # pulse repetition frequency (Hz)
    
    ##### platform location parameters
    orbit_filename  # filename to take orbit data from
    SAR_duration    # synthetic aperture duration (s)
    SAR_start_time  # SAR imaging start time (s)
    nplat           # number of platforms
    
    ##### target location parameters (volumetric grid) defined in geo (θϕh)
    
    t_θ         # deg latitude
    t_ϕ         # deg longitude
    t_h         # m  heights
    
    ###### image/scene pixel coordinate parameters
    s_θ         # deg latitude
    s_ϕ         # deg longitude
    s_h         # m  heights
    
    ##### range spread function (RSF) parameters
    Trx         # (s) - duration of RX window (may need to be increased if aperture or scene is large) TODO (adjust based on max/min range)
    τ           # (s) - pulse length
    Δt          # (s) - fast-time resolution (ADC sampling rate effect is excluded for now)
    B           # (Hz) - bandwidth
    
    ##### sync clock parameters
    sync_pri                # (s) - repetition interval of sync
    sync_processing_time    # (s) - processing time between stage 1 and stage 2 sync
    sync_signal_len         # sync signal's waveform length
    sync_fc                 # (Hz) - sync signal's center frequency
    sync_fs                 # (samp/s) - sync sensor sampling clock frequency
    sync_fbw                # (Hz) - sync receiver bandwidth
    sync_fmin               # (Hz) - minimum frequency > 0 in Hz to window clock phase error PSD
    f_osc                   # (Hz) - local oscillator frequency
    sync_clk_fs             # (samp/s) -  sample rate of clock error process
    a_coeff_dB              # oscillator quality coefficients (either 5x1 used for all platforms, or 5 x nplat to define individually
    osc_coeffs              # oscillator quality coefficients (matrix of values for each receiver)
    
end #struct

function Parameter()
    ##### planet parameters
    radius_a = 6378.137e3
    eccentricity = sqrt(0.00669437999015)
    SOL = 299792458
    
    ##### radar parameters
    mode = 1
    tx_el = 1
    fc = 1e9
    fp = 3
    
    ##### platform location parameters
    SAR_duration = 2
    SAR_start_time = 0
    nplat = NaN
    
    # grab number of platforms based on orbit data
    # orbit_dataset=Dataset(orbit_filename) # Read orbits data in NetCDF format
    # orbit_pos=orbit_dataset["position"][:,:,1] # grab first data point
    # nplat=size(orbit_pos)[2] 
    
    ##### target location parameters
    
    t_θ = 0
    t_ϕ = 0
    t_h = 0
    
    ###### image/scene pixel coordinate parameters
    s_θ = -0.0005:0.00002:0.0005
    s_ϕ = -0.001:0.00002:0.001
    s_h = 0
    
    ##### range spread function (RSF) parameters
    Trx = 300e-6
    τ = 10e-6
    Δt = 1e-8
    B = 10e6
    
    ##### sync clock parameters
    sync_pri = 10
    sync_processing_time = 0.001
    sync_signal_len = 1024
    sync_fc = 1e9
    sync_fs = 25e6
    sync_fbw = sync_fs
    sync_fmin = .01
    f_osc = 10e6
    sync_clk_fs =  20e3;
    a_coeff_dB = [-28 -40 -200 -130 -155] # [USRP E312]
    #a_coeff_dB = [-95 -90 -200 -130 -155] # [Krieger]
    
    return Parameter(radius_a,
        eccentricity,
        SOL,
        mode,
        tx_el,
        fc,
        fp,
        orbit_filename,
        SAR_duration,
        SAR_start_time,
        nplat,
        t_θ,
        t_ϕ,
        t_h,
        s_θ,
        s_ϕ,
        s_h,
        Trx,
        τ,
        Δt,
        B,
        sync_pri,
        sync_processing_time,
        sync_signal_len,
        sync_fc,
        sync_fs,
        sync_fbw,
        sync_fmin,
        f_osc,
        sync_clk_fs,
        a_coeff_dB,
        osc_coeffs
        )
    end#function
    
    ##TODO Add any validation or retrieval (get wavelength etc.) functions here
    function setClkCoeff(a_coeff_dB) #TODO test function
        if size(a_coeff_dB)[1] == 1
            osc_coeffs = repeat(a_coeff_dB,nplat)
        else
            osc_coeffs = a_coeff_dB
        end #if
    end#function
    
    function getNumPlat(new_filename...) #TODO test function
        if length(new_filename) > 0 # if function is passed a new filename
            # grab number of platforms based on orbit data
            orbit_dataset=Dataset(new_filename) # Read orbits data in NetCDF format
            orbit_pos=orbit_dataset["position"][:,:,1] # grab first data point
            nplat=size(orbit_pos)[2] 
        elseif isnan(nplat) # empty argument, but nplat not defined within structure yet
             # grab number of platforms based on orbit data
             orbit_dataset=Dataset(orbit_filename) # Read orbits data in NetCDF format
             orbit_pos=orbit_dataset["position"][:,:,1] # grab first data point
             nplat=size(orbit_pos)[2] 
         end#if
             
         return nplat # should return nplat after it has been defined as a number (not NaN)
     end#function
    
end#module
    