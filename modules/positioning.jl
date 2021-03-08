module positioning
using Random
using LinearAlgebra
using Statistics #for debugging

# TODO  more code here, will call the GPS observation model (GPSobsErr2)



##---------------------------------------------------------------------------------
# Reproduction of function "GPSobsErr2.m"
function GPSobsErr2(x_rx,x_tx,P_tx,time)
    #x_rx (6xn): the rx position and velocity vector
    #x_tx (6xn): the tx position and velocity vector
    #P_tx (1xn): the power of the signals at rx antenna
    #time : (integer number), time at which the measurement is computed
    #n : number of visible GPS sats
    
    # Assume we take inputs in units of meters and m/s rather than the internal km used here.
    # We convert to km and km/s because it's easier than changing everything else within the function
    x_rx = x_rx./1000
    x_tx = x_tx./1000
    
    
    size_x_tx = size(x_tx,2); # number of visible satellites
    
    #instantiate vars
    pseudorange = zeros(size_x_tx);
    range_thermal_noise = zeros(size_x_tx);
    pseudorange_rate = zeros(size_x_tx);
    range_rate_thermal_noise = zeros(size_x_tx);
    pseudorange_error = zeros(size_x_tx);
    phase_meas = zeros(size_x_tx);
    phase_meas_error = zeros(size_x_tx);
    
    # error budget [1]
    rx_clock_offset = 10 # [km] initial receiver clock bias
    rx_clock_drift = 0.1 # [km/sec] receiver clock drift
    receiver_noise = 0.1e-3
    
    if norm(x_rx[1:3])< 6400+1000
        iono_error = 7e-3
    else
        iono_error = 0 #15e-3; # [km] std of the residual iono delay after corrections are applied
    end#if
        
    tropo_error = 0;#0.2e-3; # [km] std of the residual tropo delay after corrections are applied
    multipath_error = 0.2e-3 # [km] multipath
    tx_ephemeris = 0.8e-3 # [km] broadcast ephemeris error
    tx_clock = 1.1e-3 # [km] broadcast clock error
    
    f_carrier = 1575.42e6 # [Hz] carrier frequency
    c_light = 299792.458
    wavelength = c_light/f_carrier;
    multipath_phase_error = 0.005e-3 #/wavelength; % [km] multipath carrier (phase)
    antenna_pcv_error =  0.01e-3 #/wavelength; % [km] antenna phase center variation
    
    # This code can be used to produce estimates of the code and carrier measurements. Range is the true value, added to the estimation errors
    # receiver_mat=repeat(x_rx[1:3]',size_x_tx); # could be used to vectorize, but won't do here
    
    # # range = sqrt(dot((x_tx[1:3,:]-receiver_mat),(x_tx[1:3,:]-receiver_mat))); # [km] range [3] -- why is this here? it just gets overwritten
    # range=zeros(size_x_tx)
    # for j = 1:size_x_tx
    #     range[j] = norm(x_tx[1:3,j]-x_rx[1:3])
    #     # range[j] = distance(x_tx[1:3,j],x_rx[1:3])
    # end#for
    # 
    # 
    # range_rate=zeros(size_x_tx)
    # for j = 1:size_x_tx
    #     # range_rate[j] = dot( (x_tx[1:3,j]-receiver_mat), (x_tx[4:6,:] .- x_rx[4:6]') ) ./range[j]  # [km/sec] range rate [4]
    #     range_rate[j] = dot( (x_tx[1:3,j] .- x_rx[1:3]), (x_tx[4:6,j] .- x_rx[4:6]) ) ./range[j]  # [km/sec] range rate [4]        
    # end#for
    

    # range_thermal_noise(i) = sqrt(A/C2N*(1+B/C2N))*c/f_code; % [km] std of the thermal noise on pseudorange [2]

    length_pseudorange=size_x_tx; # = length(range)
    (range_thermal_noise,range_rate_thermal_noise,carrier_thermal_noise_k) = thermal_noise_gps_L1(size_x_tx,P_tx)
    
    Tx_ephemeris_err = tx_ephemeris .*randn(size_x_tx)
    Tx_clock_err = tx_clock .*randn(size_x_tx)
    Iono_err = iono_error .*randn(size_x_tx)
    Tropo_err = tropo_error .*randn(length_pseudorange)
    Multipath_err = multipath_error .*randn(length_pseudorange)
    Receiver_noise_err = receiver_noise .*randn(length_pseudorange)
    Range_thermal_noise_err = range_thermal_noise .*randn(length_pseudorange)
    
    # --------------- Pseudorange ---------------%
    # pseudorange = range' + rx_clock_offset + rx_clock_drift*(t) + Tx_ephemeris_err + Tx_clock_err ...
    #     + Iono_err + Tropo_err + Multipath_err + Receiver_noise_err + Range_thermal_noise_err ; % [km] pseudorange [1]
    
    pseudorange_error = Tx_ephemeris_err .+ Tx_clock_err .+ Iono_err .+ Tropo_err .+ Multipath_err .+ Receiver_noise_err .+ Range_thermal_noise_err

    # --------------- Range rate ---------------%
    
    
    Range_rate_thermal_noise_err = range_rate_thermal_noise .* randn(length_pseudorange)
    
    #pseudorange_rate = range_rate .+ rx_clock_drift .+ Range_rate_thermal_noise_err # [km/sec] pseudorange rate [6]
    pseudorange_rate_error =  rx_clock_drift .+ Range_rate_thermal_noise_err 
    
    # --------------- Phase meas ---------------%
    
    Multipath_phase_err = multipath_phase_error .*randn(length_pseudorange)
    Antenna_pcv_err = antenna_pcv_error .* randn(length_pseudorange)
    Carrier_thermal_noise_err = carrier_thermal_noise_k .*randn(length_pseudorange)
    Ambiguity = rand(-1000:1000, length_pseudorange)
    
    # -------------- Windup --------------------%
    # %windup1 = Windup(x_rx',x_tx);
    # % windup is significant, but could be mostly compensated if you have 
    # % really good knowledge of your spacecraft attitude (1 arcsec).
    windup = Windup(x_rx',x_tx) .*c_light ./ f_carrier;

    #phase_meas  = range' .+ rx_clock_offset .+ rx_clock_drift .* time .+ Tx_ephemeris_err .+ Tx_clock_err .- Iono_err .+ Tropo_err .+ Multipath_phase_err .+ Antenna_pcv_err .+ Ambiguity.*wavelength .+ Carrier_thermal_noise_err .+ windup
    
#    %add error from phase wind up here
    
    # % phase meas with DD
    # phase_meas  = range' + Multipath_phase_err + Antenna_pcv_err + Carrier_thermal_noise_err;
    
    # complete phase measurement err
    phase_meas_error = Tx_ephemeris_err .+ Tx_clock_err .- Iono_err .+ Tropo_err .+ Multipath_phase_err .+ Antenna_pcv_err .+ Ambiguity.*wavelength .+ Carrier_thermal_noise_err .+ windup;
    
    # simplification by removing all the stuff that can be canceled out by DD
    # only valid for short separation (<10km)
    phase_meas_error_DD = Multipath_phase_err .+ Antenna_pcv_err .+ Ambiguity.*wavelength .+ Carrier_thermal_noise_err .+ windup;
    
    beatphase_error = Multipath_phase_err .+ Antenna_pcv_err .+ Carrier_thermal_noise_err .+ windup
    #and here
    
    # Range_rate_thermal_noise_err1 = Range_rate_thermal_noise_err' # unused variables??
    # Range_thermal_noise_err1 = Range_thermal_noise_err' #???
    
    # convert everything back to units of meters from km for return
    pseudorange_error       = pseudorange_error      .*1000
    pseudorange_rate_error  = pseudorange_rate_error .*1000
    phase_meas_error        = phase_meas_error       .*1000
    phase_meas_error_DD     = phase_meas_error_DD    .*1000
    beatphase_error         = beatphase_error        .*1000
    
    return (pseudorange_error,pseudorange_rate_error, phase_meas_error,phase_meas_error_DD, beatphase_error)
    
    # % References:
    # 
    # % [1] E. Kaplan, C. Hegarty; "Understanding GPS: principles and applications" 2nd edition, Artech house, 2006, pp. 302- 322;
    # % [2] V. Capuano, C. Botteron, Y. Wang, J. Tian, J. Leclere, P. A. Farine; "GNSS/INS/Star Tracker Integrated Navigation System for Earth-Moon Transfer Orbit", ION GNSS+ 2014;
    # % [3] E. Kaplan, C. Hegarty; "Understanding GPS: principles and applications" 2nd edition, Artech house, 2006, pp. 52;
    # % [4] E. Kaplan, C. Hegarty; "Understanding GPS: principles and applications" 2nd edition, Artech house, 2006, pp. 59;
    # % [5] D. Borio, N. Sokolova, G. Lachapelle; "Doppler Measurement Accuracy in Standard and High-Sensitivity GNSS Receivers", IET RADAR, SONAR & NAVIGATION, NOVEMBER 2010;
    # % [6] E. Kaplan, C. Hegarty; "Understanding GPS: principles and applications" 2nd edition, Artech house, 2006, pp. 61;
    # % [7] M.L Psiaki, S. Mohiuddin; "Modeling, Analysis, and Simulation of GPS Carrier Phase for Spacecraft Relative Navigation"

end#function


#----Support functions-------------
function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
    return dist
end#function
#-----------------------------------------------------
function dot(a,b)
    c = sum(conj(a).*b)
    return c
end#function
#-----------------------------------------------------
function thermal_noise_gps_L1(length_pseudorange,Power_k)
    #function [thermal_noise_k,rate_thermal_noise_k,carrier_thermal_noise_k] = thermal_noise_gps_L1(length_pseudorange,Power_k)
    # #% thermal_noise2.m
    # 
    # % This function estimate the thermal noise on pseudorange and pseudorange
    # % rate from the GPS measures
    # %------------------------------------------
    # %%%%%%%%%% Input %%%%%%%%%%
    # % length_pseudorange : length of matrix (nx1) GPS pseudorange measurements
    # %------------------------------------------
    # %%%%%%%%%% Output %%%%%%%%%%
    # % thermal_noise_k : matrix (nx1) estimated thermal noise on pseudoranges
    # % rate_thermal_noise_k : matrix (nx1) estimated thermal noise on pseudoranges rate
    # %------------------------------------------
    f_code = 1.023e6; # [Hz] code frequency
    f_carrier = 1575.42e6; # [Hz] carrier frequency
    
    B_fe = 13e6; # [Hz] double-sided front-end bandwidth [2]
    R_c = 1.023e6; # [chip/sec] chipping rate [2]
    T_c =1/R_c;
    D = 0.5; # [chip] early-to-late correlator spacing [2]
    B_n = 4.5; # [Hz] code loop noise bandwidth [2]
    #B_n2 = 1; # [Hz] code loop noise bandwidth [2]
    B_nPLL = 25; # [Hz] carrier loop noise bandwidth [2]
    T = 10e-3; # [sec] coherent integration time [2]
    F = 1; # = 1 for high C2N0 and 2 for low C2N
    
    
    C2N = 10 .^((Power_k.+174)/10) # [Hz] carrier to noise ratio [2]
    
    if (D >= pi * R_c/B_fe)
        thermal_noise_k         = sqrt.( ( (B_n./2 ./C2N).*D) .* (1 .+ (2 ./( T.* (2 .-D) ) ) ./ C2N) ).*299792.458./f_code # [km] std of the thermal noise on pseudorange [2]
        rate_thermal_noise_k    = (4 * F * B_n ./C2N .*(1 .+ 1 ./ (2 .* T .* C2N) ) ).^0.5 .*299792.458 ./f_carrier ./ (2*pi) ./T # [km] std of the thermal noise on the pseudorange rate [3] eq.28 (for PLL)
        carrier_thermal_noise_k = (B_nPLL ./ C2N .* (1 .+1 ./ (2 .* T.* C2N)) ).^0.5.*299792.458 ./f_carrier ./2 ./pi; # [km] std of the PLL thermal noise on the pseudorange rate [1]
        
    elseif (D <= R_c/B_fe)
        thermal_noise_k         = sqrt.((B_n ./ 2 ./C2N).*(1/(B_fe*T_c)).*(1+(1/T)./C2N))*299792.458/f_code; # [km] std of the thermal noise on pseudorange [2]
        rate_thermal_noise_k    = (B_n ./ C2N .* (1 .+ 1 ./ (2 .* T .* C2N) ) ) .^0.5 .*299792.458 ./f_carrier ./ (2 * pi) ./T; # [km] std of the thermal noise on the pseudorange rate [3]
        carrier_thermal_noise_k = (B_nPLL ./ C2N .* (1 .+1 ./ (2 .* T.* C2N)) ).^0.5.*299792.458 ./f_carrier ./2 ./pi; # [km] std of the PLL thermal noise on the pseudorange rate [1]

    else
        thermal_noise_k         = sqrt.((B_n /2 ./C2N) .* ((1 ./ (B_fe.*T_c)) .+ (B_fe*T_c./(pi-1)) .* (D.-(1 ./( B_fe.*T_c))).^2 ).* (1 .+( 2 ./(T.*(2 .-D )))./ C2N)) .*299792.458 ./f_code; # [km] std of the thermal noise on pseudorange [2]
        rate_thermal_noise_k    = (B_n2./C2N.* (1 .+1 ./(2 .*T.*C2N))).^0.5 .*299792.458 ./f_carrier/(2*pi)/T; # [km] std of the thermal noise on the pseudorange rate [3]
        carrier_thermal_noise_k = (B_nPLL ./ C2N .* (1 .+1 ./ (2 .* T.* C2N)) ).^0.5.*299792.458 ./f_carrier ./2 ./pi; # [km] std of the PLL thermal noise on the pseudorange rate [1]
    end#if
    
    
    # # change output units to METERS from km
    # thermal_noise_k = thermal_noise_k .*1000
    # rate_thermal_noise_k = rate_thermal_noise_k.*1000
    # carrier_thermal_noise_k = carrier_thermal_noise_k.*1000
    
    return (thermal_noise_k, rate_thermal_noise_k, carrier_thermal_noise_k) # return outputs
    
    
    # % References:
    # [1] J. Perello, A. Garcia; GNSS navigation on the Moon (preliminary link budgets), Technical report, 2010.
    # [2] V. Capuano, C. Botteron, Y. Wang, J. Tian, J. Leclere, P. A. Farine; "GNSS/INS/Star Tracker Integrated Navigation System for Earth-Moon Transfer Orbit", ION GNSS+ 2014;
    # [3] D. Borio, N. Sokolova, G. Lachapelle; "Doppler Measurement Accuracy in Standard and High-Sensitivity GNSS Receivers", IET RADAR, SONAR & NAVIGATION, NOVEMBER 2010;
    # [4] Psiaki and Mohiuddin "Modeling, Analysis, and Simulation of GPS Carrier Phase for Spacecraft Relative Navigation", 2007
end#function


#-----------------------------------------------------
function Windup(x_rx,x_tx)
    # % x_rx : Matrix (6X1), true position and velocity of the receiver;
    # % x_tx : Matrix (6Xn), true position and velocity of the n GNSS satellites;
    # % j = gps satellite, a = reciever
    # % output: windup [nx1]

    n = size(x_tx,2);
    windup=zeros(n)
    for i = 1 : n
        x_tx_i = x_tx[1:3,i];
        vx_tx_i = x_tx[4:6,i];
        p = x_rx[1:3] - x_tx_i;
        p = p./norm(p);

        y_j = x_tx[4:6,i]/norm(vx_tx_i);
        z_j = x_tx[1:3,i]/norm(x_tx[1:3,i]);
        x_j = cross(y_j,z_j);
    
        y_a = x_rx[4:6]/norm(x_rx[4:6]);
        z_a = -x_rx[1:3]/norm(x_rx[1:3]);
        x_a = cross(y_a,z_a);
        
        xeffr = xEffR(p,x_a,y_a);
        xefft = xEffT(p,x_j,y_j);
        
        windup[i] = 1/(2*pi)*atan.(dot(p,cross(xefft,xeffr)),dot(xefft,xeffr))
    end#for
    return windup
end#function
#-----------------------------------------------------

function xEffR(p,x,y)
 # find effective x unit vector of rx antenna
 # p: LOS unit vector between Tx and Rx
 # (x,y): unit vector antenna pair, their cross product is center of antenna's FOV
 xEff = x-p * dot(p,x) + cross(p,y)
 return xEff
end#function
#-----------------------------------------------------

function xEffT(p,x,y)
 # find effective x unit vector of rx antenna
 # p: LOS unit vector between Tx and Rx
 # (x,y): unit vector antenna pair, their cross product is center of antenna's FOV
 xEff = x-p*dot(p,x)-cross(p,y)
 return xEff
end#function




end#module