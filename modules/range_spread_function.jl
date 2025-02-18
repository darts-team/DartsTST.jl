module RSF # Range Spread Function (matched filter output)

using Parameters

function ideal_RSF(Trx, params) # ideal RSF for LFM pulse, no windowing
    @unpack pulse_length,Δt,bandwidth = params
    
    #ft=-pulse_length/2:Δt:pulse_length/2 # fast-time axis, RSF is zero outside +/- τ #Test
    #MF=(bandwidth*pulse_length)*ones.(length(ft)) # equation for LFM chirp matched filter output, 0° constant phase #Test
    ft=-pulse_length:Δt:pulse_length # fast-time axis, RSF is zero outside +/- τ
    MF=(bandwidth*pulse_length)*sinc.((1.0.-abs.(ft)/pulse_length).*ft*bandwidth).*(1.0.-abs.(ft)/pulse_length) # equation for LFM chirp matched filter output, 0° constant phase
    MF[findall(x->0<=x.<0.0001,MF)].=0.0001
    MF[findall(x->-0.0001<x.<0,MF)].=-0.0001
    t_rx=-Trx/2:Δt:Trx/2 # RX window
    Nft=length(t_rx) # number of fast-time samples
    Srx=zeros(Nft) # RX window
    #start_index=Int(round(Nft/2))-Int(round(pulse_length/2/Δt)) #Test
    start_index=Int(round(Nft/2))-Int(round(pulse_length/Δt))
    Srx[start_index:start_index+length(ft)-1]=MF # MF centered on RX window
    return Srx,MF,ft,t_rx
end

function ideal_RSF_pulsewindow(Trx, params) # ideal RSF for LFM pulse, no windowing
    @unpack pulse_length,Δt,bandwidth = params
    
    ft=-Trx/2:Δt:Trx/2 # fast-time axis, RSF is zero outside +/- τ
    MF=(bandwidth*pulse_length)*sinc.((1.0.-abs.(ft)/pulse_length).*ft*bandwidth).*(1.0.-abs.(ft)/pulse_length) # equation for LFM chirp matched filter output, 0° constant phase
    MF[findall(x->0<=x.<0.0001,MF)].=0.0001
    MF[findall(x->-0.0001<x.<0,MF)].=-0.0001
    return MF,ft
end


function non_ideal_RSF(τ,Δt,B,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing
end

end
