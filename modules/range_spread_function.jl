module RSF # Range Spread Function (matched filter output)

function ideal_RSF(τ,Δt,B,Trx) # ideal RSF for LFM pulse, no windowing
    ft=-τ:Δt:τ # fast-time axis, RSF is zero outside +/- τ
    MF=(B*τ)*sinc.((1.0.-abs.(ft)/τ).*ft*B).*(1.0.-abs.(ft)/τ) # equation for LFM chirp matched filter output, 0° constant phase
    t_rx=-Trx/2:Δt:Trx/2 # RX window
    Nft=length(t_rx) # number of fast-time samples
    Srx=zeros(1,Nft) # RX window
    start_index=Int(round(Nft/2))-Int(round(τ/Δt))
    Srx[start_index:start_index+length(ft)-1]=MF # MF centered on RX window
    return Srx,MF,ft,t_rx
end

function non_ideal_RSF(τ,Δt,B,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing
end

end
