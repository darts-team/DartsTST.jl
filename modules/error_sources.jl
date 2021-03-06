module Error_Sources

function random_noise(rawdata,SNR,enable_fast_time,mode)
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

function synchronization_errors
end

function position_errors
end

end
