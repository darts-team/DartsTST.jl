module Generate_Raw_Data
using Random

c=299792458 # speed of light (m/s)

function main(t_xyz_grid,p_xyz_grid,mode,tx_el,fc) # no RSF
    λ=c/fc # wavelength (m)
    Nt=size(t_xyz_grid)[2] # number of targets
    Np=size(p_xyz_grid)[2] # number of platforms
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata=zeros(ComplexF64,Np)
    elseif mode==3
        rawdata=zeros(ComplexF64,Np,Np)
    end
    for j=1:Nt # targets
        if mode==2;range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,tx_el]);end
        for i=1:Np # RX platform
            range_rx=distance(t_xyz_grid[:,j],p_xyz_grid[:,i])
            if mode==1 # SAR (ping-pong)
                range_tx=range_rx
                rawdata[i]=rawdata[i]+exp(-im*4*pi/λ*range_tx)
            elseif mode==2 # SIMO
                rawdata[i]=rawdata[i]+exp(-im*2*pi/λ*(range_tx+range_rx))
            elseif mode==3 # MIMO
                for k=1:Np # TX platform for MIMO
                    range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,k])
                    rawdata[i,k]=rawdata[i,k]+exp(-im*2*pi/λ*(range_tx+range_rx))
                end
            end
        end
    end
    return rawdata
end

function main_RSF(t_xyz_grid,p_xyz_grid,mode,tx_el,fc,Srx,t_rx,ref_range) # with RSF
    # TODO add descriptions of inputs and output
    λ=c/fc # wavelength (m)
    Nt=size(t_xyz_grid)[2] # number of targets
    Np=size(p_xyz_grid)[2] # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Δt=t_rx[2]-t_rx[1] # fast-time resolution
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata=zeros(ComplexF64,Np,Nft)
    elseif mode==3
        rawdata=zeros(ComplexF64,Np,Np,Nft)
    end
    ref_delay=2*ref_range/c # reference delay
    for j=1:Nt # targets
        if mode==2;range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,tx_el]);end
        for i=1:Np # RX platform
            range_rx=distance(t_xyz_grid[:,j],p_xyz_grid[:,i])
            if mode==1 # SAR (ping-pong)
                range_tx=range_rx
                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                rel_delay_ind=Int(round(rel_delay/Δt))
                if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                    Srx_shifted=[zeros(1,rel_delay_ind) Srx[1:Nft-rel_delay_ind]']
                elseif rel_delay_ind<0
                    Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft]' zeros(1,abs(rel_delay_ind))]
                end
                rawdata[i,:]=rawdata[i,:]'+exp(-im*4*pi/λ*range_tx)*Srx_shifted
            elseif mode==2 # SIMO
                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                rel_delay_ind=Int(round(rel_delay/Δt))
                if rel_delay_ind>=0
                    Srx_shifted=[zeros(1,rel_delay_ind) Srx[1:Nft-rel_delay_ind]']
                elseif rel_delay_ind<0
                    Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft]' zeros(1,abs(rel_delay_ind))]
                end
                rawdata[i,:]=rawdata[i,:]'+exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
            elseif mode==3 # MIMO
                for k=1:Np # TX platform for MIMO
                    range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,k])
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt))
                    if rel_delay_ind>=0
                        Srx_shifted=[zeros(1,rel_delay_ind) Srx[1:Nft-rel_delay_ind]']
                    elseif rel_delay_ind<0
                        Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft]' zeros(1,abs(rel_delay_ind))]
                    end
                    rawdata[i,k,:]=rawdata[i,k,:]'+exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
                end
            end
        end
    end
    return rawdata
end

function main_RSF_slowtime(t_xyz_grid,p_xyz_3D,mode,tx_el,fc,Srx,t_rx,ref_range,t_ref) # with RSF and slow-time
    # TODO add descriptions of inputs and output
    λ=c/fc # wavelength (m)
    Nt=size(t_xyz_grid,2) # number of targets
    Np=size(p_xyz_3D,2) # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D,3) # number of slow-time samples
    Δt_ft=t_rx[2]-t_rx[1] # fast-time resolution
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata=zeros(ComplexF64,Nst,Np,Nft)
    elseif mode==3
        rawdata=zeros(ComplexF64,Nst,Np,Np,Nft)
    end
    ref_delay=2*ref_range/c # reference delay
    for j=1:Nt # targets
        for s=1:Nst # slow-time (pulses)
            if mode==2;range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);end
            for i=1:Np # RX platform
                range_rx=distance(t_xyz_grid[:,j],p_xyz_3D[:,i,s])
                if mode==1 # SAR (ping-pong)
                    range_tx=range_rx
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt_ft))
                    if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                        Srx_shifted=cat(zeros(rel_delay_ind),Srx[1:Nft-rel_delay_ind],dims=1)
                    elseif rel_delay_ind<0
                        Srx_shifted=cat(Srx[1+abs(rel_delay_ind):Nft],zeros(abs(rel_delay_ind)),dims=1)
                    end
                    rawdata[s,i,:]=rawdata[s,i,:]+t_ref[j]*exp(-im*4*pi/λ*range_tx)*Srx_shifted
                elseif mode==2 # SIMO
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt_ft))
                    if rel_delay_ind>=0
                        Srx_shifted=[zeros(rel_delay_ind);Srx[1:Nft-rel_delay_ind]]
                    elseif rel_delay_ind<0
                        Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft];zeros(abs(rel_delay_ind))]
                    end
                    rawdata[s,i,:]=rawdata[s,i,:]+t_ref[j]*exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
                elseif mode==3 # MIMO
                    for k=1:Np # TX platform for MIMO
                        range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,k,s])
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=Int(round(rel_delay/Δt_ft))
                        if rel_delay_ind>=0
                            Srx_shifted=[zeros(rel_delay_ind);Srx[1:Nft-rel_delay_ind]]
                        elseif rel_delay_ind<0
                            Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft];zeros(abs(rel_delay_ind))]
                        end
                        rawdata[s,i,k,:]=rawdata[s,i,k,:]+t_ref[j]*exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
                    end
                end
            end
        end
    end
    return rawdata
end

function main_noRSF_slowtime(t_xyz_grid,p_xyz_3D,mode,tx_el,fc,t_ref) # without fast-time and slow-time
    # TODO add descriptions of inputs and output
    λ=c/fc # wavelength (m)
    Nt=size(t_xyz_grid)[2] # number of targets
    Np=size(p_xyz_3D)[2] # number of platforms
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata=zeros(ComplexF64,Nst,Np)
    elseif mode==3
        rawdata=zeros(ComplexF64,Nst,Np,Np)
    end
    for j=1:Nt # targets
        for s=1:Nst # slow-time (pulses)
            if mode==2;range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);end
            for i=1:Np # RX platform
                range_rx=distance(t_xyz_grid[:,j],p_xyz_3D[:,i,s])
                if mode==1 # SAR (ping-pong)
                    range_tx=range_rx
                    rawdata[s,i]=rawdata[s,i]+t_ref[j]*exp(-im*4*pi/λ*range_tx)
                elseif mode==2 # SIMO
                    rawdata[s,i]=rawdata[s,i]+t_ref[j]*exp(-im*2*pi/λ*(range_tx+range_rx))
                elseif mode==3 # MIMO
                    for k=1:Np # TX platform for MIMO
                        range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,k,s])
                        rawdata[s,i,k]=rawdata[s,i,k]+t_ref[j]*exp(-im*2*pi/λ*(range_tx+range_rx))
                    end
                end
            end
        end
    end
    return rawdata
end

function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
end

end
