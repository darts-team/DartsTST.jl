module Process_Raw_Data

c=299792458 # speed of light (m/s)

function main(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc) # no RSF
    Ns=size(s_xyz_grid)[2] # number of pixels in the scene
    Np=size(p_xyz_grid)[2] # number of platforms
    processed_image=zeros(ComplexF64,Ns) # intensity image vector
    λ=c/fc # wavelength (m)
    for j=1:Ns # for each pixel
        if mode==2;range_tx=distance(s_xyz_grid[:,j],p_xyz_grid[:,tx_el]);end
        for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
            range_rx=distance(s_xyz_grid[:,j],p_xyz_grid[:,i])
            if mode==1 # SAR (ping-pong)
                range_tx=range_rx
                processed_image[j]=processed_image[j]+rawdata[i]*exp(im*4*pi/λ*range_tx)
            elseif mode==2 # SIMO
                processed_image[j]=processed_image[j]+rawdata[i]*exp(im*2*pi/λ*(range_tx+range_rx))
            elseif mode==3 # MIMO
                for k=1:Np # TX platform
                    range_tx=distance(s_xyz_grid[:,j],p_xyz_grid[:,k])
                    processed_image[j]=processed_image[j]+rawdata[i,k]*exp(im*2*pi/λ*(range_tx+range_rx))
                end
            end
        end
    end
    return abs.(processed_image) # square for power?
end

function main_RSF(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc,t_rx,ref_range) # with RSF
    Ns=size(s_xyz_grid)[2] # number of pixels in the scene
    Np=size(p_xyz_grid)[2] # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Δt=t_rx[2]-t_rx[1]
    processed_image=zeros(ComplexF64,Ns) # intensity image vector
    λ=c/fc # wavelength (m)
    ref_delay=2*ref_range/c # reference delay
    for j=1:Ns # for each pixel
        if mode==2;range_tx=distance(s_xyz_grid[:,j],p_xyz_grid[:,tx_el]);end
        for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
            range_rx=distance(s_xyz_grid[:,j],p_xyz_grid[:,i])
            if mode==1 # SAR (ping-pong)
                range_tx=range_rx
                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                rel_delay_ind=Int(round(rel_delay/Δt))
                processed_image[j]=processed_image[j]+rawdata[i,Int(round(Nft/2))+rel_delay_ind]*exp(im*4*pi/λ*range_tx)
            elseif mode==2 # SIMO
                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                rel_delay_ind=Int(round(rel_delay/Δt))
                processed_image[j]=processed_image[j]+rawdata[i,Int(round(Nft/2))+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
            elseif mode==3 # MIMO
                for k=1:Np # TX platform
                    range_tx=distance(s_xyz_grid[:,j],p_xyz_grid[:,k])
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt))
                    processed_image[j]=processed_image[j]+rawdata[i,k,Int(round(Nft/2))+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                end
            end
        end
    end
    return abs.(processed_image) # square for power?
end

function main_RSF_slowtime(rawdata,s_xyz_grid,p_xyz_3D,mode,tx_el,fc,t_rx,ref_range) # with RSF and slow-time
    Ns=size(s_xyz_grid)[2] # number of pixels in the scene
    Np=size(p_xyz_3D)[2] # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    Δt=t_rx[2]-t_rx[1]
    processed_image=zeros(ComplexF64,Ns) # intensity image vector
    λ=c/fc # wavelength (m)
    ref_delay=2*ref_range/c # reference delay
    for j=1:Ns # for each pixel
        for s=1:Nst # slow-time (pulses)
            if mode==2;range_tx=distance(s_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);end
            for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                range_rx=distance(s_xyz_grid[:,j],p_xyz_3D[:,i,s])
                if mode==1 # SAR (ping-pong)
                    range_tx=range_rx
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt))
                    processed_image[j]=processed_image[j]+rawdata[s,i,Int(round(Nft/2))+rel_delay_ind]*exp(im*4*pi/λ*range_tx)
                elseif mode==2 # SIMO
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt))
                    processed_image[j]=processed_image[j]+rawdata[s,i,Int(round(Nft/2))+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                elseif mode==3 # MIMO
                    for k=1:Np # TX platform
                        range_tx=distance(s_xyz_grid[:,j],p_xyz_3D[:,k,s])
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=Int(round(rel_delay/Δt))
                        processed_image[j]=processed_image[j]+rawdata[s,i,k,Int(round(Nft/2))+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                    end
                end
            end
        end
    end
    return abs.(processed_image) # square for power?
end

function main_noRSF_slowtime(rawdata,s_xyz_grid,p_xyz_3D,mode,tx_el,fc) # without RSF and with slow-time
    Ns=size(s_xyz_grid)[2] # number of pixels in the scene
    Np=size(p_xyz_3D)[2] # number of platforms
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    processed_image=zeros(ComplexF64,Ns) # intensity image vector
    λ=c/fc # wavelength (m)
    for j=1:Ns # for each pixel
        for s=1:Nst # slow-time (pulses)
            if mode==2;range_tx=distance(s_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);end
            for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                range_rx=distance(s_xyz_grid[:,j],p_xyz_3D[:,i,s])
                if mode==1 # SAR (ping-pong)
                    range_tx=range_rx
                    processed_image[j]=processed_image[j]+rawdata[s,i]*exp(im*4*pi/λ*range_tx)
                elseif mode==2 # SIMO
                    processed_image[j]=processed_image[j]+rawdata[s,i]*exp(im*2*pi/λ*(range_tx+range_rx))
                elseif mode==3 # MIMO
                    for k=1:Np # TX platform
                        range_tx=distance(s_xyz_grid[:,j],p_xyz_3D[:,k,s])
                        processed_image[j]=processed_image[j]+rawdata[s,i,k]*exp(im*2*pi/λ*(range_tx+range_rx))
                    end
                end
            end
        end
    end
    return abs.(processed_image) # square for power?
end

function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
end

end
