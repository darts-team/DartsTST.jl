module Process_Raw_Data

c=299792458 # speed of light (m/s)

function main_RSF(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc,t_rx,ref_range) # with fast-time and tomographic processing, no slowtime
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

function main_RSF_slowtime(rawdata,s_xyz_grid,p_xyz_3D,mode,tx_el,fc,t_rx,ref_range) # with fast-time, slow-time, and tomographic processing
    Ns=size(s_xyz_grid)[2] # number of pixels in the scene
    Np=size(p_xyz_3D)[2] # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    Δt=t_rx[2]-t_rx[1]
    processed_image=zeros(ComplexF64,Ns) # intensity image vector
    λ=c/fc # wavelength (m)
    ref_delay=2*ref_range/c # reference delay

    if mode==1 # SAR (ping-pong)
        for j=1:Ns # for each pixel
            pixel_j = @view(s_xyz_grid[:,j])
            pixel_sum = 0.0im;
            for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                for s=1:Nst # slow-time (pulses)
                    range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                    range_tx=range_rx
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=round(Int,rel_delay/Δt)
                    pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(im*4*pi/λ*range_tx)
                end
            end
            processed_image[j] = pixel_sum
        end
    elseif mode==2 # SIMO
        for j=1:Ns # for each pixel
            pixel_j = @view(s_xyz_grid[:,j])
            pixel_sum = 0.0im;
            for s=1:Nst # slow-time (pulses)
                range_tx=distance(pixel_j,@view(p_xyz_3D[:,tx_el,s]))
                for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                    range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=round(Int,rel_delay/Δt)
                    pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                end
            end
            processed_image[j] = pixel_sum
        end
    elseif mode==3 # MIMO
        for j=1:Ns # for each pixel
            pixel_j = @view(s_xyz_grid[:,j])
            pixel_sum = 0.0im;
            for s=1:Nst # slow-time (pulses)
                for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                    range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                    for k=1:Np # TX platform
                        range_tx=distance(pixel_j,@view(p_xyz_3D[:,k,s]))
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=round(Int,rel_delay/Δt)
                        pixel_sum=pixel_sum+rawdata[s,i,k,round(Int,Nft/2)+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                    end
                end
            end
            processed_image[j] = pixel_sum
        end
    end
    return abs.(processed_image) # square for power?
end

function main_SAR_tomo_3D(rawdata,s_xyz_grid,Ns_1,Ns_2,Ns_3,p_xyz_3D,mode,tx_el,fc,t_rx,ref_range) # with fast-time, slow-time, and tomographic processing; pixels in 3D
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    s_xyz_3D=reshape(s_xyz_grid,3,Ns_3,Ns_2,Ns_1) # convert scene to 3D
    #s_xyz_3D=Scene.convert_scene_3xN_to_3D(s_xyz_grid,Ns_1,Ns_2,Ns_3) # convert scene to 3D
    Np=size(p_xyz_3D)[2] # number of platforms
    Δt=t_rx[2]-t_rx[1]
    processed_image=zeros(ComplexF64,Ns_1,Ns_2,Ns_3) # intensity image vector
    λ=c/fc # wavelength (m)
    ref_delay=2*ref_range/c # reference delay

    if mode==1 # SAR (ping-pong)
        for j1=1:Ns_1 # for each pixel in axis-1
            for j2=1:Ns_2 # for each pixel in axis-2
                for j3=1:Ns_3 # for each pixel in axis-3
                    pixel_j = @view(s_xyz_3D[:,j3,j2,j1])
                    pixel_sum = 0.0im;
                    for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                        for s=1:Nst # slow-time (pulses)
                            range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                            range_tx=range_rx
                            rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind=round(Int,rel_delay/Δt)
                            pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(im*4*pi/λ*range_tx)
                        end
                    end
                    processed_image[j1,j2,j3] = pixel_sum
                end
            end
        end
    elseif mode==2 # SIMO
        for j1=1:Ns_1 # for each pixel in axis-1
            for j2=1:Ns_2 # for each pixel in axis-2
                for j3=1:Ns_3 # for each pixel in axis-3
                    pixel_j = @view(s_xyz_3D[:,j3,j2,j1])
                    pixel_sum = 0.0im;
                    for s=1:Nst # slow-time (pulses)
                        range_tx=distance(pixel_j,@view(p_xyz_3D[:,tx_el,s]))
                        for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                            range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                            rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind=round(Int,rel_delay/Δt)
                            pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                        end
                    end
                    processed_image[j1,j2,j3] = pixel_sum
                end
            end
        end
    elseif mode==3 # MIMO
        for j1=1:Ns_1 # for each pixel in axis-1
            for j2=1:Ns_2 # for each pixel in axis-2
                for j3=1:Ns_3 # for each pixel in axis-3
                    pixel_j = @view(s_xyz_3D[:,j3,j2,j1])
                    pixel_sum = 0.0im;
                    for s=1:Nst # slow-time (pulses)
                        for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                            range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                            for k=1:Np # TX platform
                                range_tx=distance(pixel_j,@view(p_xyz_3D[:,k,s]))
                                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                                rel_delay_ind=round(Int,rel_delay/Δt)
                                pixel_sum=pixel_sum+rawdata[s,i,k,round(Int,Nft/2)+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                            end
                        end
                    end
                    processed_image[j1,j2,j3] = pixel_sum
                end
            end
        end
    end
    return abs.(processed_image) # square for power?
end

function SAR_processing(rawdata,s_xyz_grid,p_xyz_3D,mode,tx_el,fc,t_rx,ref_range) # slow-time processing of rawdata with fast-time

    return SAR_images_3D
end

function tomo_processing_afterSAR(SAR_images_3D,s_xyz_grid,p_xyz_3D,mode,tx_el,fc,t_rx,ref_range) # tomographic processing of slow-time processed data

    return image_3D
end

function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
end

end
