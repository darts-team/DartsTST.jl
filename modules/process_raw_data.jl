module Process_Raw_Data

using Parameters

const c=299792458 # speed of light (m/s)

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

function main_RSF_slowtime(rawdata,s_xyz_grid,p_xyz_3D, params, t_rx, ref_range) # with RSF and slow-time
    @unpack mode, tx_el, λ = params
    Ns=size(s_xyz_grid)[2] # number of pixels in the scene
    Np=size(p_xyz_3D)[2] # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    Δt=t_rx[2]-t_rx[1]
    processed_image=zeros(ComplexF64,Ns) # intensity image vector
    ref_delay=2*ref_range/c # reference delay

    if mode==1 # SAR (ping-pong)
        for j=1:Ns # for each pixel
            pixel_j = @view(s_xyz_grid[:,j])
            pixel_sum = 0.0im
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
            pixel_sum = 0.0im
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
            pixel_sum = 0.0im
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

function main_SAR_tomo_3D(rawdata,s_xyz_grid,p_xyz_3D,t_rx, ref_range, params) # with fast-time, slow-time, and tomographic processing; pixels in 3D
    @unpack Ns_1, Ns_2, Ns_3, mode, tx_el, fc = params

    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    s_xyz_3D=reshape(s_xyz_grid,3,Ns_3,Ns_2,Ns_1) # convert scene to 3D
    Np=size(p_xyz_3D)[2] # number of platforms
    Δt=t_rx[2]-t_rx[1] # fast-time resolution
    processed_image=zeros(ComplexF64,Ns_1,Ns_2,Ns_3) # 3D image array
    λ=c/fc # wavelength (m)
    ref_delay=2*ref_range/c # reference delay
    rel_delay_ind::Int = 0

    if mode==1 # SAR (ping-pong)
        for j1=1:Ns_1 # for each pixel in axis-1
            for j2=1:Ns_2 # for each pixel in axis-2
                for j3=1:Ns_3 # for each pixel in axis-3
                    pixel_j = @view(s_xyz_3D[:,j3,j2,j1])
                    pixel_sum = 0.0im
                    for i=1:Np # TX or RX platform
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
                    pixel_sum = 0.0im
                    for s=1:Nst # slow-time (pulses)
                        range_tx=distance(pixel_j,@view(p_xyz_3D[:,tx_el,s]))
                        for i=1:Np # RX platform
                            range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                            rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind=round(Int,rel_delay/Δt)
                            pixel_sum += rawdata[s,i,rel_delay_ind + round(Int,Nft/2)]*exp(im*2*pi/λ*(range_tx+range_rx))
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
                    pixel_sum = 0.0im
                    for s=1:Nst # slow-time (pulses)
                        for i=1:Np # RX platform
                            range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                            for k=1:Np # TX platform
                                range_tx=distance(pixel_j,@view(p_xyz_3D[:,k,s]))
                                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                                rel_delay_ind=round(Int,rel_delay/Δt)
                                pixel_sum += rawdata[s,i,k,round(Int,Nft/2)+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                            end
                        end
                    end
                    processed_image[j1,j2,j3] = pixel_sum
                end
            end
        end
    end
    return abs.(processed_image)
end

function main_SAR_tomo_3D_new(rawdata,s_xyz_grid,p_xyz_3D,t_rx, ref_range, params) # with fast-time, slow-time, and tomographic processing; pixels in 3D
    @unpack Ns_1, Ns_2, Ns_3, mode, tx_el, fc = params

    Nft         = length(t_rx) # number of fast-time samples
    s_xyz_3D    = reshape(s_xyz_grid,3,Ns_3,Ns_2,Ns_1) # convert scene to 3D
    Δt          = t_rx[2]-t_rx[1] # fast-time resolution
    processed_image = zeros(ComplexF64,Ns_3,Ns_2,Ns_1) # 3D image array
    λ           = c/fc # wavelength (m)
    ref_delay   = 2*ref_range/c # reference delay
    CI          = CartesianIndices(processed_image) 

    if mode==1 # SAR (ping-pong)
        for j1=1:length(CI)
            pixel_j         = @view(s_xyz_3D[:,CI[j1]])
            range_rx        = permutedims( sum((pixel_j.-p_xyz_3D).^2,dims=1).^0.5, [3,2,1]) 
            range_tx        = range_rx
            rel_delay       = (range_tx.+range_rx)./c.-ref_delay
            rel_delay_ind   = round.(Int,rel_delay./Δt)
            index3          = round(Int,Nft/2).+rel_delay_ind
            CI2             = CartesianIndices(index3[:,:,1])
            index_all       = CartesianIndex.(getindex.(CI2,1)[:],getindex.(CI2,2)[:],index3[:])
            RanF            = exp.(im.*4. *pi./λ.*range_tx)
            processed_image[CI[j1]] = sum(rawdata[index_all]  .* RanF[CI2[:]])
        end
        processed_image = permutedims(processed_image,[3,2,1])
    elseif mode==2 # SIMO
        for j1=1:length(CI)
            pixel_j         = @view(s_xyz_3D[:,CI[j1]])
            range_rx        = permutedims( sum((pixel_j.-p_xyz_3D).^2,dims=1).^0.5, [3,2,1]) 
            range_tx        = permutedims(repeat(range_rx2[:,tx_el,:],outer = [10,1,1]), [3,1,2])
            rel_delay       = (range_tx.+range_rx)./c.-ref_delay
            rel_delay_ind   = round.(Int,rel_delay./Δt)
            index3          = round(Int,Nft/2).+rel_delay_ind
            CI2             = CartesianIndices(index3[:,:,1])
            index_all       = CartesianIndex.(getindex.(CI2,1)[:],getindex.(CI2,2)[:],index3[:])
            RanF            = exp.(im.*4. *pi./λ.*range_tx)
            processed_image[CI[j1]] = sum(rawdata[index_all]  .* RanF[CI2[:]])
        end
        processed_image = permutedims(processed_image,[3,2,1])
    elseif mode==3 # MIMO
        for j1=1:Ns_1 # for each pixel in axis-1
            for j2=1:Ns_2 # for each pixel in axis-2
                for j3=1:Ns_3 # for each pixel in axis-3
                    pixel_j = @view(s_xyz_3D[:,j3,j2,j1])
                    pixel_sum = 0.0im
                    for s=1:Nst # slow-time (pulses)
                        for i=1:Np # RX platform
                            range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                            for k=1:Np # TX platform
                                range_tx=distance(pixel_j,@view(p_xyz_3D[:,k,s]))
                                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                                rel_delay_ind=round(Int,rel_delay/Δt)
                                pixel_sum += rawdata[s,i,k,round(Int,Nft/2)+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                            end
                        end
                    end
                    processed_image[j1,j2,j3] = pixel_sum
                end
            end
        end
    end
    return abs.(processed_image)
end

function SAR_processing(rawdata, s_xyz_grid, p_xyz_3D, t_rx, ref_range, params) # slow-time processing of rawdata with fast-time
    @unpack Ns_1, Ns_2, Ns_3, mode, tx_el, fc, λ = params

    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    s_xyz_3D=reshape(s_xyz_grid,3,Ns_3,Ns_2,Ns_1) # convert scene to 3D
    #Np=size(p_xyz_3D)[2] # number of platforms
    Np=size(rawdata)[2] # number of platforms
    Δt=t_rx[2]-t_rx[1] # fast-time resolution
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        SAR_images_3D=zeros(ComplexF64,Np,Ns_1,Ns_2,Ns_3) # complex SAR images array (4D)
    elseif mode==3
        SAR_images_3D=zeros(ComplexF64,Np,Np,Ns_1,Ns_2,Ns_3) # complex SAR images array (5D)
    end
    ref_delay::Float64=2*ref_range/c # reference delay
    if mode==1 # SAR (ping-pong)
        for i=1:Np # TX or RX platform
            for j1=1:Ns_1 # for each pixel in axis-1
                for j2=1:Ns_2 # for each pixel in axis-2
                    for j3=1:Ns_3 # for each pixel in axis-3
                        pixel_j = @view(s_xyz_3D[:,j3,j2,j1])
                        pixel_sum = 0.0im
                        for s=1:Nst # slow-time (pulses)
                            range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                            range_tx=range_rx
                            rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind=round(Int,rel_delay/Δt)
                            pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(im*4*pi/λ*range_tx)
                        end
                        SAR_images_3D[i,j1,j2,j3] = pixel_sum # i : Tx or Rx platform
                    end
                end
            end
        end
    elseif mode==2 # SIMO
        for i=1:Np # RX platform
            for j1=1:Ns_1 # for each pixel in axis-1
                for j2=1:Ns_2 # for each pixel in axis-2
                    for j3=1:Ns_3 # for each pixel in axis-3
                        pixel_j = @view(s_xyz_3D[:,j3,j2,j1])
                        pixel_sum = 0.0im
                        for s=1:Nst # slow-time (pulses)
                            range_tx=distance(pixel_j,@view(p_xyz_3D[:,tx_el,s]))
                            range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                            rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind=round(Int,rel_delay/Δt)
                            pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                        end
                        SAR_images_3D[i,j1,j2,j3] = pixel_sum # i : Rx platform
                    end
                end
            end
        end
    elseif mode==3 # MIMO
        for i=1:Np # RX platform
            for k=1:Np # TX platform
                for j1=1:Ns_1 # for each pixel in axis-1
                    for j2=1:Ns_2 # for each pixel in axis-2
                        for j3=1:Ns_3 # for each pixel in axis-3
                            pixel_j = @view(s_xyz_3D[:,j3,j2,j1])
                            pixel_sum = 0.0im
                            for s=1:Nst # slow-time (pulses)
                                range_tx=distance(pixel_j,@view(p_xyz_3D[:,k,s]))
                                range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                                rel_delay_ind=round(Int,rel_delay/Δt)
                                pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(im*2*pi/λ*(range_tx+range_rx))
                            end
                            SAR_images_3D[i,k,j1,j2,j3] = pixel_sum # i : Rx platform, k: Tx platform
                        end
                    end
                end
            end
        end
    end
    return SAR_images_3D
end

function tomo_processing_afterSAR(SAR_images_3D) # tomographic processing of slow-time processed data
    if ndims(SAR_images_3D)==4
        image_3D=sum(SAR_images_3D,dims=1)
        Ns_1=size(SAR_images_3D,2);Ns_2=size(SAR_images_3D,3);Ns_3=size(SAR_images_3D,4)
        image_3D=reshape(image_3D,Ns_1,Ns_2,Ns_3)
    elseif ndims(SAR_images_3D)==5
        image_3D=sum(sum(SAR_images_3D,dims=1),dims=2)
        Ns_1=size(SAR_images_3D,3);Ns_2=size(SAR_images_3D,4);Ns_3=size(SAR_images_3D,5)
        image_3D=reshape(image_3D,Ns_1,Ns_2,Ns_3)
    end
    return abs.(image_3D)
end

function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
end

end
