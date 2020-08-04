module Process_Raw_Data

c=299792458 # speed of light (m/s)

function main(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc,a,e)

    Ns=size(s_xyz_grid)[2] # number of pixels in the scene
    Np=size(p_xyz_grid)[2] # number of platforms
    processed_image=zeros(ComplexF64,Ns) # intensity image vector
    λ=c/fc # wavelength (m)

        for j=1:Ns # for each pixel
            if mode==2;range_tx=distance(s_xyz_grid[:,j],p_xyz_grid[:,tx_el]);end
            for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                range_rx=distance(s_xyz_grid[:,j],p_xyz_grid[:,i])
                if mode==1 # SAR (ping-pong)
                    range_tx=distance(s_xyz_grid[:,j],p_xyz_grid[:,i])
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

    return processed_image
end

function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
end

end
