module Generate_Raw_Data

c=299792458 # speed of light (m/s)

function main(t_xyz_grid,p_xyz_grid,mode,tx_el,fc,a,e)

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
            for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                range_rx=distance(t_xyz_grid[:,j],p_xyz_grid[:,i])
                if mode==1 # SAR (ping-pong)
                    range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,i])
                    rawdata[i]=rawdata[i]+exp(-im*4*pi/λ*range_tx)
                elseif mode==2 # SIMO
                    rawdata[i]=rawdata[i]+exp(-im*2*pi/λ*(range_tx+range_rx))
                elseif mode==3 # MIMO
                    for k=1:Np # TX platform
                        range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,k])
                        rawdata[i,k]=rawdata[i,k]+exp(-im*2*pi/λ*(range_tx+range_rx))
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
