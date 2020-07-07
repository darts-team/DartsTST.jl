module Raw_Data

include("geometry.jl")
include("scene.jl")

c=299792458 # speed of light (m/s)

function main(t_θ,t_ϕ,t_h,p_θ,p_ϕ,p_h,mode,tx_el,fc,a,e)
    λ=c/fc # wavelength (m)
    Nt=length(t_θ)*length(t_ϕ)*length(t_h) # number of targets
    Np=length(p_θ)*length(p_ϕ)*length(p_h) # number of platforms
    t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops
    p_geo_grid=Scene.form3Dgrid_for(p_θ,p_ϕ,p_h) # using 3 nested for loops
    #t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing
    #p_geo_grid=Scene.form3Dgrid_array(p_θ,p_ϕ,p_h) # using array processing
    t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,a,e)
    p_xyz_grid=Geometry.geo_to_xyz(p_geo_grid,a,e)
    if mode==1 # SAR (ping-pong)
        rawdata=zeros(ComplexF64,Np)
        for i=1:Np # TX-RX platform
            for j=1:Nt
                range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,i])
                rawdata[i]=rawdata[i]+exp(-im*4*pi/λ*range_tx)
            end
        end
    elseif mode==2 # SIMO
        rawdata=zeros(ComplexF64,Np)
        for i=1:Np # RX platform
            for j=1:Nt
                range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,tx_el])
                range_rx=distance(t_xyz_grid[:,j],p_xyz_grid[:,i])
                rawdata[i]=rawdata[i]+exp(-im*2*pi/λ*(range_tx+range_rx))
            end
        end
    elseif mode==3 # MIMO
        rawdata=zeros(ComplexF64,Np,Np)
        for i=1:Np # TX platform
            for k=1:Np # RX platform
                for j=1:Nt
                    range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,i])
                    range_rx=distance(t_xyz_grid[:,j],p_xyz_grid[:,k])
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
