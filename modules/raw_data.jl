module Raw_Data

function main(t_xyz_grid,pb,fc,h,α_b,mode,tx_element)
    λ=c/fc # wavelength (m)
    Nt=size(t_xyz_grid,2) # number of targets
    Np=size(pb) # number of platforms
    ranges=calculate_ranges(t_xyz_grid,pb,h,α_b,mode,tx_element)
    rawdata=add_phase_vectors(ranges,fc)
end

function calculate_ranges(t_xyz_grid,pb,h,α_b,mode,tx_element)
    p_xyz=
    for i=1:Np
        ranges=sqrsum???(t_xyz_grid-p_xyz(i))
    end
end

function add_phase_vectors(ranges,fc)
    if mode==1 # SAR (ping-pong)
        for i=1:Np
            rawdata(i)=sum(exp(-im*4*pi/λ*ranges(i,:)))
        end
    elseif mode==2 # SIMO
    elseif mode==3 # MIMO
    end
end

end
