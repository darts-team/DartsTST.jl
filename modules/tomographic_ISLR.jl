module tomographic_ISLR

"Inputs:
    mode: tomographic mode, 1: SAR, 2: SIMO, 3: MIMO
    plat_distr: 1D row array specifying locations of platforms projected perpendicular to look angle direction relative to the reference platform location (m)
    look_angle: perpendicular to along-track and with respect to nadir (deg)
    local_slope: local slope of terrain (deg)
    altitude: of the reference platform (m)
    max_veg_H: maximum vegetation height (m)
    freq: radar center frequency (Hz)
    scene_res: pixel resolution of the scene (m)
Output:
    PSF: tomographic point spread function (dB) (peak is normalized to 0 dB)
    ISLR: tomographic integrated side-lobe ratio along perpendicular to look angle direction inside the scene (dB)
    scene_axis: scene axis to display PSF (m) (depends on scene_size and scene_res)
Notes:
    ISLR in linear (smaller the better) (is usually negative and more negative is better)
    scene_size along perpendicular to look angle direction is calculated from max_veg_H, look_angle, local_slope
    larger max_veg_H and smaller difference between look_angle and local_slope means larger scene_size (resulting in worse ISLR)
    spherical earth is assumed
    point target is assumed to be where the look angle vector of reference platform at plat_distr=0 intersects the scene in perpendicular to look angle direction"
function main(mode,plat_distr,look_angle,local_slope,altitude,max_veg_H,freq,scene_res)
    c=299792458 # speed of light (m/s)
    earth_radius=6378137 # earth radius at equator (m)
    wl=c/freq # radar center wavelength (m)
    #slant_range=altitude/cosd(look_angle) # flat earth
    slant_range,ground_range=lookangle_to_range(look_angle,altitude,0,earth_radius) # spherical earth
    xt=0 # target location at the scene in perpendicular to look angle direction (0 means target is where the look angle vector of reference platform at 0 intersects the scene in perpendicular to look angle direction)
    N=length(plat_distr) # number of platforms
    scene_size=max_veg_H*cosd(local_slope)/sind(look_angle-local_slope) # scene extend along perpendicular to look angle direction is calculated from max_veg_H and look_angle (m)
    scene_axis=-scene_size/2:scene_res:scene_size/2 # scene size and resolution
    K=length(scene_axis) # number of pixels in the scene
    PSF=zeros(K,1) # tomographic point spread function
    if mode==1 # SAR (ping-pong)
        Pr=zeros(ComplexF64,N,1)
        for i=1:N
           Ri=((xt-plat_distr[i]).^2+slant_range^2).^0.5
           Pr[i]=sum(exp(-1im*4*pi/wl*Ri))
        end
        for i=1:K
            x=scene_axis[i]
            Ri=((x.-plat_distr).^2 .+ slant_range^2).^0.5
            PSF[i]=abs(sum(Pr.*exp.(1im*4*pi/wl*Ri)))
        end
    elseif mode==2 # SIMO (assuming first platform is TX)
        Pr=zeros(ComplexF64,N,1)
        R1=((xt-plat_distr[1]).^2+slant_range^2).^0.5
        for i=1:N
            R2=((xt-plat_distr[i]).^2+slant_range^2).^0.5
            Pr[i]=sum(exp(-1im*2*pi/wl*(R1+R2)))
        end
        for i=1:K
            x=scene_axis[i]
            R1=((x-plat_distr[1]).^2+slant_range^2).^0.5
            R2=((x.-plat_distr).^2 .+ slant_range^2).^0.5
            PSF[i]=abs(sum(Pr.*exp.(1im*2*pi/wl*(R1.+R2))))
        end
    elseif mode==3 # MIMO
        Pr=zeros(ComplexF64,N,N)
        for i=1:N
            R1=((xt-plat_distr[i]).^2+slant_range^2).^0.5
            for j=1:N
                R2=((xt-plat_distr[j]).^2+slant_range^2).^0.5
                Pr[i,j]=sum(exp(-1im*2*pi/wl*(R1+R2)))
            end
        end
        for i=1:K
            x=scene_axis[i]
            Pr1=zeros(ComplexF64,N,1)
            for j=1:N
                R1=((x-plat_distr[j]).^2+slant_range^2).^0.5
                R2=((x.-plat_distr).^2 .+ slant_range^2).^0.5
                Pr1[j]=sum(Pr[:,j].*exp.(1im*2*pi/wl*(R1.+R2)))
            end
            PSF[i]=abs(sum(Pr1))
        end
    end
    PSF=PSF/maximum(PSF)
    ISLR=calculateISLR(PSF)
    return PSF,ISLR,scene_axis
end

function calculateISLR(PSF)
    res_dB=4
    PSF_dB=20*log10.(PSF)
    #max_ind=findall(PSF_dB .==maximum(PSF_dB));max_ind=max_ind[1]
    max_value,max_ind=findmax(PSF);max_ind=max_ind[1]
    res_ind_2=findfirst(PSF_dB[max_ind].-PSF_dB[max_ind:end] .>=res_dB)+max_ind-1
    res_ind_1=findlast(PSF_dB[max_ind].-PSF_dB[1:max_ind] .>=res_dB)
    sidelobe_energy=sum(PSF[1:res_ind_1-1])+sum(PSF[res_ind_2+1:end])
    mainlobe_energy=sum(PSF[res_ind_1:res_ind_2])
    ISLR=20*log10(sidelobe_energy/mainlobe_energy)
    return ISLR
end

function lookangle_to_range(θ_l,p_h,t_h,ra)
    ra=ra.+t_h
    p_h=p_h.-t_h
    inc=asind.(sind.(θ_l).*(ra+p_h)./ra) # deg incidence angle
    α=inc-θ_l # deg planet-central angle
    rg=ra.*α*pi/180 # ground range
    rs=(((ra+p_h).*sind.(α)).^2+((ra.+p_h).*cosd.(α).-ra).^2).^0.5 # slant range
    return rs,rg
end

end
