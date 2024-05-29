


include("../modules/geometry.jl")
include("../modules/scene.jl")
include("../modules/orbits.jl")
include("../modules/antenna.jl")
include("../modules/simsetup.jl")
include("../modules/global_scale_support.jl")

using Plots
using NCDatasets
using Dates
using Statistics

#orbit_dataset               = Dataset("/Users/joshil/Documents/Orbits/From_Eric/cartwheel/ROSEL_cartwheel_orbits 2/Cartwheel 4sat/ROSE_L_Cartwheel_4sat.nc")
orbit_dataset               = Dataset("/Users/joshil/Documents/Orbits/From_Eric/helix/ROSEL_Helical_orbits/ROSE_L_Helical_4sat.nc")



t12_orbits 		        = orbit_dataset["time"][1:2] # first two time samples
orbit_time_index      = Int(1):length(orbit_dataset["time"])
 # index range for orbit times for time interval of interest
orbit_time1 		      = orbit_dataset["time"][orbit_time_index] # read in time data
orbit_pos_ECI 	      = 1e3*orbit_dataset["position"][:,:,orbit_time_index] # read in position data, 3 x Np x Nt
orbit_vel_ECI         = 1e3*orbit_dataset["velocity"][:,:,orbit_time_index] # read in velocity data, 3 x Np x Nt (used optionally in avg peg and heading calculation)
dv 				            = orbit_dataset.attrib["epoch"];
if typeof(dv) == String
  temp = dv
  dv= []
  dv = [parse(Int,temp[1:4]);parse(Int,temp[6:7]);parse(Int,temp[9:10]);parse(Int,temp[12:13]);parse(Int,temp[15:16]);parse(Int,temp[18:19])];
end
epoch 			          = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
#global dcm 		        = Orbits.eci_dcm(orbit_time1, epoch);
global dcm = orbit_dataset["dcm"][:,:,orbit_time_index] 
orbit_pos1,orbit_vel1 = Orbits.ecef_orbitpos(orbit_pos_ECI,orbit_vel_ECI,dcm)


orbit_time_all = orbit_time1
orbit_pos_all = orbit_pos1
orbit_vel_all = orbit_vel1


#global mast_plat            = 1
#flag_plat                   = 0 #descending orbit
#orbit_time_all, orbit_pos_all, orbit_vel_all, orbit_pos_geo_all = Global_Scale_Support.get_orbit_info_fromfile(orbit_dataset, mast_plat, flag_plat)

orbit_pos_geo = zeros(size(orbit_pos_all)) 

for i=1:size(orbit_pos_all)[2]
    for j=1:size(orbit_pos_all)[3]
        orbit_pos_geo[:,i,j] = Geometry.xyz_to_geo(orbit_pos_all[:,i,j])
    end
end

#=
pos = orbit_pos_all
vel=orbit_vel_all
itimes = 1
    #create geocetric TCN frame for reference satellite
    that, chat, nhat = Geometry.get_tcn(pos[:,1,itimes], vel[:,1,itimes]);

    lhat                = Geometry.tcn_lvec(that, chat, nhat, -30.0 *(pi/180), 0.0 *(pi/180))

    #iterate over platforms and compute full perp-baseline matrix
    for iplat = 1:Nplats
        for jplat = 1:Nplats
            baseline = pos[:,jplat,itimes] - pos[:,iplat,itimes]

            bnorm[iplat,jplat,itimes] = norm(baseline)
            # use the imaging plane baseline only
            baseline_nc = dot(baseline,nhat)*nhat + dot(baseline,chat)*chat + dot(baseline,that)*that
            # find baseline component perpendicular to look direction
            bperp[iplat,jplat,itimes] = norm(baseline_nc - dot(baseline_nc,lhat)*lhat)
            # find the along-track baseline component
            b_at[iplat,jplat,itimes]  = abs(dot(baseline,that))
        end
    end
=#


look_angle=30
mast_plat=1
p_mode=1
Np=4
earth_radius=6378.137e3

global Norm_baseline_max = zeros(length(orbit_time_all))
global Norm_baseline_min = zeros(length(orbit_time_all))
global Perp_baseline_max = zeros(length(orbit_time_all))
global Perp_baseline_min = zeros(length(orbit_time_all))
global Par_baseline_max = zeros(length(orbit_time_all))
global Par_baseline_min = zeros(length(orbit_time_all))
global avg_sep = zeros(length(orbit_time_all))
global res_theory_n = zeros(length(orbit_time_all))
global amb_H = zeros(length(orbit_time_all))
global amb_N = zeros(length(orbit_time_all))

global amb_H_12 = zeros(length(orbit_time_all))
global amb_H_13 = zeros(length(orbit_time_all))
global amb_H_14 = zeros(length(orbit_time_all))

global amb_H_12_ic = zeros(length(orbit_time_all))
global amb_H_13_ic = zeros(length(orbit_time_all))
global amb_H_14_ic = zeros(length(orbit_time_all))

global Norm_baseline_max2 = zeros(length(orbit_time_all))
global Norm_baseline_min2 = zeros(length(orbit_time_all))


bperp, b_at, bnorm  = Orbits.get_perp_baselines_new(orbit_pos_all[:,:,:], orbit_vel_all[:,:,:], look_angle, 0.0, "left",  mast_plat)

bperp2, b_at2, bnorm2  = Orbits.get_perp_baselines_new(orbit_pos_all[:,2:end,:], orbit_vel_all[:,2:end,:], look_angle, 0.0, "left",  mast_plat)

for i=1:length(orbit_time_all)

    global Norm_baseline_max[i] = maximum(bnorm[:,:,i]) ./ 1e3
    global Norm_baseline_min[i] = minimum(filter(!iszero,bnorm[:,:,i])) ./ 1e3

    global Norm_baseline_max2[i] = maximum(bnorm2[:,:,i]) ./ 1e3
    global Norm_baseline_min2[i] = minimum(filter(!iszero,bnorm2[:,:,i])) ./ 1e3
    
    global Perp_baseline_max[i] = maximum(bperp[:,:,i]) ./ 1e3
    global Perp_baseline_min[i] = minimum(filter(!iszero,bperp[:,:,i])) ./ 1e3
    
    global Par_baseline_max[i]  = maximum(b_at[:,:,i]) ./ 1e3
    global Par_baseline_min[i]  = minimum(filter(!iszero,b_at[:,:,i])) ./ 1e3
    
    global avg_sep[i]             = maximum(bperp[:,:,i])/(Np - 1)

    # theoretical resolution along-n
    range_s, range_g = Scene.lookangle_to_range(look_angle, mean(Geometry.xyz_to_geo(orbit_pos_all[:,mast_plat,i])[3,:]), 0, earth_radius)
    global res_theory_n[i]      = (299792458/1.24e9)*range_s/p_mode/ (Perp_baseline_max[i]*1e3)

    global amb_H[i]             =  (299792458/1.24e9)*range_s/p_mode/avg_sep[i]*sind(look_angle)
    global amb_N[i]             =  (299792458/1.24e9)*range_s/p_mode/avg_sep[i]

    global amb_H_12[i] =  (299792458/1.24e9)*range_s/p_mode/(bperp[1,2,i]/(Np - 1))*sind(look_angle)
    global amb_H_13[i] =  (299792458/1.24e9)*range_s/p_mode/(bperp[1,3,i]/(Np - 1))*sind(look_angle)
    global amb_H_14[i] =  (299792458/1.24e9)*range_s/p_mode/(bperp[1,4,i]/(Np - 1))*sind(look_angle)


    global amb_H_12_ic[i] =  (299792458/1.24e9)*range_s/4/(bperp[1,2,i])*sind(look_angle)
    global amb_H_13_ic[i] =  (299792458/1.24e9)*range_s/4/(bperp[1,3,i])*sind(look_angle)
    global amb_H_14_ic[i] =  (299792458/1.24e9)*range_s/4/(bperp[1,4,i])*sind(look_angle)


end


p1=(plot(orbit_pos_geo[1,2,:],orbit_pos_geo[3,2,:]./1e3,ylabel="Altitude [km]",xlabel="Latitude [deg]",title ="", label="Master", lw=2, color=:tan1,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,3,:],orbit_pos_geo[3,3,:]./1e3,ylabel="Altitude [km]",xlabel="Latitude [deg]",title ="", label="Co-flier 1", lw=2,color=:gold,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,4,:],orbit_pos_geo[3,4,:]./1e3,ylabel="Altitude [km]",xlabel="Latitude [deg]",title ="", label="Co-flier 2", lw=2,color=:maroon4,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,1,:],orbit_pos_geo[3,1,:]./1e3,ylabel="Altitude [km]",xlabel="Latitude [deg]",title ="", label="Co-flier 3", lw=2, xlim=(-100,100),ylim=(690,730), color=:steelblue,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm ))
display(p1)

p1=(plot(orbit_pos_geo[1,2,:],(orbit_pos_geo[3,1,:]./1e3) .- (orbit_pos_geo[3,2,:]./1e3),ylabel="Vertical Seperation [km]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm ))
display(p1)

p1=(plot(orbit_pos_geo[1,2,:],(orbit_pos_geo[3,1,:]./1e3) .- (orbit_pos_geo[3,3,:]./1e3),ylabel="Vertical Seperation [km]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm ))
display(p1)

p1=(plot(orbit_pos_geo[1,2,:],(orbit_pos_geo[3,1,:]./1e3) .- (orbit_pos_geo[3,4,:]./1e3),ylabel="Vertical Seperation [km]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm ))
display(p1)




display(plot(orbit_pos_geo[1,1,:],Perp_baseline_max,ylabel="Perpendicular baseline [km]",xlabel="Latitude [deg]",title ="", legend=:false, lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))


display(plot(orbit_pos_geo[1,1,:],Norm_baseline_max,ylabel="Baseline [km]",xlabel="Latitude [deg]",title ="", legend=:false, lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))

display(plot(orbit_pos_geo[1,1,:],Norm_baseline_min,ylabel="Baseline [km]",xlabel="Latitude [deg]",title ="", legend=:false, lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))

display(plot(orbit_pos_geo[1,1,:],Norm_baseline_min2,ylabel="Baseline [km]",xlabel="Latitude [deg]",title ="", legend=:false, lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))


display(plot(orbit_pos_geo[1,1,:],Par_baseline_max,ylabel="Along-track Baseline [km]",xlabel="Latitude [deg]",title ="", legend=:false, lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))



p1=(plot(orbit_pos_geo[1,1,:],bperp[1,2,:]./1e3,ylabel="Perpendicular baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 2", lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,1,:],bperp[1,3,:]./1e3,ylabel="Perpendicular baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 3", lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,1,:],bperp[1,4,:]./1e3,ylabel="Perpendicular baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 4", lw=2, rightmargin=5mm, leftmargin=5mm,bottommargin=5mm, 
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480), legend=:bottomright,xlim=(-90,90) ))
display(p1)

p1=(plot(orbit_pos_geo[1,1,:],bnorm[1,2,:]./1e3,ylabel=" Baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 2", lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,1,:],bnorm[1,3,:]./1e3,ylabel="Baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 3", lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,1,:],bnorm[1,4,:]./1e3,ylabel="Baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 4", lw=2, rightmargin=5mm, leftmargin=5mm,bottommargin=5mm, 
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480), legend=:bottomright,xlim=(-90,90)  ))
display(p1)



p1=(plot(orbit_pos_geo[1,1,:],b_at[1,2,:]./1e3,ylabel="AT baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 2", lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,1,:],b_at[1,3,:]./1e3,ylabel="AT baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 3", lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
p1=(plot!(orbit_pos_geo[1,1,:],b_at[1,4,:]./1e3,ylabel="AT baseline [km]",xlabel="Latitude [deg]",title ="", label="Master - Coflier 4", lw=2, rightmargin=5mm, leftmargin=5mm,bottommargin=5mm, 
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480), legend=:bottomright,xlim=(-90,90)  ))
display(p1)


p1=(plot(orbit_pos_geo[1,2,:],amb_H_12,ylabel="HOA h [m]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm ))
p1=(plot!(orbit_pos_geo[1,2,:],amb_H_13,ylabel="HOA h [m]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm ))
p1=(plot!(orbit_pos_geo[1,2,:],amb_H_14,ylabel="HOA h [m]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm,ylim=(200,400) ))
display(p1)



p1=(plot(orbit_pos_geo[1,2,:],amb_H_12_ic,ylabel="pi HOA [m]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm ))
p1=(plot!(orbit_pos_geo[1,2,:],amb_H_13_ic,ylabel="pi HOA [m]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm ))
p1=(plot!(orbit_pos_geo[1,2,:],amb_H_14_ic,ylabel="pi HOA [m]",xlabel="Latitude [deg]",title ="", label="Vertical Seperation ", lw=2, legend=:false,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(750,480),leftmargin=5mm,rightmargin=5mm,bottommargin=5mm,ylim=(15,50) ))
display(p1)

display(plot(orbit_time_all[1:720]./60,Norm_baseline_min[1:720],ylabel="Baseline [km]",xlabel="Time [min]",title ="", legend=:false, lw=2,
tickfont=font(13), xtickfont=font(13), ytickfont=font(13), guidefont=font(13), titlefontsize=13, size=(500,360) ))
