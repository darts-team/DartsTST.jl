using XLSX
using Plots
using LinearAlgebra

figsavepath_common          = "/Users/joshil/Documents/Code/Plots_custom_orbits_study/Equal_spacing_5km_v3/MIMO/5x5/"
UAV_sims_file   = figsavepath_common*"Output_stat.xlsx"
xf              = XLSX.readxlsx(UAV_sims_file);
sh              = xf["Sheet1"]

Num_iter        = 6
factor          = 0 #0.25*6

theo_n          = zeros(1,Num_iter)
sim_n_bpa       = zeros(1,Num_iter)
sim_n_beamforming = zeros(1,Num_iter)
sim_n_capon     = zeros(1,Num_iter)

theo_s          = zeros(1,Num_iter)
sim_s_bpa       = zeros(1,Num_iter)
sim_s_beamforming = zeros(1,Num_iter)
sim_s_capon     = zeros(1,Num_iter)

theo_c          = zeros(1,Num_iter)
sim_c_bpa       = zeros(1,Num_iter)
sim_c_beamforming = zeros(1,Num_iter)
sim_c_capon     = zeros(1,Num_iter)

islr_bpa       = zeros(1,Num_iter)
islr_beamforming = zeros(1,Num_iter)
islr_capon     = zeros(1,Num_iter)

for i=1:Num_iter

    idx             = (10*i)+11 - 20
    theo_n[i]       = sh[:][idx,1]
    sim_n_bpa[i]    = sh[:][idx+4,1]
    sim_n_beamforming[i] = sh[:][idx+5,1]
    sim_n_capon[i]  = sh[:][idx+6,1]

    theo_s[i]       = sh[:][idx,2] + factor
    sim_s_bpa[i]    = sh[:][idx+1,1]
    sim_s_beamforming[i] = sh[:][idx+2,1]
    sim_s_capon[i]  = sh[:][idx+3,1]

    theo_c[i]       = 6.70 #(sh[:][idx,3] / sind(30)) + factor
    sim_c_bpa[i]    = sh[:][idx+1,2]
    sim_c_beamforming[i] = sh[:][idx+2,2]
    sim_c_capon[i]  = sh[:][idx+3,2]

    islr_bpa[i]    = sh[:][idx+4,3]
    islr_beamforming[i] = sh[:][idx+5,3]
    islr_capon[i]  = sh[:][idx+6,3]
    

end

Numplat = [2; 3; 4; 5; 6; 7]
#Numplat = [3 ;5 ;7 ;9 ;11; 15; 17]

Numplat = Int64.(Numplat)

p1 = plot(Numplat, transpose(theo_n), label="Theoretical", xlabel="Number of platforms", ylabel="Resolution along tomographic axis (n axis) [m]", title="Resolution verus Number of platforms for 1 target",
linewidth=3, linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,25), size=(800,600))
p1 =plot!(Numplat, transpose(sim_n_bpa), label="Back projection",linewidth=3)
p1 =plot!(Numplat, transpose(sim_n_beamforming), label="Beamforming",linewidth=3)
p1 =(plot!(Numplat, transpose(sim_n_capon), label="CAPON",linewidth=3))

p2 = plot(Numplat, transpose(theo_s), label="Theoretical", xlabel="Number of platforms", ylabel="Resolution along along-track (S axis) [m]", title="Resolution verus Number of platforms for 1 target",
linewidth=3,linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,5),size=(800,600))
p2 = plot!(Numplat, transpose(sim_s_bpa), label="Back projection",linewidth=3)
p2 = plot!(Numplat, transpose(sim_s_beamforming), label="Beamforming",linewidth=3)
p2 = (plot!(Numplat, transpose(sim_s_capon), label="CAPON",linewidth=3))

p3 = plot(Numplat, transpose(theo_c), label="Theoretical", xlabel="Number of platforms", ylabel="Resolution along ground range (C axis) [m]", title="Resolution verus Number of platforms for 1 target",
linewidth=3,linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,15), size=(800,600))
p3 = plot!(Numplat, transpose(sim_c_bpa), label="Back projection",linewidth=3)
p3 = plot!(Numplat, transpose(sim_c_beamforming), label="Beamforming",linewidth=3)
p3 = (plot!(Numplat, transpose(sim_c_capon), label="CAPON",linewidth=3))


savefig(p1, figsavepath_common*"Resolution_N.png")
savefig(p2, figsavepath_common*"Resolution_S.png")
savefig(p3, figsavepath_common*"Resolution_C.png")

p1 = plot(Numplat, transpose(theo_n), label="Theoretical", xlabel="Number of platforms", ylabel="Resolution along tomographic axis (n axis) [m]", title="Resolution verus Number of platforms for 1 target",
linewidth=3, linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,25), size=(800,600))
p1 =plot!(Numplat, transpose(sim_n_bpa), label="Back projection",linewidth=3)

p2 = plot(Numplat, transpose(theo_s), label="Theoretical", xlabel="Number of platforms", ylabel="Resolution along along-track (S axis) [m]", title="Resolution verus Number of platforms for 1 target",
linewidth=3,linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,5),size=(800,600))
p2 = plot!(Numplat, transpose(sim_s_bpa), label="Back projection",linewidth=3)

p3 = plot(Numplat, transpose(theo_c), label="Theoretical", xlabel="Number of platforms", ylabel="Resolution along ground range (C axis) [m]", title="Resolution verus Number of platforms for 1 target",
linewidth=3,linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,15), size=(800,600))
p3 = plot!(Numplat, transpose(sim_c_bpa), label="Back projection",linewidth=3)


savefig(p1, figsavepath_common*"Resolution_N2.png")
savefig(p2, figsavepath_common*"Resolution_S2.png")
savefig(p3, figsavepath_common*"Resolution_C2.png")



p1 = plot(Numplat, transpose(islr_bpa), label="Back projection - MIMO", xlabel="Number of platforms", ylabel="Integrated Side Lobe Ratio (ISLR)", title="ISLR verus Number of platforms for 1 target",
linewidth=3, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:bottomright,legendfont=font(15), ylim=(0,25), size=(800,600))

p2 = plot(Numplat, transpose(islr_beamforming), label="Beamforming", xlabel="Number of platforms", ylabel="Integrated Side Lobe Ratio (ISLR)", title="ISLR verus Number of platforms for 1 target",
linewidth=3, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:bottomright,legendfont=font(15), ylim=(0,25),size=(800,600))

p3 = plot(Numplat, transpose(islr_capon), label="CAPON", xlabel="Number of platforms", ylabel="Integrated Side Lobe Ratio (ISLR)", title="ISLR verus Number of platforms for 1 target",
linewidth=3, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:bottomright,legendfont=font(15), ylim=(0,25), size=(800,600))


savefig(p1, figsavepath_common*"ISLR_BPA_report.png")
savefig(p2, figsavepath_common*"ISLR_Beam.png")
savefig(p3, figsavepath_common*"ISLR_Capon.png")



#For paper
mycolor3 = theme_palette(:auto)
p1 = plot(Numplat, sar_theo, label="SISO - Theoretical", xlabel="Number of platforms", ylabel="Resolution along tomographic axis (n axis) [m]", title="Resolution verus Number of platforms for 1 target",
linewidth=3, linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,25), size=(800,600), color=mycolor3[1])
p1 =plot!(Numplat, sar_sim, label="SISO - Back projection",linewidth=3, color=mycolor3[1])
p1 = plot!(Numplat, simo_theo, label="SIMO - Theoretical", xlabel="Number of platforms", ylabel="Resolution along tomographic axis (n axis) [m]", title="Resolution verus Number of platforms for 1 target",
linewidth=3, linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,25), size=(800,600), color=mycolor3[2])
p1 =plot!(Numplat, simo_sim, label="SIMO - Back projection",linewidth=3, color=mycolor3[2])
p1 = plot!(Numplat, mimo_theo, label="MIMO - Theoretical", xlabel="Number of platforms", ylabel="Resolution along tomographic axis (n axis) [m]", title="Tomographic Resolution verus Number of platforms for 1 target",
linewidth=3, linestyle=:dash, tickfont=font(15), ytickfont=font(15), guidefont=font(15),legend=:topright,legendfont=font(15), ylim=(0,25), size=(800,600), color=mycolor3[3])
p1 =plot!(Numplat, mimo_sim, label="MIMO - Back projection",linewidth=3, color=mycolor3[3])
savefig(p1, "/Users/joshil/Documents/Journal_papers/TGRS/Fig_res_plat_2.png")