using Distributed
addprocs(1) 

@everywhere include("../../modules/generate_raw_data.jl")
@everywhere include("../../modules/process_raw_data.jl")
@everywhere include("../../modules/geometry.jl")
@everywhere include("../../modules/scene.jl")
@everywhere include("../../modules/range_spread_function.jl") # as RSF
@everywhere include("../../modules/orbits.jl")
@everywhere include("../../modules/sync.jl")
@everywhere include("../../modules/error_sources.jl")
@everywhere include("../../modules/performance_metrics.jl")
@everywhere include("../../modules/antenna.jl")
@everywhere include("../../modules/simsetup.jl")
@everywhere include("../../modules/user_parameters.jl")
@everywhere include("../../modules/data_processing.jl")
@everywhere include("../../modules/data_plotting.jl")
@everywhere include("../../modules/scattering.jl")

@everywhere using NCDatasets
@everywhere using Statistics
@everywhere using Parameters
@everywhere using Dates
@everywhere using StaticArrays
@everywhere using Plots
@everywhere using LinearAlgebra
@everywhere using .UserParameters
@everywhere using TimerOutputs
@everywhere using Interpolations
@everywhere using XLSX

to = TimerOutput()

@everywhere c = 299792458 #TODO does not work without redefining c here
@everywhere earth_radius = 6378.137e3 # Earth semi-major axis at equator

# Define parametsr for signal processing - Beamforming , CAPON etc
@everywhere filt_len                    = 5 # Covariance matrix box filter size 3->5x5 box filter
#figsavepath_common          = "/Users/joshil/Documents/Code/Plots_custom_orbits_study/Equal_spacing_5km_correct/SIMO/5x5/"

#figsavepath_common          = "/Users/joshil/Documents/Code/global_scale_outputs/plots_5p_v2/temp/1/"
#figsavepath_common          = "/Users/joshil/Documents/Conferences/Postdoc_research_day/plots/"
@everywhere figsavepath_common          = "/Users/joshil/Documents/Code/Plots_custom_orbits_study/Equal_spacing_5km_v3/SIMO/test/"


@sync @distributed for loop_i=2:7 #2:7

try

p_spacing = 5
#pos_n_ip   = transpose( (range(-loop_i,loop_i,1+(2*loop_i))*p_spacing) *ones(1,1) )*1e3 

    if loop_i==2
        pos_n_ip   = [-2.5 2.5]*1e3 
    elseif loop_i==3
        pos_n_ip   = [-5 0 5]*1e3 
    elseif loop_i==4
        pos_n_ip   = [-7.5 -2.5 2.5 7.5]*1e3 
    elseif loop_i==5
        pos_n_ip   = [-10 -5 0 5 10]*1e3 
    elseif loop_i==6
        pos_n_ip   = [-12.5 -7.5 -2.5 2.5 7.5 12.5]*1e3 
    elseif loop_i==7
        pos_n_ip   = [-15 -10 -5 0 5 10 15]*1e3 
    else
        pos_n_ip=[0]
    end
# Define user parameters
#include("../inputs/predefined-input-parameters.jl") TODO gives errors
params = UserParameters.inputParameters(mode= 2, pos_n = pos_n_ip)

# Check consistency of input parameters
paramsIsValid = UserParameters.validateInputParams(params)

    # theoretical resolution
    if params.mode == 1 # SAR
        p_mode = 2
    elseif params.mode == 2 # SIMO
        p_mode = 1
    elseif params.mode == 3 # MIMO
        p_mode = 1.38
    end

# Compute orbits time, position, and velocity

orbit_time, orbit_pos, orbit_vel = Orbits.computeTimePosVel(params)


# interpolate orbit to slow time, 3 x Np x Nst, convert km to m
 p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)
 v_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_vel, params)

 orbit_pos_geo_all=  Geometry.xyz_to_geo(p_xyz[:,1,:])

# Create target/scene location
targets_loc, targets_ref, Nt = Scene.construct_targets_str(params) # Nt: number of targets, targets: structure array containing target locations and reflectivities
 s_loc_3xN  = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops
t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions
#t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, orbit_vel, params) ## calculate avg heading from platform positions/velocities

# Read number of platforms (todo: move into a struct)
 Np  = size(orbit_pos)[2] # number of platforms

#Theoretical computations
    # input max baseline along-n (Bn = Np x dn)
    spacing = params.pos_n[2]-params.pos_n[1] # spacing for equally spaced platforms
    #spacing = mean(diff(params.pos_n,dims=2)) # average spacing for unequally spaced platforms
    max_baseline_n = (maximum(params.pos_n)-minimum(params.pos_n)) + spacing
    println();println("input max baseline along-n: ",max_baseline_n)

    # theoretical resolution along-n
    range_s, range_g = Scene.lookangle_to_range(params.look_angle, params.p_t0_LLH[3], 0, earth_radius)
    res_theory_n = (c/params.fc)*range_s/p_mode/max_baseline_n
    println("theoretical resolution along-n: ",round(res_theory_n,digits=2))

    # theoretical resolution along-track
    mu = 3.986004418e14
    sc_speed = sqrt(mu./(params.p_t0_LLH[3]+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
    Lsa = sc_speed*params.SAR_duration + sc_speed/params.fp
    res_theory_s = (c/params.fc)*range_s/2/Lsa
    println("theoretical resolution along track: ",round(res_theory_s,digits=2))

    res_theory_sr = (c/ (2*params.bandwidth))
    println("theoretical resolution along slant-range: ",round(res_theory_sr,digits=2))

    mean_p_xyz = mean(mean(p_xyz,dims=3),dims=2)
    Rs1 = Geometry.distance(mean_p_xyz, t_xyz_3xN)
    Rs2 = Rs1 + res_theory_sr
    Rg1,teta1 = Scene.slantrange_to_lookangle(earth_radius,Rs1,Geometry.xyz_to_geo(mean_p_xyz)[3],0) 
    Rg2,teta2 = Scene.slantrange_to_lookangle(earth_radius,Rs2,Geometry.xyz_to_geo(mean_p_xyz)[3],0) 
    del_Rg = Rg2-Rg1
    println("Theoretical ground range resolution: ", del_Rg)




# Apply antenna pattern
if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
    Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
end

# Generate range spread function (matched filter output)
 min_range, max_range = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
 Trx = 2*(max_range-min_range)/c + 5*params.pulse_length # s duration of RX window #100 for co-flyer
 Srx, MF, ft, t_rx = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
# Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

# Generate TomoSAR raw data
 ref_range = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
 rawdata = Generate_Raw_Data.main_RSF_slowtime_surfaceBRCS(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
if params.enable_thermal_noise # adding random noise based on SNR after range (fast-time) processing
     rawdata = Error_Sources.random_noise(rawdata, params)
end

# Add phase error
 sync_osc_coeffs = repeat(params.sync_a_coeff_dB, Np)
if params.enable_sync_phase_error
     rawdata = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, t_xyz_3xN, sync_osc_coeffs, params)
end
#=
@timeit to "BPA using existing function " begin
# Process raw data to generate image
if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
    const image_3D = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
    const SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
    const image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
end
end
=#
#@timeit to "BPA using new function " begin
#const image_3D = Process_Raw_Data.main_SAR_tomo_3D_new(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
#end

#@timeit to "BPA using old function " begin
#const image_3D2 = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
#end


@timeit to "BPA using new function " begin
# Process raw data to generate image
if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
     image_3D = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
    if params.processing_mode == 1
        SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
    elseif params.processing_mode == 2
        # for co-flyer
        SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata[:,2:size(p_xyz)[2],:], s_xyz_3xN, p_xyz, t_rx, ref_range, params)
    end
    image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
end
end

norm_flag  					= 1
plot_idx 					= [Int64(ceil(length(params.s_loc_1)/2)),Int64(ceil(length(params.s_loc_2)/2)),Int64(ceil(length(params.s_loc_3)/2))]
figsavepath                 = figsavepath_common*string(loop_i)*"/"
if ~ispath(figsavepath)
    mkdir(figsavepath)
end
Data_Plotting.plot_tomo_output(image_3D, params,norm_flag, plot_idx, figsavepath, 20, "Back-projection output \n Normalized Intensity [dB]", "_orig")

image_3D = Data_Processing.average_2D_data(image_3D, filt_len)

#include("../modules/data_plotting.jl")

Data_Plotting.plot_tomo_output(image_3D, params,norm_flag, plot_idx, figsavepath, 20, "Back-projection output \n Normalized Intensity [dB]", "_avg")


# Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
if params.left_right_look == "left";PSFcutdir=-1;elseif params.left_right_look == "right";PSFcutdir=1;end
params.PSF_direction[3]=PSFcutdir*params.PSF_direction[3]
scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)

# Calculate point target performance metrics
if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
     resolutions, PSLRs, ISLRs, loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
    println("Resolutions: ",round.(resolutions,digits=8)," in scene axes units")
    println("Location Errors: ",round.(loc_errors,digits=8)," in scene axes units")
    println("PSLRs: ",round.(PSLRs,digits=2)," dB")
    println("ISLRs: ",round.(ISLRs,digits=2)," dB")
    println("PSF Peak Amplitude: ",round(maximum(20*log10.(image_3D)),digits=2)," dB")
else
    @warn "PSF related performance metrics cannot be calculated for more than 1 target."
end


# Define parametsr for signal processing - Beamforming , CAPON etc
Master_platform             = 1#Int64(ceil(length(pos_n_ip)/2)) # Master platform
#filt_len                    = 0 # Covariance matrix box filter size 3->5x5 box filter
azimuth_lim 				= [1, Int64(ceil(length(params.s_loc_1)))]
srange_lim 					= [1, Int64(ceil(length(params.s_loc_2)))]
ref_hl                      = Int64(ceil(length(params.s_loc_3)/2))
heights_t 					= params.s_loc_3 #2*minimum(params.s_loc_3):params.s_loc_3[2]-params.s_loc_3[1]:2*maximum(params.s_loc_3);#params.s_loc_3

# Modify input matrix to desired shape
input_SP_temp               = SAR_images_3D[:,:,:,ref_hl];
input_SP                    = permutedims(input_SP_temp,[2,1,3])

Ns1                         = size(input_SP)[1] # Nst
Ns2                         = size(input_SP)[2] # Np
Ns3                         = size(input_SP)[3] # Nr

Cov_mat                     = zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Ns2, Ns2);
Corr_mat                    = zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Ns2, Ns2);
Pbf                         = zeros(Ns1,Ns3,length(heights_t));
PC                          = zeros(Ns1,Ns3,length(heights_t));

s_xyz_3xN_2D_all            = zeros(3,Int64(ceil(length(params.s_loc_1))),Int64(ceil(length(params.s_loc_2))),Int64(ceil(length(params.s_loc_3))))
s_xyz_3xN_2D                = zeros(3,Int64(ceil(length(params.s_loc_1))),Int64(ceil(length(params.s_loc_2))))
l                           = 1
for i=1:length(params.s_loc_1)
    for j=1:length(params.s_loc_2)
        for k=1:length(params.s_loc_3)
            s_xyz_3xN_2D_all[:,i,j,k] = s_xyz_3xN[:,k] 
            l               = l + 1
        end
    end
end
s_xyz_3xN_2D                = s_xyz_3xN_2D_all[:,:,:,ref_hl]

#@time Cov_mat, Corr_mat           = Data_Processing.get_covariance_correlation_matrices(input_SP, azimuth_lim, srange_lim, Ns2, filt_len, 0);
#@time Cov_mat2, Corr_mat2           = Data_Processing.get_covariance_correlation_matrices_2(input_SP, azimuth_lim, srange_lim, Ns2, filt_len, 0);
@time Cov_mat3, Corr_mat3           = Data_Processing.get_covariance_correlation_matrices_new(input_SP,  Ns2, filt_len);

steering_mat                = Data_Processing.get_steering_matrix(p_xyz, s_xyz_3xN_2D, azimuth_lim, srange_lim, heights_t, Ns2, Master_platform, params.λ, params.mode, params.processing_mode);

Pbf                         = Data_Processing.tomo_beamforming(Cov_mat3, steering_mat, azimuth_lim, srange_lim, [size(Cov_mat)[1] size(Cov_mat)[2] size(heights_t)[1]])
Pbf = Pbf[:,:,end:-1:1];
Data_Plotting.plot_tomo_output(Pbf, params,norm_flag, plot_idx,figsavepath, 10, "Beamforming output \n Normalized Intensity [dB]", "_orig", heights_t)
Pbf2                        = Data_Processing.tomocoordinates_to_scenecoordinates(Pbf, heights_t, params.s_loc_2, params.s_loc_3, params.look_angle, params.left_right_look, mean(orbit_pos_geo_all[3,:]))
Pbf2[Pbf2.==0].=1e-6
Data_Plotting.plot_tomo_output(Pbf2, params,norm_flag, plot_idx,figsavepath, 10, "Beamforming output \n Normalized Intensity [dB]", "_2", heights_t)

PC                          = Data_Processing.tomo_CAPON(Cov_mat3, steering_mat, azimuth_lim, srange_lim, [size(Cov_mat)[1] size(Cov_mat)[2] size(heights_t)[1]])
PC = PC[:,:,end:-1:1];
Data_Plotting.plot_tomo_output(PC, params,norm_flag, plot_idx,figsavepath, 10, "CAPON output \n Normalized Intensity [dB]", "_orig", heights_t)
PC2                         = Data_Processing.tomocoordinates_to_scenecoordinates(PC, heights_t, params.s_loc_2, params.s_loc_3, params.look_angle, params.left_right_look, mean(orbit_pos_geo_all[3,:]))
PC2[PC2.==0].=1e-6
Data_Plotting.plot_tomo_output(PC2, params,norm_flag, plot_idx,figsavepath, 10, "CAPON output  \n Normalized Intensity [dB]", "_2", heights_t)


scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(sqrt.(Pbf2), params)
# Calculate point target performance metrics Beamforming and CAPON
if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
     Beamforming_resolutions, Beamforming_PSLRs, Beamforming_ISLRs, Beamforming_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
    println("Beamforming Resolutions: ",round.(Beamforming_resolutions,digits=8)," in scene axes units")
    println("Beamforming Location Errors: ",round.(Beamforming_loc_errors,digits=8)," in scene axes units")
    println("Beamforming PSLRs: ",round.(Beamforming_PSLRs,digits=2)," dB")
    println("Beamforming ISLRs: ",round.(Beamforming_ISLRs,digits=2)," dB")
    println("Beamforming PSF Peak Amplitude: ",round(maximum(20*log10.(Pbf2)),digits=2)," dB")
else
    @warn "Beamforming PSF related performance metrics cannot be calculated for more than 1 target."
end

scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(sqrt.(PC2), params)
# Calculate point target performance metrics Beamforming and CAPON
if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
     CAPON_resolutions, CAPON_PSLRs, CAPON_ISLRs, CAPON_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
    println("CAPON Resolutions: ",round.(CAPON_resolutions,digits=8)," in scene axes units")
    println("CAPON Location Errors: ",round.(CAPON_loc_errors,digits=8)," in scene axes units")
    println("CAPON PSLRs: ",round.(CAPON_PSLRs,digits=2)," dB")
    println("CAPON ISLRs: ",round.(CAPON_ISLRs,digits=2)," dB")
    println("CAPON PSF Peak Amplitude: ",round(maximum(20*log10.(PC2)),digits=2)," dB")
else
    @warn "CAPON PSF related performance metrics cannot be calculated for more than 1 target."
end


params = UserParameters.inputParameters(PSF_cuts = 2, PSF_direction = [0 1 -tand(params.inc_angle)]) 

# Along tomography axis
if params.left_right_look == "left";PSFcutdir=-1;elseif params.left_right_look == "right";PSFcutdir=1;end
params.PSF_direction[3]=PSFcutdir*params.PSF_direction[3]
#scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)
image_1D_1, scene_axis11, scene_axis22, scene_axis33 = Scene.obtain_1D_slice_tilted(image_3D, params.s_loc_1, params.s_loc_2, params.s_loc_3, params.PSF_direction)
scene_res=((scene_axis11[2]-scene_axis11[1])^2+(scene_axis22[2]-scene_axis22[1])^2+(scene_axis33[2]-scene_axis33[1])^2)^0.5 # scene resolution along the PSF direction


# Calculate point target performance metrics
if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
     BPA_tomo_resolutions, BPA_tomo_PSLRs, BPA_tomo_ISLRs, BPA_tomo_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
    println("Tomo_Resolutions Back projection: ",round.(BPA_tomo_resolutions,digits=8)," in scene axes units")
else
    @warn "PSF related performance metrics cannot be calculated for more than 1 target."
end

image_1D_1, scene_axis11, scene_axis22, scene_axis33 = Scene.obtain_1D_slice_tilted(sqrt.(Pbf2), params.s_loc_1, params.s_loc_2, params.s_loc_3, params.PSF_direction)
scene_res=((scene_axis11[2]-scene_axis11[1])^2+(scene_axis22[2]-scene_axis22[1])^2+(scene_axis33[2]-scene_axis33[1])^2)^0.5 # scene resolution along the PSF direction

# Calculate point target performance metrics
if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
     Beamforming_tomo_resolutions, Beamforming_tomo_PSLRs, Beamforming_tomo_ISLRs, Beamforming_tomo_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
    println("Tomo_Resolutions beamforming: ",round.(Beamforming_tomo_resolutions,digits=8)," in scene axes units")
else
    @warn "PSF related performance metrics cannot be calculated for more than 1 target."
end

image_1D_1, scene_axis11, scene_axis22, scene_axis33 = Scene.obtain_1D_slice_tilted(sqrt.(PC2), params.s_loc_1, params.s_loc_2, params.s_loc_3, params.PSF_direction)
scene_res=((scene_axis11[2]-scene_axis11[1])^2+(scene_axis22[2]-scene_axis22[1])^2+(scene_axis33[2]-scene_axis33[1])^2)^0.5 # scene resolution along the PSF direction

# Calculate point target performance metrics
if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
     CAPON_tomo_resolutions, CAPON_tomo_PSLRs, CAPON_tomo_ISLRs, CAPON_tomo_loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
    println("Tomo_Resolutions CAPON: ",round.(CAPON_tomo_resolutions,digits=8)," in scene axes units")
else
    @warn "PSF related performance metrics cannot be calculated for more than 1 target."
end

li1 = (loop_i*10) +1 #line number
#UAV_sims_file = "/Users/joshil/Documents/Code/Plots_analysis/Output_stat simo.xlsx"
UAV_sims_file = figsavepath_common*"Output_stat.xlsx"
XLSX.openxlsx(UAV_sims_file, mode = "rw") do xf2
    sheet2 = xf2[1]
    sheet2["A"*string(li1)]        = round.(res_theory_n,digits=8)
    sheet2["B"*string(li1)]        = round.(res_theory_s,digits=8)
    sheet2["C"*string(li1)]        = round.(res_theory_sr,digits=8)
    sheet2["D"*string(li1)]        = loop_i
    sheet2["A"*string(li1+1), dim=2] = round.(resolutions,digits=8)
    sheet2["D"*string(li1+1), dim=2] = round.(loc_errors,digits=8)
    sheet2["G"*string(li1+1), dim=2] = round.(PSLRs,digits=2)
    sheet2["J"*string(li1+1), dim=2] = round.(ISLRs,digits=2)
    sheet2["M"*string(li1+1)]        = round(maximum(20*log10.(image_3D)),digits=2)
    sheet2["A"*string(li1+2), dim=2] = round.(Beamforming_resolutions,digits=8)
    sheet2["D"*string(li1+2), dim=2] = round.(Beamforming_loc_errors,digits=8)
    sheet2["G"*string(li1+2), dim=2] = round.(Beamforming_PSLRs,digits=2)
    sheet2["J"*string(li1+2), dim=2] = round.(Beamforming_ISLRs,digits=2)
    sheet2["M"*string(li1+2)]        = round(maximum(10*log10.(Pbf2)),digits=2)
    sheet2["A"*string(li1+3), dim=2] = round.(CAPON_resolutions,digits=8)
    sheet2["D"*string(li1+3), dim=2] = round.(CAPON_loc_errors,digits=8)
    sheet2["G"*string(li1+3), dim=2] = round.(CAPON_PSLRs,digits=2)
    sheet2["J"*string(li1+3), dim=2] = round.(CAPON_ISLRs,digits=2)
    sheet2["M"*string(li1+3)]        = round(maximum(10*log10.(PC2)),digits=2)
    sheet2["A"*string(li1+4)]        = round.(BPA_tomo_resolutions,digits=8)
    sheet2["A"*string(li1+5)]        = round.(Beamforming_tomo_resolutions,digits=8)
    sheet2["A"*string(li1+6)]        = round.(CAPON_tomo_resolutions,digits=8)
    sheet2["B"*string(li1+4)]        = round.(BPA_tomo_PSLRs,digits=8)
    sheet2["B"*string(li1+5)]        = round.(Beamforming_tomo_PSLRs,digits=8)
    sheet2["B"*string(li1+6)]        = round.(CAPON_tomo_PSLRs,digits=8)
    sheet2["C"*string(li1+4)]        = round.(BPA_tomo_ISLRs,digits=8)
    sheet2["C"*string(li1+5)]        = round.(Beamforming_tomo_ISLRs,digits=8)
    sheet2["C"*string(li1+6)]        = round.(CAPON_tomo_ISLRs,digits=8)
    sheet2["D"*string(li1+4)]        = round.(BPA_tomo_loc_errors,digits=8)
    sheet2["D"*string(li1+5)]        = round.(Beamforming_tomo_loc_errors,digits=8)
    sheet2["D"*string(li1+6)]        = round.(CAPON_tomo_loc_errors,digits=8)
end

catch
    continue
end

end

[rmprocs(p) for p in workers()]




#=
mean_p_xyz = mean(mean(p_xyz,dims=3),dims=2)
Rs1 = Geometry.distance(mean_p_xyz, t_xyz_3xN)
Rs2 = Rs1 + (c/ (2*params.bandwidth))
Rg1,teta1 = Scene.slantrange_to_lookangle(earth_radius,Rs1,Geometry.xyz_to_geo(mean_p_xyz)[3],0) 
Rg2,teta2 = Scene.slantrange_to_lookangle(earth_radius,Rs2,Geometry.xyz_to_geo(mean_p_xyz)[3],0) 
del_Rg = Rg2-Rg1
println("Theoretical ground range resolution: ", del_Rg)
=#
