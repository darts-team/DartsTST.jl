include("../modules/generate_raw_data.jl")
include("../modules/process_raw_data.jl")
include("../modules/geometry.jl")
include("../modules/scene.jl")
include("../modules/range_spread_function.jl") # as RSF
include("../modules/orbits.jl")
include("../modules/sync.jl")
include("../modules/error_sources.jl")
include("../modules/performance_metrics.jl")
include("../modules/antenna.jl")
include("../modules/simsetup.jl")
include("../modules/user_parameters.jl")
using NCDatasets
using Statistics
using Parameters
using Dates
using StaticArrays
using .UserParameters
using JLD
using Plots
gr()

c = 299792458 #TODO does not work without redefining c here
earth_radius = 6378.137e3 # Earth semi-major axis at equator

# Define study parameters
Ntr = 5# number of trials
n_platforms = 11 # number of platforms

# for equal spacing
init_spc = 1e3 # initial spacing
spc_inc = 1e3 # spacing increment
# for unequal spacing
max_dist = 45e3 # max distance along-n among platforms

res_theory_array=zeros(Ntr)
res_measured_array=zeros(Ntr)
pos_n_all=zeros(Ntr,n_platforms)

anim = Plots.Animation()
anim2 = Plots.Animation()
for i = 1:Ntr
    # Define user parameters
    params = UserParameters.inputParameters(
    mode = 2, #1: SAR (ping-pong), 2:SIMO, 3:MIMO
    look_angle = 30, # in cross-track direction, required only if SCH coordinates, using same look angle for targets and scene (deg)
    user_defined_orbit = 2, # 1: use orbits file; 2: user defined orbits in TCN
    p_t0_LLH = [0;0;750e3], # initial lat/lon (deg) and altitude (m) of reference platform (altitude is assumed constant over slow-time if SCH option)
    PSF_cuts = 2, # 1: principal axes (SCH, LLH, XYZ based on ts_coord_sys), 2: a single cut along PSF_direction_xyz in scene coordinates relative to center of scene
    Torbit    = 30, # orbital duration (s) (should be larger than 2 x (SAR_start_time+SAR_duration) )
    dt_orbits = 1, # orbit time resolution (s)
    p_heading = 0, # heading (deg), all platforms assumed to have the same heading, 0 deg is north
    left_right_look = "right", # left or right looking geometry
    display_custom_orbit = false, #whether to show orbit on Earth sphere (for a duration of Torbit)
    display_1D_cuts = true, # whether to 1D cuts from Scene module
    display_tomograms = true, # how to display tomograms, 0: do not display, 1: display only 3 slices at the reference point, 2: display all slices in each dimension, 3: display as 3D scatter plot
    display_geometry = false, # whether to display geometry plots
    s_loc_1 = 0, # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2 = -15:0.5:15, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3 = -15:0.5:15, # m  heights if LLH or SCH, Z if XYZ
    #pos_n   = ((1:n_platforms).-(n_platforms+1)/2)*(init_spc+(i-1)*spc_inc), # relative position of each platform along n (m), 0 is the reference location, equal spacing
    pos_n   = sort(round.([-max_dist/2 -max_dist/2*rand(1,Int((n_platforms-3)/2)) 0 max_dist/2*rand(1,Int((n_platforms-3)/2)) max_dist/2]),dims=2),
    #pos_n   = sort(round.([-max_dist/2*rand(1,Int((n_platforms-1)/2)) 0 max_dist/2*rand(1,Int((n_platforms-1)/2))]),dims=2),
    res_dB = 3.89 # dB two-sided resolution relative power level
    # Np:res_dB [2:3.01 3:3.52 4:3.70 5:3.78 6:3.82 7:3.85 8:3.87 9:3.88 10:3.89] for baseline = max distance + 1 spacing (for analytical resolution formula)
    )

    # theoretical resolution
    if params.mode == 1 # SAR
        p_mode = 2
    elseif params.mode == 2 # SIMO
        p_mode = 1
    elseif params.mode == 3 # MIMO
        p_mode = 1.38
    end

    # platform distributions along-n
    println("pos_n: ",Integer.(params.pos_n))
    pos_n_all[i,:]=params.pos_n

    # input max baseline along-n (Bn = Np x dn)
    #spacing = params.pos_n[2]-params.pos_n[1] # spacing for equally spaced platforms
    spacing = mean(diff(params.pos_n,dims=2)) # average spacing for unequally spaced platforms
    max_baseline_n = (maximum(params.pos_n)-minimum(params.pos_n)) + spacing
    println();println("input max baseline along-n: ",max_baseline_n)

    # theoretical resolution along-n
    range_s, range_g = Scene.lookangle_to_range(params.look_angle, params.p_t0_LLH[3], 0, earth_radius)
    res_theory_n = (c/params.fc)*range_s/p_mode/max_baseline_n
    println("theoretical resolution along-n: ",round(res_theory_n,digits=2))
    res_theory_array[i] = res_theory_n

    # theoretical resolution along-track
    mu = 3.986004418e14
    sc_speed = sqrt(mu./(params.p_t0_LLH[3]+earth_radius)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
    Lsa = sc_speed*params.SAR_duration + sc_speed/params.fp
    res_theory_s = (c/params.fc)*range_s/2/Lsa
    println("theoretical resolution along track: ",round(res_theory_s,digits=2))

    # Check consistency of input parameters
    paramsIsValid = UserParameters.validateInputParams(params)

    # Compute orbits time, position, and velocity
    orbit_time, orbit_pos, orbit_vel = Orbits.computeTimePosVel(params)

    # measured baselines
    refind = findall(params.pos_n.==0)[1][2]
    bperp, b_at, bnorm = Orbits.get_perp_baselines(orbit_pos, orbit_vel, params.look_angle, refind)
    avg_bperp = mean(bperp, dims=3)
    max_bperp = maximum(avg_bperp)
    println("measured max baseline along-n: ",round(max_bperp+spacing,digits=0))

    # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
    p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)

    # Create target/scene location
    targets_loc, targets_ref, Nt = Scene.construct_targets_str(params) # Nt: number of targets, targets: structure array containing target locations and reflectivities
    s_loc_3xN  = Scene.form3Dgrid_for(params.s_loc_1, params.s_loc_2, params.s_loc_3) # using 3 nested for loops
    t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, params) ## calculate avg heading from platform positions
    #t_xyz_3xN, s_xyz_3xN, avg_peg = Scene.convert_target_scene_coord_to_XYZ(s_loc_3xN, targets_loc, orbit_pos, orbit_vel, params) ## calculate avg heading from platform positions/velocities

    # Read number of platforms (todo: move into a struct)
    Np  = size(orbit_pos)[2] # number of platforms

    # Apply antenna pattern
    if params.include_antenna # calculate look angle (average over platforms and slow-time positions)
        Antenna.applyAntennaPattern!(targets_ref, p_xyz, orbit_vel, params)
    end

    # Generate range spread function (matched filter output)
    min_range, max_range = Geometry.find_min_max_range(t_xyz_3xN, p_xyz)
    Trx = 2*(max_range-min_range)/c + 5*params.pulse_length # s duration of RX window
    Srx, MF, ft, t_rx = RSF.ideal_RSF(Trx, params) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
    # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(params.pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing

    # Generate TomoSAR raw data
    ref_range = Geometry.distance(mean(t_xyz_3xN, dims=2), mean(mean(p_xyz,dims=2), dims=3)) # reference range (equal to slant_range in sch?)
    rawdata = Generate_Raw_Data.main_RSF_slowtime(t_xyz_3xN, p_xyz, Srx, t_rx, ref_range, targets_ref, params) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
    if params.enable_thermal_noise # adding random noise based on SNR after range (fast-time) processing
        rawdata = Error_Sources.random_noise(rawdata, params)
    end

    # Add phase error
    sync_osc_coeffs = repeat(params.sync_a_coeff_dB, Np)
    if params.enable_sync_phase_error
        rawdata = Error_Sources.synchronization_errors!(rawdata, slow_time, p_xyz, sync_osc_coeffs, params)
    end

    # Process raw data to generate image
    if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
        image_3D = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
    elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
        SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
        image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
    end

    # Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
    if params.left_right_look == "left";PSFcutdir=-1;elseif params.left_right_look == "right";PSFcutdir=1;end
    params.PSF_direction[3]=PSFcutdir*params.PSF_direction[3]
    scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)
    Plots.frame(anim)
    # Calculate point target performance metrics
    if size(t_xyz_3xN,2) == 1 # PSF related performance metrics are calculated when there is only one point target
        resolutions, PSLRs, ISLRs, loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
        res_measured_array[i]=resolutions[1]
        println("Resolutions: ",round.(resolutions,digits=2)," in scene axes units")
        println("Location Errors: ",round.(loc_errors,digits=2)," in scene axes units")
        println("PSLRs: ",round.(PSLRs,digits=2)," dB")
        println("ISLRs: ",round.(ISLRs,digits=2)," dB")
        println("PSF Peak Amplitude: ",round(maximum(20*log10.(image_3D)),digits=2)," dB")
        println("ratio of theoretical/measured resolution along-n: ",res_theory_n/resolutions)
    else
        @warn "PSF related performance metrics cannot be calculated for more than 1 target."
    end

    # Relative Radiometric Accuracy (amplitude difference between input 3D scene and output 3D image, max normalized to 1)
    inputscene_3D = Scene.generate_input_scene_3D(targets_ref, Nt, params)
    diff_image3D, mean_diff_image, std_diff_image = Performance_Metrics.relative_radiometric_accuracy(inputscene_3D, image_3D)
    println("Relative Radiometric Accuracy: Mean: ", round(mean_diff_image, digits=2),", Std: ",round(std_diff_image, digits=2)) # mean=0 & std_dev=0 means perfect result


    # Plots (1D PSF cuts are displayed by default in the performance.metrics module)
    if params.display_geometry || params.display_RSF_rawdata || params.display_input_scene || params.display_tomograms != 0
        include("../modules/plotting.jl")
        display_geometry_coord_txt=Plotting.coordinates(params.display_geometry_coord)
        ts_coord_txt=Plotting.coordinates(params.ts_coord_sys)
        if params.display_RSF_rawdata; Plotting.plot_RSF_rawdata(ft, t_rx, MF, Srx, Np, Nst, rawdata, params); end
        if params.display_geometry
            # convert platform and target locations to desired coordinate system
            p_loc,t_loc,s_loc=Geometry.convert_platform_target_scene_coordinates(Np,Nst,Nt,p_xyz,t_xyz_3xN,targets_loc,s_xyz_3xN,s_loc_3xN,avg_peg, params)
            Plotting.plot_geometry(orbit_time,orbit_pos,p_loc,t_loc,s_loc,display_geometry_coord_txt)
        end
        if params.display_input_scene; Plotting.plot_input_scene(inputscene_3D, ts_coord_txt, params);end
        if params.display_tomograms != 0;Plotting.plot_tomogram(image_3D, ts_coord_txt, scene_axis11, scene_axis22, scene_axis33, params);end
        Plots.frame(anim2)
        if params.display_input_scene; Plotting.plot_input_scene(diff_image3D, ts_coord_txt, params);end
    end

end

@save "platform_distributions.jld" pos_n_all
display(scatter(pos_n_all./1e3,leg=false,xlabel="iteration",ylabel="platform distribution (km)",title="Platform Positions",size=(1200,800),markersize=6,titlefont=22,xguidefontsize=18,yguidefontsize=18,xtickfontsize=15,ytickfontsize=15,xticks=1:Ntr))
savefig("platform_spacings.png")

gif(anim, "anim_1Dcuts.gif", fps = 1)
gif(anim2, "anim_2Dtomogram.gif", fps = 1)

xax=(init_spc.+((1:Ntr).-1).*spc_inc)./1e3
plot(xax,[res_theory_array res_measured_array] ,xaxis=("trials"),yaxis=("resolution (m)"),labels=permutedims(["theory","simulation"]))
savefig("resolution_vs_spacing.png")
