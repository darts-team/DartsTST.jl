include("../../modules/generate_raw_data.jl")
include("../../modules/process_raw_data.jl")
include("../../modules/geometry.jl")
include("../../modules/scene.jl")
include("../../modules/range_spread_function.jl") # as RSF
include("../../modules/orbits.jl")
include("../../modules/sync.jl")
include("../../modules/error_sources.jl")
include("../../modules/performance_metrics.jl")
include("../../modules/antenna.jl")
include("../../modules/simsetup.jl")
include("../../modules/user_parameters.jl")
using NCDatasets
using Statistics
using Parameters
using Dates
using StaticArrays
using .UserParameters
c = 299792458 #TODO does not work without redefining c here
numApertures= 20 #sets number of aperture time points to look at
dtime = 120 # sec, time resolution between apertures -- Looking for total time dtime*setNumApertures to be roughly 100 min (?)


# Define user parameters
params = UserParameters.inputParameters(PSF_image_point=1,PSF_cuts=1,display_tomograms=0,user_defined_orbit=0,include_antenna=false)

# Check consistency of input parameters
paramsIsValid = UserParameters.validateInputParams(params)


mode = params.mode
fc   = params.fc
#define p_𝛿 for analytical resolution calc
if mode == 1
    p_𝛿 = 2
elseif mode == 2
    p_𝛿 = 1
elseif mode ==3
    p_𝛿 = 1.38
end

aperture_time_vec = 0 : dtime : numApertures*dtime - dtime

# initialize arrays for saving data

# Save the data for each location
#-- target position (convert to LLH)
#-- platform position(s) (convert to LLH) (or take mean platform position?)
#--vertical resolution, along+cross track res are probably good too.
#--ambiguity location? PSLR? ISLR?

platform_pos_ecef       = Array{Float64,2}(undef,3,numApertures) # average position of platforms
target_pos_ecef         = Array{Float64,2}(undef,3,numApertures)
# target_pos_llh          = Array{Float64,2}(undef,3,numApertures)
orbit_resolutions_ndB   = Array{Float64}(undef,1,numApertures)
orbit_ISLRs             = Array{Float64}(undef,1,numApertures)
orbit_PSLRs             = Array{Float64}(undef,1,numApertures)
orbit_res_analytical    = Array{Float64}(undef,1,numApertures)
orbit_perp_baseline     = Array{Float64}(undef,1,numApertures)
orbit_ranges            = Array{Float64}(undef,1,numApertures)
orbit_peaks             = Array{Float64}(undef,1,numApertures)

for naperture = 1 : numApertures
    
    ## find the orbit positions at the current aperture
    start_time = aperture_time_vec[naperture]
    
    global params = UserParameters.inputParameters(PSF_image_point=1,display_tomograms=0,user_defined_orbit=0,include_antenna=false,SAR_start_time = start_time,look_angle=30,res_dB = 4.93)
    # Compute orbits time, position, and velocity
    orbit_time, orbit_pos, orbit_vel = Orbits.computeTimePosVel(params)
    
    # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
    p_xyz, Nst, slow_time = Orbits.interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)
    v_xyz = Orbits.interp_orbit(orbit_time,orbit_vel,slow_time)
    
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
    
    # Process raw data to generate image
    if params.processing_steps === :bp3d # 1-step processing TODO do we need this option?
        image_3D = Process_Raw_Data.main_SAR_tomo_3D(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
    elseif params.processing_steps === :bp2d3d # 2-step processing, first SAR (along-track), then tomographic
        SAR_images_3D = Process_Raw_Data.SAR_processing(rawdata, s_xyz_3xN, p_xyz, t_rx, ref_range, params)
        image_3D = Process_Raw_Data.tomo_processing_afterSAR(SAR_images_3D)
    end

    # Take 1D cuts from the 3D tomogram and plot the cuts (for multiple targets cuts are taken from the center of the scene)
    scene_axis11, scene_axis22, scene_axis33, image_1D_1, image_1D_2, image_1D_3, scene_res = Scene.take_1D_cuts(image_3D, params)

    # Calculate point target performance metrics for both 3dB and 5dB points
    #3 dB was walready written in above
    resolutions_ndB,PSLRs,ISLRs,loc_errors = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)
    
    # global params = UserParameters.inputParameters(PSF_image_point=1,display_tomograms=0,user_defined_orbit=0,include_antenna=false,SAR_start_time = start_time,res_dB=5) # overwrites params
    # resolutions_5dB,PSLRs,ISLRs,loc_errors,scene_axis11,scene_axis22,scene_axis33 = Performance_Metrics.computePTPerformanceMetrics(image_1D_1, image_1D_2, image_1D_3, scene_res, params)

    
    # Save the data for each location
    #-- target position (convert to LLH)
    #-- platform position(s) (convert to LLH) (or take mean platform position?)
    #--vertical resolution, along+cross track res are probably good too.
    #--ambiguity location? PSLR? ISLR?
    target_pos_ecef[:,naperture] = t_xyz_3xN
    # llh = Geometry.ecef_to_geodetic([t_xyz_3xN[1];t_xyz_3xN[2];t_xyz_3xN[3]]) # LLH in rad,rad, scene units
    # target_pos_llh[:,naperture] = [rad2deg(llh[1]) rad2deg(llh[2]) llh[3]]
    orbit_resolutions_ndB[naperture] = resolutions_ndB
    platform_pos_ecef[:,naperture] = reshape(mean(mean(p_xyz,dims=2),dims=3),3)
    orbit_ISLRs[naperture] = ISLRs
    orbit_PSLRs[naperture] = PSLRs

    #find analytical resolution
    λ = params.λ
    # find max baseline
    baselines, b_at, bnorm=Orbits.get_perp_baselines(p_xyz,v_xyz,params.look_angle)
    time_index=1:length(slow_time)
    maxbaselines=reshape(maximum(maximum(baselines,dims=1),dims=2),length(slow_time))
    maxbaseline,maxb_time_ind=findmax(maxbaselines)

    L_n = maxbaseline
    r_0 = ref_range #TODO is this correct?
    w_n = 1 # assume no windowing
    orbit_res_analytical[naperture] = w_n * (λ * r_0) / (p_𝛿 * L_n)
    orbit_perp_baseline[naperture]  = maxbaseline
    orbit_ranges[naperture]         = r_0
    orbit_peaks[naperture]          = maximum(image_3D)

    println("$naperture out of  $numApertures")
end#numApertures