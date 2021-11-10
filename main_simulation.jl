include("modules/generate_raw_data.jl")
include("modules/process_raw_data.jl")
include("modules/geometry.jl")
include("modules/scene.jl")
#include("inputs/input_parameters_SCH_tree.jl")
#include("inputs/input_parameters_antenna_pattern_grid.jl")
#include("inputs/input_parameters_CR_nadirlooking.jl")
#include("inputs/input_parameters_CR_slantlooking.jl")
#include("inputs/input_parameters_CR_nadirlooking_tiltedcuts.jl")
include("inputs/input_parameters_CR_slantlooking_tiltedcuts.jl")
#include("inputs/input_parameters_CRs_cross.jl")
include("modules/range_spread_function.jl") # as RSF
include("modules/orbits.jl")
include("modules/sync.jl")
include("modules/error_sources.jl")
include("modules/performance_metrics.jl")
include("modules/antenna.jl")
include("modules/simsetup.jl")
using NCDatasets
using Statistics
using Dates
## PLATFORM LOCATIONS and HEADINGS
orbit_dataset=Dataset("inputs/"*orbit_filename) # Read orbits data in NetCDF format
t12_orbits=orbit_dataset["time"][1:2] # first two time samples
dt_orbits=t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
orbit_time_index=(Int(round(SAR_start_time/dt_orbits))+1:1:Int(round((SAR_start_time+SAR_duration)/dt_orbits))+1) # index range for orbit times for time interval of interest
orbit_time=orbit_dataset["time"][orbit_time_index] # read in time data
orbit_pos_ECI=1e3*orbit_dataset["position"][:,:,orbit_time_index] # read in position data, 3 x Np x Nt
orbit_vel_ECI=1e3*orbit_dataset["velocity"][:,:,orbit_time_index] # read in velocity data, 3 x Np x Nt (used optionally in avg peg and heading calculation)
try #does file have dcm already?
    global dcm=orbit_dataset["dcm"];
catch #if not generate from Orbits
    dv = orbit_dataset.attrib["epoch"];
    local epoch = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
    global dcm = Orbits.eci_dcm(orbit_time, epoch);
end
#orbit_pos=orbit_pos_ECI
#orbit_vel=orbit_vel_ECI
#orbit_pos=Orbits.ecef_orbitpos(orbit_pos_ECI,dcm)# convert ECI to ECEF
orbit_pos,orbit_vel=Orbits.ecef_orbitpos(orbit_pos_ECI,orbit_vel_ECI,dcm) # ECI to ECEF TODO velocity conversion function not ready yet
slow_time=(SAR_start_time:1/fp:SAR_start_time+SAR_duration) # create slow time axis
if length(slow_time)==1;p_xyz=orbit_pos
else;p_xyz=Orbits.interp_orbit(orbit_time,orbit_pos,slow_time);end # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
Np=size(orbit_pos)[2] # number of platforms
Nst=size(slow_time)[1] # number of slow-time samples (pulses processed)
## TARGET/SCENE LOCATIONS
targets,Nt=Scene.construct_targets_str(target_pos_mode,t_loc_1,t_loc_2,t_loc_3,t_ref) # Nt: number of targets, targets: structure array containing target locations and reflectivities
targets_loc=zeros(3,Nt);for i=1:Nt;targets_loc[:,i]=targets[i].loc;end # 3xN
s_loc_3xN=Scene.form3Dgrid_for(s_loc_1,s_loc_2,s_loc_3) # using 3 nested for loops
if ts_coord_sys=="XYZ" || ts_coord_sys=="LLH";look_angle=[];end
#t_xyz_3xN,s_xyz_3xN,avg_peg=Scene.convert_target_scene_coord_to_XYZ(ts_coord_sys,s_loc_3xN,targets_loc,orbit_pos,look_angle,earth_radius,earth_eccentricity) ## calculate avg heading from platform positions
t_xyz_3xN,s_xyz_3xN,avg_peg=Scene.convert_target_scene_coord_to_XYZ(ts_coord_sys,s_loc_3xN,targets_loc,orbit_pos,orbit_vel,look_angle,earth_radius,earth_eccentricity) # calculate avg heading from platform velocities
## TARGET REFLECTIVITIES
targets_ref=zeros(1,Nt);for i=1:Nt;targets_ref[i]=targets[i].ref;end
## ANTENNA PATTERN
if include_antenna # calculate look angle (average over platforms and slow-time positions)
    avg_p_xyz=reshape(mean(mean(p_xyz,dims=2),dims=3),3)
    avg_p_vel=reshape(mean(mean(orbit_vel,dims=2),dims=3),3)
    if ts_coord_sys=="SCH"
        look_ang=look_angle
    elseif ts_coord_sys=="XYZ" || ts_coord_sys=="LLH"
        platform_heights=zeros(Np);slant_ranges=zeros(Np)
        avg_t_xyz=mean(t_xyz_3xN,dims=2) # average target location in XYZ
        avg_rs=Geometry.distance(avg_t_xyz,avg_p_xyz) # average slant range
        for i=1:Np
            p_xyz_i=p_xyz[:,i,:] # p_xyz: 3 x Np x Nst
            p_xyz_i=reshape(p_xyz_i,3,Nst) # p_xyz: 3 x Nst
            p_LLH=Geometry.xyz_to_geo(p_xyz_i)
            platform_heights[i]=mean(p_LLH[3,:]) # average platform heights over slow-time for each platform
        end
        avg_p_h=mean(platform_heights) # average platform height over platforms and slow-time
        if avg_rs<avg_p_h;avg_rs=avg_p_h;end
        avg_rg,look_ang=Scene.slantrange_to_lookangle(earth_radius,avg_rs,avg_p_h,0) # assuming target height is 0 (negligible effect), look_ang: average look angle
    end
    vgrid = Antenna.AntGrid("inputs/darts_ant_03192021.nc") # read in vpol grid, takes time to load
    ant = SimSetup.sc_ant(vgrid); #create antenna structure, additional arguments are rotation and origin
    sc = SimSetup.spacecraft(avg_p_xyz, Float64.(avg_p_vel), ant = ant, look_angle = look_ang, side = "right"); ##create spacecraft structure; ant, look_angle, side are optional
    co_pol,cross_pol = SimSetup.interpolate_pattern(sc, t_xyz_3xN);#inteprolate pattern (cp:co-pol, xp: cross-pol), outputs are 1xNt complex vectors
    targets_ref=targets_ref.*transpose(co_pol).^2/maximum(abs.(co_pol))^2 #TODO separate TX and RX, include range effect?
    if target_pos_mode=="grid" # plotting projected pattern only for grid type target distribution
        projected_pattern_3D=Scene.convert_image_1xN_to_3D(abs.(co_pol),length(t_loc_1),length(t_loc_2),length(t_loc_3))#take magnitude and reshape to 3D
        include("modules/plotting.jl");coords_txt=Plotting.coordinates(ts_coord_sys)
        Nt_1=length(t_loc_1)
        Nt_2=length(t_loc_2)
        Nt_3=length(t_loc_3)
        using Plots;gr()
        if Nt_2>1 && Nt_1>1;display(heatmap(t_loc_2,t_loc_1, 20*log10.(projected_pattern_3D[:,:,1]), ylabel=coords_txt[1],xlabel=coords_txt[2],title = "Antenna Pattern Projected on Targets (V-copol)", fill=true,size=(1600,900)));end #, clim=(-80,40),aspect_ratio=:equal
        if Nt_3>1 && Nt_2>1;display(heatmap(t_loc_3,t_loc_2, 20*log10.(projected_pattern_3D[1,:,:]),ylabel=coords_txt[2],xlabel=coords_txt[3],title = "Antenna Pattern Projected on Targets (V-copol)", fill=false,size=(1600,900)));end #, clim=(-80,40),aspect_ratio=:equal
        if Nt_3>1 && Nt_1>1;display(heatmap(t_loc_3,t_loc_1, 20*log10.(projected_pattern_3D[:,1,:]),ylabel=coords_txt[1],xlabel=coords_txt[3],title = "Antenna Pattern Projected on Targets (V-copol)", fill=false,size=(1600,900)));end #, clim=(-80,40),aspect_ratio=:equal
    end
end
## RANGE SPREAD FUNCTION (matched filter output)
min_range,max_range=Geometry.find_min_max_range(t_xyz_3xN,p_xyz)
Trx=2*(max_range-min_range)/c+5*pulse_length # s duration of RX window
if enable_fast_time # matched filter gain is included in Srx
    Srx,MF,ft,t_rx=RSF.ideal_RSF(pulse_length,Δt,bandwidth,Trx) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
    # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing
end
## GENERATE RAW DATA
ref_range=Geometry.distance(mean(t_xyz_3xN,dims=2),mean(mean(p_xyz,dims=2),dims=3)) # reference range (equal to slant_range in sch?)
if enable_fast_time # with fastime and slowtime; matched filter gain is included in Srx
    rawdata=Generate_Raw_Data.main_RSF_slowtime(t_xyz_3xN,p_xyz,mode,tx_el,fc,Srx,t_rx,ref_range,targets_ref) # rawdata is a: 3D array of size Nst x Np x Nft (SAR/SIMO), 4D array of size Nst x Np(RX) x Np(TX) x Nft (MIMO)
else # without fastime, with slowtime; matched filter gain is included inside the function
    rawdata=Generate_Raw_Data.main_noRSF_slowtime(t_xyz_3xN,p_xyz,mode,tx_el,fc,targets_ref) # rawdata is a: 2D array of size Nst x Np (SAR/SIMO), 3D array of size Nst x Np(RX) x Np(TX) (MIMO)
end
if !enable_fast_time;SNR=SNR*pulse_length*bandwidth;end # SNR increases after matched filter
if enable_thermal_noise;rawdata=Error_Sources.random_noise(rawdata,SNR,enable_fast_time,mode);end # adding random noise based on SNR after range (fast-time) processing
## PROCESS RAW DATA TO GENERATE IMAGE
image_1xN=Process_Raw_Data.main_RSF_slowtime(rawdata,s_xyz_3xN,p_xyz,mode,tx_el,fc,t_rx,ref_range)
Ns_1=length(s_loc_1);Ns_2=length(s_loc_2);Ns_3=length(s_loc_3)
image_3D=Scene.convert_image_1xN_to_3D(image_1xN,Ns_1,Ns_2,Ns_3)
## PERFORMANCE METRICS
# PSF metrics
include("modules/performance_metrics.jl")
if size(t_xyz_3xN,2)==1 # PSF related performance metrics are calculated when there is only one point target
    target_index1=findall(t_loc_1 .==s_loc_1)
    target_index2=findall(t_loc_2 .==s_loc_2)
    target_index3=findall(t_loc_3 .==s_loc_3)
    if isempty(target_index1) || isempty(target_index2) || isempty(target_index3)
        println("PSF related performance metrics cannot be calculated since target is not inside the scene!")
        PSF_metrics=false
    else
        PSF_metrics=true
        target_location=[t_loc_1 t_loc_2 t_loc_3] # point target location
        resolutions,PSLRs,ISLRs,loc_errors,scene_axis11,scene_axis22,scene_axis33=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_loc_1,s_loc_2,s_loc_3,PSF_image_point,PSF_cuts,PSF_direction) # resolutions in each of the 3 axes
    end
else
    PSF_metrics=false
    println("PSF related performance metrics cannot be calculated since there are more than 1 targets!")
end
if PSF_metrics
    println("Resolutions: ",round.(resolutions,digits=8)," in scene axes units")
    println("Location Errors: ",round.(loc_errors,digits=8)," in scene axes units")
    println("PSLRs: ",round.(PSLRs,digits=2)," dB")
    println("ISLRs: ",round.(ISLRs,digits=2)," dB")
    println("PSF Peak Amplitude: ",round(maximum(20*log10.(image_3D)),digits=2)," dB")
end
# Relative Radiometric Accuracy (amplitude difference between input 3D scene and output 3D image, max normalized to 1)
inputscene_3D=Scene.generate_input_scene_3D(s_loc_1,s_loc_2,s_loc_3,t_loc_1,t_loc_2,t_loc_3,targets_ref,Nt,target_pos_mode)
diff_image3D,mean_diff_image,std_diff_image=Performance_Metrics.relative_radiometric_accuracy(inputscene_3D,image_3D)
println("Relative Radiometric Accuracy: Mean: ",round(mean_diff_image,digits=2),", Std: ",round(std_diff_image,digits=2)) # mean=0 & std_dev=0 means perfect result
## PLOTS (1D PSF cuts are displayed by default in the performance.metrics module)
if display_geometry || display_RSF_rawdata || display_input_scene || display_tomograms!=0
    include("modules/plotting.jl")
    display_geometry_coord_txt=Plotting.coordinates(display_geometry_coord)
    ts_coord_txt=Plotting.coordinates(ts_coord_sys)
    if display_RSF_rawdata;Plotting.plot_RSF_rawdata(enable_fast_time,mode,ft,t_rx,MF,Srx,Np,Nst,rawdata);end
    if display_geometry
        # convert platform and target locations to desired coordinate system
        p_loc,t_loc,s_loc=Geometry.convert_platform_target_scene_coordinates(Np,Nst,Nt,p_xyz,t_xyz_3xN,targets_loc,s_xyz_3xN,s_loc_3xN,avg_peg,display_geometry_coord,ts_coord_sys)
        Plotting.plot_geometry(orbit_time,orbit_pos,p_loc,t_loc,s_loc,display_geometry_coord_txt)
    end
    if display_input_scene;Plotting.plot_input_scene(inputscene_3D,s_loc_1,s_loc_2,s_loc_3,ts_coord_txt);end
    if display_tomograms!=0;Plotting.plot_tomogram(PSF_image_point,display_tomograms,image_1xN,image_3D,s_loc_1,s_loc_2,s_loc_3,s_loc_3xN,t_loc_1,t_loc_2,t_loc_3,ts_coord_txt,mode,scene_axis11,scene_axis22,scene_axis33,PSF_cuts,PSF_metrics);end
    if display_input_scene;Plotting.plot_input_scene(diff_image3D,s_loc_1,s_loc_2,s_loc_3,ts_coord_txt);end
end
