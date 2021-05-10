include("modules/generate_raw_data.jl")
include("modules/process_raw_data.jl")
include("modules/geometry.jl")
include("modules/scene.jl")
#include("inputs/input_parameters_LLH_nadirlooking.jl")
#include("inputs/input_parameters_LLH_slantlooking.jl")
include("inputs/input_parameters_SCH_lookangle.jl")
include("modules/range_spread_function.jl") # as RSF
include("modules/orbits.jl")
include("modules/sync.jl")
include("modules/error_sources.jl")

using NCDatasets
using Statistics

## PLATFORM LOCATIONS
orbit_filename="../darts-simtool/inputs/orbitOutput_082020.nc" # position in km, time in sec
# orbit_dataset=Dataset("inputs/"*orbit_filename) # Read orbits data in NetCDF format
orbit_dataset=Dataset(orbit_filename) # Read orbits data in NetCDF format

t12_orbits=orbit_dataset["time"][1:2] # first two time samples
dt_orbits=t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
orbit_time_index=(Int(round(SAR_start_time/dt_orbits))+1:1:Int(round((SAR_start_time+SAR_duration)/dt_orbits))+1) # index range for orbit times for time interval of interest
orbit_time=orbit_dataset["time"][orbit_time_index] # read in time data
orbit_pos=orbit_dataset["position"][:,:,orbit_time_index] # read in position data #TODO convert ECI to ECEF?
slow_time=(SAR_start_time:1/fp:SAR_start_time+SAR_duration) # create slow time axis
orbit_pos_interp=Orbits.interp_orbit(orbit_time,orbit_pos,slow_time) # interpolate orbit to slow time
p_xyz=1e3*orbit_pos_interp # convert km to m
Np=size(orbit_pos)[2] # number of platforms
Nst=size(slow_time)[1] # number of slow-time samples (pulses processed)
## TARGET/SCENE LOCATIONS
targets,Nt=Scene.construct_targets_str(target_pos_mode,t_loc_1,t_loc_2,t_loc_3,t_ref) # Nt: number of targets, targets: structure array containing target locations and reflectivities
targets_loc=zeros(3,Nt);for i=1:Nt;targets_loc[:,i]=targets[i].loc;end # 3xN
s_loc_3xN=Scene.form3Dgrid_for(s_loc_1,s_loc_2,s_loc_3) # using 3 nested for loops
t_xyz_3xN,s_xyz_3xN=Scene.convert_target_scene_coord_to_XYZ(ts_coord_sys,s_loc_3xN,targets_loc,p_xyz,look_angle,p_avg_heading,earth_radius,earth_eccentricity)
## TARGET REFLECTIVITIES
targets_ref=zeros(1,Nt);for i=1:Nt;targets_ref[i]=targets[i].ref;end

## RANGE SPREAD FUNCTION (matched filter output)
min_range,max_range=Geometry.find_min_max_range(t_xyz_3xN,p_xyz)
Trx=2*(max_range-min_range)/c+2*pulse_length # s duration of RX window
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

## Add Sync effects
include("inputs/input_parameters_sync.jl")
@time rawdata_sync = Error_Sources.synchronization_errors(rawdata,slow_time,orbit_pos_interp,enable_fast_time,parameters)

## PROCESS RAW DATA TO GENERATE IMAGE
#image_3xN=Process_Raw_Data.main(rawdata,s_xyz_3xN,p_xyz_grid,mode,tx_el,fc) # without fastime, without slowtime
#image_3xN=Process_Raw_Data.main_RSF(rawdata,s_xyz_3xN,p_xyz,mode,tx_el,fc,t_rx,ref_range)  # with fastime, without slowtime
if enable_fast_time # with fastime, with slowtime
    image_1xN=Process_Raw_Data.main_RSF_slowtime(rawdata_sync,s_xyz_3xN,p_xyz,mode,tx_el,fc,t_rx,ref_range)
else # without fastime, with slowtime
    image_1xN=Process_Raw_Data.main_noRSF_slowtime(rawdata_sync,s_xyz_3xN,p_xyz,mode,tx_el,fc)
end
Ns_1=length(s_loc_1);Ns_2=length(s_loc_2);Ns_3=length(s_loc_3)
image_3D=Scene.convert_image_1xN_to_3D(image_1xN,Ns_1,Ns_2,Ns_3)## PERFORMANCE METRICS

## PERFORMANCE METRICS
# PSF metrics
if size(t_xyz_3xN,2)==1 # PSF related performance metrics are calculated when there is only one point target
    target_index1=findall(t_loc_1 .==s_loc_1)
    target_index2=findall(t_loc_2 .==s_loc_2)
    target_index3=findall(t_loc_3 .==s_loc_3)
    if isempty(target_index1) || isempty(target_index2) || isempty(target_index3)
        show("PSF related performance metrics cannot be calculated since target is not inside the scene!")
        PSF_metrics=false
    else
        include("modules/performance_metrics.jl")
        PSF_metrics=true
        target_location=[t_loc_1 t_loc_2 t_loc_3] # point target location
        resolutions,PSLRs,ISLRs,loc_errors=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_loc_1,s_loc_2,s_loc_3,PSF_peak_target) # resolutions in each of the 3 axes
    end
else
    PSF_metrics=false
    show("PSF related performance metrics cannot be calculated since there are more than 1 targets!")
end
if PSF_metrics
    println("Resolutions: ",round.(resolutions,digits=8)," in scene axes units")
    println("Location Errors: ",round.(loc_errors,digits=8)," in scene axes units")
    println("PSLRs: ",round.(PSLRs,digits=2)," dB")
    println("ISLRs: ",round.(ISLRs,digits=2)," dB")
    println("PSF Peak Amplitude: ",round(maximum(20*log10.(image_3D)),digits=2)," dB")
end

## PLOTS (1D PSF cuts are displayed by default in the performance.metrics module)
if display_geometry || display_RSF_rawdata || display_tomograms!=0
    include("modules/plotting.jl")
    coords=Plotting.coordinates(ts_coord_sys)
    if display_RSF_rawdata;Plotting.plot_RSF_rawdata(enable_fast_time,mode,ft,t_rx,MF,Srx,Np,Nst,rawdata);end
    if display_geometry;Plotting.plot_geometry(orbit_time,orbit_pos,p_xyz,t_xyz_3xN,s_loc_3xN,s_xyz_3xN,coords);end
    if display_tomograms!=0;Plotting.plot_tomogram(display_tomograms,image_1xN,image_3D,s_loc_1,s_loc_2,s_loc_3,s_loc_3xN,s_xyz_3xN,coords);end
end