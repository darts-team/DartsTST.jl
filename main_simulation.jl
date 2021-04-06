include("modules/generate_raw_data.jl")
include("modules/process_raw_data.jl")
include("modules/geometry.jl")
include("modules/scene.jl")
#include("inputs/input_parameters_nadirlooking.jl")
include("inputs/input_parameters_slantlooking.jl")
include("modules/range_spread_function.jl") # as RSF
include("modules/orbits.jl")
include("modules/sync.jl")
include("modules/error_sources.jl")
include("modules/performance_metrics.jl")
using NCDatasets
using Statistics
## RANGE SPREAD FUNCTION (matched filter output)
if enable_fast_time # matched filter gain is included in Srx
    Srx,MF,ft,t_rx=RSF.ideal_RSF(pulse_length,Δt,bandwidth,Trx) # Srx: RX window with MF centered, MF: ideal matched filter output (range spread function, RSF) for LFM pulse, ft: fast-time axis for MF, t_rx: RX window
    # Srx,MF,ft,t_rx=RSF.non_ideal_RSF(pulse_length,Δt,bandwidth,Trx,SFR,window_type) # TODO non-ideal RSF for LFM pulse with system complex frequency response (SFR) and fast-time windowing
end
## PLATFORM LOCATIONS
orbit_dataset=Dataset("inputs/"*orbit_filename) # Read orbits data in NetCDF format
t12_orbits=orbit_dataset["time"][1:2] # first two time samples
dt_orbits=t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
orbit_time_index=(Int(round(SAR_start_time/dt_orbits))+1:1:Int(round((SAR_start_time+SAR_duration)/dt_orbits))+1) # index range for orbit times for time interval of interest
orbit_time=orbit_dataset["time"][orbit_time_index] # read in time data
orbit_pos=orbit_dataset["position"][:,:,orbit_time_index] # read in position data
#TODO convert ECI to ECEF
slow_time=(SAR_start_time:1/fp:SAR_start_time+SAR_duration) # create slow time axis
orbit_pos_interp=Orbits.interp_orbit(orbit_time,orbit_pos,slow_time) # interpolate orbit to slow time
p_xyz=1e3*orbit_pos_interp # convert km to m
Np=size(orbit_pos)[2] # number of platforms
Nst=size(slow_time)[1] # number of slow-time samples (pulses processed)
## TARGET LOCATIONS
t_geo_grid=Scene.form3Dgrid_for(t_θ,t_ϕ,t_h) # using 3 nested for loops
#t_geo_grid=Scene.form3Dgrid_array(t_θ,t_ϕ,t_h) # using array processing
t_xyz_grid=Geometry.geo_to_xyz(t_geo_grid,earth_radius,earth_eccentricity)
## GENERATE RAW DATA
#rawdata=Generate_Raw_Data.main(t_xyz_grid,p_xyz_grid,mode,tx_el,fc) # without RSF
ref_range=Generate_Raw_Data.distance(mean(t_xyz_grid,dims=2),mean(mean(p_xyz,dims=2),dims=3)) # reference range
#rawdata=Generate_Raw_Data.main_RSF(t_xyz_grid,p_xyz,mode,tx_el,fc,Srx,t_rx,ref_range) # with fasttime, without slowtime #TODO use structure as input
if enable_fast_time # with fastime and slowtime; matched filter gain is included in Srx
    rawdata=Generate_Raw_Data.main_RSF_slowtime(t_xyz_grid,p_xyz,mode,tx_el,fc,Srx,t_rx,ref_range)
else # without fastime, with slowtime; matched filter gain is included inside the function
    rawdata=Generate_Raw_Data.main_noRSF_slowtime(t_xyz_grid,p_xyz,mode,tx_el,fc)
end
if !enable_fast_time
    SNR=SNR*pulse_length*bandwidth # SNR increases after matched filter
end
if enable_thermal_noise
    rawdata=Error_Sources.random_noise(rawdata,SNR,enable_fast_time,mode)
end
## IMAGE SCENE
Ns_θ=length(s_θ)
Ns_ϕ=length(s_ϕ)
Ns_h=length(s_h)
s_geo_grid=Scene.form3Dgrid_for(s_θ,s_ϕ,s_h) # using 3 nested for loops
#s_geo_grid=Scene.form3Dgrid_array(s_θ,s_ϕ,s_h) # using array processing
s_xyz_grid=Geometry.geo_to_xyz(s_geo_grid,earth_radius,earth_eccentricity)
## PROCESS RAW DATA TO GENERATE IMAGE
#image_3xN=Process_Raw_Data.main(rawdata,s_xyz_grid,p_xyz_grid,mode,tx_el,fc) # without fastime, without slowtime
#image_3xN=Process_Raw_Data.main_RSF(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc,t_rx,ref_range)  # with fastime, without slowtime
if enable_fast_time # with fastime, with slowtime
    image_1xN=Process_Raw_Data.main_RSF_slowtime(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc,t_rx,ref_range)
else # without fastime, with slowtime
    image_1xN=Process_Raw_Data.main_noRSF_slowtime(rawdata,s_xyz_grid,p_xyz,mode,tx_el,fc)
end
image_3D=Scene.convert_image_1xN_to_3D(image_1xN,Ns_θ,Ns_ϕ,Ns_h)
## PERFORMANCE METRICS
# PSF metrics
if size(t_xyz_grid)[2]==1 # PSF related performance metrics are calculated when there is only one point target
    target_index1=findall(t_θ .==s_θ)
    target_index2=findall(t_ϕ .==s_ϕ)
    target_index3=findall(t_h .==s_h)
    if isempty(target_index1) || isempty(target_index2) || isempty(target_index3)
        show("PSF related performance metrics cannot be calculated since target is not inside the scene!")
        PSF_metrics=false
    else
        PSF_metrics=true
        target_location=[t_θ t_ϕ t_h] # point target location
        resolutions,PSLRs,ISLRs,loc_errors=Performance_Metrics.PSF_metrics(image_3D,res_dB,target_location,s_θ,s_ϕ,s_h,PSF_peak_target) # resolutions in each of the 3 axes
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
## PLOTS
# 1D PSF cuts are displayed by default in the performance.metrics module
if display_geometry || display_RSF_rawdata || display_tomograms!=0
    include("modules/plotting.jl")
    orbit_pos_all=reshape(p_xyz,3,Np*Nst) # platform positions in xyz; for each platform, its position at each pulse (PRI) is plotted; output loops over platforms first, then slow-time
    Plotting.main(enable_fast_time,display_geometry,display_RSF_rawdata,display_tomograms,mode,rawdata,image_3D,image_1xN,ft,t_rx,MF,Srx,bandwidth,pulse_length,orbit_time,orbit_pos,orbit_pos_all,t_xyz_grid,s_xyz_grid,s_geo_grid,p_xyz,s_θ,s_ϕ,s_h)
end
