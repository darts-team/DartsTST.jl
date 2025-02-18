module Orbits

#include("geometry.jl")

using ..Geometry

using SatelliteToolbox
using NCDatasets
using Dates
using LinearAlgebra
using Dierckx
using Parameters
using Plots

include("../modules/scene.jl")

function interpolateOrbitsToSlowTime(orbit_time, orbit_pos, params)
    @unpack SAR_start_time, SAR_duration, fp = params

    # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
    slow_time = SAR_start_time : 1/fp : SAR_start_time+SAR_duration # create slow time axis

    Nst = size(slow_time)[1] # number of slow-time samples (pulses processed)

    if Nst == 1;
        p_xyz = orbit_pos
    else
        p_xyz = Orbits.interp_orbit(orbit_time,orbit_pos,slow_time)
    end # interpolate orbit to slow time, 3 x Np x Nst, convert km to m

    return p_xyz, Nst, slow_time
end

function computeTimePosVel(params)
    @unpack SAR_start_time, dt_orbits, SAR_duration, user_defined_orbit, pos_n,
    Torbit, p_t0_LLH, p_heading, look_angle, display_custom_orbit, orbit_filename, left_right_look, orbit_file_coords = params

    ## PLATFORM LOCATIONS and HEADINGS
    #if left_right_look == "left";N_dir=1;elseif left_right_look == "right";N_dir=-1;end
    if left_right_look == "left";C_dir=-1;elseif left_right_look == "right";C_dir=1;end
    if user_defined_orbit==1 # orbits from file
        orbit_dataset=Dataset("../inputs/"*orbit_filename) # Read orbits data in NetCDF format
        t12_orbits=orbit_dataset["time"][1:2] # first two time samples
        dt_orbits=t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
        orbit_time_index=(Int(round(SAR_start_time/dt_orbits))+1:1:Int(ceil((SAR_start_time+SAR_duration)/dt_orbits))+1) # index range for orbit times for time interval of interest
        orbit_time=orbit_dataset["time"][orbit_time_index] # read in time data
        if orbit_file_coords=="ECI"
        orbit_pos_ECI=1e3*orbit_dataset["position"][:,:,orbit_time_index] # read in position data, 3 x Np x Nt
        orbit_vel_ECI=1e3*orbit_dataset["velocity"][:,:,orbit_time_index] # read in velocity data, 3 x Np x Nt (used optionally in avg peg and heading calculation)
        #try #does file have dcm already?
        #    global dcm=orbit_dataset["dcm"];
        #catch #if not generate from Orbits
            dv = orbit_dataset.attrib["epoch"];
            if typeof(dv) == String
                temp = dv
                dv= []
                dv = [parse(Int,temp[1:4]);parse(Int,temp[6:7]);parse(Int,temp[9:10]);parse(Int,temp[12:13]);parse(Int,temp[15:16]);parse(Int,temp[18:19])];
            end
            local epoch = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
            global dcm = Orbits.eci_dcm(orbit_time, epoch);
        #end
        #orbit_pos=Orbits.ecef_orbitpos(orbit_pos_ECI,dcm)# convert ECI to ECEF
        orbit_pos,orbit_vel=Orbits.ecef_orbitpos(orbit_pos_ECI,orbit_vel_ECI,dcm) # ECI to ECEF
        elseif orbit_file_coords=="ECEF"
            orbit_pos=1e3*orbit_dataset["position"][:,:,orbit_time_index] # read in position data, 3 x Np x Nt
            orbit_vel=1e3*orbit_dataset["velocity"][:,:,orbit_time_index] # read in velocity data, 3 x Np x Nt (used optionally in avg peg and heading calculation)
        else
            @error(raw"Unsupported orbit file coordinate system: Should be 'ECI' or 'ECEF'")
        end

    elseif user_defined_orbit==2 # user defined, TCN option
        #pos_T = zeros(1,length(pos_n)) # no along-track spacings
        #pos_C =   C_dir * pos_n * cosd(look_angle)
        #pos_N =  -1 * pos_n * sind(look_angle)
        #pos_C = 1 * pos_n * cosd(look_angle)
        #pos_N = N_dir * pos_n * sind(look_angle)
        #pos_TCN = [pos_T;pos_C;pos_N]
        pos_XYZ=Geometry.geo_to_xyz(p_t0_LLH)
        orbit_time_all=-Torbit/2:dt_orbits:Torbit/2
        #orbit_pos_all,orbit_vel_all=Orbits.make_orbit(pos_XYZ,p_heading,pos_TCN,orbit_time_all)
        orbit_pos_all,orbit_vel_all = Orbits.make_orbit_new(pos_XYZ, p_heading, pos_n, orbit_time_all, look_angle, 0.0, left_right_look)
    end
    if user_defined_orbit==2
        if display_custom_orbit  #plot orbit on surface sphere
            lats=-90:1:90;lons=-180:1:180;hgts=0; # background spherical grid on surface
            spherical_grid=Geometry.geo_to_xyz(Scene.form3Dgrid_for(lats,lons,hgts)); #create grid in LLH and convert to XYZ
            plotly();scatter(spherical_grid[1,:]/1e3,spherical_grid[2,:]/1e3,spherical_grid[3,:]/1e3,xlabel = "X-ECEF", ylabel="Y-ECEF", zlabel="Z-ECEF",size=(1600,900),leg=false,markersize=0.1)#display background grid in 3D
            Np=length(pos_n)
            for i=1:Np
                if i==Np
                    display(scatter!(orbit_pos_all[1,i,:]/1e3,orbit_pos_all[2,i,:]/1e3,orbit_pos_all[3,i,:]/1e3,markersize=0.5))
                else
                    scatter!(orbit_pos_all[1,i,:]/1e3,orbit_pos_all[2,i,:]/1e3,orbit_pos_all[3,i,:]/1e3,markersize=0.5)
                end
            end
        end
        t0_index = findall(orbit_time_all .== 0)[1]
        if isempty(t0_index);@warn "orbit time should include reference time of 0";end
        orbit_time_index=(Int(round(SAR_start_time/dt_orbits))+t0_index:1:Int(round((SAR_start_time+SAR_duration)/dt_orbits))+t0_index) # index range for orbit times for time interval of interest
        orbit_pos=orbit_pos_all[:,:,orbit_time_index]
        orbit_time=orbit_time_all[orbit_time_index]
        orbit_vel=orbit_vel_all[:,:,orbit_time_index]
    end

    return orbit_time, orbit_pos, orbit_vel
end

function interpolateOrbitsToSlowTime(orbit_time, orbit_pos, SAR_start_time, params)
    @unpack SAR_duration, fp = params

    # interpolate orbit to slow time, 3 x Np x Nst, convert km to m
    slow_time = SAR_start_time : 1/fp : SAR_start_time+SAR_duration # create slow time axis

    Nst = size(slow_time)[1] # number of slow-time samples (pulses processed)

    if Nst == 1;
        p_xyz = orbit_pos
    else
        p_xyz = Orbits.interp_orbit(orbit_time,orbit_pos,slow_time)
    end # interpolate orbit to slow time, 3 x Np x Nst, convert km to m

    return p_xyz, Nst, slow_time
end

function computeTimePosVel_refpoint(params, ref_lat, ref_lon)
    @unpack SAR_start_time, SAR_duration, p_t0_LLH, orbit_filename, master_p = params

    orbit_dataset       = Dataset("inputs/"*orbit_filename) # Read orbits data in NetCDF format
    t12_orbits 		    = orbit_dataset["time"][1:2] # first two time samples
    dt_orbits 		    = t12_orbits[2]-t12_orbits[1] # time resolution of orbits (s)
    orbit_time_index    = Int(1):Int(dt_orbits):length(orbit_dataset["time"])
     # index range for orbit times for time interval of interest
    orbit_time 		    = orbit_dataset["time"][orbit_time_index] # read in time data
    orbit_pos_ECI 	    = 1e3*orbit_dataset["position"][:,:,orbit_time_index] # read in position data, 3 x Np x Nt
    orbit_vel_ECI       = 1e3*orbit_dataset["velocity"][:,:,orbit_time_index] # read in velocity data, 3 x Np x Nt (used optionally in avg peg and heading calculation)
    dv 				    = orbit_dataset.attrib["epoch"];
    epoch 			    = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
    global dcm 		    = Orbits.eci_dcm(orbit_time, epoch);
    orbit_pos,orbit_vel = Orbits.ecef_orbitpos(orbit_pos_ECI,orbit_vel_ECI,dcm)

    orbit_pos_geo       = Geometry.xyz_to_geo(orbit_pos[:,master_p,:])
    #close_val_lat_lon   = findmin(abs.(orbit_pos_geo[2,:].-(p_t0_LLH[2])) + abs.(orbit_pos_geo[1,:].-(p_t0_LLH[1])))
    close_val_lat_lon   = findmin(abs.(orbit_pos_geo[2,:].-(ref_lon)) + abs.(orbit_pos_geo[1,:].-(ref_lat)))

    mid_pt_geo          = orbit_pos_geo[:,close_val_lat_lon[2]]
    mid_pt_geo2          = orbit_pos_geo[:,close_val_lat_lon[2]+1]

    ref_lat             = mid_pt_geo[1]
    ref_lon             = mid_pt_geo[2]
    ref_ht              = mid_pt_geo[3]
    SAR_start_time      = orbit_time[close_val_lat_lon[2]] .+ SAR_start_time

    p_heading           = Geometry.compute_heading(ref_lat, ref_lon, orbit_vel[:,master_p,Int64(close_val_lat_lon[2])])


    orbit_time_index2    = (Int(floor(SAR_start_time/dt_orbits))-1:1:Int(ceil((SAR_start_time+SAR_duration)/dt_orbits))+1)
    orbit_time2          = orbit_dataset["time"][orbit_time_index2] # read in time data
    #orbit_pos_ECI       = 1e3*orbit_dataset["position"][:,:,orbit_time_index] # read in position data, 3 x Np x Nt
    #orbit_vel_ECI       = 1e3*orbit_dataset["velocity"][:,:,orbit_time_index] # read in velocity data, 3 x Np x Nt (used optionally in avg peg and heading calculation)
    #dv                  = orbit_dataset.attrib["epoch"];
    #local epoch         = DateTime(dv[1], dv[2], dv[3], dv[4], dv[5], dv[6]);
    #global dcm          = Orbits.eci_dcm(orbit_time, epoch);

    #orbit_pos,orbit_vel = Orbits.ecef_orbitpos(orbit_pos_ECI,orbit_vel_ECI,dcm)

    orbit_pos2=orbit_pos[:,:,orbit_time_index2]
    orbit_vel2=orbit_vel[:,:,orbit_time_index2]


    orbit_pos_geo       = Geometry.xyz_to_geo(orbit_pos2[:,master_p,:])

    p_heading2          = Geometry.compute_heading([orbit_pos_geo[1,1],orbit_pos_geo[1,end]], [orbit_pos_geo[2,1],orbit_pos_geo[2,end]])
    ref_peg             = Geometry.PegPoint(ref_lat, ref_lon, p_heading2[1])



    return orbit_time2, orbit_pos2, orbit_vel2, SAR_start_time, p_heading, ref_peg, p_heading2
end


"convert eci orbit posoitions to ECEF based using input DCM"
function ecef_orbitpos(eci_pos, dcm)
    #eci_pos = 3 x N_plat x N_time
    #dcm = 3 x 3 x 1 or 3 x 3 x N_time
    szp = size(eci_pos)
    szdcm = size(dcm)
    @assert szp[1]==3 "POS needs to be 3 x Np x Nt"
    if ndims(eci_pos) == 3
        nplat = szp[2]
        ntimes = szp[3]
    elseif ndims(eci_pos) == 2
        nplat = 1
        ntimes = szp[2]
    end

    if ndims(dcm) == 3 #DCM per time instance in orbit
        @assert szdcm[3] == ntimes "POS and DCM need to have same no. of entries"
    elseif ndims(dcm)==2 #A scalar DCM for all times (quasi-fixed ref. frame)
            @assert szdcm[1]==3 & szdcm[2]==3 "DCM needs to be 3 x 3 matrix"
            dcm = repeat(dcm,1,1,ntimes)
    else
        error("DCM is badly shaped")
    end

    ecef_pos = zeros(3,nplat, ntimes);
    for iplat=1:nplat
        for itime=1:ntimes
            ecef_pos[:,iplat,itime] = dcm[:,:,itime]*eci_pos[:,iplat,itime];
        end
    end
    return ecef_pos
end

"convert eci orbit posoitions and velocities to ECEF based using input DCM"
function ecef_orbitpos(eci_pos, eci_vel, dcm)
    #eci_pos = 3 x N_plat x N_time
    #dcm = 3 x 3 x 1 or 3 x 3 x N_time
    szp = size(eci_pos)
    szdcm = size(dcm)
    @assert szp[1]==3 "POS needs to be 3 x Np x Nt"
    if ndims(eci_pos) == 3
        nplat = szp[2]
        ntimes = szp[3]
    elseif ndims(eci_pos) == 2
        nplat = 1
        ntimes = szp[2]
    end

    if ndims(dcm) == 3 #DCM per time instance in orbit
        @assert szdcm[3] == ntimes "POS and DCM need to have same no. of entries"
    elseif ndims(dcm)==2 #A scalar DCM for all times (quasi-fixed ref. frame)
            @assert szdcm[1]==3 & szdcm[2]==3 "DCM needs to be 3 x 3 matrix"
            dcm = repeat(dcm,1,1,ntimes)
    else
        error("DCM is badly shaped")
    end

    # Earth rotation rate.
    w = 7.292115146706979e-5

    #set aside space for ECEF position and velocity
    ecef_pos = zeros(3,nplat, ntimes);
    ecef_vel = zeros(3,nplat, ntimes);

    for iplat=1:nplat
        for itime=1:ntimes
            ecef_pos[:,iplat,itime] = dcm[:,:,itime]*eci_pos[:,iplat,itime];
            ecef_vel[:,iplat,itime] = dcm[:,:,itime]*eci_vel[:,iplat,itime] - cross([0;0;w],ecef_pos[:,iplat,itime]);
        end
    end
    return ecef_pos, ecef_vel
end


"compute DCM to convert ECI to ECEF based on epoch and IERS EOP data"
function eci_dcm(time, epoch::DateTime, eop_data)
    dcm = zeros(3,3,length(time))
    for ii = 1:length(time)
          dt = unix2datetime(datetime2unix(epoch)+time[ii]);
          dcm[:,:,ii] = convert(Array{Float64}, r_eci_to_ecef(J2000(), ITRF(), date_to_jd(dt), eop_data));
    end
    return dcm
end

"compute DCM to convert ECI to ECEF based on epoch, wrapper method that downloads EOP data"
function eci_dcm(time, epoch::DateTime)
    eop_data = fetch_iers_eop();
    return eci_dcm(time, epoch, eop_data)
end

"interpolate darts orbits"
function interp_orbit(time_old, pos, time_new)
    #pos = 3 x N_plat x N_time
    szp = size(pos)
    @assert szp[1]==3 "POS needs to be 3 x [Np x Nt] or 3 x Nt"
    if ndims(pos) == 3
        nplat = szp[2]
    elseif ndims(pos) == 2
        nplat = 1
    end

    pos_i = zeros(szp[1],nplat, length(time_new));
    for iplat=1:nplat
        for iaxis=1:szp[1]
            itp = Spline1D(time_old, pos[iaxis, iplat,:], w=ones(length(time_old)), k=3, bc = "error");
            pos_i[iaxis,iplat,:] = evaluate(itp,time_new);
        end
    end
    return pos_i
end

"helper function to compute perpendicular baselines.
Inputs include 3D position vector (3xNpxNt) [m], velocity vector (3xNpxNt) [m/s]
look angle (θ in degrees), and reference index"
function get_perp_baselines(pos, vel, θ,refind=1)
    @assert ndims(pos)==3 "POS needs to be 3 x Np x Nt"
    @assert ndims(vel)==3 "VEL needs to be 3 x Np x Nt"
    @assert size(pos)==size(vel) "POS and VEL need to have same size"
    @assert refind <= size(pos,2) "Reference index needs to be <= Np"

    #set aside some space
    Nplats = size(pos,2); #number of platforms
    Ntimes = size(pos,3); #number of time steps
    bperp  = zeros(Nplats, Nplats, Ntimes); #output perp-baseline matrix
    b_at  = zeros(Nplats, Nplats, Ntimes); #output along-track baseline matrix
    bnorm  = zeros(Nplats, Nplats, Ntimes); #output along-track baseline matrix

    for itimes = 1:Ntimes
        #create geocetric TCN frame for reference satellite
        that, chat, nhat = Geometry.get_tcn(pos[:,refind,itimes], vel[:,refind,itimes]);
        #get look vector based on TCN-frame and look angle TODO: add azimuth
        lhat = cosd(θ)*nhat + sind(θ)*chat;

        #iterate over platforms and compute full perp-baseline matrix
        for iplat = 1:Nplats
            for jplat = 1:Nplats
                baseline = pos[:,jplat,itimes] - pos[:,iplat,itimes];

                bnorm[iplat,jplat,itimes] = norm(baseline);
                # use the imaging plane baseline only
                baseline_nc = dot(baseline,nhat)*nhat + dot(baseline,chat)*chat;
                # find baseline component perpendicular to look direction
                bperp[iplat,jplat,itimes] = norm(baseline_nc - dot(baseline_nc,lhat)*lhat);
                # find the along-track baseline component
                b_at[iplat,jplat,itimes]  = abs(dot(baseline,that));
            end
        end

    end
    return bperp, b_at, bnorm

end

"""
helper function to create synthetic orbits.
Arguments
    - `position::Array{3x1}`, ECEF xyz position (m,m,m)
    - `heading::Float64`, heading (deg)
    - `baselines::Array{3,Np}`, TCN baselines for Np platforms
    - `tvec::Array{3x1}`, time vector
   # Output
    - `platf_pos`
    - `platf_vel`
"""
function make_orbit(pos, hdg, baselines, tvec, mu = 3.986004418e14)
    @assert ndims(pos)==1 "POS needs to be 3 x 1 Vector"

    # create velocity vector from heading
    pos_llh = Geometry.xyz_to_geo(pos); #position in geodetic coords
    ehat,nhat,uhat = Geometry.enu_from_geo(pos_llh[1], pos_llh[2]); #ENU basis
    vel = cosd(hdg)*nhat + sind(hdg)*ehat;

    #compute TCN frame for that position
    tcnquat = Geometry.tcn_quat(pos,vel[:,1]);

    #create platform position and velocity vectors
    platf_pos = zeros(3,size(baselines,2), length(tvec));
    platf_vel = repeat(vel, 1, size(baselines,2), length(tvec)); #repeat velocity vector for now. TODO: update if necessary

    #spacecraf speed
    sc_speed = sqrt(mu./norm(pos)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed

    platf_vel = platf_vel .* sc_speed
    #iterate over all baselines to create synthetic orbits
    for ii=1:size(baselines,2)
        # compute platform position for each platform by adding baseline (in XYZ)
        # to platform position and repeating for all slow time samples
        # Note: we use rotate_vec instead of rotate_frame because tcn_quat describes
        #       rotation from XYZ to TCN
        #platf_pos[:,ii,:] = pos .+ repeat(Geometry.rotate_vec(baselines[:,ii], tcnquat[1]),1,length(tvec));
        for jj=1:length(tvec)
            platf_pos[:,ii,jj] = pos .+ Geometry.rotate_vec(baselines[:,ii], tcnquat[1]) + platf_vel[:,ii,jj]*tvec[jj];
        end
    end

    return platf_pos, platf_vel
end

"""
Updated helper function to create synthetic orbits.
Arguments
    - `position::Array{3x1}`, ECEF xyz position (m,m,m)
    - `heading::Float64`, heading (deg)
    - `baselines::Array{3,Np}`, TCN baselines for Np platforms
    - `tvec::Array{3x1}`, time vector
    - `θ::Float`: elevation angle (rad)
    - `ϕ::Float`: azimuth angle (rad)
   # Output
    - `platf_pos`
    - `platf_vel`
"""
function make_orbit_new(pos, hdg, baseline_pos, tvec, theta, phi=0.0, left_right_look="right", mu = 3.986004418e14)
    @assert ndims(pos)==1 "POS needs to be 3 x 1 Vector"

    # create velocity vector from heading
    pos_llh             = Geometry.xyz_to_geo(pos); #position in geodetic coords
    ehat,nhat,uhat      = Geometry.enu_from_geo(pos_llh[1], pos_llh[2]); #ENU basis
    vel                 = cosd(hdg)*nhat + sind(hdg)*ehat;
    
    #spacecraf speed
    sc_speed            = sqrt(mu./norm(pos)); #sqrt(GM/R)->https://en.wikipedia.org/wiki/Orbital_speed
    platf_vel           = vel .* sc_speed;

    that, chat, nhat    = Geometry.get_tcn(pos, platf_vel[:,1]);
    if left_right_look == "right"
        lhat                = Geometry.tcn_lvec(that, chat, nhat, theta*(pi/180), phi*(pi/180))
        bhat                = cross(-that,lhat)
    elseif left_right_look == "left"
        lhat                = Geometry.tcn_lvec(that, chat, nhat, -theta*(pi/180), phi*(pi/180))
        bhat                = cross(that,lhat)
    end

    #create platform position and velocity vectors
    platf_pos           = zeros(3,size(baseline_pos,2), length(tvec));
    platf_vel           = repeat(platf_vel, 1, size(baseline_pos,2), length(tvec));

    for jj=1:length(tvec)  
        #platf_pos[:,:,jj] = ( pos .+ ( bhat .* baseline_pos) ) .+ (platf_vel[:,:,jj]*tvec[jj]);
        platf_pos[:,:,jj] = ( pos .+ ( nhat .* baseline_pos) ) .+ (platf_vel[:,:,jj]*tvec[jj]); # for platform spacing along vertical axis

    end

    return platf_pos, platf_vel
end

"Updated helper function to compute perpendicular baselines.
Inputs include 3D position vector (3xNpxNt) [m], velocity vector (3xNpxNt) [m/s]
look angle (θ in degrees), and reference index"
function get_perp_baselines_new(pos, vel, theta, phi=0.0, left_right_look="right", refind=1)
    @assert ndims(pos)==3 "POS needs to be 3 x Np x Nt"
    @assert ndims(vel)==3 "VEL needs to be 3 x Np x Nt"
    @assert size(pos)==size(vel) "POS and VEL need to have same size"
    @assert refind <= size(pos,2) "Reference index needs to be <= Np"

    #set aside some space
    Nplats          = size(pos,2); #number of platforms
    Ntimes          = size(pos,3); #number of time steps
    bperp           = zeros(Nplats, Nplats, Ntimes); #output perp-baseline matrix
    b_at            = zeros(Nplats, Nplats, Ntimes); #output along-track baseline matrix
    bnorm           = zeros(Nplats, Nplats, Ntimes); #output along-track baseline matrix

    for itimes = 1:Ntimes
        #create geocetric TCN frame for reference satellite
        that, chat, nhat = Geometry.get_tcn(pos[:,refind,itimes], vel[:,refind,itimes]);
        #get look vector based on TCN-frame and look angle TODO: add azimuth
        if left_right_look == "right"
            lhat                = Geometry.tcn_lvec(that, chat, nhat, theta*(pi/180), phi*(pi/180))
        elseif left_right_look == "left"
            lhat                = Geometry.tcn_lvec(that, chat, nhat, -theta*(pi/180), phi*(pi/180))
        end

        #iterate over platforms and compute full perp-baseline matrix
        for iplat = 1:Nplats
            for jplat = 1:Nplats
                baseline = pos[:,jplat,itimes] - pos[:,iplat,itimes];

                bnorm[iplat,jplat,itimes] = norm(baseline);
                # use the imaging plane baseline only
                baseline_nc = dot(baseline,nhat)*nhat + dot(baseline,chat)*chat;
                # find baseline component perpendicular to look direction
                bperp[iplat,jplat,itimes] = norm(baseline_nc - dot(baseline_nc,lhat)*lhat);
                # find the along-track baseline component
                b_at[iplat,jplat,itimes]  = abs(dot(baseline,that));
            end
        end

    end
    return bperp, b_at, bnorm

end

end #end of module
