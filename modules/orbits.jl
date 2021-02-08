module Orbits

using SatelliteToolbox
using NCDatasets
using Dates
using Interpolations
using LinearAlgebra

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
    #println("No. platforms: ",nplat)
    #println("No. locations: ",ntimes)

    ecef_pos = zeros(3,nplat, ntimes);
    for iplat=1:nplat
        for itime=1:ntimes
            ecef_pos[:,iplat,itime] = dcm[:,:,itime]*eci_pos[:,iplat,itime];
        end
    end
    #println("Done with it, ", size(ecef_pos))
    return ecef_pos
end

"compute DCM to convert ECI to ECEF based on epoch and IERS EOP data"
function eci_dcm(time::Array{Float32,1}, epoch::DateTime, eop_data)
    dcm = zeros(3,3,length(time))
    for ii = 1:length(time)
          dt = unix2datetime(datetime2unix(epoch)+time[ii]);
          dcm[:,:,ii] = convert(Array{Float64}, rECItoECEF(J2000(), ITRF(), DatetoJD(dt), eop_data));
      end
      return dcm
end

"compute DCM to convert ECI to ECEF based on epoch"
function eci_dcm(time::Array{Float32,1}, epoch::DateTime)
    eop_data = get_iers_eop();
    dcm = zeros(3,3,length(time))
    for ii = 1:length(time)
          dt = unix2datetime(datetime2unix(epoch)+time[ii]);
          dcm[:,:,ii] = convert(Array{Float64}, rECItoECEF(J2000(), ITRF(), DatetoJD(dt), eop_data));
      end
      return dcm
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
            #TODO: using CubicSplineInterpolation instead of Linear Interpolations
            itp = LinearInterpolation(time_old, pos[iaxis, iplat, :])
            # itp = CubicSplineInterpolation(time_old, pos[iaxis, iplat,:])
            pos_i[iaxis,iplat,:] = itp(time_new);
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


    Nplats = size(pos,2); #number of platforms
    Ntimes = size(pos,3); #number of time steps

    print
    bperp = zeros(Nplats, Nplats, Ntimes);

    for itimes = 1:Ntimes
        #create geocetric TCN frame for reference satellite
        nhat = -pos[:,refind,itimes]/norm(pos[:,refind,itimes]);
        vrad = dot(vel[:,refind,itimes], -nhat)*(-nhat); #radial velocity
        vtan = vel[:,refind,itimes]  - vrad; #tangential velocity
        that = vtan/norm(vtan); #track vector
        chat = cross(nhat, that); #cross-track vector
        #get look vector based on TCN-frame and look angle TODO: add azimuth
        lhat = cosd(θ)*nhat + sind(θ)*chat;
        #iterate over platforms and compute full perp-baseline matrix
        for iplat = 1:Nplats
            for jplat = 1:Nplats
                baseline = pos[:,jplat,itimes] - pos[:,iplat,itimes];
                # use the imaging plane baseline only
                baseline = dot(baseline,nhat)*nhat + dot(baseline,chat)*chat;
                # find baseline component perpendicular to look direction
                bperp[iplat,jplat,itimes] = norm(baseline - dot(baseline,lhat)*lhat);
            end
        end

    end
    return bperp

end



end #end of module
