module Interferometry

include("../modules/geometry.jl")
include("../modules/scene.jl")
include("../modules/orbits.jl")
include("../modules/data_processing.jl")

using Statistics

c               = 299792458
earth_radius    = 6378.137e3 # Earth semi-major axis at equator

function get_scene_geometry_values(p_xyz, v_xyz, s_xyz_3xN, N_all, ref_plat, sec_plat, p_mode, params, grid, num_targ_vol )
    #geometry computations based on scene
    slant_range_ref                 = zeros(size(s_xyz_3xN,2),1)
    look_angle_ref                  = zeros(size(s_xyz_3xN,2),1)
    incidence_angle_ref             = zeros(size(s_xyz_3xN,2),1)
    slant_range_sec                 = zeros(size(s_xyz_3xN,2),1)
    look_angle_sec                  = zeros(size(s_xyz_3xN,2),1)
    incidence_angle_sec             = zeros(size(s_xyz_3xN,2),1)
    Perp_baseline_ref               = zeros(size(s_xyz_3xN,2),1)
    Vert_wavnum_ref                 = zeros(size(s_xyz_3xN,2),1)
    local_incidence_angle_ref       = zeros(size(s_xyz_3xN,2),1)
    range_slope_angle_ref           = zeros(size(s_xyz_3xN,2),1)
    Critical_baseline_ref           = zeros(size(s_xyz_3xN,2),1)
    Correlation_theo_ref            = zeros(size(s_xyz_3xN,2),1)

    mean_plats_pos_ref              = mean(p_xyz[:,ref_plat,:], dims=2)
    mean_plats_pos_sec              = mean(p_xyz[:,sec_plat,:], dims=2)

    ti2 = 1

    for ti = 1:size(s_xyz_3xN,2)
        slant_range_ref[ti]         = Geometry.distance( mean_plats_pos_ref  , s_xyz_3xN[:,ti] )
        look_angle_ref[ti]          = Scene.slantrange_to_lookangle(earth_radius,slant_range_ref[ti],Geometry.xyz_to_geo(mean_plats_pos_ref)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3])[2]
        incidence_angle_ref[ti]     = Scene.lookangle_to_incangle(look_angle_ref[ti],Geometry.xyz_to_geo(mean_plats_pos_ref)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3],earth_radius)

        slant_range_sec[ti]         = Geometry.distance( mean_plats_pos_sec  , s_xyz_3xN[:,ti] )
        look_angle_sec[ti]          = Scene.slantrange_to_lookangle(earth_radius,slant_range_sec[ti],Geometry.xyz_to_geo(mean_plats_pos_sec)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3])[2]
        incidence_angle_sec[ti]     = Scene.lookangle_to_incangle(look_angle_sec[ti],Geometry.xyz_to_geo(mean_plats_pos_sec)[3],Geometry.xyz_to_geo(s_xyz_3xN[:,ti])[3],earth_radius)

        bs_perp, bs_at, bs_norm     = Orbits.get_perp_baselines_new(mean(p_xyz[:,:,:],dims=3), mean(v_xyz[:,:,:],dims=3), look_angle_ref[ti], 0.0, params.left_right_look, 1)
        Perp_baseline_ref[ti]       = bs_perp[1,sec_plat,1]

        Va                          = mean(p_xyz[:, ref_plat, :], dims=2) - s_xyz_3xN[:, ti]
        Vb                          = mean(p_xyz[:, sec_plat, :], dims=2) - s_xyz_3xN[:, ti]
        angle_ip                    = Data_Processing.angle_2vec(Va, Vb) * 1

        if params.mode == 1
            Vert_wavnum_ref[ti]     = (4 * pi * (angle_ip * pi / 180) ) / (params.λ * sind(look_angle_ref[ti]))
        elseif params.mode == 2
            Vert_wavnum_ref[ti]     = (2 * pi * (angle_ip * pi / 180) ) / (params.λ * sind(look_angle_ref[ti]))
        end

        plat_pt_xyz                 =  mean(mean(p_xyz,dims=2),dims=3)[:] #??????
        look_vec_xyz                = (plat_pt_xyz - s_xyz_3xN[:,ti])
        look_vec_xyz_norm           = (plat_pt_xyz - s_xyz_3xN[:,ti]) / Geometry.distance( plat_pt_xyz,s_xyz_3xN[:,ti])

        Geo_location                = Geometry.xyz_to_geo(s_xyz_3xN[:,ti])
        pegθ                        = Geo_location[1]*π/180
        pegϕ                        = Geo_location[2]*π/180
        #ENU to XYZ transformation matrix
        Menu_xyz                    = [-sin(pegϕ) -sin(pegθ)*cos(pegϕ) cos(pegθ)*cos(pegϕ);
                                    cos(pegϕ) -sin(pegθ)*sin(pegϕ) cos(pegθ)*sin(pegϕ);
                                    0            cos(pegθ)             sin(pegθ)]
        #XYZ to ENU transformation matrix
        Mxyz_enu                    = [-sin(pegϕ)           cos(pegϕ)             0;
                                    -sin(pegθ)*cos(pegϕ) -sin(pegθ)*sin(pegϕ)  cos(pegθ)  ;
                                    cos(pegθ)*cos(pegϕ)   cos(pegθ)*sin(pegϕ)  sin(pegθ)]

        if grid == "Scene"
            Nxyz                        = Menu_xyz * [0,0,1];
        elseif grid == "Target"
            Nxyz                        = Menu_xyz * N_all[ti2,:];
        end
        look_vec_enu                = Mxyz_enu * look_vec_xyz

        local_incidence_angle_ref[ti] = Data_Processing.angle_2vec(look_vec_xyz_norm, Nxyz)

        look_direction_norm         = sqrt( look_vec_enu[1] * look_vec_enu[1]  + look_vec_enu[2] * look_vec_enu[2])

        if grid == "Scene"
            range_slope_angle_ref[ti]   = atand((0 * (look_vec_enu[1] / look_direction_norm)) + (0 * (look_vec_enu[2] / look_direction_norm)))    
        elseif grid == "Target"
            range_slope_angle_ref[ti]   = atand((N_all[ti2,1] * (look_vec_enu[1] / look_direction_norm)) + (N_all[ti2,2] * (look_vec_enu[2] / look_direction_norm)))    
        end            
        
        Critical_baseline_ref[ti]   = params.λ * ((2*params.bandwidth)/c) * slant_range_ref[ti] * tand(local_incidence_angle_ref[ti] - range_slope_angle_ref[ti]) / p_mode
        Correlation_theo_ref[ti]    = 1 - (Perp_baseline_ref[ti] ./ (Critical_baseline_ref[ti]))

        if (mod(ti,num_targ_vol)==0)
            ti2 = ti2 + 1
        end

    end
    return slant_range_ref, look_angle_ref, incidence_angle_ref, slant_range_sec, look_angle_sec, incidence_angle_sec, Perp_baseline_ref, Vert_wavnum_ref, local_incidence_angle_ref, range_slope_angle_ref, Critical_baseline_ref, Correlation_theo_ref
end



end