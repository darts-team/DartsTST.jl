## scattering module
module Scattering
using Parameters, LinearAlgebra, Random
include("../modules/geometry.jl")

"""
Calculate surface BRCS
# Usage
    - brcs = get_surface_brcs(params,tx_ecef,rx_ecef,tgt_ecef,patch_area)

# Arguments
    - `params::Parameters`: input parameters struct
    - `tx_ecef::3x1-element Array`: Transmitter position in ECEF (m)
    - `rx_ecef::3x1-element Array`: Receiver position in ECEF (m)
    - `tgt_ecef::3x1-element Array`: Target position in ECEF (m)

# Return
- `brcs_complex::ComplexF64`: BRCS value (linear magnitude, unitless)
"""
function get_surface_brcs(params,tx_ecef,rx_ecef,tgt_ecef)
    @unpack λ, polarization, s_loc_1, s_loc_2, σ, l, θᵥ, fc = params

    # this function should probably be called from scene.jl to calculate the brcs for a given target:
    # in there, geometry of Tx to Tgt and Tgt to Rx needs be solved to find inc angles and scattering angles
    # generate raw data has to be modified to use target reflectivites that are based on the brcs, and unfortunately there are 
    # changing values/sizes for all modes (with MIMO having different matrix size).


    # ------in this function-----
    # determine geometry to select backscatter codes or bistatic codes (backscattering faster)
    # - for geometry: calculate position of tx and rx in ENU frame local to the target position
    # check constraints to determine SPM or KA usage
    # determine brcs value from selected value. Report in linear units?
    # allow for profile of reflectivity. profile should be in a magnitude - phase format so that phase values can be altered by some phasing function to set based on profile data
    patch_area = (s_loc_1[2]-s_loc_1[1]) * (s_loc_2[2]-s_loc_2[1])
    
    
    #determine scattering geometry given Rx/tx/target locations
    tx_llh = Geometry.xyz_to_geo(tx_ecef)
    rx_llh = Geometry.xyz_to_geo(rx_ecef)
    tgt_llh = Geometry.xyz_to_geo(tgt_ecef)

    #next 3 lines for testing geometry
    # tx_llh = [0; 0; 500e3]; rx_llh = [0;20; 500e3]; tgt_llh = [0;10;0] # bistatic - specular
    # tx_llh = [0; 0; 500e3]; rx_llh = [0; 0; 500e3]; tgt_llh = [0;10;0] # monostatic
    # tx_llh = [0; 0; 500e3]; rx_llh = [.1; 0; 500e3]; tgt_llh = [0;10;0] # quasi-monostatic (?)

    tx_enu = Geometry.llh_to_enu_new_org(tx_llh[1], tx_llh[2], tx_llh[3], tgt_llh)
    rx_enu = Geometry.llh_to_enu_new_org(rx_llh[1], rx_llh[2], rx_llh[3], tgt_llh)
    tx_sph = Geometry.enu_to_sph(tx_enu[1], tx_enu[2], tx_enu[3])
    rx_sph = Geometry.enu_to_sph(rx_enu[1], rx_enu[2], rx_enu[3])
    #correct transmitter azimuth for change of origin direction (origin calculated as target location)
    tx_sph[3] = tx_sph[3] + 180; tx_sph[3] = tx_sph[3]%360 # add 180 degrees and mod360 to enure 0 < azimuth < 360

    # change elevation angles to incidence angles
    tx_incidence = 90 - tx_sph[2]
    rx_incidence = 90 - rx_sph[2]
    
    # check for "monostatic" bounds -- TODO: How to define what we approximate as monostatic? 
    # i.e. what scattering error tolerance is acceptable --> angular difference from direct backscatter that meets tolerance
    #for now we'll say 2 degrees in both azimuth and incidence...
    tol = 2.0
    if (abs(tx_incidence - rx_incidence) < tol) & ( abs(  tx_sph[3] - rx_sph[3] + 180 ) < tol)
        monostatic = true
    else
        monostatic = false
    end


    # ulaby_terrain_switch will calculate the backscatter RCS value at the incident angle
    if ulaby_terrain_switch == 1       
        if monostatic == false
            @warn raw"Suggest only to use Ulaby Terrain for monostatic-like geometries."
        end
        σᵘ_vh,  σᵘ_hv,  σᵘ_vv,  σᵘ_hh = Ulaby_book_terrain_backscatter_values(rx_incidence,fc,terrain)
        rx_unit = rx_enu/norm(rx_enu)
        
        patch_normal = [0; 0; 1] # using this a placeholder; future work might include a DEM, in which the surface normal may be tilted
        patch_area_eff = patch_area * dot(patch_normal,rx_unit) # take dot product between surface normal [0 0 1] in ENU and unit vector to receiver
        brcs = σʳ*patch_area_eff
    
    #TODO fix this version where the values are scaled to Ulaby at nadir. But how? Would need to calculate at a level above and fit an RCS curve that fits through the Ulaby value. This only calculates at this level
    
        # elseif ulaby_terrain_switch == 1 # this will use the standard scattering calculations but scale to the nearest Ulaby terrain value
        
    #     σᵘ_vh,  σᵘ_hv,  σᵘ_vv,  σᵘ_hh = Ulaby_book_terrain_backscatter_values(rx_incidence,fc,terrain)

    #     if monostatic #monostatic geometry
    #         # Auto-select surface scattering mechanism based on surface properties
    #         m = 2 * σ / l
    #         if λ*σ < l^2/2.76
    #             #use GO
    #             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA_backscatter(tx_incidence,σ,l,θᵥ)
                
    #         elseif σ < λ/21 & m < 0.3
    #             #use SPM
    #             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_backscatter(λ,tx_incidence,l,σ,θᵥ)
    #         else
    #             error(raw"Scattering surface not supported")
    #         end
    #     else # generalized bistatic geometry
    #         # Auto-select surface scattering mechanism based on surface properties
    #         m = 2 * σ / l
    #         if λ*σ < l^2/2.76
    #             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],σ,l,θᵥ)
                
    #         elseif (σ < λ/21) & (m < 0.3)
    #             #use SPM
    #             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_tsang(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],l,σ,θᵥ)
    #         else
    #             error(raw"Scattering surface not supported")
    #         end
    #     end#if monostatic

    #     rx_unit = rx_enu/norm(rx_enu)
    #     patch_normal = [0; 0; 1] # using this a placeholder; future work might include a DEM, in which the surface normal may be tilted
    #     patch_area_eff = patch_area * dot(patch_normal,rx_unit) # take dot product between surface normal [0 0 1] in ENU and unit vector to receiver
    #     #select polarization for BRCS calc #sadly theres no "switch:case" statement in Julia
    #     if polarization == 1 #vh
    #         brcs = σʳ_vh * patch_area_eff
    #     elseif polarization == 2 # hv
    #         brcs = σʳ_hv * patch_area_eff
    #     elseif polarization == 3 # vv
    #         brcs = σʳ_vv * patch_area_eff
    #     elseif polarization == 4 # hh
    #         brcs = σʳ_hh * patch_area_eff
    #     end    


    else # use standard scattering calculations
        if monostatic #monostatic geometry
            # Auto-select surface scattering mechanism based on surface properties
            m = 2 * σ / l
            if λ*σ < l^2/2.76
                #use GO
                σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA_backscatter(tx_incidence,σ,l,θᵥ)
                
            elseif σ < λ/21 & m < 0.3
                #use SPM
                σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_backscatter(λ,tx_incidence,l,σ,θᵥ)
            else
                error(raw"Scattering surface not supported")
            end
        else # generalized bistatic geometry
            # Auto-select surface scattering mechanism based on surface properties
            m = 2 * σ / l
            if λ*σ < l^2/2.76
                σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],σ,l,θᵥ)
                
            elseif (σ < λ/21) & (m < 0.3)
                #use SPM
                σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_tsang(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],l,σ,θᵥ)
            else
                error(raw"Scattering surface not supported")
            end
        end#if monostatic

        rx_unit = rx_enu/norm(rx_enu)
        
        patch_normal = [0; 0; 1] # using this a placeholder; future work might include a DEM, in which the surface normal may be tilted
        patch_area_eff = patch_area * dot(patch_normal,rx_unit) # take dot product between surface normal [0 0 1] in ENU and unit vector to receiver
        # deciding not to take dot product of transmitter and surface normal; thought experiment is that for monostatic case, you don't have cos^2(θ) effect 

        
        #select polarization for BRCS calc #sadly theres no "switch:case" statement in Julia
        if polarization == 1 #vh
            brcs = σʳ_vh * patch_area_eff
        elseif polarization == 2 # hv
            brcs = σʳ_hv * patch_area_eff
        elseif polarization == 3 # vv
            brcs = σʳ_vv * patch_area_eff
        elseif polarization == 4 # hh
            brcs = σʳ_hh * patch_area_eff
        end    


        # check for nadir-like geometry and calculate nadir specular component. use same tolerance as backscatter for area around nadir
        # seeing if elevation is directly up for tx and rx
        if (abs(90.0 - tx_incidence) < tol) & (abs(90.0 - rx_incidence) < tol)
            brcs_coh = RCS_coherent(σ,θᵥ,λ,tx_incidence,ϕ,rx_incidence,ϕₛ)
            brcs = brcs .+ brcs_coh[polarization]
        end
    end

    phase = 0 # keep a fixed phase value
    brcs_complex = brcs*cos(phase)+ (1im*brcs*sin(phase))

    #TODO return a scattering matrix instead. We're already calculating all 4 components
    return brcs_complex
end#function

"""
Calculate surface BRCS
# Usage
    - brcs = get_brcs(params,tx_ecef,rx_ecef,tgt_ecef,patch_area)

# Arguments
    - `params::Parameters`: input parameters struct
    - `tx_ecef::3x1-element Array`: Transmitter position in ECEF (m)
    - `rx_ecef::3x1-element Array`: Receiver position in ECEF (m)
    - `tgt_ecef::3x1-element Array`: Target position in ECEF (m)
    - `patch_area::Float32`: Area of surface patch (m^2)

# Return
- `brcs_complex::ComplexF64`: BRCS value (linear magnitude, unitless)
"""
function get_surface_brcs(params,tx_ecef,rx_ecef,tgt_ecef,patch_area)
    @unpack λ, polarization, σ, l, θᵥ = params

    # this function should probably be called from scene.jl to calculate the brcs for a given target:
    # in there, geometry of Tx to Tgt and Tgt to Rx needs be solved to find inc angles and scattering angles
    # generate raw data has to be modified to use target reflectivites that are based on the brcs, and unfortunately there are 
    # changing values/sizes for all modes (with MIMO having different matrix size).


    # ------in this function-----
    # determine geometry to select backscatter codes or bistatic codes (backscattering faster)
    # - for geometry: calculate position of tx and rx in ENU frame local to the target position
    # check constraints to determine SPM or KA usage
    # determine brcs value from selected value. Report in linear units?
    # allow for profile of reflectivity. profile should be in a magnitude - phase format so that phase values can be altered by some phasing function to set based on profile data
    
    
    
    #determine scattering geometry given Rx/tx/target locations
    tx_llh = Geometry.xyz_to_geo(tx_ecef)
    rx_llh = Geometry.xyz_to_geo(rx_ecef)
    tgt_llh = Geometry.xyz_to_geo(tgt_ecef)

    #next 3 lines for testing geometry
    # tx_llh = [0; 0; 500e3]; rx_llh = [0;20; 500e3]; tgt_llh = [0;10;0] # bistatic - specular
    # tx_llh = [0; 0; 500e3]; rx_llh= [0; 0; 500e3]; tgt_llh = [0;10;0] # monostatic
    # tx_llh = [0; 0; 500e3]; rx_llh = [.1; 0; 500e3]; tgt_llh = [0;10;0] # quasi-monostatic (?)

    tx_enu = Geometry.llh_to_enu_new_org(tx_llh[1], tx_llh[2], tx_llh[3], tgt_llh)
    rx_enu = Geometry.llh_to_enu_new_org(rx_llh[1], rx_llh[2], rx_llh[3], tgt_llh)
    tx_sph = Geometry.enu_to_sph(tx_enu[1], tx_enu[2], tx_enu[3])
    rx_sph = Geometry.enu_to_sph(rx_enu[1], rx_enu[2], rx_enu[3])
    #correct transmitter azimuth for change of origin direction (origin calculated as target location)
    tx_sph[3] = tx_sph[3] + 180; tx_sph[3] = tx_sph[3]%360 # add 180 degrees and mod360 to enure 0 < azimuth < 360

    # change elevation angles to incidence angles
    tx_incidence = 90 - tx_sph[2]
    rx_incidence = 90 - rx_sph[2]

    # check for "monostatic" bounds -- TODO: How to define what we approximate as monostatic? 
    # i.e. what scattering error tolerance is acceptable --> angular difference from direct backscatter that meets tolerance
    #for now we'll say 2 degrees in both azimuth and incidence...
    tol = 2.0
    if (abs(tx_incidence - rx_incidence) < tol) & ( abs(  tx_sph[3] - rx_sph[3] + 180 ) < tol)
        monostatic = true
    else
        monostatic = false
    end



    if monostatic #monostatic geometry
        # Auto-select surface scattering mechanism based on surface properties
        m = 2 * σ / l
        if λ*σ < l^2/2.76
            #use GO
            σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA_backscatter(tx_incidence,σ,l,θᵥ)
            
        elseif σ < λ/21 & m < 0.3
            #use SPM
            σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_backscatter(λ,tx_incidence,l,σ,θᵥ)
        else
            error(raw"Scattering surface not supported")
        end
    else # generalized bistatic geometry
         # Auto-select surface scattering mechanism based on surface properties
         m = 2 * σ / l
         if λ*σ < l^2/2.76
             #use GO
             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],σ,l,θᵥ)
             
         elseif (σ < λ/21) & (m < 0.3)
             #use SPM
             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_tsang(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],l,σ,θᵥ)
         else
             error(raw"Scattering surface not supported")
         end
    end#if monostatic

    patch_normal = [0; 0; 1] # using this a placeholder; future work might include a DEM, in which the surface normal may be tilted
    patch_area_eff = patch_area * dot(patch_normal,rx_unit) # take dot product between surface normal [0 0 1] in ENU and unit vector to receiver
    # deciding not to take dot product of transmitter and surface normal; thought experiment is that for monostatic case, you don't have cos^2(θ) effect 

    
    #select polarization for BRCS calc #sadly theres no "switch:case" statement in Julia
    if polarization == 1 #vh
        brcs = σʳ_vh * patch_area_eff
    elseif polarization == 2 # hv
        brcs = σʳ_hv * patch_area_eff
    elseif polarization == 3 # vv
        brcs = σʳ_vv * patch_area_eff
    elseif polarization == 4 # hh
        brcs = σʳ_hh * patch_area_eff
    end    


    # check for nadir-like geometry and calculate nadir specular component. use same tolerance as backscatter for area around nadir
    # seeing if elevation is directly up for tx and rx
    if abs(90 - tx_incidence) < tol & abs(90 - rx_incidence) < tol
        brcs_coh = RCS_coherent(σ,θᵥ,λ,tx_incidence,ϕ,rx_incidence,ϕₛ)
        brcs = brcs + brcs_coh[polarization]
    end

    phase = 0 # keep a fixed phase value
    brcs_complex = brcs*cos(phase)+ (1im*brcs*sin(phase))



    return brcs_complex
end#function

"""
Calculate surface BRCS
# Usage
    - brcs = get_surface_brcs(params,tx_ecef,rx_ecef,tgt_ecef,patch_area)

# Arguments
    - `params::Parameters`: input parameters struct
    - `tx_ecef::3x1-element Array`: Transmitter position in ECEF (m)
    - `rx_ecef::3x1-element Array`: Receiver position in ECEF (m)
    - `tgt_ecef::3x1-element Array`: Target position in ECEF (m)

# Return
- `brcs_complex::ComplexF64`: BRCS value (linear magnitude, unitless)
"""
function get_surface_brcs(params,tx_ecef,rx_ecef,tgt_ecef)
    @unpack λ, polarization, s_loc_1, s_loc_2, σ, l, θᵥ, fc = params

    # this function should probably be called from scene.jl to calculate the brcs for a given target:
    # in there, geometry of Tx to Tgt and Tgt to Rx needs be solved to find inc angles and scattering angles
    # generate raw data has to be modified to use target reflectivites that are based on the brcs, and unfortunately there are 
    # changing values/sizes for all modes (with MIMO having different matrix size).


    # ------in this function-----
    # determine geometry to select backscatter codes or bistatic codes (backscattering faster)
    # - for geometry: calculate position of tx and rx in ENU frame local to the target position
    # check constraints to determine SPM or KA usage
    # determine brcs value from selected value. Report in linear units?
    # allow for profile of reflectivity. profile should be in a magnitude - phase format so that phase values can be altered by some phasing function to set based on profile data
    patch_area = (s_loc_1[2]-s_loc_1[1]) * (s_loc_2[2]-s_loc_2[1])
    
    
    #determine scattering geometry given Rx/tx/target locations
    tx_llh = Geometry.xyz_to_geo(tx_ecef)
    rx_llh = Geometry.xyz_to_geo(rx_ecef)
    tgt_llh = Geometry.xyz_to_geo(tgt_ecef)

    #next 3 lines for testing geometry
    # tx_llh = [0; 0; 500e3]; rx_llh = [0;20; 500e3]; tgt_llh = [0;10;0] # bistatic - specular
    # tx_llh = [0; 0; 500e3]; rx_llh = [0; 0; 500e3]; tgt_llh = [0;10;0] # monostatic
    # tx_llh = [0; 0; 500e3]; rx_llh = [.1; 0; 500e3]; tgt_llh = [0;10;0] # quasi-monostatic (?)

    tx_enu = Geometry.llh_to_enu_new_org(tx_llh[1], tx_llh[2], tx_llh[3], tgt_llh)
    rx_enu = Geometry.llh_to_enu_new_org(rx_llh[1], rx_llh[2], rx_llh[3], tgt_llh)
    tx_sph = Geometry.enu_to_sph(tx_enu[1], tx_enu[2], tx_enu[3])
    rx_sph = Geometry.enu_to_sph(rx_enu[1], rx_enu[2], rx_enu[3])
    #correct transmitter azimuth for change of origin direction (origin calculated as target location)
    tx_sph[3] = tx_sph[3] + 180; tx_sph[3] = tx_sph[3]%360 # add 180 degrees and mod360 to enure 0 < azimuth < 360

    # change elevation angles to incidence angles
    tx_incidence = 90 - tx_sph[2]
    rx_incidence = 90 - rx_sph[2]
    
    # check for "monostatic" bounds -- TODO: How to define what we approximate as monostatic? 
    # i.e. what scattering error tolerance is acceptable --> angular difference from direct backscatter that meets tolerance
    #for now we'll say 2 degrees in both azimuth and incidence...
    tol = 2.0
    if (abs(tx_incidence - rx_incidence) < tol) & ( abs(  tx_sph[3] - rx_sph[3] + 180 ) < tol)
        monostatic = true
    else
        monostatic = false
    end


    # ulaby_terrain_switch will calculate the backscatter RCS value at the incident angle
    if ulaby_terrain_switch == 1       
        if monostatic == false
            @warn raw"Suggest only to use Ulaby Terrain for monostatic-like geometries."
        end
        σᵘ_vh,  σᵘ_hv,  σᵘ_vv,  σᵘ_hh = Ulaby_book_terrain_backscatter_values(rx_incidence,fc,terrain)
        rx_unit = rx_enu/norm(rx_enu)
        
        patch_normal = [0; 0; 1] # using this a placeholder; future work might include a DEM, in which the surface normal may be tilted
        patch_area_eff = patch_area * dot(patch_normal,rx_unit) # take dot product between surface normal [0 0 1] in ENU and unit vector to receiver
        brcs = σʳ*patch_area_eff
    
    #TODO fix this version where the values are scaled to Ulaby at nadir. But how? Would need to calculate at a level above and fit an RCS curve that fits through the Ulaby value. This only calculates at this level
    
        # elseif ulaby_terrain_switch == 1 # this will use the standard scattering calculations but scale to the nearest Ulaby terrain value
        
    #     σᵘ_vh,  σᵘ_hv,  σᵘ_vv,  σᵘ_hh = Ulaby_book_terrain_backscatter_values(rx_incidence,fc,terrain)

    #     if monostatic #monostatic geometry
    #         # Auto-select surface scattering mechanism based on surface properties
    #         m = 2 * σ / l
    #         if λ*σ < l^2/2.76
    #             #use GO
    #             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA_backscatter(tx_incidence,σ,l,θᵥ)
                
    #         elseif σ < λ/21 & m < 0.3
    #             #use SPM
    #             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_backscatter(λ,tx_incidence,l,σ,θᵥ)
    #         else
    #             error(raw"Scattering surface not supported")
    #         end
    #     else # generalized bistatic geometry
    #         # Auto-select surface scattering mechanism based on surface properties
    #         m = 2 * σ / l
    #         if λ*σ < l^2/2.76
    #             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],σ,l,θᵥ)
                
    #         elseif (σ < λ/21) & (m < 0.3)
    #             #use SPM
    #             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_tsang(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],l,σ,θᵥ)
    #         else
    #             error(raw"Scattering surface not supported")
    #         end
    #     end#if monostatic

    #     rx_unit = rx_enu/norm(rx_enu)
    #     patch_normal = [0; 0; 1] # using this a placeholder; future work might include a DEM, in which the surface normal may be tilted
    #     patch_area_eff = patch_area * dot(patch_normal,rx_unit) # take dot product between surface normal [0 0 1] in ENU and unit vector to receiver
    #     #select polarization for BRCS calc #sadly theres no "switch:case" statement in Julia
    #     if polarization == 1 #vh
    #         brcs = σʳ_vh * patch_area_eff
    #     elseif polarization == 2 # hv
    #         brcs = σʳ_hv * patch_area_eff
    #     elseif polarization == 3 # vv
    #         brcs = σʳ_vv * patch_area_eff
    #     elseif polarization == 4 # hh
    #         brcs = σʳ_hh * patch_area_eff
    #     end    


    else # use standard scattering calculations
        if monostatic #monostatic geometry
            # Auto-select surface scattering mechanism based on surface properties
            m = 2 * σ / l
            if λ*σ < l^2/2.76
                #use GO
                σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA_backscatter(tx_incidence,σ,l,θᵥ)
                
            elseif σ < λ/21 & m < 0.3
                #use SPM
                σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_backscatter(λ,tx_incidence,l,σ,θᵥ)
            else
                error(raw"Scattering surface not supported")
            end
        else # generalized bistatic geometry
            # Auto-select surface scattering mechanism based on surface properties
            m = 2 * σ / l
            if λ*σ < l^2/2.76
                σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],σ,l,θᵥ)
                
            elseif (σ < λ/21) & (m < 0.3)
                #use SPM
                σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_tsang(λ,tx_incidence,tx_sph[3],rx_incidence,rx_sph[3],l,σ,θᵥ)
            else
                error(raw"Scattering surface not supported")
            end
        end#if monostatic

        rx_unit = rx_enu/norm(rx_enu)
        
        patch_normal = [0; 0; 1] # using this a placeholder; future work might include a DEM, in which the surface normal may be tilted
        patch_area_eff = patch_area * dot(patch_normal,rx_unit) # take dot product between surface normal [0 0 1] in ENU and unit vector to receiver
        # deciding not to take dot product of transmitter and surface normal; thought experiment is that for monostatic case, you don't have cos^2(θ) effect 

        
        #select polarization for BRCS calc #sadly theres no "switch:case" statement in Julia
        if polarization == 1 #vh
            brcs = σʳ_vh * patch_area_eff
        elseif polarization == 2 # hv
            brcs = σʳ_hv * patch_area_eff
        elseif polarization == 3 # vv
            brcs = σʳ_vv * patch_area_eff
        elseif polarization == 4 # hh
            brcs = σʳ_hh * patch_area_eff
        end    


        # check for nadir-like geometry and calculate nadir specular component. use same tolerance as backscatter for area around nadir
        # seeing if elevation is directly up for tx and rx
        if (abs(90.0 - tx_incidence) < tol) & (abs(90.0 - rx_incidence) < tol)
            brcs_coh = RCS_coherent(σ,θᵥ,λ,tx_incidence,ϕ,rx_incidence,ϕₛ)
            brcs = brcs + brcs_coh[polarization]
        end
    end

    phase = 0 # keep a fixed phase value
    brcs_complex = brcs*cos(phase)+ (1im*brcs*sin(phase))

    #TODO return a scattering matrix instead. We're already calculating all 4 components
    return brcs_complex
end#function

# KA-GO brcs calculation from eq. 12.43 in Microwave Remote Sensing Vol. 2, by Ulaby, Moore, Fung   (page 935)
# assumptions: stationary phase approximation- derived from Kirchhoff approximation (tangent-plane). Rough surface scattering where local diffraction,
# shadowing, multiple scattering is ignored.

# Equation: σʳ_pq = [ k q |U_pq| ]² / [ 2 (q_z)⁴ σ² |ρ''(0)| ] exp[- ( (q_x)² + (q_y)²) ) / (2 (q_z)² σ²  )]
# where:
# k: wavenumber (in medium 1, assume free space)
# q_x = k (sinθₛ cosϕₛ - sinθ cosϕ ) 
# q_y = k (sinθₛ sinϕₛ - sinθ sinϕ ) 
# q_z = k (cosθₛ + cosθ ) 
# surface heights are Gaussian distributed with p(z) = (2π σ^2)^(-1/2) exp( -z^2 / (2σ^2) ), σ = stdev of surface heights
# U_pq: polarization factors. U_hh, U_vv, U_hv, U_vh
    # U_vh = -R_prp0(1+cosθcosθₛ)sin(ϕₛ-ϕ) - [R_prp0 sinθcosθₛ+ R_prp1(1+cosθcosθₛ)]* sin(ϕₛ-ϕ)(Z_x cosϕ + Z_y sinϕ) ----- (Z_x,Z_y are local surface slopes)
    # U_hv = R_pll0(1+cosθcosθₛ)sin(ϕₛ-ϕ) + [R_pll0 sinθcosθₛ+ R_pll1(1+cosθcosθₛ)]* sin(ϕₛ-ϕ)(Z_x cosϕ + Z_y sinϕ)
    # U_hh = -R_prp0(cosθ+cosθₛ)cos(ϕₛ-ϕ) + {R_prp0[sinθₛ - sinθcos(ϕₛ-ϕ)]-R_prp1 (cosθ+cosθₛ) * cos(ϕₛ-ϕ)} (Z_x cosϕ + Z_y sinϕ)
    # U_vv = -R_pll0 (cosθ+cosθₛ)cos(ϕₛ-ϕ)+ {R_pll0[sinθₛ - sinθcos(ϕₛ-ϕ)] + R_pll1 (cosθ+cosθₛ) * cos(ϕₛ-ϕ)}(Z_x cosϕ + Z_y sinϕ)
# R_pll is the parallel reflection coeff, R_prp is the perpendicular. 0 and 1 are the first/second order Taylor series terms. Assumes surface slope < 0.25
# R_pll0 = (η₁cosθ - η₂cosθₜ)/(η₁cosθ + η₂cosθₜ)
# R_pll1 = - [η₁sinθ - η₂sinθₜ - R_pll0 (η₁sinθ + η₂sinθₜ) ] / (η₁cosθ + η₂cosθₜ)
# R_prp0 = (η₂cosθ - η₁cosθₜ)/(η₂cosθ + η₁cosθₜ)
# R_prp1 = -R_prp0 (η₂sinθ + η₁sinθₜ)/(η₂cosθ + η₁cosθₜ)
#-----------------------------------------------------------------------------------------------------------------------------------------

"""
# this function calculates the normalized bistatic radar crossection (nbrcs) of a rough surface using KA-GO for vh,hv,vv,hh polatizations based on incident and scattering angles, and surface rougness

# Arguments
- `λ::Float64`: slow time vector
- `θ::Float64: incidence angle of incident wave (deg)
- `ϕ::Float64`: azimuth angle of incident wave (deg)
- `θₛ::Float64`: incidence angle of scattering direction (deg)
- `ϕₛ::Float64`: azimuth angle of scattering direction (deg)
- `σ::Float64`: standard deviation of surface heights (m)
- `l::Float64`: surface correlation length (m)
- `θᵥ::Float64: soil moisture volume fraction [0,1]

"""
function BRCS_KA(λ,θ,ϕ,θₛ,ϕₛ,σ,l,θᵥ)
    # TODO Replace surface input values with input_parameters struct?
    k = 2*pi / λ
    f = 299792458 / λ  # approx center freq derived from free space λ

    q_x,q_y,q_z = get_q_vec(k,θ,ϕ,θₛ,ϕₛ)
    q = sqrt(q_x^2 + q_y^2 + q_z^2)
    σ² = σ^2 # create sigma squared notation

    ϵ₂ = soil_dielectric(θᵥ)
    μ₂ = 12.566e-7 # free space permeability
    
    #mean squared slope = σ² * abs(p_dp0)
    s = 2 * σ / l
    #----------- using mean slope values ------------
    mss = s^2/2 #since m = s/sqrt(2) and mss = m^2, where m = rms slope
    Z_x = mss; Z_y = mss;
    U_vh,U_hv,U_hh,U_vv = get_pol_factors_U_pq(ϵ₂,μ₂,θ,θₛ,ϕ,ϕₛ,f,Z_x,Z_y) # Here we define slopes Z_x and Z_y as the mss values, assuming 2D gaussian shape(?)
    σʳ_vh  = (k*q* abs(U_vh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_hv  = (k*q* abs(U_hv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_vv  = (k*q* abs(U_vv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_hh  = (k*q* abs(U_hh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    #----------------------------------------------------


    #below is code for calculating a mean value of the scattering using averages over a distribution of surface slopes
    # #local surface normal is nhat defined as =-xhat*Z_x - yhat*Z_y + zhat
    # # let Z_x, Z_y be Gaussian random numbers with  zero-mean and variance = mss
    # numInstances = 100
    # Z_x = randn(numInstances).*sqrt(σ²/2/pi)
    # Z_y = randn(numInstances).*sqrt(σ²/2/pi)
    
    # # now, we'll need to calculate average brcs based on instances of Z_x, Z_y
    # σʳ_vh = zeros(1,numInstances)
    # σʳ_hv = zeros(1,numInstances)
    # σʳ_vv = zeros(1,numInstances)
    # σʳ_hh = zeros(1,numInstances)
    # for i = 1 : numInstances
    #     U_vh,U_hv,U_hh,U_vv = get_pol_factors_U_pq(ϵ₂,μ₂,θ,θₛ,ϕ,ϕₛ,f,Z_x[i],Z_y[i]) # no idea how to define Z right now

    #     σʳ_vh[i]  = (k*q* abs(U_vh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    #     σʳ_hv[i]  = (k*q* abs(U_hv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    #     σʳ_vv[i]  = (k*q* abs(U_vv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    #     σʳ_hh[i]  = (k*q* abs(U_hh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    # end
    # #find mean BRCS values
    # σʳ_vh = sum(σʳ_vh)/numInstances
    # σʳ_hv = sum(σʳ_hv)/numInstances
    # σʳ_vv = sum(σʳ_vv)/numInstances
    # σʳ_hh = sum(σʳ_hh)/numInstances
    
    return  σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh
end#function

"""
# this function calculates the normalized bistatic radar crossection (nbrcs) of a rough surface using KA-GO for vh,hv,vv,hh polatizations based on incident and scattering angles, and surface rougness
# overloaded function call with mean slope instead of std of heights and correlation length

# Arguments
- `λ::Float64`: wavelength
- `θ::Float64: incidence angle of incident wave (deg)
- `ϕ::Float64`: azimuth angle of incident wave (deg)
- `θₛ::Float64`: incidence angle of scattering direction (deg)
- `ϕₛ::Float64`: azimuth angle of scattering direction (deg)
- `s::Float64`: mean slope (m), is square root of mss (mean squared slope)
- `θᵥ::Float64: soil moisture volume fraction [0,1]

"""
function BRCS_KA(λ,θ,ϕ,θₛ,ϕₛ,s,θᵥ)
    # TODO Replace surface input values with input_parameters struct?
    k = 2*pi / λ
    f = 299792458 / λ  # approx center freq derived from free space λ

    q_x,q_y,q_z = get_q_vec(k,θ,ϕ,θₛ,ϕₛ)
    q = sqrt(q_x^2 + q_y^2 + q_z^2)
    

    ϵ₂ = soil_dielectric(θᵥ)
    μ₂ = 12.566e-7 # free space permeability
    
    #local surface normal is nhat defined as =-xhat*Z_x - yhat*Z_y + zhat
    # let Z_x = 0, Z_y = 0?
    m = s/sqrt(2)
    Z_x = m; Z_y = m; # these are the slopes (?)

    U_vh,U_hv,U_hh,U_vv = get_pol_factors_U_pq(ϵ₂,μ₂,θ,θₛ,ϕ,ϕₛ,f,Z_x,Z_y) # no idea how to define Z right now
    
    # p_dp0 is p''(0), the second order derivative of surface height pdf evaluated at zero ??
    σʳ_vh  = (k*q* abs(U_vh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_hv  = (k*q* abs(U_hv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_vv  = (k*q* abs(U_vv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_hh  = (k*q* abs(U_hh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    #----------------------------------------------------

    return  σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh
end#function

"""
# this function calculates the normalized bistatic radar crossection (nbrcs) of a rough surface using KA-GO for vh,hv,vv,hh polatizations in backscattering direction only

# Arguments
- `λ::Float64`: wavelength
- `θ::Float64: incidence angle of incident wave (deg)
- `σ::Float64`: standard deviation of surface heights (m)
- `l::Float64`: surface correlation length (m)
- `θᵥ::Float64: soil moisture volume fraction [0,1]

"""
function BRCS_KA_backscatter(θ,σ,l,θᵥ)
    # TODO Replace surface input values with input_parameters struct?
   
    s = 2 * σ / l
    ϵ₁ = 8.854e-12 
    ϵ₂ = soil_dielectric(θᵥ)

    Γ  = abs((sqrt(ϵ₂)-sqrt(ϵ₁))/(sqrt(ϵ₂)+sqrt(ϵ₁)))^2 # at normal incidence
    σʳ_vv = Γ .* exp.(-(tand.(θ)./s).^2) ./ (cosd.(θ).^4*s.^2)
    σʳ_hh = σʳ_vv 
    σʳ_vh = 0
    σʳ_hv = 0


    return  σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh
end

"""
# this function calculates the normalized bistatic radar crossection (nbrcs) of a the coherent component in the nadir direction

# Arguments
- `λ::Float64`: wavelength
- `θ::Float64: incidence angle of incident wave (deg)
- `θᵥ::Float64: soil moisture volume fraction [0,1]
- `σ::Float64`: standard deviation of surface heights (m)
- `patch_area::Float64: resolution cell area (m^2)

"""
# function RCS_coherent(σ,θ,θᵥ,λ,patch_area) 
#     #going to assume a homogeneous infinite surface approximation for coherent reflection component. Only exists in angular area around specular direction, 
#     # RCS is going to be approximately area of patch divided by First Fresnel zone size * total reflected power (infinite surface assumption)
#     ϵ₁ = 8.854e-12 
#     ϵ₂ = soil_dielectric(θᵥ)
#     ϵᵣ =  ϵ₂/ϵ₁ # assumes free space into soil reflection
#     Γ  = abs((sqrt(ϵ₂)-sqrt(ϵ₁))/(sqrt(ϵ₂)+sqrt(ϵ₁)))^2 # at normal incidence

#     σ₀ = Γ * 4*pi*patch_area/λ^2 * ϵᵣ * exp(-16*pi^2*(cosd(θ)*σ/λ)^2*ϵᵣ) # this equation comes from Ilgin's REASON document. used for nadir (but added incidence angle)

#     σᶜ_vh = 0;  σᶜ_hv = 0;
#     σᶜ_vv = σ₀
#     σᶜ_hh  = σ₀
# end#function

"""
# this function calculates the  bistatic radar crossection (brcs) of a the coherent component in the forward direction

# Arguments
- `σ::Float64`: standard deviation of surface heights (m)
- `θᵥ::Float64: soil moisture volume fraction [0,1]
- `λ::Float64`: wavelength
- `θ::Float64: incidence angle of incident wave (deg)
- `ϕ::Float64: azimuth angle of incident wave (deg)
- `θₛ::Float64: incidence angle of scattered wave (deg)
- `ϕₛ::Float64: azimuth angle of scattered wave (deg)

"""

function RCS_coherent(σ,θᵥ,λ,θ,ϕ,θₛ,ϕₛ) #Incidence angle dependent version here
    #going to assume a homogeneous infinite surface approximation for coherent reflection component. Only exists in angular area around specular direction, 
    # RCS is going to be approximately area of patch divided by First Fresnel zone size * total reflected power (infinite surface assumption)
    f = 299792458 / λ
    ϵ₁ = 8.854e-12 
    μ₂ = 12.566e-7 # free space permeability
    ϵ₂ = soil_dielectric(θᵥ)
    ϵᵣ =  ϵ₂/ϵ₁ # assumes free space into soil reflection
    Γ  = abs((sqrt(ϵ₂)-sqrt(ϵ₁))/(sqrt(ϵ₂)+sqrt(ϵ₁)))^2 # at normal incidence
    k = 2*pi/λ
    q_x,q_y,q_z = get_q_vec(k,θ,ϕ,θₛ,ϕₛ)
    
    # tol = 2  # angular tolerance
    # # Check if the scattered signal direction is in the specular direction
    # if theta == 0 && theta_s == 0 # checks if nadir (then azimuth doesn't matter) or if in specular direction
    #     σ₀ = abs(Γ)^2 * k * dirac(q_x) * dirac(q_y) * exp(-q_z^2 * σ^2)   
    #     σᶜ_vv = σ₀
    #     σᶜ_hh  = σ₀
    # elseif abs(θₛ - θ) < tol && (abs(ϕₛ - ϕₛ - π) < tol) # if off nadir, checks for forward scattering direction
    #         σ₀ = abs(Γ)^2 * k * dirac(q_x) * dirac(q_y) * exp(-q_z^2 * σ^2)   
    #         σᶜ_vv = σ₀
    #         σᶜ_hh  = σ₀
    #     else
    #         σᶜ_vv = 0
    #         σᶜ_hh  = 0
    #     end
    # end
    # σ₀ = abs(Γ)^2 * k * dirac(q_x) * dirac(q_y) * exp.(-q_z.^2 .* σ.^2)   
    R_pll0, R_pll1, R_prp0, R_prp1 = get_Fresnel_coeffs(ϵ₂,μ₂,θ,f)
    a₀_hh = -R_prp0*(cosd(θ)+cosd(θₛ))*cosd(ϕ-ϕₛ)
    a₀_vv = -R_pll0*(cosd(θ)+cosd(θₛ))*cosd(ϕ-ϕₛ)
    σᶜ_vv = π * k^2 * abs(a₀_vv)^2 *  dirac(q_x) * dirac(q_y) * exp.(-q_z.^2 .* σ.^2)   
    σᶜ_hh = π * k^2 * abs(a₀_hh)^2 *  dirac(q_x) * dirac(q_y) * exp.(-q_z.^2 .* σ.^2)   
    σᶜ_vh = 0.0;  
    σᶜ_hv = 0.0;
        
     return σᶜ_vh,  σᶜ_hv,  σᶜ_vv,  σᶜ_hh
end#function

# this function calculates the scattering vector components q based on incidence and reflection angles
function get_q_vec(k,θ,ϕ,θₛ,ϕₛ)
    q_x = k*(sind(θₛ)*cosd(ϕₛ) - sind(θ)*cosd(ϕ)) 
    q_y = k*(sind(θₛ)*sind(ϕₛ) - sind(θ)*sind(ϕ)) #note can write qₓ but q_y in same subscript is forbidden in Julia....
    q_z = k*(cosd(θₛ) + cosd(θ))
    return q_x,q_y,q_z
end#function

function dirac(q) # simple dirac-delta function
    if q == 0
        return 1.0
    else 
        return 0.0
    end
end

"""
# this function calculates the polarization factors U_pq for pq = {vh hv vv hh}

# Arguments
- `ϵ₂::Float64`: permeability of medium 2 (ground)
- `μ₂::Float64:  permittivity of medium 2 (ground)
- `θ::Float64`:  incidence angle (deg)
- `f::Float64`:  frequency
- `Z_x::Float64`: slope in x-direction
- `Z_y::Float64`: slope in y-direction

"""
function get_pol_factors_U_pq(ϵ₂,μ₂,θ,θₛ,ϕ,ϕₛ,f,Z_x,Z_y)
    R_pll0, R_pll1, R_prp0, R_prp1 = get_Fresnel_coeffs(ϵ₂,μ₂,θ,f)
    U_vh = -R_prp0* (1+cosd(θ)*cosd(θₛ)) *sind(ϕₛ-ϕ) - (R_prp0 * sind(θ)*cosd(θₛ) + R_prp1*(1+cosd(θ)*cosd(θₛ)))* sind(ϕₛ-ϕ)*(Z_x*cosd(ϕ) + Z_y*sind(ϕ))
    U_hv = -R_pll0* (1+cosd(θ)*cosd(θₛ)) *sind(ϕₛ-ϕ) + (R_pll0 * sind(θ)*cosd(θₛ)+ R_pll1*(1+cosd(θ)*cosd(θₛ)))* sind(ϕₛ-ϕ)*(Z_x*cosd(ϕ) + Z_y*sind(ϕ))
    U_hh = -R_prp0* (cosd(θ)+cosd(θₛ))*cosd(ϕₛ-ϕ) + (R_prp0*(sind(θₛ) - sind(θ)*cos(ϕₛ-ϕ))- R_prp1*(cosd(θ)+cosd(θₛ)) * cosd(ϕₛ-ϕ)) * (Z_x*cosd(ϕ) + Z_y*sind(ϕ))
    U_vv = -R_pll0* (cosd(θ)+cosd(θₛ))*cosd(ϕₛ-ϕ)+ (R_pll0*(sind(θₛ) - sind(θ)*cos(ϕₛ-ϕ)) + R_pll1*(cosd(θ)+cosd(θₛ)) * cosd(ϕₛ-ϕ)) * (Z_x*cosd(ϕ) + Z_y*sind(ϕ))

    return U_vh,U_hv,U_hh,U_vv
end#function

#this function finds the zeroth and first order Taylor series Fresnel reflection coeffs
function get_Fresnel_coeffs(ϵ₂,μ₂,θ,f)
    # find instrinsic impedances
    # assumes free space in upper half-space
    μ₁ = 12.566e-7 
    ϵ₁ = 8.854e-12 
    η₁ = sqrt(μ₁/ϵ₁)
    k₁ = 2* π * f * sqrt(ϵ₁*μ₁)

    η₂ = sqrt(μ₂/ϵ₂)
    k₂ = 2* π * f * sqrt(ϵ₂*μ₂) 
    θₜ = k₁*sind(θ)/k₂ # transmission angle derived from Snell's law

    R_pll0 = (η₁*cosd(θ) - η₂*cosd(θₜ)) / (η₁*cosd(θ) + η₂*cosd(θₜ))
    R_pll1 = - (η₁*sind(θ) - η₂*sind(θₜ) - R_pll0 * (η₁*sind(θ)+ η₂*sind(θₜ)) ) / (η₁*cosd(θ) + η₂*cosd(θₜ))
    R_prp0 = (η₂*cosd(θ) - η₁*cosd(θₜ))/(η₂*cosd(θ) + η₁*cosd(θₜ))
    R_prp1 = -R_prp0 * (η₂*sind(θ) + η₁*sind(θₜ)) / (η₂*cosd(θ) + η₁*cosd(θₜ))
    return R_pll0, R_pll1, R_prp0, R_prp1
end#function

"""
# this function calculates dielectric of soil for a given water content

# Arguments
- `θᵥ::Float64`: soil water content [0,1]
- `f::Float64`: frequency

"""
function soil_dielectric(θᵥ)
    #using the Topp 1980 model modified with Wang, Schmugge, Williams 1978
    # K = K' + j [K'' + (σᵪ/ω*ϵₒ)]
    #where K is complex dielectric, K' is the real part, K'' is the imaginary part, and ϵₒ = free space permittivity. Assumes magnetic properties of free space
    # this is an empirical model derived from measurements from 20 MHz to 1 GHz. Take care using for freq >>1 GHz

    # Combining two models for real (K') and imaginary parts in total, so K = K' + K''
    # ω = f*2*π
    ϵₒ = 8.854e-12 
  
    # real part from Topp
    K_prime = 3.03 + 9.3*θᵥ + 146 * θᵥ^2 - 76.7*θᵥ^3
    
    #Imaginary part from Wang, Schmugge, Williams 1978 -----------
    # ϵₛ = 75 
    # ϵᵢ = 5 
    # σᵢ = 1
    # τ = 1e-11
    # K_double_prime = (ω*τ*(ϵₛ - ϵᵢ)/(1+ω^2*τ^2) + σᵢ/(ω*ϵₒ)) / 5
    #------------

    #Imaginary part taken from Wang + Schmugge measurements, curve fitting polynomial through data. Imaginary part mostly insensitive to frequency
    K_double_prime =  16.42*θᵥ^2 + 1.616*θᵥ + .1379
    
    K = Complex(K_prime * ϵₒ, K_double_prime * ϵₒ)
    return K
end#function


"""
# this function calculates bistatic radar crossection using Small Perturbation Method 

# Arguments
- `λ::Float64`: wavelength
- `θᵢ::Float64: incidence angle of incident wave (deg)
- `ϕᵢ::Float64`: azimuth angle of incident wave (deg)
- `θₛ::Float64`: incidence angle of scattering direction (deg)
- `ϕₛ::Float64`: azimuth angle of scattering direction (deg)
- `l::Float64`: correlation length
- `σ::Float64`: stdev of surface heights (m)
- `θᵥ::Float64: soil moisture volume fraction [0,1]

"""
function BRCS_SPM_tsang(λ,θᵢ,ϕᵢ,θₛ,ϕₛ,l,σ,θᵥ)

    ϕₛ = ϕₛ-180#testing flipped axis

    k = 2*pi / λ
    f = 299792458 / λ  # approx center freq derived from free space λ
    σ² = σ^2 # create sigma squared notation
    # l = 3 # correlation length 

    μ₁ = 12.566e-7 
    μ₂ = μ₁
    ϵ₁ = 8.854e-12 
    η₁ = sqrt(μ₁/ϵ₁)
    k = 2* π * f * sqrt(ϵ₁*μ₁)

    ϵ₂ = soil_dielectric(θᵥ)
    η₂ = sqrt(μ₂/ϵ₂)
    k_1 = 2* π * f * sqrt(ϵ₂*μ₂)  
    k_1zi = sqrt(k_1^2 - k^2 * sind(θᵢ)^2)
    k_1z =k_1 *cosd(θₛ) # this one is a guess, needs to be confirmed
    k_zi = k * cosd(θᵢ)

    ## WHAT IS k_z????
    k_z = (μ₂/μ₁)*cosd(θₛ)

    

    k_dp = k.^2 * ( sind(θₛ).^2 .+ sind(θᵢ).^2 .- 2 .*sind(θₛ).*sind(θᵢ).*cosd(ϕₛ-ϕᵢ) )
    
    pre_polarized_field = (4 .* k.^2 .* σ² .* l.^2 .* cosd(θₛ).^2 * cosd(θᵢ).^2) / cosd(θᵢ) .* exp(-1/4 * k_dp ^.2 * l.^2)
    # f_ba = 4 different polarization factors: hh,vv,hv,vh
    f_hh = abs( (k_1.^2-k.^2) / ( (k_z + k_1z) * (k_zi + k_1zi) ) ).^2 * cosd(ϕₛ-ϕᵢ).^2
    f_vh = abs( (k_1.^2-k.^2)* k * k_1z / ( (k_1.^2 * k_z + k.^2 * k_1z) * (k_zi + k_1zi) ) ).^2 * sind(ϕₛ-ϕᵢ).^2
    f_hv = abs( (k_1.^2-k.^2)* k * k_1zi / ( (k_z + k_1z) * (k_1.^2 * k_zi + k.^2 * k_1zi) ) ).^2 * sind(ϕₛ-ϕᵢ).^2
    f_vv = abs( (k_1.^2-k.^2) / ((k_1.^2*k_z + k^2*k_1z)*(k_1.^2*k_zi + k.^2*k_1zi)) * (k_1.^2 *k.^2 .*sind(θₛ).*sind(θᵢ) -k.^2 * k_1z * k_1zi .*cosd(ϕₛ-ϕᵢ)) ).^2
    
    σʳ_vh = pre_polarized_field .* f_vh
    σʳ_vv = pre_polarized_field .* f_vv
    σʳ_hv = pre_polarized_field .* f_hv
    σʳ_hh = pre_polarized_field .* f_hh

    return σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh
end

function BRCS_SPM_backscatter(λ,θᵢ,l,σ,θᵥ)
    k = 2*pi / λ
    f = 299792458 / λ  # approx center freq derived from free space λ
    σ² = σ^2 # create sigma squared notation
    # l = 3 # correlation length 

    μ₁ = 12.566e-7 
    ϵ₂ = soil_dielectric(θᵥ)
    ϵ₁ = 8.854e-12 
    # η₁ = sqrt(μ₁/ϵ₁)
    k = 2* π * f * sqrt(ϵ₁*μ₁)
    # η₂ = sqrt(μ₂/ϵ₂)
    ϵᵣ = ϵ₂/ϵ₁

    Γ  = abs((sqrt(ϵ₂)-sqrt(ϵ₁))/(sqrt(ϵ₂)+sqrt(ϵ₁)))^2 # at normal incidence
    σʳ_hh = Γ * 4 * k^4 * (l*σ)^2 * cosd(θᵢ)^4 * exp(-1*(k*l*sind(θᵢ))^2)
    σʳ_vv = 4 * k^4 * (l*σ)^2 * cosd(θᵢ)^4 * exp(-1*(k*l*sind(θᵢ))^2) * abs( (ϵᵣ-1)*(sind(θᵢ)^2-ϵᵣ*(1+sind(θᵢ)^2)) / (ϵᵣ*cosd(θᵢ) + sqrt(ϵᵣ-sind(θᵢ)^2))^2     )^2
    σʳ_vh = 0
    σʳ_hv = 0
    return σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh
end#function



"""
# this function calculates the emperically-derived backscatter values from the Ulaby-Dobson book "Handbook of Radar Scattering Statistics for Terrain"

# Arguments
- `θ::Float64: incidence angle (deg)
- `band::String`: frequency band 
- `terrain::String`: terrain type (from set of available)

# Output
- `σ ::Float64: Backscatter RCS for HH, HV, or VV (dB). Book gives HV = VH
"""
function Ulaby_book_terrain_backscatter_values(θ,fc,terrain)
    if fc >= 1e9 && fc < 2e9
        band = "L"
    elseif fc >= 2e9 && fc < 4e9
        band = "S"
    elseif fc >= 4e9 && fc < 8e9
        band = "C"
    else
        error(raw"Frequency out of bounds")
    end

    θ_rad = pi/180 .* θ # convert to radians

    if terrain == "soil"   
        if band == "L" # HH, HV, VV coeffs in tables
            P₁ = [2442 -30.2 -94.36] # recalculated HH coefficients because book values gave wrong curve
            P₂ = [-3.771 15.261 99.0] 
            P₃ = [-2.193 3.56 0.365]
            P₄ = [-2442 -0.424 -3.398]
            P₅ = [ 0.1711 0.0 5.0]
            P₆ = [-0.07967 0.0 -1.739]
        elseif band == "S"
            P₁ = [-91.2 -46.467 -97.016]
            P₂ = [99.0 31.788 99.0]
            P₃ = [0.433 2.189 0.270]
            P₄ = [5.063 -17.99 -2.056]
            P₅ = [2.941 1.34 5.0]
            P₆ = [-3.142 1.583 -1.754]
        elseif band == "C"
            P₁ = [-24.855 -26.7 -24.951]
            P₂ = [26.351 15.055 28.742]
            P₃ = [1.146 1.816 1.043]
            P₄ = [0.204 -0.499 -1.681]
            P₅ = [0.0 0.0 0.0]
            P₆ = [0.0 0.0 0.0]
        end#band

    elseif terrain == "grass" 
        if band == "L"
            P₁ = [-29.325 -40.166 -28.022]
            P₂ = [37.55 26.833 36.59]
            P₃ = [2.332 2.029 2.53]
            P₄ = [-2.615 -1.473 -1.53]
            P₅ = [5.0 3.738 5.0]
            P₆ = [-1.616 -1.324 -1.513]
        elseif band == "S"
            P₁ = [-20.361 -29.035 -21.198]
            P₂ = [25.727 18.055 26.694]
            P₃ = [2.979 2.8 2.828]
            P₄ = [-1.13 -1.556 -0.612]
            P₅ = [5.0 4.554 5.0]
            P₆ = [-1.916 -0.464 -2.079]
        elseif band == "C"
            P₁ = [-15.75 -23.109 -93.606]
            P₂ = [17.931 13.591 99]
            P₃ = [2.369 1.508 0.22]
            P₄ = [-1.502 -.757 -5.509]
            P₅ = [4.592 4.491 -2.964]
            P₆ = [-3.142 -3.142 1.287]
        end#band

    elseif terrain == "short_veg" 
        if band == "L"
            P₁ = [-27.265 -41.6 -24.614]
            P₂ = [32.39 22.872 27.398]
            P₃ = [2.133 0.689 2.265]
            P₄ = [1.438 -1.238 -1.080]
            P₅ = [-3.847 0.0 5.0]
            P₆ = [3.142 0.0 -1.999]
        elseif band == "S"
            P₁ = [-20.779 -99.0 -20.367]
            P₂ = [21.867 85.852 21.499]
            P₃ = [2.434 0.179 2.151]
            P₄ = [0.347 3.687 -1.069]
            P₅ = [-0.013 2.121 5.0]
            P₆ = [-0.393 -3.142 -1.95]
        elseif band == "C"
            P₁ = [-87.727 -99.0 -88.593]
            P₂ = [99.0 93.293 99.0]
            P₃ = [0.322 0.181 0.326]
            P₄ = [10.188 5.359 9.574]
            P₅ = [-1.747 1.948 1.969]
            P₆ = [3.143 -3.142 -3.142]
        end#band
    elseif terrain == "shrubs" 
        if band == "L"
            P₁ = [-26.688 -99.0 -81.371]
            P₂ = [29.454 99.0 99.0]
            P₃ = [1.814 0.086 0.567]
            P₄ = [0.873 -21.298 16.2]
            P₅ = [4.135 0.0 -1.948]
            P₆ = [-3.142 0.0 3.142]
        elseif band == "S"
            P₁ = [-21.202 -89.222 -20.566]
            P₂ = [99.0 91.002 99.0]
            P₃ = [0.270 0.156 0.294]
            P₄ = [6.98 3.948 8.107]
            P₅ = [-5.0 -0.355 5.0]
            P₆ = [-3.141 0.526 -1.983]
        elseif band == "C"
            P₁ = [-91.95 -99.0 -91.133]
            P₂ = [99.0 91.003 99.0]
            P₃ = [0.270 0.156 0.294]
            P₄ = [6.980 3.948 8.107]
            P₅ = [1.922 2.239 2.112]
            P₆ = [-3.142 -3.142 -3.142]
        end#band    

    elseif terrain == "dry_snow" 
        if band == "L"
            P₁ = [-74.019 -91.341 -77.032]
            P₂ = [99.0 99.0 99.0]
            P₃ = [1.592 1.202 1.415]
            P₄ = [-30.0 30.0 -30.0]
            P₅ = [1.928 1.790 1.720]
            P₆ = [0.905 -2.304 0.997]
        elseif band == "S"
            P₁ = [-47.055 -54.29 -40.652]
            P₂ = [30.164 13.292 18.826]
            P₃ = [5.788 10.0 9.211]
            P₄ = [30.0 -30.0 30.0]
            P₅ = [1.188 -0.715 0.690]
            P₆ = [-0.629 3.142 0.214]
        elseif band == "C"
            P₁ = [-42.864 -25.543 -19.765]
            P₂ = [20.762 16.64 19.83]
            P₃ = [10.0 10.0 10.0]
            P₄ = [30.0 -2.959 7.089]
            P₅ = [0.763 3.116 1.54]
            P₆ = [-0.147 2.085 -0.012]
        end#band

    elseif terrain == "wet_snow" 
        if band == "L"
            P₁ = [-73.069 -90.98 -75.156]
            P₂ = [95.221 99.0 99.0]
            P₃ = [1.548 1.129 1.446]
            P₄ = [30.0 30.0 30.0]
            P₅ = [1.795 1.827 1.793]
            P₆ = [-2.126 -2.308 -2.179]
        elseif band == "S"
            P₁ = [-45.772 -42.940 -39.328]
            P₂ = [25.16 9.935 18.594]
            P₃ = [5.942 15.0 8.046]
            P₄ = [30.0 30.0 30.0]
            P₅ = [0.929 0.438 0.666]
            P₆ = [-0.284 0.712 0.269]
        elseif band == "C"
            P₁ = [-31.91 -24.622 4.288]
            P₂ = [17.749 15.102 15.642]
            P₃ = [11.854 15.0 15.0]
            P₄ = [30.0 -3.401 30.0]
            P₅ = [0.421 2.431 0.535]
            P₆ = [0.74 3.142 1.994]
        end#band
    else
        @error(raw"Terrain type not supported")
    end#terrain
    
    σᵘ_hh = P₁[1] .+ P₂[1] .* exp.(-P₃[1].*θ_rad) .+ P₄[1] .* cos.(P₅[1] .* θ_rad .+ P₆[1])
    σᵘ_hv = P₁[2] .+ P₂[2] .* exp.(-P₃[2].*θ_rad) .+ P₄[2] .* cos.(P₅[2] .* θ_rad .+ P₆[2])
    σᵘ_vv = P₁[3] .+ P₂[3] .* exp.(-P₃[3].*θ_rad) .+ P₄[3] .* cos.(P₅[3] .* θ_rad .+ P₆[3])
    σᵘ_vh = σᵘ_hv #book only gives vh results, assuming hv is similar
    
    # convert from dB to linear units - use real component, something funky with Ulaby coefficients and equation. Checked to ensure db->linear->dB conversion correct here
    σᵘ_hh = 10 .^(real.(σᵘ_hh) ./10)
    σᵘ_hv = 10 .^(real.(σᵘ_hv) ./10)
    σᵘ_vv = 10 .^(real.(σᵘ_vv) ./10)
    σᵘ_vh = 10 .^(real.(σᵘ_vh) ./10)

    return  σᵘ_vh,  σᵘ_hv,  σᵘ_vv,  σᵘ_hh
end#function


end#module