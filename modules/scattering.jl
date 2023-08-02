## scattering module
module scattering
using Parameters, LinearAlgebra, Random


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
- `elev::N-element Array`: elev angle (deg)
- `phi::N-element Array`: azimuth angle (deg)
- `range::N-element Array`: range from origin to spherical point
"""
function get_surface_brcs(params,tx_ecef,rx_ecef,tgt_ecef,patch_area)
    @unpack λ, polarization = params

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
    # tx_llh = [0; 0; 500e3]; tx_llh; tgt_llh = [0;10;0] # monostatic
    # tx_llh = [0; 0; 500e3]; rx_llh = [.1; 0; 500e3]; tgt_llh = [0;10;0] # quasi-monostatic (?)

    tx_enu = Geometry.llh_to_enu_new_org(tx_llh[1], tx_llh[2], tx_llh[3], tgt_llh)
    rx_enu = Geometry.llh_to_enu_new_org(rx_llh[1], rx_llh[2], rx_llh[3], tgt_llh)
    tx_sph = Geometry.enu_to_sph(tx_enu[1], tx_enu[2], tx_enu[3])
    rx_sph = Geometry.enu_to_sph(rx_enu[1], rx_enu[2], rx_enu[3])
    #correct transmitter azimuth for change of origin direction (origin calculated as target location)
    tx_sph[3] = tx_sph[3] + 180; tx_sph[3] = tx_sph[3]%360 # add 180 degrees and mod360 to enure 0 < azimuth < 360

    # check for "monostatic" bounds -- TODO: How to define what we approxiamte as monostatic? 
    # i.e. what scattering error tolerance is acceptable --> angular difference from direct backscatter that meets tolerance
    #for now we'll say 2 degrees in both azimuth and incidence...
    tol = 2
    if (abs(tx_sph[2] - rx_sph[2]) < tol) & ( abs(  tx_sph[3] - rx_sph[3] + 180 ) < tol)

        monostatic = true
    else
        monostatic = false
    end



    if monostatic #monostatic geometry
        # Auto-select surface scattering mechanism based on surface properties
        m = 2 * σ / l
        if λ*σ < l^2/2.76
            #use GO
            σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA_backscatter(θ,σ,l,θᵥ)
            
        elseif σ < λ/21 & m < 0.3
            #use SPM
            σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_backscatter(λ,θᵢ,l,σ,θᵥ)
            error(raw"Scattering surface not supported")
    else # generalized bistatic geometry
         # Auto-select surface scattering mechanism based on surface properties
         m = 2 * σ / l
         if λ*σ < l^2/2.76
             #use GO
             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_KA(λ,θ,ϕ,θₛ,ϕₛ,σ,l,θᵥ)
             
         elseif σ < λ/21 & m < 0.3
             #use SPM
             σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh = BRCS_SPM_tsang(λ,θᵢ,ϕᵢ,θₛ,ϕₛ,l,σ,θᵥ)
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


    # add a random phase because we are looking at non-coherent scattering. Need to consider phase stability over motion?
    # s= 2*pi; phase = s*rand()
    phase = 0
    brcs_complex = brcs*cos(phase)+ (1im*brcs*sin(phase))
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
    #mss = 0.01 # ?? Not sure what value this should be
    # mss = (σ²/2/pi)
    s = 2 * σ / l

    #local surface normal is nhat defined as =-xhat*Z_x - yhat*Z_y + zhat
    # let Z_x, Z_y be Gaussian random numbers with  zero-mean and variance = mss
    numInstances = 100
    Z_x = randn(numInstances).*sqrt(σ²/2/pi)
    Z_y = randn(numInstances).*sqrt(σ²/2/pi)
    
    # now, we'll need to calculate average brcs based on instanes of Z_x, Z_y
    σʳ_vh = zeros(1,numInstances)
    σʳ_hv = zeros(1,numInstances)
    σʳ_vv = zeros(1,numInstances)
    σʳ_hh = zeros(1,numInstances)
    for i = 1 : numInstances
        U_vh,U_hv,U_hh,U_vv = get_pol_factors_U_pq(ϵ₂,μ₂,θ,f,Z_x[i],Z_y[i]) # no idea how to define Z right now

        σʳ_vh[i]  = (k*q* abs(U_vh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
        σʳ_hv[i]  = (k*q* abs(U_hv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
        σʳ_vv[i]  = (k*q* abs(U_vv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
        σʳ_hh[i]  = (k*q* abs(U_hh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    end
    #find mean BRCS values
    σʳ_vh = sum(σʳ_vh)/numInstances
    σʳ_hv = sum(σʳ_hv)/numInstances
    σʳ_vv = sum(σʳ_vv)/numInstances
    σʳ_hh = sum(σʳ_hh)/numInstances
    
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

    U_vh,U_hv,U_hh,U_vv = get_pol_factors_U_pq(ϵ₂,μ₂,θ,f,Z_x,Z_y) # no idea how to define Z right now
    
    # p_dp0 is p''(0), the second order derivative of surface height pdf evaluated at zero ??
    σʳ_vh  = (k*q* abs(U_vh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_hv  = (k*q* abs(U_hv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_vv  = (k*q* abs(U_vv) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))
    σʳ_hh  = (k*q* abs(U_hh) )^2 / ( (q_z)^4 * s.^2 ) * exp(-((q_x)^2 + (q_y)^2) / ((q_z)^2 * s.^2 ))

    return  σʳ_vh,  σʳ_hv,  σʳ_vv,  σʳ_hh
end#function

"""
# this function calculates the normalized bistatic radar crossection (nbrcs) of a rough surface using KA-GO for vh,hv,vv,hh polatizations in backscattering direction only

# Arguments
- `λ::Float64`: slow time vector
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

# this function calculates the scattering vector components q based on incidence and reflection angles
function get_q_vec(k,θ,ϕ,θₛ,ϕₛ)
    q_x = k*(sind(θₛ)*cosd(ϕₛ) - sind(θ)*cos(ϕ)) 
    q_y = k*(sind(θₛ)*sind(ϕₛ) - sind(θ)*sin(ϕ)) #note can write qₓ but q_y in same subscript is forbidden in Julia....
    q_z = k*(cos(θₛ) + cosd(θ))
    return q_x,q_y,q_z
end#function


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
function get_pol_factors_U_pq(ϵ₂,μ₂,θ,f,Z_x,Z_y)
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
    ω = f*2*π
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

end#module