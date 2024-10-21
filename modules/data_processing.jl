module Data_Processing

#using ..Waveform

using FFTW
using LinearAlgebra
using Statistics
using Interpolations
using Plots
using Parameters
using FIRInterp
using ImageFiltering


##
c = 299792458 # speed of light, used for computations in this module

##
"""
This function performs range compression on received signals for different simulation modes.
It utilizes the range_compression function from the waveform.jl module.

Arguments:
- rawdata: The received signal data
- RCRF_F: Range compression filter in frequency domain
- params_densesim: Parameters for dense simulation
- Np: Number of platforms
- Nst: Number of slow times
- Ncr: Number of cross range

Returns:
- RC_signal: Range compressed signal data
"""
function generate_rangecompressed_received_data(rawdata, RCRF_F, params_densesim, Np, Nst, Ncr)
    # Unpack necessary parameters from params_densesim
    @unpack waveform_domain_process, mode = params_densesim 

    # Initialize RC_signal based on the mode
    if mode == 1 || mode == 2 # SAR (ping-pong) or SIMO
        RC_signal = zeros(ComplexF64, Nst, Np, Ncr)
    elseif mode == 3 # MIMO
        RC_signal = zeros(ComplexF64, Nst, Np, Np, Ncr)
    end
    
    # Perform range compression based on different modes
    if mode == 1 || mode == 2
        # Loop through platforms and slow times for SAR or SIMO mode
        for idx_p = 1:Np
            for idx_st = 1:Nst
                # Apply range compression to the received signal data
                RC_signal[idx_st, idx_p, :] = transpose(Waveform.range_compression(fft(rawdata[idx_st, idx_p, :]), RCRF_F, waveform_domain_process))
            end
        end
    elseif mode == 3
        # Loop through platforms and slow times for MIMO mode 
        for idx_p1 = 1:Np
            for idx_p2 = 1:Np
                for idx_st = 1:Nst
                    # Apply range compression to the received signal data
                    RC_signal[idx_st, idx_p1, idx_p2, :] = transpose(Waveform.range_compression(fft(rawdata[idx_st, idx_p1, idx_p2, :]), RCRF_F, waveform_domain_process))
                end
            end
        end
    end

    return RC_signal  # Return the range-compressed signal data
end


"""
Function for performing range Doppler correction of the received range compressed signal for all platforms
Uses range_doppler_algorithm function from this module

Arguments:
- RC_signal: Range compressed signal data
- Np: Number of platforms
- Nst: Number of slow times
- Ncr: Number of cross range
- Veff: Effective velocity array for platforms
- fDopp: Doppler frequencies
- slrng: Slant range array
- mid_idx: Mid indices for platforms
- max_range_time: Maximum range time
- ref_range: Reference range
- params_densesim: Parameters for dense simulation

Returns:
- RC_RD_signal: Range Doppler corrected signal data
"""
function generate_rangeDoppler_received_data(RC_signal, Np, Nst, Ncr, Veff, fDopp, slrng, mid_idx, max_range_time, ref_range, params_densesim)
    # Unpack necessary parameters from params_densesim
    @unpack fs, fp, fc 	= params_densesim
    
    # Initialize RC_RD_signal
    RC_RD_signal 		= zeros(ComplexF64, Nst, Np, Ncr)
    
    # Loop through platforms to perform range Doppler correction
    for idx_p = 1:Np
        # Apply range Doppler correction using range_doppler_algorithm method 1 function
        RC_RD_signal[:, idx_p, :] = range_doppler_algorithm_method1(RC_signal[:, idx_p, :], Ncr, fc, fs, fp,transpose(Veff[idx_p, :]), Nst, c,  slrng[idx_p, mid_idx[idx_p]], max_range_time, fDopp[idx_p, :])
        # Alternative (commented) method using range_doppler_algorithm method 2 function
        # RC_RD_signal[:, idx_p, :] = range_doppler_algorithm_method2(RC_signal[:, idx_p, :],Ncr,fc,fs,fp,transpose(Veff[idx_p, :]),Nst,c,ref_range,fDopp[idx_p, :])
    end
    
    return RC_RD_signal  # Return the range Doppler corrected signal data
end

##
"""
Function for range Doppler algorithm, used or range curve correction
Note: Two implementations 
1. based on slant ranges (mean in azimuth) 
2. based on fixed range (mean in azimuth and range)

Arguments
 - signal:  Input range compressed signal
 - Ncr: Number of samples along slant range
 - Na: Number of samples along azimuth
 - fc:  Carrier freqency
 - fs:  Sampling freqency
 - PRF:  Pulse repetition freqency
 - vel: Effective velocity
 - ref_range: Reference range
 - fDoppler:, Doppler frequency
 - slrng: Slant ranges
 - max_range_time: Maximum range time
Output
 - srd2: Range corrected data
"""
function range_doppler_algorithm_method1(signal, Ncr, fc, fs, PRF, vel, Na, c, slrng, max_range_time, fDoppler)
	fr   	   		= fftfreq(Ncr,fs).* ones(1,Na);
	fa   	   		= fDoppler;

	slant_range 	= ((1:Ncr).*(c/(fs)) ) .- (max_range_time/2 * c) .+ slrng
	detaR      		= (c^2 .* (slant_range) .* (ones(Ncr,1).*transpose(fa)).^2 ./ fc.^2 ./ 8 ./ vel.^2);
	# another implementation
	#fD 			= sqrt.( 1 .- ( ((c/fc)^2 .* (ones(N,1).*transpose(fa)).^2) ./ (4 .* (vel).^2) ))
	#detaR 			= ref_range .* (1 .- fD)./fD

	srd2 			= zeros(ComplexF64,Na,Ncr)
	for i=1:Na
		srd2[i,:] 	= ifft(fft(signal[i,:]).*exp.(1im .* 4 .* pi .* fr[:,i] .* ((detaR[:,i])/c)))
	end
	return srd2
end

function range_doppler_algorithm_method2(signal, N, fc, fs, PRF, vel, Na, c, ref_range,fDoppler)
	fr   	   		= fftfreq(N,fs).* ones(1,Na);
	fa   	   		= fDoppler;

	detaR      		= (c^2 .* (ref_range) .* (ones(N,1).*transpose(fa)).^2 ./ fc.^2 ./ 8 ./ vel.^2);
	# another implementation
	#fD 			= sqrt.( 1 .- ( ((c/fc)^2 .* (ones(N,1).*transpose(fa)).^2) ./ (4 .* (vel).^2) ))
	#detaR 			= ref_range .* (1 .- fD)./fD

	srd2 			= zeros(ComplexF64,Na,N)
	for i=1:Na
		srd2[i,:] 	= ifft(fft(signal[i,:]).*exp.(1im .* 4 .* pi .* fr[:,i] .* ((detaR[:,i])/c)))
	end
	return srd2
end

"""
Perform azimuth compression of the received signal for all platforms and range samples.

Arguments:
- RC_RD_signal: Range Doppler corrected signal data
- Np: Number of platforms
- Nst: Number of slow times
- Ncr: Number of cross range
- Veff: Effective velocities for platforms
- fDopp: Doppler frequencies 
- ref_fix_range: Reference fixed range
- params_densesim: Parameters for dense simulation

Returns:
- AzC_signal: Azimuth compressed signal data
"""
function generate_azimuthcompressed_received_data(RC_RD_signal, Np, Nst, Ncr, Veff, fDopp, ref_fix_range, params_densesim)
    # Unpack necessary parameters from params_densesim
    @unpack λ, fs, fp, fc 	= params_densesim
    
    # Initialize AzC_signal to store azimuth compressed signal data
    AzC_signal 				= zeros(ComplexF64, Nst, Np, Ncr)
    
    # Loop over each platform
    for idx_p = 1:Np
        # Calculate fD based on parameters
        fD 					= sqrt.(1 .- ((λ^2 .* (fDopp[idx_p, :]).^2) ./ (4 .* (Veff[idx_p, 1]).^2)))
        
        # Calculate azimuth reference for compression
        azi_ref 			= fft(conj(exp.(-1im .* 4 .* pi .* ((1 ./ fD) .* ref_fix_range) .* fc ./ c)))
        
        # Loop over each compression point (range samples)
        for idx_cr = 1:Ncr
            # Perform azimuth compression
            AzC_signal[:, idx_p, idx_cr] = fftshift(ifft(fft(RC_RD_signal[:, idx_p, idx_cr]) .* azi_ref))
        end
    end
    
    return AzC_signal  # Return azimuth compressed signal data
end

"""
Calculates covariance and correlation matrices for the given signal.

Arguments:
- signal: Input signal data (size: A x B x C)
- Np: Number of platforms
- filt_len: Filter length for averaging

Returns:
- Rxx: Covariance matrices (size: A x C x Np x Np)
- Corr: Correlation matrices (size: A x C x Np x Np)
"""
function get_covariance_correlation_matrices_new(signal, Np, filt_len)

	A, B, C 					= size(signal)
	Rxxp                     	= zeros(ComplexF64, A, C, Np);
	Rxx                     	= zeros(ComplexF64, A, C, Np, Np);
    Corr                     	= zeros(ComplexF64, A, C, Np, Np);
	
	# Define a kernel for filtering
	kernel 						= ones(filt_len,filt_len)
	weight 						= kernel ./ sum(kernel)

	# Calculate Rxxp for each platform
	for idx_p=1:Np
		signal_pow 				= (abs.(signal[:,idx_p,:])).^2
		Rxxp[:,:,idx_p] 		= imfilter(signal_pow, weight)
	end

	# Calculate Rxx and Corr matrices
	for idx_p1=1:Np
		for idx_p2=1:Np
			sam_cov_mat 		= signal[:,idx_p1,:] .* conj(signal[:,idx_p2,:])
			#sam_cov_mat_real 	= imfilter(real(sam_cov_mat), weight)
			#sam_cov_mat_imag 	= imfilter(imag(sam_cov_mat), weight)
			#Rxx[:,:,idx_p1,idx_p2] = sam_cov_mat_real .+ (1im .* sam_cov_mat_imag)
			Rxx[:,:,idx_p1,idx_p2] = imfilter(sam_cov_mat, weight)

			# Calculate correlation matrix
			Corr[:,:,idx_p1,idx_p2] = Rxx[:,:,idx_p1,idx_p2] ./ ((Rxxp[:,:,idx_p1] .* Rxxp[:,:,idx_p2]).^ 0.5)
		end
	end

	return Rxx, Corr # Return covariance and correlation matrice

end


"""
Calculates the angle between two vectors a and b in degrees using the dot product.

Arguments:
- a: First input vector
- b: Second input vector

Returns:
- angle: Angle between vectors a and b in degrees
"""
function angle_2vec(a, b)
    # Calculate the dot product of vectors a and b
    dot_product = a⋅b #dot(a, b)
    
    # Calculate the product of the norms of vectors a and b
    product_norms = norm(a) * norm(b)
    
    # Ensure the value passed to acosd is within the valid range [-1, 1]
    cosine_angle = clamp(dot_product / product_norms, -1, 1)
    
    # Calculate the angle in degrees using the inverse cosine function (acosd)
    angle = acosd(cosine_angle)
    
    return angle  # Return the angle between vectors a and b in degrees
end


"""
Generates the steering matrix for beamforming based on the given parameters.

Arguments:
- p_xyz: Platform coordinates (3 x Np x Nst)
- s_xyz_3xN_2D: Target coordinates (3 x azimuth x range)
- azimuth_lim: Azimuth limits for processing
- srange_lim: Range limits for processing
- heights_z: Height values for processing
- Np: Number of platforms
- Master_platform: Index of the master platform
- λ: Wavelength of the signal
- mode: Mode of the signal processing
- processing_mode: Processing mode for steering matrix calculation

Returns:
- aT: Steering matrix (3D Complex matrix: azimuth x range x height x platform)
"""
function get_steering_matrix(p_xyz, s_xyz_3xN_2D, azimuth_lim, srange_lim, heights_z, Np, Master_platform, λ, mode, processing_mode)
    # Initialize the steering matrix
    aT = zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), length(heights_z), Np)
    
    # Iterate based on processing mode
    if processing_mode == 1
        # Processing Mode 1
        for idx_st = azimuth_lim[1]:azimuth_lim[2]
            for idx_r = srange_lim[1]:srange_lim[2]
                for idx_z = 1:length(heights_z)
                    for idx_p = 1:Np
                        # Calculate Va and Vb
                        Va = mean(p_xyz[:, idx_p, :], dims=2) - s_xyz_3xN_2D[:, idx_st, idx_r]
                        Vb = mean(p_xyz[:, Master_platform, :], dims=2) - s_xyz_3xN_2D[:, idx_st, idx_r]
                        
                        # Calculate angle between Va and Vb
                        if cross(Va[:], Vb[:])[3] < 0.0
                            angle_ip = angle_2vec(Va, Vb) * 1
                        else
                            angle_ip = angle_2vec(Va, Vb) * -1
                        end
                        
                        # Calculate steering vector (kaz)
                        if mode == 1
                            kaz = (4 * pi * (angle_ip * pi / 180) * heights_z[idx_z]) / λ
                        elseif mode == 2
                            kaz = (2 * pi * (angle_ip * pi / 180) * heights_z[idx_z]) / λ
                        end
                        
                        # Compute steering vector and assign to aT
                        aT[idx_st, idx_r, idx_z, idx_p] = exp.(1im * kaz)
                    end
                end
            end
        end
    elseif processing_mode == 2
        # Processing Mode 2
        for idx_st = azimuth_lim[1]:azimuth_lim[2]
            for idx_r = srange_lim[1]:srange_lim[2]
                for idx_z = 1:length(heights_z)
                    for idx_p = 1:Np
                        # Calculate Va and Vb
                        Va = mean(p_xyz[:, idx_p + 1, :], dims=2) - s_xyz_3xN_2D[:, idx_st, idx_r]
                        Vb = mean(p_xyz[:, Master_platform, :], dims=2) - s_xyz_3xN_2D[:, idx_st, idx_r]
                        
                        # Calculate angle between Va and Vb
                        if cross(Va[:], Vb[:])[3] < 0.0
                            angle_ip = angle_2vec(Va, Vb) * 1
                        else
                            angle_ip = angle_2vec(Va, Vb) * -1
                        end
                        
                        # Calculate steering vector (kaz)
                        if mode == 1
                            kaz = (4 * pi * (angle_ip * pi / 180) * heights_z[idx_z]) / λ
                        elseif mode == 2
                            kaz = (2 * pi * (angle_ip * pi / 180) * heights_z[idx_z]) / λ
                        end
                        
                        # Compute steering vector and assign to aT
                        aT[idx_st, idx_r, idx_z, idx_p] = exp.(1im * kaz)
                    end
                end
            end
        end
    end
    
    return aT  # Return the generated steering matrix
end


"""
Performs tomographic beamforming to generate a beamformed power matrix.

Arguments:
- Cov_mat: Covariance matrix (4D Complex matrix: azimuth x range x height x platform)
- steering_mat: Steering matrix (4D Complex matrix: azimuth x range x height x platform)
- azimuth_lim: Azimuth limits for processing
- srange_lim: Range limits for processing
- size_op: Size of the output matrix (3-element array)

Returns:
- Pbf: Beamformed power matrix (3D matrix: azimuth x range x height)
"""
function tomo_beamforming(Cov_mat, steering_mat, azimuth_lim, srange_lim, size_op)
    # Initialize the beamformed power matrix
    Pbf = zeros(size_op[1], size_op[2], size_op[3])
    
    # Perform tomographic beamforming
    for idx_st = azimuth_lim[1]:azimuth_lim[2]
        for idx_r = srange_lim[1]:srange_lim[2]
            a = transpose(steering_mat[idx_st, idx_r, :, :])  # Transpose of steering vector
            
            # Calculate beamformed power using conventional beamformer formula
            Pbf[idx_st, idx_r, :] = abs.(diag(conj(transpose(a)) * Cov_mat[idx_st, idx_r, :, :] * a))
        end
    end
    
    return Pbf  # Return the beamformed power matrix
end


"""
Performs tomographic beamforming using the CAPON method to generate a power matrix.

Arguments:
- Cov_mat: Covariance matrix (4D Complex matrix: azimuth x range x height x platform)
- steering_mat: Steering matrix (4D Complex matrix: azimuth x range x height x platform)
- azimuth_lim: Azimuth limits for processing
- srange_lim: Range limits for processing
- size_op: Size of the output matrix (3-element array)

Returns:
- PC: Power matrix generated using CAPON method (3D matrix: azimuth x range x height)
"""
function tomo_CAPON(Cov_mat, steering_mat, azimuth_lim, srange_lim, size_op)
    # Initialize the CAPON power matrix
    PC = zeros(size_op[1], size_op[2], size_op[3])
    
    # Calculate the size of the data matrix
    Data_mat_size = size(Cov_mat)[3]
    
    # Perform tomographic beamforming using CAPON method
    for idx_st = azimuth_lim[1]:azimuth_lim[2]
        for idx_r = srange_lim[1]:srange_lim[2]
            a = transpose(steering_mat[idx_st, idx_r, :, :])  # Transpose of steering vector
            #PC[idx_st,idx_r,:] 		=  abs.(diag(1 ./ (conj(transpose(a)) * inv(Cov_mat[idx_st,idx_r,:,:]) * (a))) )# power of conventional beamformer

            # Calculate the diagonal loading factor
            diagL = std(diag(Cov_mat[idx_st, idx_r, :, :]))
            
            try
                # Calculate the CAPON beamformed power
                PC[idx_st, idx_r, :] = abs.(diag(1 ./ (conj(transpose(a)) * inv(Cov_mat[idx_st, idx_r, :, :] .+ (diagL .* Matrix{Float64}(I, Data_mat_size, Data_mat_size))) * (a))))
            catch
                continue  # If an error occurs, skip and continue
            end
        end
    end
    
    return PC  # Return the CAPON power matrix
end


"""
Converts tomographic coordinates to scene coordinates based on the provided parameters.

Arguments:
- signal: Input signal data (3D matrix: frequency x range x platform)
- heights_t: Target heights
- scene_axis2: Scene axis 2
- scene_axis3: Scene axis 3
- look_angle: Look angle for scene
- left_right_look: Direction of look ("left" or "right")
- plat_height: Platform height

Returns:
- op_signal: Transformed signal in scene coordinates (3D matrix: frequency x range x platform)
"""
function tomocoordinates_to_scenecoordinates(signal, heights_t, scene_axis2, scene_axis3, look_angle, left_right_look, plat_height)
    test_mat_h 		= zeros(size(signal, 2), size(signal, 3))
    test_mat_c 		= zeros(size(signal, 2), size(signal, 3))
    op_signal 		= zeros(size(signal, 1), size(signal, 2), size(signal, 3))

    earth_radius 	= 6378.137e3
    l_inc_angle 	= asind.(sind.(look_angle) * (earth_radius + plat_height) / earth_radius)
    heights_z 		= heights_t .* sind(l_inc_angle)

    # Compute test_mat_c
    for j = 1:size(signal, 2)
        grange 		= heights_z ./ tand(l_inc_angle)
        test_mat_c[j, :] = scene_axis2[j] .+ grange
    end
    test_mat_h = transpose(repeat(heights_z, 1, size(signal, 2)))

    # Interpolate signal for scene axis 2
    for i = 1:size(signal, 1)
        for j = 1:size(signal, 3)
            itp 	= linear_interpolation(test_mat_c[:, j], signal[i, :, j], extrapolation_bc=0)
            op_signal[i, :, j] = itp.(scene_axis2)
        end
    end

    # Interpolate signal for scene axis 3
    for i = 1:size(signal, 1)
        for j = 1:size(signal, 2)
            itp 	= linear_interpolation(test_mat_h[j, :], op_signal[i, j, :], extrapolation_bc=0)
            op_signal[i, j, :] = itp.(heights_t)
        end
    end

    # Reverse signal if look direction is right
    if left_right_look == "right"
        op_signal = op_signal[:, end:-1:1, :]
    end

    return op_signal  # Return transformed signal in scene coordinates
end


"""
Averages 2D data (signal) using a specified filter size.

Arguments:
- signal_ip: Input signal data (3D matrix: x-dimension x y-dimension x platform)
- filt_len: Filter size for averaging

Returns:
- Avg_signal_ip: Averaged signal data (3D matrix: x-dimension x y-dimension x platform)
"""
function average_2D_data(signal_ip, filt_len)
    Avg_signal_ip 		= zeros(size(signal_ip))

    kernel 				= ones(filt_len, filt_len)
    weight 				= kernel ./ sum(kernel)

    for idx_p = 1:size(signal_ip)[3]
        Avg_signal_ip[:, :, idx_p] = imfilter(signal_ip[:, :, idx_p], weight)
    end

    return Avg_signal_ip  # Return the averaged signal data
end


##
"""
Another unwrappimg function
"""
function unwrap!(x, period = 2π)
    y = convert(eltype(x), period)
    v = first(x)
    @inbounds for k = eachindex(x)
        x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
    end
end

##
"""
Function for range migration algorithm (In development) ***

"""
function range_migration_algorithm(signal, Ncr, Nst, fc, fs, fDoppler, vel, c, ref_range )

	frequencyRange = LinRange(params_densesim.fc-(params_densesim.bandwidth/2), params_densesim.fc+(params_densesim.bandwidth/2),Ncr)
	#frequencyRange 	= params_densesim.fc .+ fftfreq(Ncr,params_densesim.fs)
	krange 			= 2 .* (2 .* pi .* frequencyRange) ./ c
	krange 			= repeat(krange,1,Nst)
	kaz 			= 2 .* pi .* (fDopp[1,:]) ./ Veff[1,:]
	kaz 			= repeat(kaz,1,Ncr)
	kazimuth 		= transpose(kaz)
	kx 				= krange.^2 - kazimuth.^2
	kx 				= sqrt.(kx .* (kx .> 0))
	kFilt 			= exp.(1im .* ((-slrng[1,mid_idx[1]] .* krange) + (slrng[1,mid_idx[1]] .* kx)) )
	sdata 			= fftshift(fft(RC_signal[:,1,:]))

	sdata_Filt 			= transpose(sdata) .* kFilt
	ifsmPol 		= transpose(ifft((sdata_Filt)))
	display(heatmap((abs.(ifsmPol[:,2500:2700]))))

	stolt_sdata_Filt = sdata_Filt
	for i=1:Nst
		interp_fn_s = LinearInterpolation((kx[:,i]),(sdata_Filt[:,i]),extrapolation_bc = Line())
		stolt_sdata_Filt[:,i] = interp_fn_s.(krange[:,i])
	end

	azcompresseddata = transpose(ifft(stolt_sdata_Filt))
	display(heatmap((abs.(azcompresseddata[:,2500:2700]))))
end

##
"""
Function for w-k (range migration) algorthm for bistatic case (In development) ***

"""
function wk_bistatic_algorithm(signal, Ncr, Nst, fc, fs, fDoppler, c, vel, Rtx, Rrx)
	#Step 1 -  2D FFT
	sdata 		= fft(signal)

	#Step 2 - Bulk Compression
	trt 		  	= 0
	kcr 			= 40e6/1e-8
	fr   	   	= fftfreq(Ncr,fs)
	fa   	   	= fDoppler
	H_ref = zeros(ComplexF64,Nst,Ncr)
	for i = 1:Nst
		H_ref[i,:] = exp.((1im .* 2 .* pi ./ c .* (fc .+ fr) .* (Rtx[i] +Rrx[i])) .+ (1im .* pi .* (fr.^2) ./ kcr) .+ (1im .* 2 .* pi .* fa[i] .* trt))
	end
	sdata2 = sdata .*H_ref
	test1 = (ifft((sdata2)))
	plotly()
	display(heatmap(abs.(test1[:,2000:3000])))

	#Step 3 - Interpolation along teh frequency axis
	frequencyRange	= fc .+ fftfreq(Ncr,fs) #LinRange(fc-fs/2,fc+fs/2,Ncr)
	kr 		= 2 .*(2 .* pi .* frequencyRange) ./ c
	kr = repeat(kr,1,Nst)
	kaz 	= 2 .* pi .*  fDoppler./ vel
	kaz = transpose(repeat(kaz,1,Ncr))
	w_k 	= c / 2 * sqrt.((kr.^2) .- (kaz.^2))
	kw 		= w_k * 2 / c
	interp_op 		= sdata2
	for i = 1:Nst
		interp_fn 			= LinearInterpolation(kw[:,i], sdata2[:,i], extrapolation_bc = Line())
		interp_op[:,i] 	    = interp_fn.(kr[:,1])
	end
	replace_nan(v) 	= map(x -> isnan(x) ? ones(x) .* 1e-30 : x, v)
	interp_op 		= map(replace_nan, interp_op)

	#Step 4 - Two dimentional IFFT
	azcompresseddata = ifft(ifft(interp_op,2),1)
	return azcompresseddata
end


## Image co-registration functions

"""
Function for co-registering the azimuth compresed received data with respect to the master image/platform
Using Linear Interpolation

"""
function coregister_azimuthcompressed_received_data_interp1(AzC_signal, Master_platform, Np, Nst, Ncr, p_xyz, lookvec, slrng, ref_fix_range, max_range_time, fs, azimuth_lim)

	AzC_signal_Itp 			= zeros(ComplexF64,Nst,Np,Ncr)
	slant_range 			= zeros(Nst, Np, Ncr)
	del_d 					= zeros(Ncr)
	del_val 				= zeros(Ncr)
	fast_t 					= collect((1:Ncr).*(1/(fs)))

	for idx_p = 1:Np
		for idx_st = azimuth_lim[1]:azimuth_lim[2]
			#slant_range[idx_st,idx_p,:] = ((1:Ncr).*(c/(fs)) ) .- (max_range_time/2 * c) .+ slrng[idx_p,idx_st]
			slant_range[idx_st,idx_p,:] = ((1:Ncr).*(c/(fs)) ) .- (max_range_time/2 * c) .+ ref_fix_range
			p_dist_diff 				= (p_xyz[:,idx_p,idx_st].-p_xyz[:,Master_platform,idx_st])
			for idx_ft = 1:Ncr
				del_d[idx_ft] 			= dot(p_dist_diff, (lookvec[:,idx_p,idx_st]./slant_range[idx_st,idx_p,idx_ft]))
			end
			interp_fn 					= LinearInterpolation(fast_t, (AzC_signal[idx_st,idx_p,:]), extrapolation_bc = Line());
			AzC_signal_Itp[idx_st,idx_p,:] 	= interp_fn.(fast_t .+ (abs.(del_d) ./ c))
		end
	end
	return AzC_signal_Itp, slant_range
end

##
"""
Function for co-registering the azimuth compresed received data with respect to the master image/platform
Using Sinc Interpolation

"""
function coregister_azimuthcompressed_received_data_interp2(AzC_signal, Master_platform, Np, Nst, Ncr, p_xyz, lookvec, slrng, ref_fix_range, max_range_time, fs, azimuth_lim)

	AzC_signal_Itp 			= zeros(ComplexF64,Nst,Np,Ncr)
	slant_range 			= zeros(Nst, Np, Ncr)
	del_d 					= zeros(Ncr)
	del_val 				= zeros(Ncr)

	itp 					= KnabInterp(5, 0.6);
	for idx_p = 1:Np
		for idx_st = azimuth_lim[1]:azimuth_lim[2]
			#slant_range[idx_st,idx_p,:] = ((1:Ncr).*(c/(fs)) ) .- (max_range_time/2 * c) .+ slrng[idx_p,idx_st]
			slant_range[idx_st,idx_p,:] = ((1:Ncr).*(c/(fs)) ) .- (max_range_time/2 * c) .+ ref_fix_range
			p_dist_diff 				= (p_xyz[:,idx_p,idx_st].-p_xyz[:,Master_platform,idx_st])
			for idx_ft = 1:Ncr
				del_d[idx_ft] 			= dot(p_dist_diff, (lookvec[:,idx_p,idx_st]./slant_range[idx_st,idx_p,idx_ft]))
				del_val[idx_ft] = (abs.(del_d[idx_ft]) ./ c) .* (fs)
				AzC_signal_Itp[idx_st,idx_p,idx_ft] = interp(itp, AzC_signal[idx_st,idx_p,:], idx_ft .+ del_val[idx_ft])
			end
		end
	end
	return AzC_signal_Itp, slant_range

end


"""
Coregisters azimuth compressed received data using correlation method.
Function for co-registering the azimuth compresed received data with respect to the master image/platform
Using Correlation - subpixel image registration
Algorithm based on:
Guizar-Sicairos, M., Thurman, S.T. and Fienup, J.R., 2008. Efficient subpixel image registration algorithms. Optics letters, 33(2), pp.156-158.
	
Arguments:
- Slc1: First input signal (azimuth compressed)
- Slc2: Second input signal (azimuth compressed)
- Upsamp_fac: Upsampling factor (0: no registration, 1: single pixel registration, >1: refinement)

Returns:
- output: Array containing [error, diffphase, row_shift, col_shift]
- Reg_Slc2: Registered version of Slc2
"""
function coregister_azimuthcompressed_received_data_corr_method(Slc1, Slc2, Upsamp_fac)
    FSlc1 = fft(Slc1)
    FSlc2 = fft(Slc2)
    
    nr, nc = size(FSlc2)
    Nr = ifftshift(-floor(nr/2):ceil(nr/2)-1)
    Nc = ifftshift(-floor(nc/2):ceil(nc/2)-1)
    
    # Initialize variables
    CCmax = 0
    Row_shift = 0
    Col_shift = 0

    if Upsamp_fac == 0
        # Simple computation of error and phase difference without registration
        CCmax = sum(FSlc1[:] .* conj(FSlc2[:]))
    elseif Upsamp_fac == 1
        # Single pixel registration
        CC = ifft(FSlc1 .* conj(FSlc2))
        CCabs = abs.(CC)
        temp_s = indexin(maximum(CCabs[:]), CCabs)
        row_shift = temp_s[1][1]
        col_shift = temp_s[1][2]
        CCmax = CC[row_shift, col_shift] * nr * nc
        row_shift = Nr[row_shift]
        col_shift = Nc[col_shift]
    elseif Upsamp_fac > 1
        # Upsampling factors > 1 (Refinement)
        CC = ifft(FTpad(FSlc1 .* conj(FSlc2), [2 * nr, 2 * nc]))
        CCabs = abs.(CC)
        temp_s = indexin(maximum(CCabs[:]), CCabs)
        row_shift = temp_s[1][1]
        col_shift = temp_s[1][2]
        CCmax = CC[row_shift, col_shift] * nr * nc
        
        # Change shifts to represent relative shifts and not indices
        Nr2 = ifftshift(-floor(nr):ceil(nr)-1)
        Nc2 = ifftshift(-floor(nc):ceil(nc)-1)
        row_shift = Nr2[row_shift] / 2
        col_shift = Nc2[col_shift] / 2
        
        # Refine estimate with matrix multiply DFT
        if Upsamp_fac > 2
            # DFT computation
            row_shift = round(row_shift * Upsamp_fac) / Upsamp_fac
            col_shift = round(col_shift * Upsamp_fac) / Upsamp_fac
            dftshift = Int64(ceil(Upsamp_fac * 1.5) / 2)
            CC = conj(dftups(FSlc2 .* conj(FSlc1), ceil(Upsamp_fac * 1.5), ceil(Upsamp_fac * 1.5), Upsamp_fac, 
                    dftshift - row_shift * Upsamp_fac, dftshift - col_shift * Upsamp_fac))
            CCabs = abs.(CC)
            temp_s = indexin(maximum(CCabs[:]), CCabs)
            rloc = temp_s[1][1]
            cloc = temp_s[1][2]
            CCmax = CC[rloc, cloc]
            rloc = rloc - dftshift - 1
            cloc = cloc - dftshift - 1
            row_shift = row_shift + rloc / Upsamp_fac
            col_shift = col_shift + cloc / Upsamp_fac
        end
        
        # Adjust shifts if it's only one row or column
        if nr == 1
            row_shift = 0
        end
        if nc == 1
            col_shift = 0
        end
    end

    # Compute error and phase difference
    rg00 = sum(abs.(FSlc1[:]) .^ 2)
    rf00 = sum(abs.(FSlc2[:]) .^ 2)
    error = 1.0 - abs(CCmax) ^ 2 / (rg00 * rf00)
    error = sqrt(abs(error))
    diffphase = angle(CCmax)
    output = [error, diffphase, row_shift, col_shift]

    # Compute registered version of Slc2
    if Upsamp_fac == 0
        Greg = FSlc2 * exp(1im * diffphase)
    elseif Upsamp_fac > 0
        Nc2 = transpose(Nc) .* ones(length(Nr))
        Nr2 = transpose(ones(length(Nc))) .* Nr
        Greg = FSlc2 .* exp.(1im * 2 * pi * (-row_shift * Nr2 / nr - col_shift * Nc2 / nc))
        Greg = Greg * exp.(1im * diffphase)
    end
    Reg_Slc2 = ifft(Greg)

    return output, Reg_Slc2  # Return computed output and registered Slc2
end

"""
Pad or crop the Fourier transformed image to a specified output size.

Arguments:
- imFT: Input Fourier transformed image
- op_size: Output size of the Fourier transformed image

Returns:
- imFTout: Padded or cropped Fourier transformed image
"""
function FTpad(imFT, op_size)

	Nout 			= op_size
	Nin 			= [size(imFT)[1] ;size(imFT)[2]]
	imFT 			= fftshift(imFT)
	center 			= floor.([size(imFT)[1] ;size(imFT)[2]] / 2).+1

	imFTout 		= ComplexF64.(zeros(op_size[1], op_size[2]))
	centerout 		= floor.([size(imFTout)[1] ;size(imFTout)[2]]/2).+1

	cenout_cen 		= Int.(centerout - center)
	imFTout[max(cenout_cen[1]+1,1):min(cenout_cen[1]+Nin[1],Nout[1]) , max(cenout_cen[2]+1,1):min(cenout_cen[2]+Nin[2],Nout[2])] .= 
		imFT[max(-cenout_cen[1]+1,1):min(-cenout_cen[1]+Nout[1],Nin[1]) , max(-cenout_cen[2]+1,1):min(-cenout_cen[2]+Nout[2],Nin[2])];

	imFTout 		= ifftshift(imFTout) * Nout[1] * Nout[2] / (Nin[1] * Nin[2])

	return imFTout
end


"""
2D Upsampling using Discrete Fourier Transform (DFT) matrix products.

Arguments:
- inp: Input data to upsample
- nor: Number of rows in the output
- noc: Number of columns in the output
- Upsamp_fac: Upsampling factor
- roff: Row offset for upsampling
- coff: Column offset for upsampling

Returns:
- out: Upsampled data using DFT matrix products
"""
function dftups(inp, nor, noc, Upsamp_fac, roff, coff)

	nr, nc 			= size(inp)
	#Compute kernels and obtain DFT by matrix products
	kernc 			= exp((-1im * 2 * pi / (nc*Upsamp_fac) ) * (transpose(ifftshift(0:nc-1)) - floor(nc/2)) * ( (0:noc-1) - coff) )
	kernr 			= exp((-1im * 2 * pi / (nr*Upsamp_fac) ) * (transpose(0:nor-1) - roff) * (ifftshift([0:nr-1]) - floor(nr/2)) )

	out 			= kernr * inp * kernc 

	return out
end


"""
Estimate the location of a target from the spectrum. assuming there's just one.

Arguments:
- x: Input signal

Returns:
- tx: Estimated x-coordinate of the target
- ty: Estimated y-coordinate of the target
"""
function measure_location_from_spectrum(x)
    X = fft(x)
    # Estimate location of target, assuming there's just one.
    tx, ty = estimate_frequency(X)
    # Scale to pixels
    tx *= -size(X, 2) / (2 * pi)
    ty *= -size(X, 1) / (2 * pi)
    # Ensure positivity by modulo operation
    tx %= size(x, 2)
    ty %= size(x, 1)
    return tx, ty
end


"""
Oversample the input signal in the frequency domain by zero-padding high frequencies.

Arguments:
- x: Input signal (2D array)
- nov: Oversampling factor
- baseband: Boolean indicating whether the signal is in baseband (default: false)
- return_slopes: Boolean to return frequency shifts (default: false)

Returns:
- y: Oversampled signal in time domain
- fx: Frequency shift in x-direction (if return_slopes=true)
- fy: Frequency shift in y-direction (if return_slopes=true)
"""
function oversample(x, nov, baseband=false, return_slopes=false)
    m, n = size(x)
    if m != n
        return error("m != n")
    end

    if baseband == false
        # Shift the data to baseband
        fx, fy = estimate_frequency(x)
        x = shift_frequency(x, -fx, -fy)
    end

    X = fft(x)

    # Zero-pad high frequencies in the spectrum.
    Y = convert(typeof(X), zeros((n * nov, n * nov)))
    n2 = Int64.(floor(n / 2))
    Y[1:n2, 1:n2] = X[1:n2, 1:n2]
    Y[end - n2 + 1:end, end - n2 + 1:end] = X[end - n2 + 1:end, end - n2 + 1:end]
    Y[1:n2, end - n2 + 1:end] = X[1:n2, end - n2 + 1:end]
    Y[end - n2 + 1:end, 1:n2] = X[end - n2 + 1:end, 1:n2]

    # Split Nyquist bins symmetrically.
    if n % 2 != 0
        return error("n % 2 != 0")
    end
    # ... (continuing the processing of Nyquist bins)

    # Back to time domain.
    y = ifft(Y)
    # Account for scaling of different-sized DFTs.
    y *= nov^2

    if baseband == false
        # Put the phase back on
        y = shift_frequency(y, fx / nov, fy / nov)
    end

    # if return_slopes == true
    if return_slopes == true
        return (y, fx, fy)
    end
    return y
end


"""
Estimate frequencies in x and y directions from the given signal.

Arguments:
- z: Input signal (2D array)

Returns:
- Estimated frequency in the x-direction (fx)
- Estimated frequency in the y-direction (fy)
"""
function estimate_frequency(z)
    cx = sum(z[:, 2:end] .* conj(z[:, 1:end-1]))
    cy = sum(z[2:end, :] .* conj(z[1:end-1, :]))
    return angle.([cx, cy])
end		

"""
Shifts the frequencies of a complex matrix 'z' along the x and y directions by the specified amounts 'fx' and 'fy'.
This function utilizes coordinate matrices 'x' and 'y' to apply frequency shifts accordingly.

Arguments:
- z: Input signal (2D complex matrix)
- fx: Frequency shift in the x-direction
- fy: Frequency shift in the y-direction

Returns:
- Signal 'z' with shifted frequencies in the x and y directions
"""
function shift_frequency(z, fx, fy)
    # Create coordinate matrices 'x' and 'y' for x and y directions
    x = Int.(transpose(1:size(z, 2)) .* ones(size(z, 1)))
    y = Int.(transpose(ones(size(z, 2))) .* (1:size(z, 1)))
    
    # Apply frequency shift along the x direction
    z *= exp(1im * fx * x)
    
    # Apply frequency shift along the y direction
    z *= exp(1im * fy * y)
    
    return z  # Return the complex matrix 'z' with shifted frequencies
end


end
