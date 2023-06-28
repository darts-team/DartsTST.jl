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
c = 299792458

##
"""
Function for performing range compression of the received signal for all platforms and slowtime
Uses range_compression function from the waveform.jl module

"""
function generate_rangecompressed_received_data(rawdata, RCRF_F, params_densesim, Np, Nst, Ncr)

	@unpack waveform_domain_process, mode = params_densesim 

	if mode==1 || mode==2 # SAR (ping-pong) or SIMO
	    RC_signal 			= zeros(ComplexF64,Nst,Np,Ncr)
	elseif mode==3
	    RC_signal 			= zeros(ComplexF64,Nst,Np,Np,Ncr)
	end
	
	if mode==1 || mode==2
		for idx_p = 1:Np
			for idx_st = 1:Nst
				RC_signal[idx_st,idx_p,:] = transpose(Waveform.range_compression(fft(rawdata[idx_st,idx_p,:]) , RCRF_F, waveform_domain_process))
			end
		end
	elseif mode==3
		for idx_p1 = 1:Np
			for idx_p2 = 1:Np
				for idx_st = 1:Nst
					RC_signal[idx_st,idx_p1,idx_p2,:] = transpose(Waveform.range_compression(fft(rawdata[idx_st,idx_p1,idx_p2,:]) , RCRF_F, waveform_domain_process))
				end
			end
		end
	end

	return RC_signal
end

##
"""
Function for performing range Doppler correction of the received range compressed signal for all platforms
Uses range_doppler_algorithm function in this module

"""
function generate_rangeDoppler_received_data(RC_signal, Np, Nst, Ncr, Veff, fDopp, slrng, mid_idx, max_range_time, ref_range, params_densesim)
	@unpack fs, fp, fc = params_densesim
	RC_RD_signal 			= zeros(ComplexF64,Nst,Np,Ncr)
	for idx_p = 1:Np
		RC_RD_signal[:,idx_p,:] = range_doppler_algorithm(RC_signal[:,idx_p,:], Ncr, fc, fs, fp, transpose(Veff[idx_p,:] ), Nst, c, slrng[idx_p,mid_idx[idx_p]], max_range_time,fDopp[idx_p,:])
		#RC_RD_signal[:,idx_p,:] = range_doppler_algorithm1(RC_signal[:,idx_p,:], Ncr, fc, fs, fp, transpose(Veff[idx_p,:] ), Nst, c, ref_range,fDopp[idx_p,:])
	end
	return RC_RD_signal
end

##
"""
Function for performing azimuth compression of the received signal for all platforms and range samples

"""
function generate_azimuthcompressed_received_data(RC_RD_signal, Np, Nst, Ncr, Veff, fDopp, ref_fix_range, params_densesim)
	@unpack λ, fs, fp, fc = params_densesim

	AzC_signal 		= zeros(ComplexF64,Nst,Np,Ncr)
	for idx_p = 1:Np
		fD 					= sqrt.( 1 .- ( (λ^2 .* (fDopp[idx_p,:]).^2) ./ (4 .* (Veff[idx_p,1]).^2) ))
		azi_ref 			= fft(conj(exp.(-1im .* 4 .* pi .* ((1 ./fD) .* (ref_fix_range)) .* fc ./ c)))
		for idx_cr = 1:Ncr
			AzC_signal[:,idx_p,idx_cr] = fftshift(ifft(fft(RC_RD_signal[:,idx_p,idx_cr]) .* (azi_ref)))
		end
	end
	return AzC_signal
end

##
"""
Function for range Doppler algorithm, used or range curve correction
Note: Two implementations 1. based on slant ranges (mean in azimuth) 2. based on fixed range (mean in azimuth and range)

# Arguments
 - `signal`,  Input range compressed signal
 - `Ncr`, Number of samples along slant range
 - `Na`, Number of samples along azimuth
 - `fc`,  Carrier freqency
 - `fs`,  Sampling freqency
 - `PRF`,  Pulse repetition freqency
 - `vel`, Effective velocity
 - `c`, Speed of light
 - `ref_range`, Reference range
 - `fDoppler`, Doppler frequency
 - `slrng`, Slant ranges
 - `max_range_time`, Maximum range time
# Output
 - `srd2`, Range corrected data
"""
function range_doppler_algorithm(signal, Ncr, fc, fs, PRF, vel, Na, c, slrng, max_range_time, fDoppler)
	fr   	   = fftfreq(Ncr,fs).* ones(1,Na);
	fa   	   = fDoppler;

	slant_range = ((1:Ncr).*(c/(fs)) ) .- (max_range_time/2 * c) .+ slrng
	detaR      =(c^2 .* (slant_range) .* (ones(Ncr,1).*transpose(fa)).^2 ./ fc.^2 ./ 8 ./ vel.^2);
	# another implementation
	#fD 			= sqrt.( 1 .- ( ((c/fc)^2 .* (ones(N,1).*transpose(fa)).^2) ./ (4 .* (vel).^2) ))
	#detaR 		= ref_range .* (1 .- fD)./fD

	srd2 		= zeros(ComplexF64,Na,Ncr)
	for i=1:Na
		srd2[i,:] = ifft(fft(signal[i,:]).*exp.(1im .* 4 .* pi .* fr[:,i] .* ((detaR[:,i])/c)))
	end
	return srd2
end

function range_doppler_algorithm1(signal, N, fc, fs, PRF, vel, Na, c, ref_range,fDoppler)
	fr   	   = fftfreq(N,fs).* ones(1,Na);
	fa   	   = fDoppler;

	detaR      =(c^2 .* (ref_range) .* (ones(N,1).*transpose(fa)).^2 ./ fc.^2 ./ 8 ./ vel.^2);
	# another implementation
	#fD 			= sqrt.( 1 .- ( ((c/fc)^2 .* (ones(N,1).*transpose(fa)).^2) ./ (4 .* (vel).^2) ))
	#detaR 		= ref_range .* (1 .- fD)./fD

	srd2 		= zeros(ComplexF64,Na,N)
	for i=1:Na
		srd2[i,:] = ifft(fft(signal[i,:]).*exp.(1im .* 4 .* pi .* fr[:,i] .* ((detaR[:,i])/c)))
	end
	return srd2
end

##
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

function get_covariance_correlation_matrices(signal, azimuth_lim, srange_lim, Np, filt_len, noise_flag=0)

    Rxx                     	= 	 zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Np, Np);
    Corr                     	= 	 zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Np, Np);

	for idx_st=azimuth_lim[1]:azimuth_lim[2]
		for idx_r=srange_lim[1]:srange_lim[2]

			if noise_flag == 1
				SNR = 30
				sz 			= (2^(-0.5))*sqrt(10^(-SNR/10)).*(randn(size(signal))+ 1im .* randn(size(signal))) #complex Gaussian noise
				x 			= signal + sz
			else
				x 			= signal#without noise
			end

			if filt_len == 0
				Rxx[idx_st, idx_r, :,:] 		= (x[idx_st,:,idx_r] .* conj(transpose(x[idx_st,:,idx_r]))) ./ 1 
				Corr[idx_st, idx_r, :,:] 		= Rxx[idx_st, idx_r, :,:] ./ sqrt(Rxx[idx_st, idx_r, :,:] .* Rxx[idx_st, idx_r, :,:] )
			else

				if (idx_st-filt_len > azimuth_lim[1])
					lims1 = idx_st-filt_len
				else
					lims1 = azimuth_lim[1];
				end
				if (idx_st+filt_len < azimuth_lim[2])
					lims2 = idx_st+filt_len
				else
					lims2 = azimuth_lim[2];
				end

				if (idx_r-filt_len > srange_lim[1])
					limr1 = idx_r-filt_len
				else
					limr1 = srange_lim[1];
				end
				if (idx_r+filt_len < srange_lim[2])
					limr2 = idx_r+filt_len
				else
					limr2 = srange_lim[2];
				end

				for l1=1:Np
					for l2=1:Np
						Rxx[idx_st, idx_r,l1,l2] = cov(x[lims1:lims2,l1,limr1:limr2][:],x[lims1:lims2,l2,limr1:limr2][:]) 
						Rxx1 =  cov(x[lims1:lims2,l1,limr1:limr2][:],x[lims1:lims2,l1,limr1:limr2][:]) 
						Rxx2 =  cov(x[lims1:lims2,l2,limr1:limr2][:],x[lims1:lims2,l2,limr1:limr2][:]) 

						Corr[idx_st, idx_r, l1,l2] = Rxx[idx_st, idx_r,l1,l2] ./ sqrt(Rxx1 .* Rxx2)

					end
				end

			end

		end
	end

	return Rxx, Corr

end

function get_covariance_correlation_matrices_2(signal, azimuth_lim, srange_lim, Np, filt_len, noise_flag=0)

    Rxx                     	= 	 zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Np, Np);
	Rxx9                     	= 	 zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Np, Np);

    Corr                     	= 	 zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), Np, Np);

	kernel 						= ones(filt_len,filt_len)
	weight 						= kernel ./ sum(kernel)

	for idx_st=azimuth_lim[1]:azimuth_lim[2]
		for idx_r=srange_lim[1]:srange_lim[2]

			if noise_flag == 1
				SNR = 30
				sz 			= (2^(-0.5))*sqrt(10^(-SNR/10)).*(randn(size(signal))+ 1im .* randn(size(signal))) #complex Gaussian noise
				x 			= signal + sz
			else
				x 			= signal#without noise
			end

			if filt_len == 0
				Rxx[idx_st, idx_r, :,:] 		= (x[idx_st,:,idx_r] .* conj(transpose(x[idx_st,:,idx_r]))) ./ 1 
				Corr[idx_st, idx_r, :,:] 		= Rxx[idx_st, idx_r, :,:] ./ sqrt(Rxx[idx_st, idx_r, :,:] .* Rxx[idx_st, idx_r, :,:] )
			else

				if (idx_st-filt_len > azimuth_lim[1])
					lims1 = idx_st-filt_len
				else
					lims1 = azimuth_lim[1];
				end
				if (idx_st+filt_len < azimuth_lim[2])
					lims2 = idx_st+filt_len
				else
					lims2 = azimuth_lim[2];
				end

				if (idx_r-filt_len > srange_lim[1])
					limr1 = idx_r-filt_len
				else
					limr1 = srange_lim[1];
				end
				if (idx_r+filt_len < srange_lim[2])
					limr2 = idx_r+filt_len
				else
					limr2 = srange_lim[2];
				end

				Rxx9[idx_st, idx_r, :,:] 		= (x[idx_st,:,idx_r] .* conj(transpose(x[idx_st,:,idx_r]))) ./ 1 
				for l1=1:Np
					for l2=1:Np
						Rxx[idx_st, idx_r,l1,l2] = mean(Rxx9[lims1:lims2,limr1:limr2,l1,l2]) 
						#Rxx[idx_st, idx_r,l1,l2] = mean(x[lims1:lims2,l1,limr1:limr2]) .* conj(transpose(mean(x[lims1:lims2,l2,limr1:limr2])))  ./ 1
						#Rxx[idx_st, idx_r,l1,l2] = cov(x[lims1:lims2,l1,limr1:limr2][:],x[lims1:lims2,l2,limr1:limr2][:]) 

						Rxx1 =  mean(Rxx9[lims1:lims2,limr1:limr2,l1,l1]) 
						Rxx2 =  mean(Rxx9[lims1:lims2,limr1:limr2,l2,l2]) 

						#Corr[idx_st, idx_r, l1,l2] = Rxx[idx_st, idx_r,l1,l2] ./ sqrt(Rxx1 .* Rxx2)

					end
				end

			end

		end
	end

	for l1=1:Np
		for l2=1:Np
			Corr[:, :,l1,l2] = imfilter(Rxx9[:,:,l1,l2], weight)

		end
	end

	return Rxx, Corr

end


function get_covariance_correlation_matrices_new(signal, Np, filt_len)

	A, B, C 					= size(signal)
	Rxxp                     	= zeros(ComplexF64, A, C, Np);
	Rxx                     	= zeros(ComplexF64, A, C, Np, Np);
    Corr                     	= zeros(ComplexF64, A, C, Np, Np);

	kernel 						= ones(filt_len,filt_len)
	weight 						= kernel ./ sum(kernel)

	for idx_p=1:Np
		signal_pow 				= (abs.(signal[:,idx_p,:])).^2
		Rxxp[:,:,idx_p] 		= imfilter(signal_pow, weight)
	end

	for idx_p1=1:Np
		for idx_p2=1:Np
			sam_cov_mat 		= signal[:,idx_p1,:] .* conj(signal[:,idx_p2,:])
			#sam_cov_mat_real 	= imfilter(real(sam_cov_mat), weight)
			#sam_cov_mat_imag 	= imfilter(imag(sam_cov_mat), weight)
			#Rxx[:,:,idx_p1,idx_p2] = sam_cov_mat_real .+ (1im .* sam_cov_mat_imag)
			Rxx[:,:,idx_p1,idx_p2] = imfilter(sam_cov_mat, weight)

			Corr[:,:,idx_p1,idx_p2] = Rxx[:,:,idx_p1,idx_p2] ./ ((Rxxp[:,:,idx_p1] .* Rxxp[:,:,idx_p2]).^ 0.5)
		end
	end

	return Rxx, Corr

end


function angle_2vec(a, b)
    return acosd(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
end

function get_steering_matrix(p_xyz, s_xyz_3xN_2D, azimuth_lim, srange_lim, heights_z, Np, Master_platform, λ, mode, processing_mode)

	aT          = zeros(ComplexF64, length(azimuth_lim[1]:azimuth_lim[2]), length(srange_lim[1]:srange_lim[2]), length(heights_z), Np)

	if processing_mode==1 
		for idx_st=azimuth_lim[1]:azimuth_lim[2]
			for idx_r=srange_lim[1]:srange_lim[2]

				for idx_z = 1:length(heights_z)
					for idx_p=1:Np

						Va = mean(p_xyz[:,idx_p,:],dims=2) - s_xyz_3xN_2D[:,idx_st,idx_r]
						Vb = mean(p_xyz[:,Master_platform,:],dims=2) - s_xyz_3xN_2D[:,idx_st,idx_r]
						
						if ((Va[2]) >= (Vb[2]))
							angle_ip =  angle_2vec(Va,Vb) * 1 #-1
						else
							angle_ip =   angle_2vec(Va,Vb) * 1  
						end

						if mode == 1
							kaz = (4 * pi * ((angle_ip)*pi/180) * heights_z[idx_z]) / (λ )#* sind(look_angle))
						elseif mode == 2
							kaz = (2 * pi * ((angle_ip)*pi/180) * heights_z[idx_z]) / (λ )#* sind(look_angle))
						end
						aT[idx_st, idx_r, idx_z,idx_p] = exp.(1im * kaz) #steering vector
					end
				end
			end
		end
	elseif processing_mode==2
		for idx_st=azimuth_lim[1]:azimuth_lim[2]
			for idx_r=srange_lim[1]:srange_lim[2]
	
				for idx_z = 1:length(heights_z)
					for idx_p=1:Np
	
						Va = mean(p_xyz[:,idx_p+1,:],dims=2) - s_xyz_3xN_2D[:,idx_st,idx_r]
						Vb = mean(p_xyz[:,Master_platform,:],dims=2) - s_xyz_3xN_2D[:,idx_st,idx_r]
						
						if ((Va[2]) >= (Vb[2]))
							angle_ip =  angle_2vec(Va,Vb) * 1 #-1
						else
							angle_ip =   angle_2vec(Va,Vb) * 1  
						end
	
						if mode == 1
							kaz = (4 * pi * ((angle_ip)*pi/180) * heights_z[idx_z]) / (λ )#* sind(look_angle))
						elseif mode == 2
							kaz = (2 * pi * ((angle_ip)*pi/180) * heights_z[idx_z]) / (λ )#* sind(look_angle))
						end
						aT[idx_st, idx_r, idx_z,idx_p] = exp.(1im * kaz) #steering vector
					end
				end
			end
		end
	end
	
	return aT
end

function tomo_beamforming(Cov_mat, steering_mat, azimuth_lim, srange_lim, size_op)

	Pbf                         = zeros(size_op[1],size_op[2],size_op[3])
    for idx_st=azimuth_lim[1]:azimuth_lim[2]
		for idx_r=srange_lim[1]:srange_lim[2]
			a 							= transpose(steering_mat[idx_st, idx_r, :,:]) #transpose of steering vector
			Pbf[idx_st,idx_r,:] 		=  abs.(diag(conj(transpose(a)) * Cov_mat[idx_st,idx_r,:,:] * (a)) )# power of conventional beamformer

		end
	end
	return Pbf
end

function tomo_CAPON(Cov_mat, steering_mat, azimuth_lim, srange_lim, size_op)

	PC                         = zeros(size_op[1],size_op[2],size_op[3])
	Data_mat_size 			   = size(Cov_mat)[3]
    for idx_st=azimuth_lim[1]:azimuth_lim[2]
		for idx_r=srange_lim[1]:srange_lim[2]
			a 							= transpose(steering_mat[idx_st, idx_r, :,:]) #transpose of steering vector
			#PC[idx_st,idx_r,:] 		=  abs.(diag(1 ./ (conj(transpose(a)) * inv(Cov_mat[idx_st,idx_r,:,:]) * (a))) )# power of conventional beamformer
			diagL = std(diag(Cov_mat[idx_st,idx_r,:,:]))
			try 
				PC[idx_st,idx_r,:] 		=  abs.(diag(1 ./ (conj(transpose(a)) * inv(Cov_mat[idx_st,idx_r,:,:] .+ (diagL .* Matrix{Float64}(I, Data_mat_size, Data_mat_size))) * (a))) )# power of conventional beamformer
			catch
				continue
			end
		end
	end
	return PC
end

function tomocoordinates_to_scenecoordinates(signal, heights_t, scene_axis2, scene_axis3, look_angle, left_right_look, plat_height)


    test_mat_h 		= zeros(size(signal,2),size(signal,3))
    test_mat_c 		= zeros(size(signal,2),size(signal,3))
    op_signal       = zeros(size(signal,1),size(signal,2),size(signal,3))

	earth_radius    = 6378.137e3
	l_inc_angle  	= asind.(sind.(look_angle).*(earth_radius+plat_height)./earth_radius)
    heights_z 		= heights_t .*sind(l_inc_angle)

    for j=1:size(signal,2)
        grange          = heights_z ./tand(l_inc_angle)
        test_mat_c[j,:] = scene_axis2[j] .+ (grange)
    end
    test_mat_h = transpose(repeat(heights_z,1,size(signal,2)))

    for i=1:size(signal,1)
        for j=1:size(signal,3)
            itp = linear_interpolation(test_mat_c[:,j], signal[i,:,j],extrapolation_bc=0) 
            op_signal[i,:,j] = itp.(scene_axis2) 
        end 
    end     

    for i=1:size(signal,1)
        for j=1:size(signal,2)
            itp = linear_interpolation(test_mat_h[j,:], op_signal[i,j,:],extrapolation_bc=0) 
            op_signal[i,j,:] = itp.(heights_t)#scene_axis3) 
        end 
    end

    if left_right_look == "right"
        op_signal = op_signal[:,end:-1:1,:];
    end

	return op_signal

end


function average_2D_data(signal_ip, filt_len)

	Avg_signal_ip               = zeros( size(signal_ip));

	kernel 						= ones(filt_len,filt_len)
	weight 						= kernel ./ sum(kernel)

	for idx_p=1:size(signal_ip)[3]
		Avg_signal_ip[:,:,idx_p] 		= imfilter(signal_ip[:,:,idx_p], weight)
	end

	return Avg_signal_ip

end

##
"""
SAR beamforming algorithm for Tomography

# Arguments
 - `signal`,  Input range+azimuth compressed signal after image co-registration
 - `SNR`, Signal to noise ratio
 - `slant_range`, Slant range
 - `z_values`,  Height values in z direction
 - `params_densesim`, Input parameter structure
# Output
 - `Pbf`, Beamforming output
"""
function beamforming_algo(signal, SNR, slant_range, z_values, params_densesim, Master_platform, azimuth_lim, srange_lim)
	@unpack λ, pos_n, look_angle = params_densesim
	Nst = size(signal)[1]
	Np = size(signal)[2]
	Nr = size(signal)[3]
	Pbf = zeros(Nst,Nr,length(z_values))
	aT = zeros(ComplexF64, length(z_values), Np)
	bn = pos_n[Master_platform] .- pos_n # Perpendicular baselines TODO

	for idx_st=azimuth_lim[1]:azimuth_lim[2]
		for idx_r=srange_lim[1]:srange_lim[2]
			aT = zeros(ComplexF64, length(z_values), Np)
			sz 			= (2^(-0.5))*sqrt(10^(-SNR/10)).*(randn(Np)+ 1im .* randn(Np)) #complex Gaussian noise
			x 			= signal[idx_st,:,idx_r] + sz; #adding noise
			#x 			= signal[idx_st,:,idx_r] #without noise
			Rxx 		= (x .* conj(transpose(x))) ./ Np #Nst Nst #data covariance matrix

			for idx_z = 1:length(z_values)
				for idx_p=1:Np
					kaz = (4 * pi * ((slant_range[idx_st,Master_platform,idx_r].-slant_range[idx_st,idx_p,idx_r])) * z_values[idx_z]) / (λ  * sind(look_angle))
					#kaz = (4 * pi * (bn[idx_p]) * z_values[idx_z]) / (λ * (slant_range[idx_st,Master_platform,idx_r]) * sind(look_angle))
					#kaz = (4 * pi * abs.(bn[idx_p]) * z_values[idx_z]) / (λ * (slant_range[idx_st,Master_platform,idx_r]) )
					aT[idx_z,idx_p] = exp.(1im * kaz) #steering vector
				end
			end

			a 							= transpose(aT) #transpose of steering vector
			Pbf[idx_st,idx_r,:] 		=  abs.(diag(conj(transpose(a)) * Rxx * (a)) )# power of conventional beamformer

		end
	end
	return Pbf
end

##
"""
SAR CAPON algorithm for Tomography (In development) ***

# Arguments
 - `signal`,  Input range+azimuth compressed signal after image co-registration
 - `SNR`, Signal to noise ratio
 - `slant_range`, Slant range
 - `z_values`,  Height values in z direction
 - `params_densesim`, Input parameter structure
# Output
 - `PC`, CAPON output
"""
function CAPON_algo(signal, SNR, slant_range, z_values, params_densesim,  Master_platform, azimuth_lim, srange_lim)
	@unpack λ, pos_n, look_angle = params_densesim
	Nst = size(signal)[1]
	Np = size(signal)[2]
	Nr = size(signal)[3]
	PC = zeros(Nst,Nr,length(z_values))
	aT = zeros(ComplexF64, length(z_values), Np)
	bn = pos_n[Master_platform] .- pos_n # Perpendicular baselines TODO

	for idx_st=azimuth_lim[1]:azimuth_lim[2]
		for idx_r=srange_lim[1]:srange_lim[2]
			aT = zeros(ComplexF64, length(z_values), Np)
			sz 			= (2^(-0.5))*sqrt(10^(-SNR/10)).*(randn(Np)+ 1im .* randn(Np)) #complex Gaussian noise
			x 			= signal[idx_st,:,idx_r] + sz; #adding noise
			#x 			= signal[idx_st,:,idx_r] ; #without noise
			Rxx 		= (x .* conj(transpose(x))) ./ Np #Nst #data covariance matrix
			iRxx 		= inv(Rxx) #inverse of covariance matrix

			for idx_z = 1:length(z_values)
				for idx_p=1:Np
					kaz = (4 * pi * (bn[idx_p]) * z_values[idx_z]) / (λ * (slant_range[idx_st,Master_platform,idx_r]) * sind(look_angle))
					#kaz = (4 * pi * abs.(bn[idx_p]) * z_values[idx_z]) / (λ * slant_range[idx_st,Master_platform,idx_r] )
					aT[idx_z,idx_p] = exp.(1im * kaz) #steering vector
				end
			end

			a 			= transpose(aT) #transpose of steering vector
			#wC 		= (iRxx .* a) / (conj(transpose(a)) .* iRxx .* a); #weight vector for Caponbeamforming
			#PC[idx_st,idx_r,:] 	= diag(abs.(conj(transpose((wC)) * Rxx * wC)) )
			PC[idx_st,idx_r,:] 	= abs.(diag( 1 ./ ((conj(transpose(a)) * iRxx * (a))))); #spatial power spectrum of Capon beamformer
		end
	end
	return PC
end

##
"""
Backup code for Back Projection Algorithm Outout, this function is duplicate from process_raw_data.jl module (used for faster processing time)
"""
function main_RSF_slowtime(rawdata,s_xyz_grid,p_xyz_3D,mode,tx_el,fc,t_rx,ref_range) # with RSF and slow-time
	c =299792458
    Ns=size(s_xyz_grid)[2] # number of pixels in the scene
    Np=size(p_xyz_3D)[2] # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D)[3] # number of slow-time samples
    Δt=t_rx[2]-t_rx[1]
    processed_image=zeros(ComplexF64,Ns) # intensity image vector
    λ=c/fc # wavelength (m)
    ref_delay=2*ref_range/c # reference delay

    if mode==1 # SAR (ping-pong)
        for j=1:Ns # for each pixel
            pixel_j = @view(s_xyz_grid[:,j])
            pixel_sum = 0.0im;
            for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                for s=1:Nst # slow-time (pulses)
                    range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                    range_tx=range_rx
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=round(Int,rel_delay/Δt)
                    pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(1im*4*pi/λ*range_tx)
                end
            end
            processed_image[j] = pixel_sum
        end
    elseif mode==2 # SIMO
        for j=1:Ns # for each pixel
            pixel_j = @view(s_xyz_grid[:,j])
            pixel_sum = 0.0im;
            for s=1:Nst # slow-time (pulses)
                range_tx=distance(pixel_j,@view(p_xyz_3D[:,tx_el,s]))
                for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                    range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=round(Int,rel_delay/Δt)
                    pixel_sum=pixel_sum+rawdata[s,i,round(Int,Nft/2)+rel_delay_ind]*exp(-im*2*pi/λ*(range_tx+range_rx))
                end
            end
            processed_image[j] = pixel_sum
        end
    elseif mode==3 # MIMO
        for j=1:Ns # for each pixel
            pixel_j = @view(s_xyz_grid[:,j])
            pixel_sum = 0.0im;
            for s=1:Nst # slow-time (pulses)
                for i=1:Np # TX or RX platform for SAR, RX platform for SIMO, RX platform for MIMO
                    range_rx=distance(pixel_j,@view(p_xyz_3D[:,i,s]))
                    for k=1:Np # TX platform
                        range_tx=distance(pixel_j,@view(p_xyz_3D[:,k,s]))
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=round(Int,rel_delay/Δt)
                        pixel_sum=pixel_sum+rawdata[s,i,k,round(Int,Nft/2)+rel_delay_ind]*exp(-im*2*pi/λ*(range_tx+range_rx))
                    end
                end
            end
            processed_image[j] = pixel_sum
        end
    end
    return abs.(processed_image) # square for power?
end

function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
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

"""
Function for co-registering the azimuth compresed received data with respect to the master image/platform
Using Correlation - subpixel image registration
Algorithm based on:
Guizar-Sicairos, M., Thurman, S.T. and Fienup, J.R., 2008. Efficient subpixel image registration algorithms. Optics letters, 33(2), pp.156-158.

"""
function coregister_azimuthcompressed_received_data_corr_method(Slc1, Slc2, Upsamp_fac)

	FSlc1 			= fft(Slc1)
	FSlc2 			= fft(Slc2)

	nr, nc 			= size(FSlc2)
	Nr 				= ifftshift(-floor(nr/2):ceil(nr/2)-1)
	Nc 				= ifftshift(-floor(nc/2):ceil(nc/2)-1)

	if(Upsamp_fac == 0)
		#Simple computation of error and phase difference without registration
		CCmax 		= sum(FSlc1[:] .* conj(FSlc2[:]))
		Row_shift 	= 0
		Col_shift 	= 0

	elseif(Upsamp_fac == 1)
		#Single pixel registration
		CC 			= ifft(FSlc1 .* conj(FSlc2))
		CCabs 		= abs.(CC)
		temp_s 		= indexin(maximum(CCabs[:]), CCabs)
		row_shift 	= temp_s[1][1]
		col_shift 	= temp_s[1][2]
		CCmax 		= CC[row_shift,col_shift] * nr * nc
		row_shift 	= Nr[row_shift]
		col_shift 	= Nc[col_shift]

	elseif(Upsamp_fac > 1)
		# Start with Upsamp_fac == 2
		CC 			= ifft(FTpad(FSlc1 .* conj(FSlc2), [2*nr, 2*nc]))
		CCabs 		= abs.(CC)
		temp_s 		= indexin(maximum(CCabs[:]), CCabs)
		row_shift 	= temp_s[1][1]
		col_shift 	= temp_s[1][2]
		CCmax 		= CC[row_shift,col_shift] * nr * nc
		# Now change the shifts so that they represent relative shifts and not indices
		Nr2 		= ifftshift(-floor(nr):ceil(nr)-1)
		Nc2 		= ifftshift(-floor(nc):ceil(nc)-1)
		row_shift 	= Nr2[row_shift]/2
		col_shift 	= Nc2[col_shift]/2
		# If upsampling > 2, then refine estimate with matrix multiply DFT 
		if(Upsamp_fac > 2)
			# DFT computation
			# Initial shift estimate in upsampled grid
			row_shift 	= round(row_shift * Upsamp_fac) / Upsamp_fac
			col_shift 	= round(col_shift * Upsamp_fac) / Upsamp_fac
			dftshift 	= Int64(ceil(Upsamp_fac * 1.5)/2)
			CC 			= conj(dftups(FSlc2 .* conj(FSlc1), ceil(Upsamp_fac * 1.5), ceil(Upsamp_fac * 1.5), Upsamp_fac, 
						dftshift-row_shift * Upsamp_fac, dftshift-col_shift * Upsamp_fac))
			CCabs 		= abs.(CC)
			temp_s 		= indexin(maximum(CCabs[:]), CCabs)
			rloc 		= temp_s[1][1]
			cloc 		= temp_s[1][2]
			CCmax 		= CC[rloc,cloc] 
			rloc 		= rloc - dftshift -1
			cloc 		= cloc - dftshift -1
			row_shift 	= row_shift + rloc/Upsamp_fac
			col_shift 	= col_shift + cloc/Upsamp_fac

		end
		# If its only one row or one column the shift along that dimension has no effect. Set to zero
		if(nr == 1)
			row_shift = 0
		end

		if(nc == 1)
			col_shift = 0
		end
	end

	rg00 			= sum(abs.(FSlc1[:]) .^2 )
	rf00 			= sum(abs.(FSlc2[:]) .^2 )
	error 			= 1.0 - abs(CCmax).^2 / (rg00*rf00)
	error 			= sqrt(abs(error))
	diffphase 		= angle.(CCmax)

	output 			= [error, diffphase, row_shift, col_shift]

	#Compute registered version of Slc1
	if(Upsamp_fac == 0)
		Greg 		= FSlc2 * exp(1im * diffphase)
	elseif(Upsamp_fac > 0)
		Nc2			= transpose(Nc) .* ones(length(Nr))
		Nr2 		= transpose(ones(length(Nc))) .* (Nr) 
		Greg 		= FSlc2 .* exp.(1im * 2 * pi * (-row_shift*Nr2/nr - col_shift*Nc2/nc))
		Greg 		= Greg * exp.(1im * diffphase)
	end

	Reg_Slc2 = ifft(Greg)

	return output, Reg_Slc2

end

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

function dftups(inp, nor, noc, Upsamp_fac, roff, coff)

	nr, nc 			= size(inp)
	#Compute kernels and obtain DFT by matrix products
	kernc 			= exp((-1im * 2 * pi / (nc*Upsamp_fac) ) * (transpose(ifftshift(0:nc-1)) - floor(nc/2)) * ( (0:noc-1) - coff) )
	kernr 			= exp((-1im * 2 * pi / (nr*Upsamp_fac) ) * (transpose(0:nor-1) - roff) * (ifftshift([0:nr-1]) - floor(nr/2)) )

	out 			= kernr * inp * kernc 

	return out
end



"""
Estimate location of a target, assuming there's just one.
Parameters
"""
function measure_location_from_spectrum(x)

    X = fft(x)
    # Estimate location of target, assuming there's just one.
    tx, ty = estimate_frequency(X)
    # scale to pixels
    tx *= -size(X,2) / (2 * pi)
    ty *= -size(X,1) / (2 * pi)
    # ensure positive
    tx %= size(x,2)
    ty %= size(x,1)
    return tx, ty

end


function oversample(x, nov, baseband=false, return_slopes=false)
	m, n = size(x)
	if m != n
		return error("m != n")
	end

	if baseband==false
		# shift the data to baseband
		fx, fy = estimate_frequency(x)
		x = shift_frequency(x, -fx, -fy)
	end

	X = fft(x)

	# Zero-pad high frequencies in the spectrum. # n2 or n2+!
	Y = convert(typeof(X),zeros((n * nov, n * nov)))
	n2 = Int64.(floor(n / 2))
	#Y=X
	Y[1:n2, 1:n2] = X[1:n2, 1:n2]
	Y[end-n2+1:end, end-n2+1:end] = X[end-n2+1:end, end-n2+1:end]
	Y[1:n2, end-n2+1:end] = X[1:n2, end-n2+1:end]
	Y[end-n2+1:end, 1:n2] = X[end-n2+1:end, 1:n2]
	# Split Nyquist bins symmetrically.
	if n%2 != 0
		return error("n%2 != 0")
	end
	Y[1:n2, n2] = Y[1:n2, end-n2+1] = 0.5 * X[1:n2, n2]
	Y[end-n2+1:end, n2] = Y[end-n2+1:end, end-n2+1] = 0.5 * X[end-n2+1:end, n2]
	Y[n2, 1:n2] = Y[end-n2+1, 1:n2] = 0.5 * X[n2, 1:n2]
	Y[n2, end-n2+1:end] = Y[end-n2+1, end-n2+1:end] = 0.5 * X[n2, end-n2+1:end]
	Y[n2, n2] = Y[n2, end-n2+1] = Y[end-n2+1, n2] = Y[end-n2+1, end-n2+1] = 0.25 * X[n2, n2]
	# Back to time domain.
	y = ifft(Y)
	# NOTE account for scaling of different-sized DFTs.
	y *= nov^2

	

	if baseband==false
		# put the phase back on
		y = shift_frequency(y, fx / nov, fy / nov)
	end
	
	#y =0 .* similar(x)

	if return_slopes==true
		return (y, fx, fy)
	end
	return y

end
		

function estimate_frequency(z)
    cx = sum(z[:, 2:end] .* conj(z[:, 1:end-1]))
    cy = sum(z[2:end, :] .* conj(z[1:end-1, :]))
    return angle.([cx, cy])
end


function shift_frequency(z, fx, fy)
	x = Int.(transpose(1:size(z,2)) .* ones(size(z,1)))
	y = Int.(transpose(ones(size(z,2))) .* (1:size(z,1)))
    z *= exp(1im * fx * x)
    z *= exp(1im * fy * y)
    return z
end


end
