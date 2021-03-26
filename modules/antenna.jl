module Antenna

function pattern_sinc(sidelobe_dB,theta, hpbw)
   #     % Nino Majurec, May 23, 2019
   # % Function simulates radiation pattern for a conical horn antenna.
   # % Function does not include gain, only the angular variation of the gain.
   # % Standard spherical coordinate system is assumed:
   # % Theta: angle between ray and z-axis (boresight along z-axis) 0 - pi
   # % Phi:   angle between ray projected on xy-plane and x-axis, 0 - 2pi
   # % Note:  Conical horn is symmetric along z-axis, so no Phi dependence here.
   # %
   # % Parameters:
   # % Sidelobe_dB - desired sidelobe suppresion in dB with respect to main
   # %               lobe. Use positive number here, e.g. 25 dB.
   # % Theta       - angle in radians with respect to z-axis. Pattern has unity
   # %               gain along z-axis (Theta = 0).
   # % Theta works as a single number or vector. Sidelobe_dB should be *only*
   # % single number, not a vector.
   # % Output:     - transmission coefficient for the E-field in linear units.
   # %               If plotted in log scale, it should go with 20*log10(),
   # %               to look like standard textbook patterns.
   # %
   # % List below contains precalculated good and realistic antenna patterns.
   # % I would suggest using B = 7.4466, C = 0.1337, which give pattern
   # % closest to the one we assumed in previous report.

   # B=   1.4304; C = 0.699106; % 111.50 deg HPBW
   # B=   2.4591; C = 0.406654; %  64.86 deg HPBW
   # B=   3.4710; C = 0.288102; %  45.95 deg HPBW
   # B=   4.4775; C = 0.223340; %  35.62 deg HPBW
   # B=   5.4816; C = 0.182429; %  29.10 deg HPBW
   # B=   6.4845; C = 0.154215; %  24.60 deg HPBW
   # B=   7.4866; C = 0.133573; %  21.30 deg HPBW
   # B=   8.4882; C = 0.117812; %  18.79 deg HPBW
   # B=   9.4894; C = 0.105382; %  16.81 deg HPBW
   # B=  10.4904; C = 0.095326; %  15.20 deg HPBW
   # B=  11.4913; C = 0.087023; %  13.88 deg HPBW
   # B=  12.4920; C = 0.080052; %  12.77 deg HPBW
   # B=  13.4926; C = 0.074116; %  11.82 deg HPBW
   # B=  14.4931; C = 0.068999; %  11.00 deg HPBW
   # B=  15.4936; C = 0.064544; %  10.29 deg HPBW
   # B=  16.4940; C = 0.060629; %   9.67 deg HPBW
   # B=  17.4943; C = 0.057162; %   9.12 deg HPBW
   # B=  18.4947; C = 0.054071; %   8.62 deg HPBW
   # B=  19.4950; C = 0.051296; %   8.18 deg HPBW
   # B=  20.4952; C = 0.048793; %   7.78 deg HPBW
   # B=  21.4954; C = 0.046523; %   7.42 deg HPBW
   # B=  22.4957; C = 0.044454; %   7.09 deg HPBW
   # B=  23.4959; C = 0.042562; %   6.79 deg HPBW
   # B=  24.4960; C = 0.040824; %   6.51 deg HPBW
   # B=  25.4962; C = 0.039223; %   6.26 deg HPBW
   # B=  26.4964; C = 0.037742; %   6.02 deg HPBW
   # B=  27.4965; C = 0.036369; %   5.80 deg HPBW
   # B=  28.4966; C = 0.035093; %   5.60 deg HPBW
   # B=  29.4968; C = 0.033903; %   5.41 deg HPBW
   # B=  30.4969; C = 0.032791; %   5.23 deg HPBW
   # B=  31.4970; C = 0.031750; %   5.06 deg HPBW
   # B=  32.4971; C = 0.030773; %   4.91 deg HPBW
   # B=  33.4972; C = 0.029854; %   4.76 deg HPBW
   # B=  34.4973; C = 0.028989; %   4.62 deg HPBW
   # B=  35.4974; C = 0.028172; %   4.49 deg HPBW
   # B=  36.4975; C = 0.027400; %   4.37 deg HPBW
   # B=  37.4975; C = 0.026669; %   4.25 deg HPBW
   # B=  38.4976; C = 0.025977; %   4.14 deg HPBW
   # B=  39.4977; C = 0.025319; %   4.04 deg HPBW
   # B=  40.4978; C = 0.024694; %   3.94 deg HPBW
   # B=  41.4978; C = 0.024099; %   3.84 deg HPBW
   # B=  42.4979; C = 0.023532; %   3.75 deg HPBW
   # B=  43.4979; C = 0.022991; %   3.67 deg HPBW
   # B=  44.4980; C = 0.022474; %   3.58 deg HPBW
   # B=  45.4981; C = 0.021980; %   3.51 deg HPBW
   # B=  46.4981; C = 0.021507; %   3.43 deg HPBW
   # B=  47.4982; C = 0.021054; %   3.36 deg HPBW
   # B=  48.4982; C = 0.020620; %   3.29 deg HPBW
   # B=  49.4983; C = 0.020204; %   3.22 deg HPBW
   # B=  50.4983; C = 0.019804; %   3.16 deg HPBW
   # B=  51.4983; C = 0.019419; %   3.10 deg HPBW
   # B=  52.4984; C = 0.019049; %   3.04 deg HPBW
   # B =  53.4984; C = 0.018693; %   2.98 deg HPBW
   # B=  55.4985; C = 0.018020; %   2.87 deg HPBW
   # B=  60.4987; C = 0.016530; %   2.64 deg HPBW
   # B=  65.4988; C = 0.015268; %   2.44 deg HPBW
   # B=  70.4990; C = 0.014186; %   2.26 deg HPBW
   # B=  75.4991; C = 0.013246; %   2.11 deg HPBW
   # B=  79.4992; C = 0.012580; %   2.01 deg HPBW
   # B=  80.4992; C = 0.012423; %   1.98 deg HPBW
   # B=  85.4993; C = 0.011697; %   1.87 deg HPBW
   # B=  90.4994; C = 0.011051; %   1.76 deg HPBW
   # B=  95.4995; C = 0.010472; %   1.67 deg HPBW
   # B= 100.4995; C = 0.009951; %   1.59 deg HPBW
   # B= 105.4996; C = 0.009480; %   1.51 deg HPBW
   # B= 110.4997; C = 0.009051; %   1.44 deg HPBW
   # B= 115.4998; C = 0.008659; %   1.38 deg HPBW
   # B= 120.4998; C = 0.008300; %   1.32 deg HPBW
   # B= 130.4999; C = 0.007664; %   1.22 deg HPBW
   # B= 140.5000; C = 0.007118; %   1.14 deg HPBW
   # B= 150.5001; C = 0.006646; %   1.06 deg HPBW
   # B= 158.5002; C = 0.006310; %   1.01 deg HPBW
   # B= 159.5002; C = 0.006271; %   1.00 deg HPBW
   # B= 160.5002; C = 0.006232; %   0.99 deg HPBW


   bArr = [1.4304,2.4591,3.471,4.4775,5.4816,6.4845,7.4866,8.4882,9.4894,10.4904,11.4913,12.492,13.4926,14.4931,15.4936,16.494,17.4943,18.4947,19.495,20.4952,21.4954,22.4957,23.4959,24.496,25.4962,26.4964,27.4965,28.4966,29.4968,30.4969,31.497,32.4971,33.4972,34.4973,35.4974,36.4975,37.4975,38.4976,39.4977,40.4978,41.4978,42.4979,43.4979,44.498,45.4981,46.4981,47.4982,48.4982,49.4983,50.4983,51.4983,52.4984,53.4984,55.4985,60.4987,65.4988,70.499,75.4991,79.4992,80.4992,85.4993,90.4994,95.4995,100.4995,10.4996,110.4997,115.4998,120.4998,130.4999,140.5,150.5001,158.5002,159.5002,160.5002]

   cArr = [0.699106,0.406654,0.288102,0.22334,0.182429,0.154215,0.133573,0.117812,0.105382,0.095326,0.087023,0.080052,0.074116,0.068999,0.064544,0.060629,0.057162,0.054071,0.051296,0.048793,0.046523,0.044454,0.042562,0.040824,0.039223,0.037742,0.036369,0.035093,0.033903,0.032791,0.03175,0.030773,0.029854,0.028989,0.028172,0.0274,0.026669,0.025977,0.025319,0.024694,0.024099,0.023532,0.022991,0.022474,0.02198,0.021507,0.021054,0.02062,0.020204,0.019804,0.019419,0.019049,0.018693,0.01802,0.01653,0.015268,0.014186,0.013246,0.01258,0.012423,0.011697,0.011051,0.010472,0.009951,0.00948,0.009051,0.008659,0.0083,0.007664,0.007118,0.006646,0.00631,0.006271,0.006232]

   hpbwArr = [111.5,64.86,45.95,35.62,29.1,24.6,21.3,18.79,16.81,15.2,13.88,12.77,11.82,11,10.29,9.67,9.12,8.62,8.18,7.78,7.42,7.09,6.79,6.51,6.26,6.02,5.8,5.6,5.41,5.23,5.06,4.91,4.76,4.62,4.49,4.37,4.25,4.14,4.04,3.94,3.84,3.75,3.67,3.58,3.51,3.43,3.36,3.29,3.22,3.16,3.1,3.04,2.98,2.87,2.64,2.44,2.26,2.11,2.01,1.98,1.87,1.76,1.67,1.59,1.51,1.44,1.38,1.32,1.22,1.14,1.06,1.01,1,0.99]

   dict = Dict(hpbwArr[i] => (bArr[i],cArr[i])  for i in 1:size(hpbwArr,1))
   if hpbw in hpbwArr
       b = dict[hpbw][1]
       c = dict[hpbw][2]
   else
       error(" HPBW doesn't match")
   end

   sls  = 10^((-sidelobe_dB + 13.26)/20);
   a    = sinc.(b/pi.*theta);
   aa_sl            = findall(abs.(theta/pi) .>= c)
   wind_flat        = ones(size(theta));
   wind_flat[aa_sl] .= sls
   wind_coef        = wind_flat;
   res = a.* wind_coef;
   return res

end

function xyz_to_sphr(xyz)
   sphr =[0.,0.,0.]
   sphr[1] = norm(xyz) #r
   sphr[2] = acos(xyz[3]/norm(xyz)) # theta
   sphr[3] = atan(xyz[2],xyz[1]) #phi
   return sphr
end

"This function converts cartesian XYZ vector to Az,El
using the TICRA Az over El convention

Arguments
   - vec::`Array{Float64,}` cartesian vector in Antenna Frame
 "
function xyz_to_azel(vec::Array{Float32,})
   @assert size(vec,1) == 3 "XYZ vector must by 3xN"
   # normalize vector
   for ii=1:size(vec,2)
      vec[:,ii] = vec[:,ii]/norm(vec[:,ii])
   end
   az = -asin.(vec[1,:])*180/π;
   el = atan.(vec[2,:], vec[3,:])*180/π;
   return az, el
end

function compute_gain(v, hpbw, slb_dB )
   v_sphr = xyz_to_sphr.(v)
   th = [v_sphr[i][2] for i in 1:size(v_sphr,1)]
   res = pattern_sinc(slb_dB,th, hpbw)
   return res
end

function interpolate_gain(v, gain_in, th_in )
   v_sphr = xyz_to_sphr.(v)
   th = [v_sphr[i][2] for i in 1:size(v_sphr,1)]

   itp = LinearInterpolation(th_in,gain_in, extrapolation_bc = Line());
   gain =  itp(th)
   return gain
end

end
