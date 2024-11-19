module Generate_Raw_Data

#using ..Scattering
include("../modules/scattering.jl")

using Distributed
using Random
using Parameters
using SharedArrays
using Interpolations

const c=299792458 # speed of light (m/s)

function distance(xyz1,xyz2)
    dist=((xyz1[1]-xyz2[1]).^2+(xyz1[2]-xyz2[2]).^2+(xyz1[3]-xyz2[3]).^2).^0.5
end

#Approximate within timesampling method
function main_gen_rawdata_method1_distributed(t_xyz_grid,p_xyz_3D,S_rsf::Vector{Float64},t_rsf,t_rx,ref_range,t_ref, params, num_proc) # with RSF and slow-time

    #set up distributed processes
    procs = addprocs(num_proc) 
    @everywhere procs begin
        Base.MainInclude.eval(include(@__FILE__))
        Base.MainInclude.eval(using .Generate_Raw_Data)
    end
        
    @unpack mode, tx_el, λ = params
    Nt              =   size(t_xyz_grid,2) # number of targets
    Np              =   size(p_xyz_3D,2) # number of platforms
    Nft             =   length(t_rx) # number of fast-time samples
    Nst             =   size(p_xyz_3D,3) # number of slow-time samples
    Δt_ft           =   t_rx[2]-t_rx[1] # fast-time resolution
    Nrsf            =   Int32(length(t_rsf))
    rawdata_const   =   2*pi/params.λ 
    ind_offset      =   Int32(floor(Nft/2)) + Int(-floor(Nrsf/2))
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata     =   SharedArray(zeros(ComplexF32, Nst,Np,Nft) )
        temp_sum    =   zeros(ComplexF32,Nft)
    elseif mode==3
        rawdata     =   SharedArray(zeros(ComplexF32,Nst,Np,Np,Nft) )
        temp_sum    =   zeros(ComplexF32,Np,Nft)
    end
    Srx_shifted     =   zeros(Nft)
    Srx_shifted2    =  zeros(Nft)

    ref_delay       =   2*ref_range/c # reference delay

    @sync @distributed for s = 1:Nst # slow-time (pulses)
        for i = 1:Np # RX platform
            temp_sum .= 0.0;
            for j = 1:Nt # targets
                if t_ref[j] != 0
                    if mode == 2;
                        range_tx = distance(t_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);
                    end
                    range_rx=distance(t_xyz_grid[:,j],p_xyz_3D[:,i,s])
                    if mode==1 # SAR (ping-pong)
                        range_tx        =   range_rx
                        rel_delay       =   (range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind   =   Int32(round(rel_delay/Δt_ft))
                        start_idx       =   ind_offset + rel_delay_ind # Start index for the RSF 
                        end_idx         =   start_idx + Nrsf - 1 # End index for the RSF          
                        temp_sum[start_idx:end_idx] .+=  (t_ref[j] .* exp(-1im*rawdata_const*(range_tx+range_rx)) .* S_rsf)
                    elseif mode==2 # SIMO
                        rel_delay       =   (range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind   =   Int32(round(rel_delay/Δt_ft))
                        start_idx       =   ind_offset + rel_delay_ind # Start index for the RSF 
                        end_idx         =   start_idx + Nrsf - 1 # End index for the RSF          
                        temp_sum[start_idx:end_idx] .+=  (t_ref[j] .* exp(-1im*rawdata_const*(range_tx+range_rx)) .* S_rsf)
                    elseif mode==3 # MIMO
                        for k=1:Np # TX platform for MIMO
                            range_tx        =   distance(t_xyz_grid[:,j],p_xyz_3D[:,k,s])
                            rel_delay       =   (range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind   =   Int32(round(rel_delay/Δt_ft))
                            start_idx       =   ind_offset + rel_delay_ind # Start index for the RSF 
                            end_idx         =   start_idx + Nrsf - 1 # End index for the RSF          
                            temp_sum[k,start_idx:end_idx] .+=  (t_ref[j] .* exp(-1im*rawdata_const*(range_tx+range_rx)) .* S_rsf)
                        end
                    end

                end
            end
            if mode==1 || mode==2 # SAR (ping-pong) or SIMO
                rawdata[s,i,:].=temp_sum
            elseif mode==3
                rawdata[s,i,:,:].=temp_sum
            end
        end
    end

    [rmprocs(p) for p in workers()]

    return rawdata
end


#Interpolation method
function main_gen_rawdata_method2_distributed(t_xyz_grid,p_xyz_3D,S_rsf::Vector{Float64},t_rsf,t_rx,ref_range,t_ref, params, num_proc) # with RSF and slow-time

    #set up distributed processes
    procs = addprocs(num_proc) 
    @everywhere procs begin
        Base.MainInclude.eval(include(@__FILE__))
        Base.MainInclude.eval(using .Generate_Raw_Data)
    end
        
    @unpack mode, tx_el, λ = params
    Nt              =   size(t_xyz_grid,2) # number of targets
    Np              =   size(p_xyz_3D,2) # number of platforms
    Nft             =   length(t_rx) # number of fast-time samples
    Nst             =   size(p_xyz_3D,3) # number of slow-time samples
    Δt_ft           =   t_rx[2]-t_rx[1] # fast-time resolution
    Nrsf            =   Int32(length(t_rsf))
    rawdata_const   =   2*pi/params.λ 
    ind_offset      =   Int32(floor(Nft/2)) + Int(-floor(Nrsf/2))
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata     =   SharedArray(zeros(ComplexF32, Nst,Np,Nft) )
        temp_sum    =   zeros(ComplexF32,Nft)
    elseif mode==3
        rawdata     =   SharedArray(zeros(ComplexF32,Nst,Np,Np,Nft) )
        temp_sum    =   zeros(ComplexF32,Np,Nft)
    end
    Srx_shifted     =   zeros(Nft)
    Srx_shifted2    =  zeros(Nft)

    ref_delay       =   2*ref_range/c # reference delay

    S_rsf_itp       = copy(S_rsf)
    interp_fn       = extrapolate(interpolate(S_rsf, BSpline(Linear())),0)

    @sync @distributed for s = 1:Nst # slow-time (pulses)
        for i = 1:Np # RX platform
            temp_sum .= 0.0;
            for j = 1:Nt # targets
                if t_ref[j] != 0
                    if mode == 2;
                        range_tx = distance(t_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);
                    end
                    range_rx=distance(t_xyz_grid[:,j],p_xyz_3D[:,i,s])
                    if mode==1 # SAR (ping-pong)
                        range_tx        =   range_rx
                        rel_delay       =   (range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind   =   Int32(round(rel_delay/Δt_ft))
                        error_delay     =   rel_delay .- (rel_delay_ind*Δt_ft)
                        if error_delay != 0
                            S_rsf_itp   =   interp_fn.((1:length(S_rsf)) .+ (error_delay / Δt_ft))
                        else
                            S_rsf_itp   =   S_rsf
                        end 
                        start_idx       =   ind_offset + rel_delay_ind # Start index for the RSF 
                        end_idx         =   start_idx + Nrsf - 1 # End index for the RSF          
                        temp_sum[start_idx:end_idx] .+=  (t_ref[j] .* exp(-1im*rawdata_const*(range_tx+range_rx)) .* S_rsf_itp)
                    elseif mode==2 # SIMO
                        rel_delay       =   (range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind   =   Int32(round(rel_delay/Δt_ft))
                        error_delay     =   rel_delay .- (rel_delay_ind*Δt_ft)
                        if error_delay != 0
                            S_rsf_itp   =   interp_fn.((1:length(S_rsf)) .+ (error_delay / Δt_ft))
                        else
                            S_rsf_itp   =   S_rsf
                        end 
                        start_idx       =   ind_offset + rel_delay_ind # Start index for the RSF 
                        end_idx         =   start_idx + Nrsf - 1 # End index for the RSF          
                        temp_sum[start_idx:end_idx] .+=  (t_ref[j] .* exp(-1im*rawdata_const*(range_tx+range_rx)) .* S_rsf_itp)
                    elseif mode==3 # MIMO
                        for k=1:Np # TX platform for MIMO
                            range_tx        =   distance(t_xyz_grid[:,j],p_xyz_3D[:,k,s])
                            rel_delay       =   (range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind   =   Int32(round(rel_delay/Δt_ft))
                            error_delay     =   rel_delay .- (rel_delay_ind*Δt_ft)
                            if error_delay != 0
                                S_rsf_itp   =   interp_fn.((1:length(S_rsf)) .+ (error_delay / Δt_ft))
                            else
                                S_rsf_itp   =   S_rsf
                            end     
                            start_idx       =   ind_offset + rel_delay_ind # Start index for the RSF 
                            end_idx         =   start_idx + Nrsf - 1 # End index for the RSF          
                            temp_sum[k,start_idx:end_idx] .+=  (t_ref[j] .* exp(-1im*rawdata_const*(range_tx+range_rx)) .* S_rsf_itp)
                        end
                    end

                end
            end
            if mode==1 || mode==2 # SAR (ping-pong) or SIMO
                rawdata[s,i,:].=temp_sum
            elseif mode==3
                rawdata[s,i,:,:].=temp_sum
            end
        end
    end

    [rmprocs(p) for p in workers()]

    return rawdata
end


#________________________________________________________________________
## Previous functions
function main_RSF(t_xyz_grid,p_xyz_grid,mode,tx_el,fc,Srx,t_rx,ref_range) # with fast-time and tomographic axes
    # TODO add descriptions of inputs and output
    λ=c/fc # wavelength (m)
    Nt=size(t_xyz_grid)[2] # number of targets
    Np=size(p_xyz_grid)[2] # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Δt=t_rx[2]-t_rx[1] # fast-time resolution
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata=zeros(ComplexF64,Np,Nft)
    elseif mode==3
        rawdata=zeros(ComplexF64,Np,Np,Nft)
    end
    ref_delay=2*ref_range/c # reference delay
    for j=1:Nt # targets
        if mode==2;range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,tx_el]);end
        for i=1:Np # RX platform
            range_rx=distance(t_xyz_grid[:,j],p_xyz_grid[:,i])
            if mode==1 # SAR (ping-pong)
                range_tx=range_rx
                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                rel_delay_ind=Int(round(rel_delay/Δt))
                if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                    Srx_shifted=[zeros(1,rel_delay_ind) Srx[1:Nft-rel_delay_ind]']
                elseif rel_delay_ind<0
                    Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft]' zeros(1,abs(rel_delay_ind))]
                end
                rawdata[i,:]=rawdata[i,:]'+exp(-im*4*pi/λ*range_tx)*Srx_shifted
            elseif mode==2 # SIMO
                rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                rel_delay_ind=Int(round(rel_delay/Δt))
                if rel_delay_ind>=0
                    Srx_shifted=[zeros(1,rel_delay_ind) Srx[1:Nft-rel_delay_ind]']
                elseif rel_delay_ind<0
                    Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft]' zeros(1,abs(rel_delay_ind))]
                end
                rawdata[i,:]=rawdata[i,:]'+exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
            elseif mode==3 # MIMO
                for k=1:Np # TX platform for MIMO
                    range_tx=distance(t_xyz_grid[:,j],p_xyz_grid[:,k])
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt))
                    if rel_delay_ind>=0
                        Srx_shifted=[zeros(1,rel_delay_ind) Srx[1:Nft-rel_delay_ind]']
                    elseif rel_delay_ind<0
                        Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft]' zeros(1,abs(rel_delay_ind))]
                    end
                    rawdata[i,k,:]=rawdata[i,k,:]'+exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
                end
            end
        end
    end
    return rawdata
end

function main_RSF_slowtime(t_xyz_grid,p_xyz_3D,Srx,t_rx,ref_range,t_ref, params) # with RSF and slow-time
    # TODO add descriptions of inputs and output
    @unpack mode, tx_el, λ = params
    Nt=size(t_xyz_grid,2) # number of targets
    Np=size(p_xyz_3D,2) # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D,3) # number of slow-time samples
    Δt_ft=t_rx[2]-t_rx[1] # fast-time resolution
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata=zeros(ComplexF64,Nst,Np,Nft)
    elseif mode==3
        rawdata=zeros(ComplexF64,Nst,Np,Np,Nft)
    end
    ref_delay=2*ref_range/c # reference delay
    for j=1:Nt # targets
        if t_ref[j]!=0
            for s=1:Nst # slow-time (pulses)
                if mode==2;range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);end
                for i=1:Np # RX platform
                    range_rx=distance(t_xyz_grid[:,j],p_xyz_3D[:,i,s])
                    if mode==1 # SAR (ping-pong)
                        range_tx=range_rx
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=Int(round(rel_delay/Δt_ft))
                        if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                            Srx_shifted=cat(zeros(rel_delay_ind),Srx[1:Nft-rel_delay_ind],dims=1)
                        elseif rel_delay_ind<0
                            Srx_shifted=cat(Srx[1+abs(rel_delay_ind):Nft],zeros(abs(rel_delay_ind)),dims=1)
                        end
                        rawdata[s,i,:]=rawdata[s,i,:]+t_ref[j]*exp(-im*4*pi/λ*range_tx)*Srx_shifted
                    elseif mode==2 # SIMO
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=Int(round(rel_delay/Δt_ft))
                        if rel_delay_ind>=0
                            Srx_shifted=[zeros(rel_delay_ind);Srx[1:Nft-rel_delay_ind]]
                        elseif rel_delay_ind<0
                            Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft];zeros(abs(rel_delay_ind))]
                        end
                        rawdata[s,i,:]=rawdata[s,i,:]+t_ref[j]*exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
                    elseif mode==3 # MIMO
                        for k=1:Np # TX platform for MIMO
                            range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,k,s])
                            rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind=Int(round(rel_delay/Δt_ft))
                            if rel_delay_ind>=0
                                Srx_shifted=[zeros(rel_delay_ind);Srx[1:Nft-rel_delay_ind]]
                            elseif rel_delay_ind<0
                                Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft];zeros(abs(rel_delay_ind))]
                            end
                            rawdata[s,i,k,:]=rawdata[s,i,k,:]+t_ref[j]*exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
                        end
                    end
                end
            end
        end
    end
    return rawdata
end

function main_RSF_slowtime_perf_opt(t_xyz_grid,p_xyz_3D,Srx::Vector{Float64},t_rx,ref_range,t_ref, params) # with RSF and slow-time
    # TODO add descriptions of inputs and output
    @unpack mode, tx_el, λ = params
    Nt=size(t_xyz_grid,2) # number of targets
    Np=size(p_xyz_3D,2) # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D,3) # number of slow-time samples
    Δt_ft=t_rx[2]-t_rx[1] # fast-time resolution
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata=zeros(ComplexF64,Nst,Np,Nft)
        temp_sum=zeros(ComplexF64,Nft)
    elseif mode==3
        rawdata=zeros(ComplexF64,Nst,Np,Np,Nft)
        temp_sum=zeros(ComplexF64,Np,Nft)
    end
    Srx_shifted = zeros(Nft)
    Srx_shifted2 = zeros(Nft)

    ref_delay=2*ref_range/c # reference delay

    for s=1:Nst # slow-time (pulses)
        for i=1:Np # RX platform
            temp_sum.=0.0;
            for j=1:Nt # targets
                if t_ref[j]!=0
                    if mode==2;range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);end
                    range_rx=distance(t_xyz_grid[:,j],p_xyz_3D[:,i,s])
                    if mode==1 # SAR (ping-pong)
                        range_tx=range_rx
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind::Int64=Int(round(rel_delay/Δt_ft))
                        if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                            circshift!(Srx_shifted, Srx, rel_delay_ind)
                            Srx_shifted[1:rel_delay_ind] .= zeros(rel_delay_ind);
                        elseif rel_delay_ind<0
                            circshift!(Srx_shifted, Srx, rel_delay_ind)
                            Srx_shifted[Nft-abs(rel_delay_ind)+1:end] .= zeros(abs(rel_delay_ind))  
                        end
                        temp_sum .= temp_sum .+  (t_ref[j].*exp(-im*4*pi/λ*range_tx).*Srx_shifted)

                    elseif mode==2 # SIMO
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=Int(round(rel_delay/Δt_ft))
                        if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                            circshift!(Srx_shifted, Srx, rel_delay_ind)
                            Srx_shifted[1:rel_delay_ind] .= zeros(rel_delay_ind);
                        elseif rel_delay_ind<0
                            circshift!(Srx_shifted, Srx, rel_delay_ind)
                            Srx_shifted[Nft-abs(rel_delay_ind)+1:end] .= zeros(abs(rel_delay_ind))  
                        end
                        temp_sum .= temp_sum .+  (t_ref[j].*exp(-im*2*pi/λ*(range_tx+range_rx)).*Srx_shifted)

                    elseif mode==3 # MIMO
                        for k=1:Np # TX platform for MIMO
                            range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,k,s])
                            rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                            rel_delay_ind=Int(round(rel_delay/Δt_ft))
                            if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                                circshift!(Srx_shifted, Srx, rel_delay_ind)
                                Srx_shifted[1:rel_delay_ind] .= zeros(rel_delay_ind);
                            elseif rel_delay_ind<0
                                circshift!(Srx_shifted, Srx, rel_delay_ind)
                                Srx_shifted[Nft-abs(rel_delay_ind)+1:end] .= zeros(abs(rel_delay_ind))  
                            end
                            temp_sum[k,:] .= temp_sum[k,:] .+  (t_ref[j].*exp(-im*2*pi/λ*(range_tx+range_rx)).*Srx_shifted)
                        end
                    end

                end
            end
            if mode==1 || mode==2 # SAR (ping-pong) or SIMO
                rawdata[s,i,:].=temp_sum
            elseif mode==3
                rawdata[s,i,:,:].=temp_sum
            end
        end
    end
    return rawdata
end

function main_RSF_slowtime_surfaceBRCS(t_xyz_grid,p_xyz_3D,Srx,t_rx,ref_range,t_ref, params) # with RSF, slow-time, and surface BRCS calculation
    # TODO add descriptions of inputs and output
    @unpack mode, tx_el, λ = params
    Nt=size(t_xyz_grid,2) # number of targets
    Np=size(p_xyz_3D,2) # number of platforms
    Nft=length(t_rx) # number of fast-time samples
    Nst=size(p_xyz_3D,3) # number of slow-time samples
    Δt_ft=t_rx[2]-t_rx[1] # fast-time resolution
    if mode==1 || mode==2 # SAR (ping-pong) or SIMO
        rawdata=zeros(ComplexF64,Nst,Np,Nft)
    elseif mode==3
        rawdata=zeros(ComplexF64,Nst,Np,Np,Nft)
    end
    ref_delay=2*ref_range/c # reference delay
    for j=1:Nt # targets
        for s=1:Nst # slow-time (pulses)
            if mode==2;range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,tx_el,s]);end
            for i=1:Np # RX platform
                range_rx=distance(t_xyz_grid[:,j],p_xyz_3D[:,i,s])
                if mode==1 # SAR (ping-pong) 
                    t_ref = Scattering.get_surface_brcs(params,p_xyz_3D[:,i,s],p_xyz_3D[:,i,s],t_xyz_grid[:,j])
                    range_tx=range_rx
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt_ft))
                    if rel_delay_ind>=0 #TODO if rel_delay_ind>=Nft Srx_shifted becomes a larger array which causes issues (also for SIMO and MIMO)
                        Srx_shifted=cat(zeros(rel_delay_ind),Srx[1:Nft-rel_delay_ind],dims=1)
                    elseif rel_delay_ind<0
                        Srx_shifted=cat(Srx[1+abs(rel_delay_ind):Nft],zeros(abs(rel_delay_ind)),dims=1)
                    end
                    rawdata[s,i,:]=rawdata[s,i,:]+t_ref*exp(-im*4*pi/λ*range_tx)*Srx_shifted
                elseif mode==2 # SIMO
                    t_ref = Scattering.get_surface_brcs(params,p_xyz_3D[:,tx_el,s],p_xyz_3D[:,i,s],t_xyz_grid[:,j])
                    rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                    rel_delay_ind=Int(round(rel_delay/Δt_ft))
                    if rel_delay_ind>=0
                        Srx_shifted=[zeros(rel_delay_ind);Srx[1:Nft-rel_delay_ind]]
                    elseif rel_delay_ind<0
                        Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft];zeros(abs(rel_delay_ind))]
                    end
                    rawdata[s,i,:]=rawdata[s,i,:]+t_ref*exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
                elseif mode==3 # MIMO
                    for k=1:Np # TX platform for MIMO
                        t_ref = Scattering.get_surface_brcs(params,p_xyz_3D[:,k,s],p_xyz_3D[:,i,s],t_xyz_grid[:,j])
                        range_tx=distance(t_xyz_grid[:,j],p_xyz_3D[:,k,s])
                        rel_delay=(range_tx+range_rx)/c-ref_delay # relative delay wrt reference delay (positive means right-shift of RSF)
                        rel_delay_ind=Int(round(rel_delay/Δt_ft))
                        if rel_delay_ind>=0
                            Srx_shifted=[zeros(rel_delay_ind);Srx[1:Nft-rel_delay_ind]]
                        elseif rel_delay_ind<0
                            Srx_shifted=[Srx[1+abs(rel_delay_ind):Nft];zeros(abs(rel_delay_ind))]
                        end
                        rawdata[s,i,k,:]=rawdata[s,i,k,:]+t_ref*exp(-im*2*pi/λ*(range_tx+range_rx))*Srx_shifted
                    end
                end
            end
        end
    end
    return rawdata
end

end
