module Performance_Metrics
using Plots

#TODO add function definitions, comments, define input types
function PSF_metrics(image_3D,res_dB,target_location,scene_axis1,scene_axis2,scene_axis3,PSF_peak_target)
    image_1D_1,image_1D_2,image_1D_3=obtain_1D_slices(image_3D,target_location,scene_axis1,scene_axis2,scene_axis3,PSF_peak_target)
    if length(image_1D_1)>1
        scene_res1=scene_axis1[2]-scene_axis1[1] # scene resolution along the 1st axis
        res_1,res_ind_1,res_ind_2=resolution_1D(image_1D_1,scene_res1,res_dB)
        PSLR_1,ISLR_1=sidelobe_1D(image_1D_1,1,res_ind_1,res_ind_2)
        loc_error_1=location_error(image_1D_1,target_location[1],scene_axis1)
    else;res_1=NaN;PSLR_1=NaN;ISLR_1=NaN;loc_error_1=NaN;end
    if length(image_1D_2)>1
        scene_res2=scene_axis2[2]-scene_axis2[1] # scene resolution along the 2nd axis
        res_2,res_ind_1,res_ind_2=resolution_1D(image_1D_2,scene_res2,res_dB)
        PSLR_2,ISLR_2=sidelobe_1D(image_1D_2,2,res_ind_1,res_ind_2)
        loc_error_2=location_error(image_1D_2,target_location[2],scene_axis2)
    else;res_2=NaN;PSLR_2=NaN;ISLR_2=NaN;loc_error_2=NaN;end
    if length(image_1D_3)>1
        scene_res3=scene_axis3[2]-scene_axis3[1] # scene resolution along the 3rd axis
        res_3,res_ind_1,res_ind_2=resolution_1D(image_1D_3,scene_res3,res_dB)
        PSLR_3,ISLR_3=sidelobe_1D(image_1D_3,3,res_ind_1,res_ind_2)
        loc_error_3=location_error(image_1D_3,target_location[3],scene_axis3)
    else;res_3=NaN;PSLR_3=NaN;ISLR_3=NaN;loc_error_3=NaN;end
    resolutions=[res_1,res_2,res_3]
    PSLRs=[PSLR_1,PSLR_2,PSLR_3]
    ISLRs=[ISLR_1,ISLR_2,ISLR_3]
    loc_errors=[loc_error_1,loc_error_2,loc_error_3]
    return resolutions,PSLRs,ISLRs,loc_errors
end

function obtain_1D_slices(image_3D,target_location,scene_axis1,scene_axis2,scene_axis3,PSF_peak_target)
    if PSF_peak_target==1 # 1D slices from PSF peak
        max_ind=findall(image_3D .==maximum(image_3D))
        slice_index1=max_ind[1][1]
        slice_index2=max_ind[1][2]
        slice_index3=max_ind[1][3]
    elseif PSF_peak_target==2 # PSF slices from actual target location
        slice_index1=findall(target_location[1] .==scene_axis1)
        slice_index2=findall(target_location[2] .==scene_axis2)
        slice_index3=findall(target_location[3] .==scene_axis3)
    end
    if PSF_peak_target==2 && (isempty(slice_index1) || isempty(slice_index2) || isempty(slice_index3))
        image_1D_1=NaN;image_1D_2=NaN;image_1D_3=NaN
    else
        gr() # or plotly()
        if length(scene_axis1)>1
            image_slice=image_3D[:,slice_index2,slice_index3]
            image_1D_1=zeros(Float64,length(image_slice))
            image_1D_1[:]=image_slice
            display(plot(scene_axis1,20*log10.(image_1D_1),xaxis=("scene axis 1 in scene units"),ylabel=("amplitude (dB)"),size=(1600,900),leg=false)) # plot the PSF along axis 1
        else;image_1D_1=NaN;println("PSF metrics along 1st dimension cannot be calculated since image has no 1st dimension.");end
        if length(scene_axis2)>1
            image_slice=image_3D[slice_index1,:,slice_index3]
            image_1D_2=zeros(Float64,length(image_slice))
            image_1D_2[:]=image_slice
            display(plot(scene_axis2,20*log10.(image_1D_2),xaxis=("scene axis 2 in scene units"),ylabel=("amplitude (dB)"),size=(1600,900),leg=false)) # plot the PSF along axis 2
        else;image_1D_2=NaN;println("PSF metrics along 2nd dimension cannot be calculated since image has no 2nd dimension.");end
        if length(scene_axis3)>1
            image_slice=image_3D[slice_index1,slice_index2,:]
            image_1D_3=zeros(Float64,length(image_slice))
            image_1D_3[:]=image_slice
            display(plot(scene_axis3,20*log10.(image_1D_3),xaxis=("scene axis 3 in scene units"),ylabel=("amplitude (dB)"),size=(1600,900),leg=false)) # plot the PSF along axis 3
        else;image_1D_3=NaN;println("PSF metrics along 3rd dimension cannot be calculated since image has no 3rd dimension.");end
    end
    return image_1D_1,image_1D_2,image_1D_3
end

function location_error(image_1D,target_location_1D,scene_axis)
    target_ind=findall(target_location_1D .==scene_axis);target_ind=target_ind[1]
    max_ind=findall(image_1D .==maximum(image_1D));max_ind=max_ind[1]
    scene_res=scene_axis[2]-scene_axis[1] # scene resolution along the axis
    loc_error=(max_ind-target_ind)*scene_res
end

function resolution_1D(image_1D,scene_res,res_dB) # image1D in linear scale (not dB)
    image_1D=20*log10.(image_1D/maximum(image_1D))
    max_ind=findall(image_1D .==maximum(image_1D));max_ind=max_ind[1]
    res_ind_2=findfirst(image_1D[max_ind].-image_1D[max_ind:end] .>=res_dB)+max_ind-1 #TODO warn if scene extent is smaller than resolution
    res_ind_1=findlast(image_1D[max_ind].-image_1D[1:max_ind] .>=res_dB)
    if isempty(res_ind_1) & isempty(res_ind_2);println("Scene is smaller than resolution on both sides! Can't calculate resolution.");res=NaN
    elseif isempty(res_ind_1) | isempty(res_ind_2)
        println("Scene is smaller than resolution on one side. Calculating resolution from one side only.")
        if isempty(res_ind_1);res=2*(res_ind_2-max_ind)*scene_res;end
        if isempty(res_ind_2);res=2*(max_ind-res_ind_1)*scene_res;end
    else; res=(res_ind_2-res_ind_1)*scene_res;end
    return res,res_ind_1,res_ind_2
end

function sidelobe_1D(image_1D,axis_no,res_ind_1,res_ind_2)
    peak_value,peak_ind=findmax(image_1D)
    SLpeaks_ind=findpeaks(image_1D)
    if length(SLpeaks_ind)>1
        PSLR=20*log10(peak_value/image_1D[SLpeaks_ind[2]])
    else;println("No sidelobes found inside the scene for dimension ",axis_no,"! Need to increase resolution and/or scene size.");PSLR=NaN;end
    if isempty(res_ind_1) & isempty(res_ind_2);println("Scene is smaller than resolution on both sides! Can't calculate ISLR.");ISLR=NaN
    elseif isempty(res_ind_1) | isempty(res_ind_2)
        println("Scene is smaller than resolution on one side. Calculating PSLR from one side only.")
        if isempty(res_ind_1)
            sidelobe_energy=sum(image_1D[res_ind_2+1:end])
            mainlobe_energy=sum(image_1D[SLpeaks_ind[1]:res_ind_2])
        elseif isempty(res_ind_2)
            sidelobe_energy=sum(image_1D[1:res_ind_1-1])
            mainlobe_energy=sum(image_1D[res_ind_1:SLpeaks_ind[1]])
        end
        ISLR=20*log10(sidelobe_energy/mainlobe_energy)
    else
        sidelobe_energy=sum(image_1D[1:res_ind_1-1])+sum(image_1D[res_ind_2+1:end])
        mainlobe_energy=sum(image_1D[res_ind_1:res_ind_2])
        ISLR=20*log10(sidelobe_energy/mainlobe_energy)
    end
    return PSLR,ISLR
end


"""
`findpeaks(y::Array{T},
x::Array{S}=collect(1:length(y))
;min_height::T=minimum(y), min_prom::T=minimum(y),
min_dist::S=0, threshold::T=0 ) where {T<:Real,S}`\n
Returns indices of local maxima (sorted from highest peaks to lowest)
in 1D array of real numbers. Similar to MATLAB's findpeaks().\n
*Arguments*:\n
`y` -- data\n
*Optional*:\n
`x` -- x-data\n
*Keyword*:\n
`min_height` -- minimal peak height\n
`min_prom` -- minimal peak prominence\n
`min_dist` -- minimal peak distance (keeping highest peaks)\n
`threshold` -- minimal difference (absolute value) between
 peak and neighboring points\n
"""
#=The Findpeaks.jl package is licensed under the MIT "Expat" License:
> Copyright (c) 2018: tungli.
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.=#
function findpeaks(
                   y :: AbstractVector{T},
                   x :: AbstractVector{S} = collect(1:length(y))
                   ;
                   min_height :: T = minimum(y),
                   min_prom :: T = zero(y[1]),
                   min_dist :: S = zero(x[1]),
                   threshold :: T = zero(y[1]),
                  ) where {T <: Real, S}

    dy = diff(y)

    peaks = in_threshold(dy, threshold)

    yP = y[peaks]
    peaks = with_prominence(y, peaks, min_prom)

    #minimal height refinement
    peaks = peaks[y[peaks] .> min_height]
    yP = y[peaks]

    peaks = with_distance(peaks, x, y, min_dist)

    peaks
end

"""
Select peaks that are inside threshold.
"""
function in_threshold(
                      dy :: AbstractVector{T},
                      threshold :: T,
                     ) where {T <: Real}

    peaks = 1:length(dy) |> collect

    k = 0
    for i = 2:length(dy)
        if dy[i] <= -threshold && dy[i-1] >= threshold
            k += 1
            peaks[k] = i
        end
    end
    peaks[1:k]
end

"""
Select peaks that have a given prominence
"""
function with_prominence(
                         y :: AbstractVector{T},
                         peaks :: AbstractVector{Int},
                         min_prom::T,
                        ) where {T <: Real}

    #minimal prominence refinement
    peaks[prominence(y, peaks) .> min_prom]
end

"""
Calculate peaks' prominences
"""
function prominence(y::AbstractVector{T}, peaks::AbstractVector{Int}) where {T <: Real}
    yP = y[peaks]
    proms = zero(yP)

    for (i, p) in enumerate(peaks)
        lP, rP = 1, length(y)
        for j = (i-1):-1:1
            if yP[j] > yP[i]
                lP = peaks[j]
                break
            end
        end
        ml = minimum(y[lP:p])
        for j = (i+1):length(yP)
            if yP[j] > yP[i]
                rP = peaks[j]
                break
            end
        end
        mr = minimum(y[p:rP])
        ref = max(mr,ml)
        proms[i] = yP[i] - ref
    end

    proms
end

"""
Select only peaks that are further apart than `min_dist`
"""
function with_distance(
                       peaks :: AbstractVector{Int},
                       x :: AbstractVector{S},
                       y :: AbstractVector{T},
                       min_dist::S,
                      ) where {T <: Real, S}

    peaks2del = zeros(Bool, length(peaks))
    inds = sortperm(y[peaks], rev=true)
    permute!(peaks, inds)
    for i = 1:length(peaks)
        for j = 1:(i-1)
            if abs(x[peaks[i]] - x[peaks[j]]) <= min_dist
                if !peaks2del[j]
                    peaks2del[i] = true
                end
            end
        end
    end

    peaks[.!peaks2del]
end

end