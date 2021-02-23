module Performance_Metrics

function resolution(image_3D,res_dB,target_location,scene_axis1,scene_axis2,scene_axis3)
    scene_res1=0;scene_res2=0;scene_res3=0
    target_index1=findall(target_location[1] .==scene_axis1)
    target_index2=findall(target_location[2] .==scene_axis2)
    target_index3=findall(target_location[3] .==scene_axis3)
    if isempty(target_index1) || isempty(target_index2) || isempty(target_index3)
        println("Target is not inside the scene!")
        resolutions=[NaN,NaN,NaN]
    else
        # resolution in dimension 1
        if length(scene_axis1)>1
            scene_res1=scene_axis1[2]-scene_axis1[1] # scene resolution along the 1st axis
            image_slice=image_3D[:,target_index2,target_index3]
            image_1D=zeros(Float64,length(image_slice))
            image_1D[:]=image_slice
            res_1=resolution_1D(image_1D,scene_res1,res_dB)
        else;res_1=NaN;println("Resolution along 1st dimension cannot be calculated since image is 2D along the other two dimensions.");end
        # resolution in dimension 2
        if length(scene_axis2)>1
            scene_res2=scene_axis2[2]-scene_axis2[1] # scene resolution along the 2nd axis
            image_slice=image_3D[target_index1,:,target_index3]
            image_1D=zeros(Float64,length(image_slice))
            image_1D[:]=image_slice
            res_2=resolution_1D(image_1D,scene_res2,res_dB)
        else;res_2=NaN;println("Resolution along 2nd dimension cannot be calculated since image is 2D along the other two dimensions.");end
        # resolution in dimension 3
        if length(scene_axis3)>1
            scene_res3=scene_axis3[2]-scene_axis3[1] # scene resolution along the 3rd axis
            image_slice=image_3D[target_index1,target_index2,:]
            image_1D=zeros(Float64,length(image_slice))
            image_1D[:]=image_slice
            res_3=resolution_1D(image_1D,scene_res3,res_dB)
        else;res_3=NaN;println("Resolution along 3rd dimension cannot be calculated since image is 2D along the other two dimensions.");end
        resolutions=[res_1,res_2,res_3]
    end
    return resolutions
end

function resolution_1D(image_1D,scene_res,res_dB) # image1D in linear scale (not dB)
    image_1D=20*log10.(image_1D/maximum(image_1D))
    max_ind=findall(image_1D .==maximum(image_1D))
    max_ind=max_ind[1]
    res_ind_2=findfirst(image_1D[max_ind].-image_1D[max_ind:end] .>=res_dB) #TODO warn if scene extent is smaller than resolution
    res_ind_1=findlast(image_1D[max_ind].-image_1D[1:max_ind] .>=res_dB)
    if isempty(res_ind_1) & isempty(res_ind_2);println("Scene is smaller than resolution on each side!");res=NaN
    elseif isempty(res_ind_1) | isempty(res_ind_2)
        println("Scene is smaller than resolution on one side. Calculating resolution from one side only.")
        if isempty(res_ind_1);res=2*(res_ind_2-1)*scene_res;end
        if isempty(res_ind_2);res=2*(max_ind-res_ind_1)*scene_res;end
    else; res=(max_ind+res_ind_2-1-res_ind_1)*scene_res;end
    return res
end

function sidelobes()
    PSLR=1
    ISLR=1
    return PSLR, ISLR
end

end
