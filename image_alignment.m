function [OutputData1, OutputData2] = image_alignment(DataInput1,DataInput2)
    % to align objects from two retrieved results
    fixed  = angle(DataInput1.object_info.real_space);
    moving = angle(DataInput2.object_info.real_space);
    [optimizer,metric] = imregconfig('multimodal');
    tform = imregtform(moving,fixed,'translation',optimizer,metric);
    %tform = imregcorr(moving,fixed);
    
    aligned_img1 = DataInput1.object_info.real_space;
    aligned_img2 = imwarp(DataInput2.object_info.real_space,tform,'OutputView',imref2d(size(fixed)));
    
    figure
    imshowpair(angle(aligned_img1),angle(aligned_img2))
    
    OutputData1 = DataInput1;
    OutputData2 = DataInput2;
    
    OutputData1.object_info.real_space = aligned_img1;
    OutputData2.object_info.real_space = aligned_img2;
    OutputData2.object_info.real_space_xaxis = OutputData1.object_info.real_space_xaxis;
    OutputData2.object_info.real_space_yaxis = OutputData1.object_info.real_space_yaxis;