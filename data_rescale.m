function data = data_rescale(data,scaling_factor)
    
    data= single(data);
    
    if scaling_factor == 1
        return
    elseif scaling_factor == 2
        ones_2_by_2 = ones(2);
        data = kron(data,ones_2_by_2);
        data = conv2(data,ones_2_by_2);
        [data_row_size,data_col_size] = size(data);
        data = data(2:data_row_size-1,2:data_col_size-1)/sum(sum(ones_2_by_2));
    else
        [original_clip_size,~] = size(data);
        resized_size = (original_clip_size-1)*scaling_factor + 1;
        data = imresize(data,[resized_size resized_size],'nearest');
    end
    
end