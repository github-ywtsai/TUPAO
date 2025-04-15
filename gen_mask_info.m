function mask_info = gen_mask_info(init_cond)
    master_fp = init_cond.master_fp;
        
    
    % fprintf('\nLoading mask from detector...\t')
    mask_info.pixel_mask = single(transpose(h5read(master_fp,'/entry/instrument/detector/detectorSpecific/pixel_mask')));
    mask_info.pixel_mask(mask_info.pixel_mask ~= 0) = 1;
    %mask_info.pixel_mask = zeros(init_cond.y_pixels_in_detector,init_cond.x_pixels_in_detector);
    % fprintf('Done.\n')
    
    % fprintf('\nLoading manual mask...\t')
    mask_info.manual_mask = single(get_manual_mask(init_cond));
    
    mask_temp = or(mask_info.pixel_mask, mask_info.manual_mask);
    
    mask_info.mask = single(mask_temp);
    
    rawdata_clip_size = init_cond.rawdata_clip_size;
    data_clip_row_start_idx = init_cond.beam_center_y - (rawdata_clip_size-1)/2;
    data_clip_row_end_idx = init_cond.beam_center_y + (rawdata_clip_size-1)/2;
    data_clip_col_start_idx = init_cond.beam_center_x - (rawdata_clip_size-1)/2;
    data_clip_col_end_idx = init_cond.beam_center_x + (rawdata_clip_size-1)/2;
    mask_info.effective_mask = mask_info.mask(data_clip_row_start_idx:data_clip_row_end_idx,data_clip_col_start_idx:data_clip_col_end_idx);
    if init_cond.probe_extending_factor ~= 1
        mask_info.effective_mask = data_rescale(mask_info.effective_mask,init_cond.probe_extending_factor);
        mask_info.effective_mask(mask_info.effective_mask~=0) = 1;
    end
    mask_info.effective_mask = logical(mask_info.effective_mask);
    mask_info.effective_active_area = sum(sum(~mask_info.effective_mask));

    %fprintf('Done.\n')

end

function manual_mask = get_manual_mask(init_cond)

    if ~isempty(init_cond.manual_mask_fp)
        manual_mask_temp = zeros(init_cond.y_pixels_in_detector,init_cond.x_pixels_in_detector,size(init_cond.manual_mask_fp,1));
        for idx = 1:size(init_cond.manual_mask_fp,1)
            manual_mask_fp = init_cond.manual_mask_fp{idx};

            temp = importdata(manual_mask_fp,',');
            textdata = temp.textdata;
            temp = temp.data;
            % new imageJ -> col_1: x, col_2: y, old imageJ -> col_1: idx, col_2: x, col_3: y, col_4: value
            % subs for accumarray is [row,col;row,col;....], so y in col_1 and x in col_2. % imageJ x and y start from 0, so +1 for matlab
            if strcmpi(textdata{1},'X') && strcmpi(textdata{2},'Y') && strcmpi(textdata{3},'Value')
                temp = fliplr(temp(:,1:2)+1);
            else
                temp = fliplr(temp(:,2:3)+1);
            end
            manual_mask_temp(:,:,idx) = accumarray(temp,1,[init_cond.y_pixels_in_detector,init_cond.x_pixels_in_detector]);
        end
        
        manual_mask = sum(manual_mask_temp,3);
        manual_mask(manual_mask~=0) = 1;
        
        if init_cond.y_pixels_in_detector == 4371 && init_cond.x_pixels_in_detector == 4150
            temp = load('/home/tsai.yw/gpfs_folder/raw data/holder_blocked_area_on_detector.mat');
            holder_blocked_area_on_detector = temp.holder_blocked_area_on_detector;
            manual_mask = manual_mask + holder_blocked_area_on_detector;
            manual_mask(manual_mask~=0) = 1;
        end
    else
        manual_mask = zeros(init_cond.y_pixels_in_detector,init_cond.x_pixels_in_detector);
    end
    
end