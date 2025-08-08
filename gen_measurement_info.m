function measurement_info = gen_measurement_info(init_cond,mask_info)

    bad_data_sn = [];

    effective_mask = mask_info.effective_mask;

    fprintf('Data preparation in progressing...\n');
    n_of_data = init_cond.n_of_data;
    measured_amp = cell(n_of_data,1);
    bad_data_mask = cell(n_of_data,1);
    individual_mask = cell(n_of_data,1);
    individual_mask_active_area = zeros(n_of_data,1);
    measured_amp_max = zeros(n_of_data,1);

    rawdata_clip_size = init_cond.rawdata_clip_size;
    data_clip_row_start_idx = init_cond.beam_center_y - (rawdata_clip_size-1)/2;
    data_clip_row_end_idx = init_cond.beam_center_y + (rawdata_clip_size-1)/2;
    data_clip_col_start_idx = init_cond.beam_center_x - (rawdata_clip_size-1)/2;
    data_clip_col_end_idx = init_cond.beam_center_x + (rawdata_clip_size-1)/2;
    data_clip_row_range = [data_clip_row_start_idx,data_clip_row_end_idx];
    data_clip_col_range = [data_clip_col_start_idx,data_clip_col_end_idx];
    
    master_info = EigerDataFunc.ReadEigerHDF5Master(init_cond.master_fp);
    for seq_num = 1:n_of_data
        if mod(seq_num,20) == 0 || seq_num == n_of_data
            fprintf('Data Sheet %d/%d...\n',seq_num,n_of_data); % show the progress
        end
        frame_sn = init_cond.img_count(init_cond.seq_num(seq_num)); % for the frame lose when suspension function enalbe in bluesky
        data = single(EigerDataFunc.ReadEigerHDF5Data(master_info,frame_sn,data_clip_col_range,data_clip_row_range));
        
        if init_cond.probe_extending_factor ~= 1
            data = data_rescale(data,init_cond.probe_extending_factor);
        end
        data = data.*~effective_mask;

        % check saturated data

        if(sum(data >= Inf,'all')~= 0)
            fprintf('!!!!!     Saturation point detected on #%d    !!!!!\n',seq_num);
            fprintf('Value = %d',max(data,[],'all'));
            bad_data_sn = [bad_data_sn seq_num];
            bad_data_mask{seq_num} = data >= 1E7;
            data = data.*~bad_data_mask{seq_num};
            individual_mask{seq_num} = or(effective_mask,bad_data_mask{seq_num});
        else
            individual_mask{seq_num} = effective_mask;
        end
        individual_mask_active_area(seq_num) = sum(~individual_mask{seq_num},'all');

        measured_amp{seq_num} = sqrt(data);
        measured_amp_max(seq_num,1) = max(measured_amp{seq_num},[],'all');
    end
    
    % convert cell to matrix for speedup when GPU applied
    measured_amp = cat(3,measured_amp{:});
    individual_mask = cat(3,individual_mask{:});
    
    measurement_info.bad_data_sn = bad_data_sn;
    measurement_info.bad_data_mask = bad_data_mask;
    measurement_info.individual_mask = individual_mask;
    measurement_info.individual_mask_active_area = individual_mask_active_area;
    measurement_info.measured_amp_max = measured_amp_max;
    measurement_info.measured_amp = measured_amp;
    fprintf('Done.\n');
end
