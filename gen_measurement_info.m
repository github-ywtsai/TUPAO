function measurement_info = gen_measurement_info(init_cond,mask_info)

    %fprintf('Analyzing data structure...\t') %20250807, replaced by EigerDatafunc, remove after 20260101

    bad_data_sn = [];

    master_fp = init_cond.master_fp;
    data_file_path = fileparts(master_fp);
    effective_mask = mask_info.effective_mask;

    %{
    %20250807, replaced by EigerDatafunc, remove after 20260101
    info = h5info(master_fp);
    [n_of_link,~] = size(info.Groups.Groups(1).Links);

    % checking n_of_link, issues occur using data from TPS 23A.
    for link_idx = 1:n_of_link
        fn = fullfile(data_file_path, info.Groups.Groups(1).Links(link_idx).Value{1});
        if ~exist(fn,'file')
            n_of_link = link_idx -1;
            break
        end
    end

    fprintf('Done.\n');
    %20250807, replaced by EigerDatafunc, remove after 20260101
    %}
    
    %link_file_list = cell(n_of_link,1); %20250807, replaced by EigerDatafunc, remove after 20260101
    %link_object_list = cell(n_of_link,1); %20250807, replaced by EigerDatafunc, remove after 20260101
    % link_index = cell(n_of_link,1); %20250807, replaced by EigerDatafunc, remove after 20260101
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
    
    %{
    %20250807, replaced by EigerDatafunc, remove after 20260101
    for i = 1:n_of_link
        link_file_list{i} = fullfile(data_file_path, info.Groups.Groups(1).Links(i).Value{1});
        link_object_list{i} = info.Groups.Groups(1).Links(i).Value{2};
        image_nr_low = h5readatt(link_file_list{i},link_object_list{i},'image_nr_low');
        image_nr_high = h5readatt(link_file_list{i},link_object_list{i},'image_nr_high');
        link_index{i} = double([image_nr_low,image_nr_high]); % start/ end index of the linked data file.
    end
    
    %% plot data in sequence
    %GUIh.PlotArea.tabgp.SelectedTab =  GUIh.PlotArea.Data_tab;
    %GUIh.PlotArea.DataSn_editfield.ValueDisplayFormat = sprintf('%%d/%d',n_of_data);
    for i = 1:n_of_data
        if mod(i,20) == 0 || i == n_of_data
            fprintf('Data Sheet %d/%d...\n',i,n_of_data);
        end
        for j = 1 : n_of_link
            if i >= link_index{j}(1) && i<= link_index{j}(2)
                target_file_sn = j;
                target_fn = [link_file_list{j}];
                break
            end
        end
        data = single( transpose( h5read( target_fn, link_object_list{target_file_sn},[data_clip_col_start_idx,data_clip_row_start_idx, i - link_index{target_file_sn}(1) + 1],[rawdata_clip_size,rawdata_clip_size,1]) ) );
        %20250807, replaced by EigerDatafunc, remove after 20260101
        %}
    master_info = EigerDataFunc.ReadEigerHDF5Master(init_cond.master_fp);
    for seq_num = 1:n_of_data
        if mod(seq_num,20) == 0 || seq_num == n_of_data
            fprintf('Data Sheet %d/%d...\n',seq_num,n_of_data); % show the progress
        end
        frame_sn = init_cond.img_count(init_cond.seq_num(seq_num)); % for the frame lose when suspension function enalbe in bluesky
        data = single(EigerDataFunc.ReadEigerHDF5Data(master_info,frame_sn,[data_clip_col_start_idx,data_clip_col_end_idx],[data_clip_row_start_idx,data_clip_row_end_idx]));
        
        effective_clip_size = init_cond.effective_clip_size;
        if init_cond.probe_extending_factor ~= 1
            data = data_rescale(data,init_cond.probe_extending_factor);
        end
        data = data.*~effective_mask;

        % check saturated data
        %GUIh.Alarm.mfig.UserData.BadDataStop = 0;
        if(sum(data >= Inf,'all')~= 0)
            %GUIh.Alarm.mfig.UserData.BadDataStop = 1;
            %mfig_xstart = GUIh.MF.mfig.Position(1);
            %mfig_ystart = GUIh.MF.mfig.Position(2);
            %mfig_xsize = GUIh.MF.mfig.Position(3);
            %mfig_ysize = GUIh.MF.mfig.Position(4);
            %GUIh.Alarm.mfig.Position = [mfig_xstart + mfig_xsize/2, mfig_ystart + mfig_ysize/2, GUIh.Alarm.mfig.Position(3), GUIh.Alarm.mfig.Position(4)];
            % GUIh.Alarm.mfig.Visible = 'On';
            fprintf('!!!!!     Saturation point detected on #%d    !!!!!\n',seq_num);
            fprintf('Value = %d',max(data,[],'all'));
            bad_data_sn = [bad_data_sn seq_num];
            % data(data >= GUIh.IterCtrl.PrepDataDeadPThres_editfield.Value) = 0;
            bad_data_mask{seq_num} = data >= 1E7;
            data = data.*~bad_data_mask{seq_num};
            individual_mask{seq_num} = or(effective_mask,bad_data_mask{seq_num});
        else
            individual_mask{seq_num} = effective_mask;
        end
        individual_mask_active_area(seq_num) = sum(~individual_mask{seq_num},'all');

        measured_amp{seq_num} = sqrt(data);
        measured_amp_max(seq_num,1) = max(measured_amp{seq_num},[],'all');
        %GUIh.PlotArea.DataSn_editfield.Value = i;
        %imagesc(GUIh.PlotArea.Data_axes,measured_amp{i},[0,measured_amp_max(i,1)]);
        %title(GUIh.PlotArea.Data_axes,sprintf('Raw data max. = %.1f',max(data,[],'all')))
        %axis(GUIh.PlotArea.Data_axes,'image')
        %drawnow
    end
    
    % convert cell to matrix for speedup when GPU applied
    measured_amp_temp = zeros(effective_clip_size,effective_clip_size,init_cond.n_of_data);
    individual_mask_temp = false(effective_clip_size,effective_clip_size,init_cond.n_of_data);
    for SN = 1:numel(measured_amp)
        measured_amp_temp(:,:,SN) = measured_amp{SN};
        individual_mask_temp(:,:,SN) = individual_mask{SN};
    end
    measured_amp = measured_amp_temp;
    individual_mask = individual_mask_temp;
    
    %% save to template file measured_amp.mat
%     measured_amp_ff = init_cond.results_path;
%     SNTemp = str2double(datestr(now,'yymmddHHMM'));
%     % check file exsit or not
%     while true
%         fn = sprintf('SN%d_MeaAmp.mat',SNTemp);
%         fp = fullfile(measured_amp_ff,fn);
%         if exist(fp,'file')
%             SNTemp = SNTemp + 1;
%         else
%             break;
%         end
%     end
    
%    measurement_info.SN = SNTemp;
%    measured_amp_fn = sprintf('SN%d_MeaAmp.mat',measurement_info.SN);
%    measured_amp_fp = fullfile(measured_amp_ff, measured_amp_fn);
%    fprintf('Saving results to %s...\t',measured_amp_fp);
%    measurement_info.measured_amp_ff = measured_amp_ff;
%    measurement_info.measured_amp_fn = measured_amp_fn;
%    measurement_info.measured_amp_fp = measured_amp_fp;
    measurement_info.bad_data_sn = bad_data_sn;
    measurement_info.bad_data_mask = bad_data_mask;
    measurement_info.individual_mask = individual_mask;
    measurement_info.individual_mask_active_area = individual_mask_active_area;
    measurement_info.measured_amp_max = measured_amp_max;
    measurement_info.measured_amp = measured_amp;
%    save(measured_amp_fp,'measured_amp','-v7.3')
%    fprintf('Done.\n');
end
