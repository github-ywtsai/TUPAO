function measurement_info = gen_measurement_info(init_cond,mask_info)

    fprintf('Analyzing data structure...\t')

    bad_data_sn = [];

    master_fp = init_cond.master_fp;
    data_file_path = fileparts(master_fp);
    effective_mask = mask_info.effective_mask;

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
    
    link_file_list = cell(n_of_link,1);
    link_object_list = cell(n_of_link,1);
    link_index = cell(n_of_link,1);
    fprintf('Data preparation in progressing...\t');
    n_of_data = init_cond.n_of_data;
    measured_amp = cell(n_of_data,1);
    bad_data_mask = cell(n_of_data,1);
    individual_mask = cell(n_of_data,1);
    individual_mask_active_area = zeros(n_of_data,1);
    measured_amp_max = zeros(n_of_data,1);

    rawdata_clip_size = init_cond.rawdata_clip_size;
    data_clip_row_start_idx = init_cond.beam_center_y - (rawdata_clip_size-1)/2;
    %data_clip_row_end_idx = init_cond.beam_center_y + (rawdata_clip_size-1)/2;
    data_clip_col_start_idx = init_cond.beam_center_x - (rawdata_clip_size-1)/2;
    %data_clip_col_end_idx = init_cond.beam_center_x + (rawdata_clip_size-1)/2;
    
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
        fprintf('Data Sheet %d/%d...\t',i,n_of_data);
        for j = 1 : n_of_link
            if i >= link_index{j}(1) && i<= link_index{j}(2)
                target_file_sn = j;
                target_fn = [link_file_list{j}];
                break
            end
        end
        data = single( transpose( h5read( target_fn, link_object_list{target_file_sn},[data_clip_col_start_idx,data_clip_row_start_idx, i - link_index{target_file_sn}(1) + 1],[rawdata_clip_size,rawdata_clip_size,1]) ) );
        
        
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
            fprintf('!!!!!     Saturation point detected on #%d    !!!!!\n',i);
            fprintf('Value = %d',max(data,[],'all'));
            bad_data_sn = [bad_data_sn i];
            % data(data >= GUIh.IterCtrl.PrepDataDeadPThres_editfield.Value) = 0;
            bad_data_mask{i} = data >= 1E7;
            data = data.*~bad_data_mask{i};
            individual_mask{i} = or(effective_mask,bad_data_mask{i});
        else
            individual_mask{i} = effective_mask;
        end
        individual_mask_active_area(i) = sum(~individual_mask{i},'all');
        
        measured_amp{i} = sqrt(data);
        measured_amp_max(i,1) = max(measured_amp{i},[],'all');
        %GUIh.PlotArea.DataSn_editfield.Value = i;
        %imagesc(GUIh.PlotArea.Data_axes,measured_amp{i},[0,measured_amp_max(i,1)]);
        %title(GUIh.PlotArea.Data_axes,sprintf('Raw data max. = %.1f',max(data,[],'all')))
        %axis(GUIh.PlotArea.Data_axes,'image')
        %drawnow
        
        fprintf('Done.\n');
    end
    
    % convert cell to matrix for speedup when GPU applied
    measured_amp_temp = zeros(effective_clip_size,effective_clip_size,numel(measured_amp));
    individual_mask_temp = false(effective_clip_size,effective_clip_size,numel(individual_mask));
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
