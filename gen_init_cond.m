function init_cond = gen_init_cond(config_tables)    
    init_cond = get_config(config_tables);
  
    master_fp = init_cond.master_fp;
    init_cond.data_precision = 'single'; % single or double
    init_cond.bit_depth_image = double(h5read(master_fp,'/entry/instrument/detector/bit_depth_image'));
    init_cond.bad_point_value = power(2,init_cond.bit_depth_image)-1;
    init_cond.x_pixels_in_detector = double(h5read(master_fp,'/entry/instrument/detector/detectorSpecific/x_pixels_in_detector'));
    init_cond.y_pixels_in_detector = double(h5read(master_fp,'/entry/instrument/detector/detectorSpecific/y_pixels_in_detector'));
    init_cond.count_time = double(h5read(master_fp,'/entry/instrument/detector/count_time'));

    %% modify parameters when probe_extending_factor ~= 1
    if init_cond.probe_extending_factor ~= 1
        init_cond.effective_clip_size = (init_cond.effective_clip_size -1)* init_cond.probe_extending_factor + 1;
        init_cond.effective_x_pixel_size = init_cond.x_pixel_size * init_cond.rawdata_clip_size/init_cond.effective_clip_size;
        init_cond.effective_y_pixel_size = init_cond.y_pixel_size * init_cond.rawdata_clip_size/init_cond.effective_clip_size;      
    end

    init_cond.pixel_res = init_cond.wavelength * init_cond.detector_distance / init_cond.effective_clip_size / init_cond.effective_x_pixel_size;
    init_cond.CDI_window = init_cond.pixel_res*init_cond.effective_clip_size;
    fprintf('Pixel resolution = %.1f [nm]\n', init_cond.pixel_res*1E9);
    fprintf('CDI window = %.1f [um]\n',init_cond.CDI_window*1E6);

    
end

% Note for get exposure position
% --The defination of the global coordination--
% +y axis along the beam
% +z axis is the up direction from the floor
% +x created by +y and +z by the right hand rule
%-- Stage direction and the sign of the exp_pos --
% If the +x and +z directions of stage along the same directions of the global
% coordination, add - sign on both read back value of x and z.
% Ex:
% exp_pos = -[exp_pos_rbv_z_x(:,1) - exp_pos_rbv_z_cen,exp_pos_rbv_z_x(:,2) - exp_pos_rbv_x_cen];
% If the +z direction and -x direction of stage along the same directions of the global
% coordination, add - sign on z axis only.
% Ex:
% exp_pos = [(exp_pos_rbv_z_x(:,1) - exp_pos_rbv_z_cen)*-1,exp_pos_rbv_z_x(:,2) - exp_pos_rbv_x_cen];

function output = get_exp_cond_bluesky(exp_cond_record_fp)
    table_temp = readtable(exp_cond_record_fp,'Delimiter',{' ',','},'VariableNamingRule','preserve');
    VariableNames = table_temp.Properties.VariableNames;
    MasterFPattern_idx = find(cellfun(@(x) ~isempty(regexpi(x,'eig\d*m_file_file_write_name_pattern')),VariableNames));
    xidx = find(cellfun(@(X)strcmpi(X,'_cisamf_x'),VariableNames));
    zidx = find(cellfun(@(X)strcmpi(X,'_cisamf_z'),VariableNames));
    MasterFPattern = table_temp.Properties.VariableNames{MasterFPattern_idx};
    xaxis_variablename = table_temp.Properties.VariableNames{xidx};
    zaxis_variablename = table_temp.Properties.VariableNames{zidx};
    cmd = sprintf('masterfilepattern = table_temp.%s;',MasterFPattern);
    eval(cmd);
    cmd = sprintf('xpos = table_temp.(''%s'');',xaxis_variablename);
    eval(cmd);
    cmd = sprintf('zpos = table_temp.(''%s'');',zaxis_variablename);
    eval(cmd);
    MasterFN = sprintf('%s_master.h5',masterfilepattern{1});
    exp_pos_rbv_z_x = [zpos,xpos];
    exp_pos_rbv_z_x = single(exp_pos_rbv_z_x);% in um
    [temp,~] = max(exp_pos_rbv_z_x);
    exp_pos_rbv_z_max = temp(1);
    exp_pos_rbv_x_max = temp(2);
    [temp,~] = min(exp_pos_rbv_z_x);
    exp_pos_rbv_z_min = temp(1);
    exp_pos_rbv_x_min = temp(2);

    exp_pos_rbv_z_cen = (exp_pos_rbv_z_max + exp_pos_rbv_z_min)/2;
    exp_pos_rbv_x_cen = (exp_pos_rbv_x_max + exp_pos_rbv_x_min)/2;
    
    % in TPS 25A2 and SP8 12XU
    % add - on z an + on x
    z_direction_modification = -1; % for TPS 25A2
    x_direction_modification = 1; % for TPS 25A2
    %{
    % old version, using relative exp_pos
    rel_exp_pos_rbv_z_x = [exp_pos_rbv_z_x(:,1) - exp_pos_rbv_z_cen,exp_pos_rbv_z_x(:,2) - exp_pos_rbv_x_cen];
    n_exp_pos = size(rel_exp_pos_rbv_z_x,1);
    exp_pos = rel_exp_pos_rbv_z_x.* [z_direction_modification*ones(n_exp_pos,1),x_direction_modification*ones(n_exp_pos,1)];
    exp_pos_cen = [z_direction_modification*exp_pos_rbv_z_cen,x_direction_modification*exp_pos_rbv_x_cen];
    %}
    % new version, using absolute position
    exp_pos = [z_direction_modification*exp_pos_rbv_z_x(:,1),x_direction_modification*exp_pos_rbv_z_x(:,2)];
    exp_pos_cen = [z_direction_modification*exp_pos_rbv_z_cen,x_direction_modification*exp_pos_rbv_x_cen];
    output.exp_pos = exp_pos*1E-6; % convert from um to m
    output.exp_pos_cen = exp_pos_cen*1E-6; % convert from um to m
    [output.n_of_data, ~] = size(exp_pos);
    output.MasterFN = MasterFN;
    output.z_direction_modification = z_direction_modification;
    output.x_direction_modification = x_direction_modification;
    
end

function init_cond = get_config(config_tables)
    init_cond.config_tables = config_tables;
    init_cond.projectFF = init_cond.config_tables.config_init_cond_table{'project_folder_path','Value'}{1};
    config_init_cond_tables = config_tables.config_init_cond_table;
    
    % get  file folder
    DataFF = config_init_cond_tables{'folder_path','Value'}{1};
    % get expeirmentcal condition file name and file path
    ScanFN = config_init_cond_tables{'scan_file_name','Value'}{1};
    buffer = dir(fullfile(DataFF,ScanFN));
    init_cond.exp_cond_record_fp = fullfile(buffer.folder,buffer.name);
    
    % get master fn and exposure position
    temp = get_exp_cond_bluesky(init_cond.exp_cond_record_fp);
    buffer = dir(fullfile(DataFF,temp.MasterFN));
    init_cond.master_fp = fullfile(buffer.folder,buffer.name);
    init_cond.exp_pos = temp.exp_pos;
    init_cond.exp_pos_cen = temp.exp_pos_cen;
    init_cond.n_of_data = temp.n_of_data;
    init_cond.x_direction_modification = temp.x_direction_modification;
    init_cond.z_direction_modification = temp.z_direction_modification;
    

    % get mask file name and file path
    for Idx = 1:4
        ManualMaskFfL{Idx} = DataFF;
        ManualMaskFnL{Idx} = config_init_cond_tables{sprintf('mask_file_%d',Idx),'Value'}{1};
        if strcmpi(ManualMaskFnL{Idx},'None')
            ManualMaskFnL{Idx} = [];
        end
    end
    ManualMaskActIdx = find(~cellfun(@isempty,ManualMaskFnL));
    if isempty(ManualMaskActIdx)
        init_cond.manual_mask_fp = [];
    else
        [~,ManualMaskActN] = size(ManualMaskActIdx);
        for loopidx = 1:ManualMaskActN
            init_cond.manual_mask_fp{loopidx} = fullfile(ManualMaskFfL{ManualMaskActIdx(loopidx)}, ManualMaskFnL{ManualMaskActIdx(loopidx)});
        end
    end

    %% get beam center x and y
    beam_center_x = config_init_cond_tables{'beam_center_X','Value'}{1};
    if strcmpi(beam_center_x,'Auto')
        beam_center_x = autoload_ExpStat(init_cond.master_fp,'beam_center_x');
    end
    init_cond.beam_center_x = round(beam_center_x);
    

    beam_center_y = config_init_cond_tables{'beam_center_Y','Value'}{1};
    if strcmpi(beam_center_y,'Auto')
        beam_center_y = autoload_ExpStat(init_cond.master_fp,'beam_center_y');
    end
    init_cond.beam_center_y = round(beam_center_y);
    
    %% get wavelength
    wavelength = config_init_cond_tables{'wavelength_A','Value'}{1};
    if strcmpi(wavelength,'Auto')
        wavelength = autoload_ExpStat(init_cond.master_fp,'wavelength'); % [m]
    else
        wavelength = wavelength*1E-10; % read in [A] and store in [m]
    end
    init_cond.wavelength = wavelength; % save as [m]
    
    
    %% get sample to detector distance and detector pixel sizes
    detector_distance = config_init_cond_tables{'detector_distance_m','Value'}{1};
    if strcmpi(detector_distance,'Auto')
        detector_distance = autoload_ExpStat(init_cond.master_fp,'detector_distance'); % [m]
    end
    init_cond.detector_distance = detector_distance; % [m]
    
    x_pixel_size = config_init_cond_tables{'pixel_size_X_um','Value'}{1};
    if strcmpi(x_pixel_size,'Auto')
        x_pixel_size = autoload_ExpStat(init_cond.master_fp,'x_pixel_size'); % [m]
    else
        x_pixel_size = x_pixel_size*1E-6; % input in [um] and convert to [m]
    end
    init_cond.x_pixel_size = x_pixel_size; % [m]
    
    y_pixel_size = config_init_cond_tables{'pixel_size_Y_um','Value'}{1};
    if strcmpi(y_pixel_size,'Auto')
        y_pixel_size = autoload_ExpStat(init_cond.master_fp,'y_pixel_size'); % [m]
    else
        y_pixel_size = y_pixel_size*1E-6; % input in [um] and convert to [m]
    end
    init_cond.y_pixel_size = y_pixel_size; % [m]
    
    init_cond.effective_x_pixel_size = init_cond.x_pixel_size;
    init_cond.effective_y_pixel_size = init_cond.y_pixel_size;
    
    
    %% get data clip size N*N
    init_cond.rawdata_clip_size = config_init_cond_tables{'data_clipping_size','Value'}{1};
    if mod(init_cond.rawdata_clip_size,2) == 0
        init_cond.rawdata_clip_size = init_cond.rawdata_clip_size + 1;
    end
    init_cond.effective_clip_size = init_cond.rawdata_clip_size;
    
    
    %% get probe range extend condition
    init_cond.probe_extending_factor = config_init_cond_tables{'probe_extend_factor','Value'}{1};

    %% get random seed
    if strcmpi(config_init_cond_tables{'random_seed','Value'}{1},'None')
        rng('shuffle')
        init_cond.rand_seed = round(rand*1E8);
    else
        init_cond.rand_seed = config_init_cond_tables{'random_seed','Value'}{1};
        rng(init_cond.rand_seed)
    end
    
    %% determine core
    init_cond.core = config_init_cond_tables{'core_algorithm','Value'}{1};
end

function ReturnValue = autoload_ExpStat(master_fp,Target)
    if ~exist(master_fp,'file')
        fprintf('!! Master file does not exist. !!\n')
        return
    end
    switch Target
        % Beam center
        case 'beam_center_x'
            ReturnValue = double(h5read(master_fp,'/entry/instrument/detector/beam_center_x'));
        case 'beam_center_y'
            ReturnValue = double(h5read(master_fp,'/entry/instrument/detector/beam_center_y'));
        % Pixel size
        case 'x_pixel_size'
            ReturnValue = double(h5read(master_fp,'/entry/instrument/detector/x_pixel_size')); % read in meter
        case 'y_pixel_size'
            ReturnValue = double(h5read(master_fp,'/entry/instrument/detector/y_pixel_size')); % read in meter
        % sample to detctor distance
        case 'detector_distance'
            ReturnValue = double(h5read(master_fp,'/entry/instrument/detector/detector_distance')); % read in meter
        % wavelength and E
        case 'wavelength'
            ReturnValue = double(h5read(master_fp,'/entry/instrument/beam/incident_wavelength'))*1E-10; % read as [A], save as [m]
    end
end
