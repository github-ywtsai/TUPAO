function init_cond = gen_init_cond(ConfigFP)  
    init_cond = get_config(ConfigFP);
  
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

    
    %fprintf('Loading probe-positions...\t')
    temp = get_exp_pos(init_cond);
    init_cond.exp_pos = temp.exp_pos;
    init_cond.n_of_data = temp.n_of_data;
    %fprintf('Done.\n')
    
    %fprintf('Done.\n\n')
end

function output = get_exp_pos(init_cond)

    pos_record_fp = init_cond.pos_record_fp;
    fid = fopen(pos_record_fp);

    [~] = fgetl(fid);
    temp = fgets(fid);

    exp_pos_rbv_z_x = []; % z and than x in mm
    exp_pos_rbv_y = [];
    while ischar(temp)
        temp = strsplit(temp,',');
        NCol = numel(temp);
        if NCol == 5 || NCol == 3
            exp_pos_rbv_z_x = [exp_pos_rbv_z_x;[str2double(temp{2}),str2double(temp{3})]];
        elseif NCol == 4
            exp_pos_rbv_z_x = [exp_pos_rbv_z_x;[str2double(temp{2}),str2double(temp{3})]];
            exp_pos_rbv_y = [exp_pos_rbv_y;str2double(temp{4})];
        end
        temp = fgets(fid);
    end
    exp_pos_rbv_z_x = single(exp_pos_rbv_z_x);
    fclose(fid);
    
    [temp,~] = max(exp_pos_rbv_z_x);
    exp_pos_rbv_z_max = temp(1);
    exp_pos_rbv_x_max = temp(2);
    [temp,~] = min(exp_pos_rbv_z_x);
    exp_pos_rbv_z_min = temp(1);
    exp_pos_rbv_x_min = temp(2);

    exp_pos_rbv_z_cen = (exp_pos_rbv_z_max + exp_pos_rbv_z_min)/2;
    exp_pos_rbv_x_cen = (exp_pos_rbv_x_max + exp_pos_rbv_x_min)/2;

    exp_pos = -[exp_pos_rbv_z_x(:,1) - exp_pos_rbv_z_cen,exp_pos_rbv_z_x(:,2) - exp_pos_rbv_x_cen];
    output.exp_pos = exp_pos * 1E-3; % covert unit from [mm] to [m]
    [output.n_of_data, ~] = size(exp_pos);
    
end

function init_cond = get_config(ConfigFP)
    Temp = readcell(ConfigFP);
    Temp(1,:) = []; % remove header
    Value = Temp(:,1); Discription = Temp(:,2);

    %% get master file name and path
    % get master file folder
    MasterFF = Value{1};
    % get master file name
    MasterFN = Value{2};   
    init_cond.master_fp = fullfile(MasterFF,MasterFN);
    
    %% get manual mask file name and path
    % get mask file folder
    ManualMaskFF = Value{3};
    % get mask file name
    for Idx = 1:4
        ManualMaskFfL{Idx} = ManualMaskFF;
        ManualMaskFnL{Idx} = Value{3+Idx};
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
    
    
    %% get scan position file name and path
    % get master file folder
    ScanFP = Value{8};
    % get master file name
    ScanFN = Value{9};
    
    init_cond.pos_record_fp = fullfile(ScanFP,ScanFN);
    
    %% get prepared data files and iteration files path
    init_cond.results_path = fullfile(Value{10});
    if ~exist(init_cond.results_path,'dir')
        mkdir(init_cond.results_path);
    end

    %% get beam center x and y
    beam_center_x = Value{11};
    if strcmpi(beam_center_x,'Auto')
        beam_center_x = autoload_ExpStat(init_cond.master_fp,'beam_center_x');
    end
    init_cond.beam_center_x = round(beam_center_x);
    

    beam_center_y = Value{12};
    if strcmpi(beam_center_y,'Auto')
        beam_center_y = autoload_ExpStat(init_cond.master_fp,'beam_center_y');
    end
    init_cond.beam_center_y = round(beam_center_y);
    
    %% get wavelength
    wavelength = Value{13};
    if strcmpi(wavelength,'Auto')
        wavelength = autoload_ExpStat(init_cond.master_fp,'wavelength'); % [m]
    else
        wavelength = wavelength*1E-10; % read in [A] and store in [m]
    end
    init_cond.wavelength = wavelength; % save as [m]
    
    
    %% get sample to detector distance and detector pixel sizes
    detector_distance = Value{14};
    if strcmpi(detector_distance,'Auto')
        detector_distance = autoload_ExpStat(init_cond.master_fp,'detector_distance'); % [m]
    end
    init_cond.detector_distance = detector_distance; % [m]
    
    x_pixel_size = Value{15};
    if strcmpi(x_pixel_size,'Auto')
        x_pixel_size = autoload_ExpStat(init_cond.master_fp,'x_pixel_size'); % [m]
    else
        x_pixel_size = x_pixel_size*1E-6; % input in [um] and convert to [m]
    end
    init_cond.x_pixel_size = x_pixel_size; % [m]
    
    y_pixel_size = Value{16};
    if strcmpi(y_pixel_size,'Auto')
        y_pixel_size = autoload_ExpStat(init_cond.master_fp,'y_pixel_size'); % [m]
    else
        y_pixel_size = y_pixel_size*1E-6; % input in [um] and convert to [m]
    end
    init_cond.y_pixel_size = y_pixel_size; % [m]
    
    init_cond.effective_x_pixel_size = init_cond.x_pixel_size;
    init_cond.effective_y_pixel_size = init_cond.y_pixel_size;
    
    
    %% get data clip size N*N
    init_cond.rawdata_clip_size = Value{17};
    if mod(init_cond.rawdata_clip_size,2) == 0
        init_cond.rawdata_clip_size = init_cond.rawdata_clip_size + 1;
    end
    init_cond.effective_clip_size = init_cond.rawdata_clip_size;
    
    
    %% get using GPU setting
    GPUDeviceNumber = Value{18};
    if strcmpi(GPUDeviceNumber,'None')
        init_cond.using_GPU = 0;
    else
        init_cond.using_GPU = GPUDeviceNumber;
    end
    
    
    %% get probe range extend condition
    init_cond.probe_extending_factor = Value{19};

    %% get random seed
    if strcmpi(Value{20},'None')
        rng('shuffle')
        init_cond.rand_seed = round(rand*1E8);
    else
        init_cond.rand_seed = Value{20};
        rng(init_cond.rand_seed)
    end
    
    %% determine core
    init_cond.core = Value{21};
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
