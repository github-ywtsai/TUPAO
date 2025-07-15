setenv('HDF5_PLUGIN_PATH','/blsw/opt/areaDetector/root/usr/lib/h5plugin');


%% define project folder and preset configurations
project_folder = '/home/tsai.yw/gpfs_folder/TUPAO/projects/20250701_8405-8429';
config_tables_temp = load(fullfile(project_folder,'config_tables_temp.mat'));
start_SN = 8417-7;
end_SN = 8417+7;
dataSN_list = start_SN:end_SN;
init_cond_list = cell(length(dataSN_list), 1);
mask_info_list = cell(length(dataSN_list), 1);
measurement_info_list = cell(length(dataSN_list), 1);

%% start merge data
for sn = 1:length(dataSN_list)
    scan_filename = sprintf('250701_ptycho-scan_id-%d-primary.csv',dataSN_list(sn));
    config_tables = config_tables_temp;
    config_tables.config_init_cond_table{'scan_file_name','Value'} = {scan_filename};
    init_cond = gen_init_cond(config_tables);
    mask_info = gen_mask_info(init_cond);
    measurement_info = gen_measurement_info(init_cond,mask_info);
    
    init_cond_list{sn} = init_cond;
    mask_info_list{sn} = mask_info;
    measurement_info_list{sn} = measurement_info;
end

% merge init_cond
init_cond = init_cond_list{1};
%init_cond.config_tables = 'merged dataset';
init_cond.exp_cond_record_fp = 'merged dataset';
init_cond.master_fp = 'merged dataset';
all_exp_pos = cell(length(dataSN_list), 1);
for sn = 1:length(dataSN_list)
    all_exp_pos{sn} = init_cond_list{sn}.exp_pos + init_cond_list{sn}.exp_pos_cen;
end
all_exp_pos = vertcat(all_exp_pos{:});
exp_pos_z_cen = (max(all_exp_pos(:,1)) + min(all_exp_pos(:,1)))/2;
exp_pos_x_cen = (max(all_exp_pos(:,2)) + min(all_exp_pos(:,2)))/2;
init_cond.exp_pos_cen = [exp_pos_z_cen,exp_pos_x_cen];
init_cond.exp_pos = all_exp_pos - init_cond.exp_pos_cen;
[init_cond.n_of_data,~] = size(init_cond.exp_pos);

% merge measurement_info
all_bad_data_mask = cell(length(dataSN_list), 1);
all_individual_mask = cell(length(dataSN_list), 1);
all_individual_mask_active_area = cell(length(dataSN_list), 1);
all_measured_amp_max = cell(length(dataSN_list), 1);
all_measured_amp = cell(length(dataSN_list), 1);
for sn = 1:length(dataSN_list)
    all_bad_data_mask{sn} = measurement_info_list{sn}.bad_data_mask;
    all_individual_mask{sn} = measurement_info_list{sn}.individual_mask;
    all_individual_mask_active_area{sn} = measurement_info_list{sn}.individual_mask_active_area;
    all_measured_amp_max{sn} = measurement_info_list{sn}.measured_amp_max;
    all_measured_amp{sn} = measurement_info_list{sn}.measured_amp;
end

measurement_info.bad_data_sn = []; % skip bad data function didn't work in merge mode
measurement_info.bad_data_mask = vertcat(all_bad_data_mask{:});
measurement_info.individual_mask = cat(3,all_individual_mask{:});
measurement_info.individual_mask_active_area = vertcat(all_individual_mask_active_area{:});
measurement_info.measured_amp_max = vertcat(all_measured_amp_max{:});
measurement_info.measured_amp = cat(3,all_measured_amp{:});

% using merged data to create probe and object
object_info = gen_object_info(init_cond);
probe_info = gen_probe_info(init_cond);
%probe_info = normalize_probe(probe_info,measurement_info);
iteration_para = gen_iteration_para(init_cond,measurement_info,object_info,probe_info);

clearvars -except init_cond mask_info measurement_info object_info probe_info iteration_para

[object_info, probe_info, iteration_para] =  ptycho(init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para);
