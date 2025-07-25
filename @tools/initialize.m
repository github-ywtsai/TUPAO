function ptycho_package = initialize(projectFF)

setenv('HDF5_PLUGIN_PATH','/blsw/opt/areaDetector/root/usr/lib/h5plugin');

tools.create_config_tables_from_files(projectFF);
config_tables = load(fullfile(projectFF,'config_tables.mat'));

init_cond = gen_init_cond(config_tables);
mask_info = gen_mask_info(init_cond);
measurement_info = gen_measurement_info(init_cond,mask_info);
object_info = gen_object_info(init_cond);
probe_info = gen_probe_info(init_cond);
%probe_info = normalize_probe(probe_info,measurement_info);
iteration_para = gen_iteration_para(init_cond,measurement_info,object_info,probe_info);

ptycho_package.init_cond = init_cond;
ptycho_package.mask_info = mask_info;
ptycho_package.measurement_info = measurement_info;
ptycho_package.object_info = object_info;
ptycho_package.probe_info = probe_info;
ptycho_package.iteration_para = iteration_para;


%% assign varialbes to workspace directly
% assignin('base','init_cond',init_cond);
% assignin('base','mask_info',mask_info);
% assignin('base','measurement_info',measurement_info);
% assignin('base','object_info',object_info);
% assignin('base','probe_info',probe_info);
% assignin('base','iteration_para',iteration_para);