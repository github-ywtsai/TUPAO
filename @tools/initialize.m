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

ptycho_package.init_cond = init_cond;
ptycho_package.mask_info = mask_info;
ptycho_package.measurement_info = measurement_info;
ptycho_package.object_info = object_info;
ptycho_package.probe_info = probe_info;

% arange GPU
idle_GPU_index = tools.find_idle_GPU();
gpuDevice(idle_GPU_index);
fprintf('Auto arrange GPU %d...\n',idle_GPU_index);