function ptycho_package = initialize(projectFF,options)
% usage:
% auto arange GPU (default) : ptycho_package = initialize(projectFF,'arange_GPU',True)
% don't arange GPU : ptycho_package = initialize(projectFF,'arange_GPU',False)
arguments
    projectFF
    
    options.arange_GPU = True
end
arange_GPU = options.arange_GPU;

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
if arange_GPU
    idle_GPU_index = tools.find_idle_GPU();
    gpuDevice(idle_GPU_index);
    fprintf('Auto arrange GPU %d...\n',idle_GPU_index);
else
    fprintf('GPU did not be arranged automatically.\n');
    fprintf('GPU must be arraged manully using command: gpuDevice(GPU#).\n');
end