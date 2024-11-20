setenv('HDF5_PLUGIN_PATH','/blsw/opt/areaDetector/root/usr/lib/h5plugin');
init_cond = gen_init_cond('/home/tsai.yw/gpfs_folder/TUPAO/projects/20241108_188/config_init_cond.txt');
mask_info = gen_mask_info(init_cond);
measurement_info = gen_measurement_info(init_cond,mask_info);
object_info = gen_object_info(init_cond);
probe_info = gen_probe_info(init_cond);
probe_info = normalize_probe(probe_info,measurement_info);
iteration_para = gen_iteration_para(init_cond,measurement_info,object_info,probe_info);


[object_info, probe_info, iteration_para] =  ptycho(init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para);