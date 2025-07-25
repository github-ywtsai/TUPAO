function run_single_data_script(projectFF)

ptycho_package = tools.initialize(projectFF);

init_cond = ptycho_package.init_cond;
mask_info = ptycho_package.mask_info;
measurement_info = ptycho_package.measurement_info;
object_info = ptycho_package.object_info;
probe_info = ptycho_package.probe_info;
iteration_para = ptycho_package.iteration_para;

[object_info, probe_info, iteration_para] =  ptycho(init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para);

assignin('base','init_cond',init_cond);
assignin('base','mask_info',mask_info);
assignin('base','measurement_info',measurement_info);
assignin('base','object_info',object_info);
assignin('base','probe_info',probe_info);
assignin('base','iteration_para',iteration_para);