function ptycho_package(filename)
init_cond = evalin('base','init_cond');
mask_info= evalin('base','mask_info');
measurement_info= evalin('base','measurement_info');
object_info= evalin('base','object_info');
probe_info= evalin('base','probe_info');
iteration_para= evalin('base','iteration_para');
filename = sprintf('%s.mat',filename);
FP = fullfile(init_cond.projectFF,filename);
save(FP,'init_cond','mask_info','measurement_info','object_info','probe_info','iteration_para','-v7.3')

