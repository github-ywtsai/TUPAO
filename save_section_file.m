function save_section_file(fn)

init_cond = evalin('base','init_cond');
iteration_para = evalin('base','iteration_para');
mask_info = evalin('base','mask_info');
%measurement_info = evalin('base','measurement_info');
object_info = evalin('base','object_info');
probe_info = evalin('base','probe_info');

fn = [fn,'.mat'];
ff = init_cond.projectFF;
fp = fullfile(ff,fn);

if exist(fp,'file')
    fprintf('File %s existed.\n',fp);
    fprintf('Skip saving process.\n')
    return
end

%save(fp,'init_cond','iteration_para','mask_info','measurement_info','object_info','probe_info','-v7.3');
save(fp,'init_cond','iteration_para','mask_info','object_info','probe_info','-v7.3'); % don't save measurement_info

