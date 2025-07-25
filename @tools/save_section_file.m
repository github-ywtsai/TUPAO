function save_section_file(fn,ptycho_package)

init_cond = ptycho_package.init_cond;
iteration_para = ptycho_package.iteration_para;
mask_info = ptycho_package.mask_info;
object_info = ptycho_package.object_info;
probe_info = ptycho_package.object_info;

fn = [fn,'.mat'];
ff = init_cond.projectFF;
fp = fullfile(ff,fn);

if exist(fp,'file')
    fprintf('File %s existed.\n',fp);
    fprintf('Skip saving process.\n')
    return
end

save(fp,'init_cond','iteration_para','mask_info','object_info','probe_info','-v7.3'); % don't save measurement_info

