function sectionfile_info = gen_sectionfile_info(init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para)
    tic
       

    section_fn = sprintf('SN%04d_Run%d.mat',measurement_info.SN,iteration_para.FinishedRun);
    section_fp = fullfile(init_cond.results_path, section_fn);
    sectionfile_info.SN = measurement_info.SN;
    sectionfile_info.CreateTime{1} = datetime;
    sectionfile_info.FileNmae{1} = section_fn;
    
    fprintf('Saving iteration results to %s ...\t',section_fn);
    save(section_fp,'init_cond','mask_info','measurement_info','object_info','probe_info','iteration_para','sectionfile_info','-v7.3');
    ElapsedT = toc;
    fprintf('Done(%.1f sec).\n',ElapsedT);
end