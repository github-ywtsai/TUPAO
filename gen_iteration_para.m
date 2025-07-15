function iteration_para = gen_iteration_para(init_cond,measurement_info,object_info,probe_info)
    config_iteration_table = init_cond.config_tables.config_iteration_table;
    iteration_para.FinishedRun = 0;
    iteration_para.chi2 = single(zeros(500000,1));
    
    % parameters for condition check in iterations
    iteration_para.scan_num = 0;
    iteration_para.pos_correction_status = 0;            
    iteration_para.pause_sign = 0;
    iteration_para.stop_sign = 0;
    
    iteration_para.max_iteration_num = config_iteration_table{'max_iteration','Value'}{1};
    iteration_para.saveing_section_file_peroid = config_iteration_table{'save_results_period','Value'}{1};
    

    iteration_para.alpha = config_iteration_table{'alpha_object_update','Value'}{1};
    iteration_para.alpha_current = 0;

    iteration_para.beta_start_pt = config_iteration_table{'probe_update_start','Value'}{1};
    iteration_para.beta = config_iteration_table{'beta_probe_update','Value'}{1};
    iteration_para.beta_current = 0;
    

    iteration_para.real_space_constraint_start_pt = config_iteration_table{'real_space_constraint_start','Value'}{1};
    iteration_para.real_space_constraint_period = 1;
    iteration_para.real_space_constraint_factor = config_iteration_table{'real_space_constraint_factor','Value'}{1};
    
    iteration_para.interesting_table = 1:init_cond.n_of_data;
    %{
    % skip data with deadpoint
    [~,n_of_bad_data] = size(measurement_info.bad_data_sn);
    if n_of_bad_data ~= 0
        for i = 1:n_of_bad_data
            iteration_para.interesting_table(iteration_para.interesting_table == init_cond.bad_data_sn(i)) = [];
        end
    end
    %}
    

    iteration_para.pos_corr_start_pt = config_iteration_table{'pos_corr_start','Value'}{1};
    iteration_para.pos_corr_period = config_iteration_table{'pos_corr_period','Value'}{1};
    iteration_para.pos_corr_extend_pts = config_iteration_table{'pos_corr_points','Value'}{1};
    iteration_para.pos_corr_extend_range = config_iteration_table{'pos_corr_range_nm','Value'}{1};
    

    iteration_para.probe_deny_area_factor = config_iteration_table{'probe_deny_area_factor','Value'}{1};
    iteration_para.probe_deny_reducing_ratio = config_iteration_table{'probe_deny_reduce_ratio','Value'}{1};
    iteration_para.probe_deny_mask = gen_probe_deny_mask(init_cond,iteration_para.probe_deny_area_factor,iteration_para.probe_deny_reducing_ratio);
    
    iteration_para.draw_results = config_iteration_table{'draw_results_period','Value'}{1};
end

function probe_deny_mask = gen_probe_deny_mask(init_cond,probe_deny_area_factor,probe_deny_reducing_ratio)
    if probe_deny_area_factor > 1
        probe_deny_area_factor = 1;
    elseif probe_deny_area_factor < 0
        probe_deny_area_factor = 0;
    end
    
    probe_deny_mask = ones(init_cond.effective_clip_size)*(1-probe_deny_reducing_ratio);
    cen_idx = (init_cond.effective_clip_size + 1)/2;
    range = round((init_cond.effective_clip_size - 1)/2*(1-probe_deny_area_factor));
    allow_idx = (cen_idx-range) : (cen_idx+range);
    probe_deny_mask(allow_idx,allow_idx) = 1;


end