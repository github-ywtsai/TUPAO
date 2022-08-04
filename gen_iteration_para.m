function iteration_para = gen_iteration_para(init_cond,measurement_info,object_info,probe_info,iteration_configFP)
    iteration_para.FinishedRun = 0;
    iteration_para.chi2 = single(zeros(500000,1));
    
    % parameters for condition check in iterations
    iteration_para.scan_num = 0;
    iteration_para.pos_correction_status = 0;            
    iteration_para.pause_sign = 0;
    iteration_para.stop_sign = 0;
    
    
    Temp = readcell(iteration_configFP);
    Temp(1,:) = []; % remove header
    Value = Temp(:,1); Discription = Temp(:,2);
    
    iteration_para.max_iteration_num = Value{1};
    iteration_para.saveing_section_file_peroid = Value{2};
    

    iteration_para.alpha = Value{3};
    

    iteration_para.beta_start_pt = Value{4};
    iteration_para.beta = Value{5};
    iteration_para.beta_current = 0;
    

    iteration_para.real_space_constraint_start_pt = Value{6};
    iteration_para.real_space_constraint_period = 1;
    iteration_para.real_space_constraint_factor = Value{7};
    
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
    

    iteration_para.pos_corr_start_pt = Value{8};
    iteration_para.pos_corr_period = Value{9};
    iteration_para.pos_corr_extend_pts = Value{10};
    iteration_para.pos_corr_extend_range = round(Value{11}/(object_info.x_res*1E9));
    

    iteration_para.probe_deny_area_factor = Value{12};
    iteration_para.probe_deny_reducing_ratio = Value{13};
    iteration_para.probe_deny_mask = gen_probe_deny_mask(init_cond,iteration_para.probe_deny_area_factor,iteration_para.probe_deny_reducing_ratio);
    
    iteration_para.draw_results = Value{14};
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