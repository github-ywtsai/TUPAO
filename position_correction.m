function updated_pos_correct_pixel = position_correction(measured_amp,init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para)
    
    % new version, only using 1st mode of probe for position collection.
    probe_info_temp = probe_info;
    if probe_info_temp.Mp > 1
        for probe_sn = 2:probe_info_temp.Mp
            probe_info_temp.real_space(:,:,probe_sn) = probe_info_temp.real_space(:,:,probe_sn)*1E-6;
        end
        probe_info_temp = normalize_probe(probe_info_temp,measurement_info);
    end
    probe = probe_info_temp.real_space;
    
    % olde version, using all probe for position correction
    % probe = probe_info.real_space;
    
    object = object_info.real_space;
    updated_pos_correct_pixel = zeros(init_cond.n_of_data,2);
    individual_mask = measurement_info.individual_mask;
    individual_mask_active_area = measurement_info.individual_mask_active_area;
    measured_amp_max = measurement_info.measured_amp_max;
    

    
    % create position extend matrix - mesh
    %{
    [row_extend_pixel,col_extend_pixel,~] = find(ones(iteration_para.pos_corr_extend_pts*2+1));
    % row_extend_pixel = reshape(row_extend_pixel,iteration_para.pos_corr_extend_pts*2+1,iteration_para.pos_corr_extend_pts*2+1); % change size to 2N + 1 matrix
    % col_extend_pixel = reshape(col_extend_pixel,iteration_para.pos_corr_extend_pts*2+1,iteration_para.pos_corr_extend_pts*2+1); % change size to 2N + 1 matrix
    row_extend_pixel = row_extend_pixel - (iteration_para.pos_corr_extend_pts + 1); % shift center to symmetry point
    col_extend_pixel = col_extend_pixel - (iteration_para.pos_corr_extend_pts + 1); % shift center to symmetry point
    row_extend_pixel  = row_extend_pixel * iteration_para.pos_corr_row_extend_res; % change res
    col_extend_pixel  = col_extend_pixel * iteration_para.pos_corr_col_extend_res; % change res
    [n_of_extended_pos,~] = size(row_extend_pixel); n_of_extended_pos = gpuArray(single(n_of_extended_pos));
    %}
    
    n_of_extended_pos = iteration_para.pos_corr_extend_pts;
    for data_sn = iteration_para.interesting_table
        % create position extend matrix - rand
        r = rand(n_of_extended_pos,1) * iteration_para.pos_corr_extend_range;
        th = rand(n_of_extended_pos,1)* 2 * pi;
        row_extend_pixel = round(r.*sin(th));
        col_extend_pixel = round(r.*cos(th));
        
        data = measured_amp(:,:,data_sn);
        mask = individual_mask(:,:,data_sn);
        active_area = individual_mask_active_area(data_sn);
        
        
        exp_cen_row_idx = object_info.exp_pos_idx(data_sn,1) + probe_info.pos_correct_pixel(data_sn,1) + row_extend_pixel;
        exp_cen_col_idx = object_info.exp_pos_idx(data_sn,2) + probe_info.pos_correct_pixel(data_sn,2) + col_extend_pixel;
        
        row_start_idx = exp_cen_row_idx - (init_cond.effective_clip_size - 1)/2;
        row_end_idx = exp_cen_row_idx + (init_cond.effective_clip_size - 1)/2;
        col_start_idx = exp_cen_col_idx - (init_cond.effective_clip_size - 1)/2;
        col_end_idx = exp_cen_col_idx + (init_cond.effective_clip_size - 1)/2;
        
        
        chi2_for_PC = zeros(1,n_of_extended_pos,'single');
        
        if init_cond.using_GPU
            chi2_for_PC = gpuArray(chi2_for_PC);
        end
        
        for extended_pos_sn = 1 : n_of_extended_pos
            clip_object = object(row_start_idx(extended_pos_sn):row_end_idx(extended_pos_sn),col_start_idx(extended_pos_sn):col_end_idx(extended_pos_sn));   

            % formula (3) S(12)
            psi = probe .* clip_object;
            % psi(:,:,k)
            % k for probe(Mp)
            % psi for real space in (S12)
            
            % formula (4)
            Psi = fftshift(fft2(ifftshift(psi)));
            % Psi(:,:,k)           
            
            % formula (5), (S11)
            Psi_amp_flat = sqrt(sum(abs(Psi).^2,3)); % flat matrix.
            % multi-probe results are non-interference. using sqrt(|a|^2 + |b|^2 + ....)
            
            % calculate chi^2
            chi2_for_PC(1,extended_pos_sn) = sum((Psi_amp_flat -data).^2.*~mask,'all')/active_area;
        end
        [~,min_chi2_idx] = min(chi2_for_PC);
        updated_pos_correct_pixel(data_sn,:) = probe_info.pos_correct_pixel(data_sn,:) + [row_extend_pixel(min_chi2_idx) col_extend_pixel(min_chi2_idx)];
        
    end
end