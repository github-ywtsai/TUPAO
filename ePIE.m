function [updated_object,updated_probe,chi2_sum] = ePIE(measured_amp,init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para)
    
    
    object = object_info.real_space;
    probe = probe_info.real_space;
    chi2_temp = zeros(1,init_cond.n_of_data);
    individual_mask = measurement_info.individual_mask;
    individual_mask_active_area = measurement_info.individual_mask_active_area;

    % info. for upstream probe constrain
    wavelength = init_cond.wavelength;
    probe_x_axis = probe_info.real_space_xaxis;
    probe_y_axis = probe_info.real_space_yaxis;
    propagating_dist = probe_info.ProbeConf.ApertureDist;
    upstream_ROI = probe_info.ProbeConf.upstream_ROI;
    
    if init_cond.using_GPU
        %object = gpuArray(object);
        % probe = gpuArray(probe);
        chi2_temp = gpuArray(chi2_temp);
        %upstream_ROI = gpuArray(upstream_ROI);
        %measured_amp = gpuArray(measured_amp);
        %individual_mask = gpuArray(individual_mask);
        %individual_mask_active_area = gpuArray(individual_mask_active_area);
    end
    
    if iteration_para.real_space_constraint_status == 1
            real_space_constraint_factor = iteration_para.real_space_constraint_factor;
            if real_space_constraint_factor > 0
                object(real(object)<0) = object(real(object)<0)*(1-real_space_constraint_factor);
            else
                real_space_constraint_factor = abs(real_space_constraint_factor);
                object(imag(object)<0) = object(imag(object)<0)*(1-real_space_constraint_factor);
            end
    end
    
    interesting_table = iteration_para.interesting_table;
    %if round(rand())
    %    interesting_table = fliplr(interesting_table);
    %end
    
    for data_sn = interesting_table
        data = measured_amp(:,:,data_sn);
        mask = individual_mask(:,:,data_sn);
        active_area = individual_mask_active_area(data_sn);
        
        exp_cen_row_idx = object_info.exp_pos_idx(data_sn,1) + probe_info.pos_correct_pixel(data_sn,1);
        exp_cen_col_idx = object_info.exp_pos_idx(data_sn,2) + probe_info.pos_correct_pixel(data_sn,2);
        row_start_idx = exp_cen_row_idx - (init_cond.effective_clip_size - 1)/2;
        row_end_idx = exp_cen_row_idx + (init_cond.effective_clip_size - 1)/2;
        col_start_idx = exp_cen_col_idx - (init_cond.effective_clip_size - 1)/2;
        col_end_idx = exp_cen_col_idx + (init_cond.effective_clip_size - 1)/2;   
        clip_object = gpuArray(object(row_start_idx:row_end_idx,col_start_idx:col_end_idx));

        % formula (3) S(12)
        psi = probe .* clip_object;
        % psi(:,:,k)
        % k for probe(Mp)
        % psi for real space in (S12)
        
        % formula (4)
        Psi = fftshift(fft2(ifftshift(psi))); % fixed 20220810
        % Psi(:,:,k)        
        
        % formula (5), (S11)
        Psi_amp_flat = sqrt(sum(abs(Psi).^2,3)); % flat matrix.
        % multi-probe results are non-interference. using sqrt(|a|^2 + |b|^2 + ....)
        Psi_amp_flat_non_zero_mask = Psi_amp_flat ~= 0;
        Psi_amp_flat(~Psi_amp_flat_non_zero_mask) = 1E30; % invoid the NaN resutls in Psi_amp_flat equals to 0;
        Psi_p = Psi.*mask + data .* Psi./Psi_amp_flat .* ~mask.*single(Psi_amp_flat_non_zero_mask);
        % calculate chi^2
        chi2_temp(1,data_sn) = sum(sum( (Psi_amp_flat -data).^2.*~mask))/active_area;
        clear Psi Psi_amp_flat_non_zero_mask Psi_amp_flat data
        psi_p = fftshift(ifft2(ifftshift(Psi_p))); % fixed 20220810
        diff_psi_p_psi = psi_p - psi;
        clear psi_p psi
        % Psi_amp_flat, Psi_p and psi_p are the matrix in [clip_size,clip_size,Mp]
        
        % update object
        % formula (6) (S21)
        first_term = iteration_para.alpha / max(max( sum(abs(probe).^2,3)));
        second_term = sum(conj(probe).* diff_psi_p_psi,3); 
        object(row_start_idx:row_end_idx,col_start_idx:col_end_idx) = gather(clip_object + first_term * second_term);
        
        %update probe
        % formula (7) (S22)
        first_term = iteration_para.beta_current / max(max( sum(abs(clip_object).^2,3)));
        second_term = conj(clip_object).* diff_psi_p_psi;
        probe = probe + first_term * second_term;
        
        if probe_info.ProbeConf.UpStreamConstrain
            [upstream_probe,upstream_x_axis,upstream_y_axis] = propagate_probe(-propagating_dist,probe,wavelength,probe_x_axis,probe_y_axis);
            upstream_probe = upstream_probe.*upstream_ROI;
            [probe,~,~] = propagate_probe(propagating_dist,upstream_probe,wavelength,upstream_x_axis,upstream_y_axis);
        end
        
        clear first_term second_term diff_psi_p_psi
        
        
    end
    
    chi2_sum = gather(sum(chi2_temp));
    updated_object = gather(object);
    updated_probe = gather(probe);
    
    
end