function ptycho_package = PIE(ptycho_package,options)
arguments
    ptycho_package
    
    options.core = 'rPIE'
    options.alpha = 0.9 % for object update
    options.beta = 0.5 % for probe update
    options.iteration_num = 5
    options.real_space_constraint = 'imaginary'
    options.real_space_constraint_factor = 0.9
    % ====== real_space_constraint & real_space_constraint_factor ======
    % real_space_constraint will be applied on imaginary or real part(using 'imaginary' or 'real' or 'none' to switch).
    % The rule is: the imaginary part or real part must be > 0 on object.
    % While pixels didn't suitable the rule, the value of these pixels will
    % be reduced by the real_space_constraint_factor.
    % ex: new_value = old_value * real_space_constraint_factor
end
core = options.core;
alpha = options.alpha;
beta = options.beta;
iteration_num = options.iteration_num;
real_space_constraint = options.real_space_constraint;
real_space_constraint_factor = options.real_space_constraint_factor;

    % reference: https://stackoverflow.com/questions/32166879/do-i-need-to-call-fftshift-before-calling-fft-or-ifft
    % when should apply ifftshift and fftshift?
    % simple answer:
    % spectrum = fftshift(fft2(ifftshift(myimage))
    % myimage = fftshift(ifft2(ifftshift(spectrum))
    
    object = ptycho_package.object_info.real_space;
    probe = ptycho_package.probe_info.real_space;
    %chi2_temp = zeros(1,init_cond.n_of_data);
    individual_mask = ptycho_package.measurement_info.individual_mask;
    individual_mask_active_area = ptycho_package.measurement_info.individual_mask_active_area;

    %gpu_chi2_temp = gpuArray(chi2_temp); % put chi2_temp into gpu


    switch real_space_constraint
        case 'real'
            object(real(object)<0) = object(real(object)<0)*real_space_constraint_factor;
        case 'imaginary'
            object(imag(object)<0) = object(imag(object)<0)*real_space_constraint_factor;
    end
    
    %{
    % didn't work well yet
    % info. for upstream probe constrain
    wavelength = init_cond.wavelength;
    probe_x_axis = probe_info.real_space_xaxis;
    probe_y_axis = probe_info.real_space_yaxis;
    propagating_dist = probe_info.ProbeConf.ApertureDist;
    upstream_ROI = probe_info.ProbeConf.upstream_ROI;
    if probe_info.ProbeConf.UpStreamConstrain
            [upstream_probe,upstream_x_axis,upstream_y_axis] = propagate_probe(-propagating_dist,probe,wavelength,probe_x_axis,probe_y_axis);
            upstream_probe = upstream_probe.*upstream_ROI;
            [probe,~,~] = propagate_probe(propagating_dist,upstream_probe,wavelength,upstream_x_axis,upstream_y_axis);
    end
    %}
    
    interesting_table = iteration_para.interesting_table;

    % assign varible to gpu
    idle_GPU_index = tools.find_idle_GPU();
    gpuDevice(idle_GPU_index);
    fprintf('Auto arrange GPU %d...\n',idle_GPU_index);
    gpu_measured_amp = gpuArray(ptycho_package.measurement_info.measured_amp);
    gpu_probe = gpuArray(probe);
    gpu_object = gpuArray(object);
    gpu_individual_mask = gpuArray(individual_mask);
    gpu_individual_mask_active_area = gpuArray(individual_mask_active_area);
    
    for data_sn = interesting_table
        data = gpu_measured_amp(:,:,data_sn);
        mask = gpu_individual_mask(:,:,data_sn);
        active_area = gpu_individual_mask_active_area(data_sn);
        
        exp_cen_row_idx = ptycho_package.object_info.exp_pos_idx(data_sn,1) + ptycho_package.probe_info.pos_correct_pixel(data_sn,1);
        exp_cen_col_idx = ptycho_package.object_info.exp_pos_idx(data_sn,2) + ptycho_package.probe_info.pos_correct_pixel(data_sn,2);
        row_start_idx = exp_cen_row_idx - (ptycho_package.init_cond.effective_clip_size - 1)/2;
        row_end_idx = exp_cen_row_idx + (ptycho_package.init_cond.effective_clip_size - 1)/2;
        col_start_idx = exp_cen_col_idx - (ptycho_package.init_cond.effective_clip_size - 1)/2;
        col_end_idx = exp_cen_col_idx + (ptycho_package.init_cond.effective_clip_size - 1)/2;   
        gpu_clip_object = gpu_object(row_start_idx:row_end_idx,col_start_idx:col_end_idx);

        % formula (3) S(12)
        psi = gpu_probe .* gpu_clip_object;
        % psi(:,:,k)
        % k for probe(Mp)
        % psi for real space in (S12)
        
        % formula (4)
        Psi = fftshift(fft2(ifftshift(psi))); % fixed 20220809
        % Psi(:,:,k)        
        
        % formula (5), (S11)
        Psi_amp_flat = sqrt(sum(abs(Psi).^2,3)); % flat matrix.
        % multi-probe results are non-interference. using sqrt(|a|^2 + |b|^2 + ....)
        Psi_amp_flat_non_zero_mask = Psi_amp_flat ~= 0;
        Psi_amp_flat(~Psi_amp_flat_non_zero_mask) = 1E30; % invoid the NaN resutls in Psi_amp_flat equals to 0;
        Psi_p = Psi.*mask + data .* Psi./Psi_amp_flat .* ~mask.*single(Psi_amp_flat_non_zero_mask);
        % calculate chi^2
        %gpu_chi2_temp(1,data_sn) = sum(sum( (Psi_amp_flat -data).^2.*~mask))/active_area;
        clear Psi Psi_amp_flat_non_zero_mask Psi_amp_flat data
        psi_p = fftshift(ifft2(ifftshift(Psi_p))); % fixed 20220809
        diff_psi_p_psi = psi_p - psi;
        clear psi_p psi
        % Psi_amp_flat, Psi_p and psi_p are the matrix in [clip_size,clip_size,Mp]
        
        switch core
            case 'rPIE'
                % update object
                % rPIE
                PMax = max(sum(abs(gpu_probe).^2,3),[],'all');
                upper_term = sum(conj(gpu_probe).* diff_psi_p_psi,3);
                lower_term = (1-alpha)*sum(abs(gpu_probe).^2,3) + alpha*PMax;
                gpu_object(row_start_idx:row_end_idx,col_start_idx:col_end_idx) = gpu_clip_object + upper_term./lower_term;
                %update probe
                % formula (7) (S22)
                oMax = max( sum(abs(gpu_clip_object).^2,3),[],'all');
                upper_term = conj(gpu_clip_object).* diff_psi_p_psi;
                lower_term = (1-beta)*abs(gpu_clip_object).^2 + beta*oMax;
                if beta ~= 0
                    gpu_probe = gpu_probe + upper_term./lower_term;
                end
            case 'ePIE'
                % update object
                % formula (6) (S21)
                first_term = alpha / max(max( sum(abs(gpu_probe).^2,3)));
                second_term = sum(conj(gpu_probe).* diff_psi_p_psi,3); 
                gpu_object(row_start_idx:row_end_idx,col_start_idx:col_end_idx) = gpu_clip_object + first_term * second_term;

                %update probe
                % formula (7) (S22)
                first_term = beta / max(max( sum(abs(gpu_clip_object).^2,3)));
                second_term = conj(gpu_clip_object).* diff_psi_p_psi;
                gpu_probe = gpu_probe + first_term * second_term;
        end
                
        
        clear first_term second_term diff_psi_p_psi

    end
    
    % gather variables from gpu
    ptycho_package.object_info.real_space = gather(gpu_object);
    ptycho_package.probe_info.real_space = gather(gpu_probe);

    %chi2_sum = gather(sum(gpu_chi2_temp));    
end