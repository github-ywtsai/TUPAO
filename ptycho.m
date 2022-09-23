function ptycho(measured_amp,init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para,sectionfile_info)    
    
    %% prepare plot axes
    if iteration_para.draw_results
        fig = figure('position',[100,100,1300,500]);
        axes_obj = axes(fig,'position',[0.0500 0.0500 0.25 0.9000]);
        axes_probe =  axes(fig,'position',[0.3500 0.0500 0.25 0.9000]);
        axes_chi2 =  axes(fig,'position',[0.6500 0.0500 0.25 0.9000]);
    end
    SectionFilePrefix = sprintf('SN%d',sectionfile_info.SN);
    
    %% arm GPU
    if init_cond.using_GPU
    fprintf('Arrange GPU %d...',init_cond.using_GPU);
    gpuDevice(init_cond.using_GPU);
    measured_amp = gpuArray(measured_amp);
    object_info.real_space = gpuArray(object_info.real_space);
    probe_info.real_space = gpuArray(probe_info.real_space);
    probe_info.ProbeConf.upstream_ROI = gpuArray(probe_info.ProbeConf.upstream_ROI);
    measurement_info.individual_mask = gpuArray(measurement_info.individual_mask);
    measurement_info.individual_mask_active_area = gpuArray(measurement_info.individual_mask_active_area);
    iteration_para.chi2 = gpuArray(iteration_para.chi2);
    fprintf('\tDone.\n');
    end

    %% iteration part    
    while iteration_para.FinishedRun + 1 <= iteration_para.max_iteration_num
        CurrentRun = iteration_para.FinishedRun + 1;
        
        %% check beta and position correction and itnesity normlization start point
        if CurrentRun >= iteration_para.beta_start_pt
            iteration_para.beta_current = iteration_para.beta;
        else
             iteration_para.beta_current = 0;
        end
        
        if mod(CurrentRun-iteration_para.real_space_constraint_start_pt,iteration_para.real_space_constraint_period) == 0 && CurrentRun >= iteration_para.real_space_constraint_start_pt && iteration_para.real_space_constraint_start_pt ~= 0
            iteration_para.real_space_constraint_status = 1;
        else
            iteration_para.real_space_constraint_status = 0;
        end

        if mod(CurrentRun-iteration_para.pos_corr_start_pt,iteration_para.pos_corr_period) == 0 && CurrentRun >= iteration_para.pos_corr_start_pt && iteration_para.pos_corr_start_pt ~= 0
            iteration_para.pos_correction_status = 1;
        else
            iteration_para.pos_correction_status = 0;
        end
        
        %% check beta and position correction and itnesity normlization start point
        
        
        %PC parts
        if iteration_para.pos_correction_status == 1
            tic;
            fprintf('%s_Run %d: position correcting...',SectionFilePrefix,CurrentRun);
            updated_pos_correct_pixel = position_correction(measured_amp,init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para);
            probe_info.pos_correct_pixel = updated_pos_correct_pixel;
            ElapsedT = toc;
            fprintf('\tDone(%.1f sec).\n',ElapsedT);
        end
        

        %% calculating parts
        fprintf('%s_Run %d in progressing...',SectionFilePrefix,CurrentRun);
        tic;       
        % ePIE or rPIE
        if strcmpi(init_cond.core,'ePIE') 
            [updated_object,updated_probe,chi2_sum] = ePIE(measured_amp,init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para);
        elseif strcmpi(init_cond.core,'rPIE') 
            [updated_object,updated_probe,chi2_sum] = rPIE(measured_amp,init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para);
        end
        [~,n_of_interesting_data] = size(iteration_para.interesting_table);
        iteration_para.chi2(CurrentRun) = chi2_sum/n_of_interesting_data;
        object_info.real_space = updated_object;
        probe_info.real_space = updated_probe;
        iteration_para.FinishedRun = CurrentRun; % finished run
        
        % DM
        %{
        [updated_object,updated_probe,chi2_sum] = DM(measured_amp,init_cond,mask_info,measurement_info,object_info,probe_info,iteration_para);
        [~,n_of_interesting_data] = size(iteration_para.interesting_table);
        iteration_para.chi2(CurrentRun) = gather(chi2_sum)/n_of_interesting_data;
        object_info.real_space = updated_object;
        probe_info.real_space = updated_probe;
        iteration_para.FinishedRun = CurrentRun;
        %}
        
        ElapsedT = toc;
        
        %% calculating parts
        
        clear updated_object updated_probe
        fprintf('\tDone(%.1f sec).\n',ElapsedT);
        
        %% plot parts
        if iteration_para.draw_results
            imagesc(axes_obj,object_info.real_space_xaxis,object_info.real_space_yaxis,angle(object_info.real_space));axes_obj.DataAspectRatio = [1,1,1];
            imagesc(axes_probe,probe_info.real_space_xaxis,probe_info.real_space_yaxis,abs(probe_info.real_space(:,:,1)));axes_probe.DataAspectRatio = [1,1,1];
            loglog(axes_chi2, iteration_para.chi2);
            drawnow
        end
        
        %% save to section file
        if  and(mod(CurrentRun,iteration_para.saveing_section_file_peroid) == 0,iteration_para.saveing_section_file_peroid ~= 0)
            SectionFN = sprintf('%s_Run%d.mat',SectionFilePrefix,CurrentRun);
            SectionFP = fullfile(init_cond.results_path,SectionFN);
            fprintf('Saving section file %s...\t',SectionFN)
            % get data from GPU
            if init_cond.using_GPU
                object_info.real_space = gather(object_info.real_space);
                probe_info.real_space = gather(probe_info.real_space);
                probe_info.ProbeConf.upstream_ROI = gather(probe_info.ProbeConf.upstream_ROI);
                measurement_info.individual_mask = gather(measurement_info.individual_mask);
                measurement_info.individual_mask_active_area = gather(measurement_info.individual_mask_active_area);
                iteration_para.chi2 = gather(iteration_para.chi2);
            end
            save(SectionFP,'init_cond','mask_info','measurement_info','object_info','probe_info','iteration_para','sectionfile_info')
            % move data to GPU
            if init_cond.using_GPU
                object_info.real_space = gpuArray(object_info.real_space);
                probe_info.real_space = gpuArray(probe_info.real_space);
                probe_info.ProbeConf.upstream_ROI = gpuArray(probe_info.ProbeConf.upstream_ROI);
                measurement_info.individual_mask = gpuArray(measurement_info.individual_mask);
                measurement_info.individual_mask_active_area = gpuArray(measurement_info.individual_mask_active_area);
                iteration_para.chi2 = gpuArray(iteration_para.chi2);
            end
            fprintf('Done.\n')
        end
        
    end
    
    
    
end