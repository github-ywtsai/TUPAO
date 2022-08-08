function ptyshow()
    current_ff = pwd;
    [fn, fp] = uigetfile([current_ff '/*.mat']);
    section_fp = [fp fn];
    if fn == 0
        return
    end
    temp = load(section_fp,'object_info');
    object_info = temp.object_info;
    temp = load(section_fp,'probe_info');
    probe_info = temp.probe_info;
    temp = load(section_fp,'iteration_para');
    iteration_para = temp.iteration_para;
    temp = load(section_fp,'init_cond');
    init_cond = temp.init_cond;

    % main figure
    handles.main_fig = figure('Units','Pixels','Position',[100 100 1280 720],'NumberTitle','OFF',...
        'Name',section_fp);
    handles.pix_res_label = uicontrol(handles.main_fig,'style','text','Units','Pixels',...
        'Position',[10 690 200 20],'String',sprintf('Pixel res.: %.1f [nm]',object_info.x_res*1E9));
    handles.tabgroup = uitabgroup(handles.main_fig,'Units','Pixels','Position',[10 10 1260 680]);
    handles.tab_obj= uitab(handles.tabgroup,'Title','Object'); 
    handles.tab_probe_amp = uitab(handles.tabgroup,'Title','Probe (amplitude)');
    handles.tab_probe_int = uitab(handles.tabgroup,'Title','Probe (intensity)');
    handles.tab_PC = uitab(handles.tabgroup,'Title','Position Correction');
    handles.tab_chi2 = uitab(handles.tabgroup,'Title','Chi^2');
    
    % Object part
    handles.slider_obj_real_H = uicontrol(handles.tab_obj,'Units','Pixels','Position',[50 150 350 20],'Style','Slider','Callback',@(src,eventdata)change_show_range(src,eventdata,'real_H'));
    handles.slider_obj_real_L = uicontrol(handles.tab_obj,'Units','Pixels','Position',[50 100 350 20],'Style','Slider','Callback',@(src,eventdata)change_show_range(src,eventdata,'real_L'));
    handles.slider_obj_abs_H = uicontrol(handles.tab_obj,'Units','Pixels','Position',[470 150 350 20],'Style','Slider','Callback',@(src,eventdata)change_show_range(src,eventdata,'abs_H'));
    handles.slider_obj_abs_L = uicontrol(handles.tab_obj,'Units','Pixels','Position',[470 100 350 20],'Style','Slider','Callback',@(src,eventdata)change_show_range(src,eventdata,'abs_L'));
    handles.slider_obj_phase_H = uicontrol(handles.tab_obj,'Units','Pixels','Position',[890 150 350 20],'Style','Slider','Callback',@(src,eventdata)change_show_range(src,eventdata,'phase_H'));
    handles.slider_obj_phase_L = uicontrol(handles.tab_obj,'Units','Pixels','Position',[890 100 350 20],'Style','Slider','Callback',@(src,eventdata)change_show_range(src,eventdata,'phase_L'));
    handles.axes_obj_real = axes(handles.tab_obj,'Units','Pixels','Position',[50 250 350 350]);
    handles.axes_obj_abs = axes(handles.tab_obj,'Units','Pixels','Position',[470 250 350 350]);
    handles.axes_obj_phase = axes(handles.tab_obj,'Units','Pixels','Position',[890 250 350 350]);
    
    handles.popupmenu_obj_real_colormap = uicontrol(handles.tab_obj,'units','pixels','position',[50 200 150 20],'style','popupmenu',...
        'String',{'parula','jet','gray','bone','i-parula','i-jet','i-gray','i-bone'},'Value',1,'Callback',{@change_colormap,'obj_real'});
    handles.popupmenu_obj_abs_colormap = uicontrol(handles.tab_obj,'units','pixels','position',[470 200 150 20],'style','popupmenu',...
        'String',{'parula','jet','gray','bone','i-parula','i-jet','i-gray','i-bone'},'Value',1,'Callback',{@change_colormap,'obj_abs'});
    handles.popupmenu_obj_phase_colormap = uicontrol(handles.tab_obj,'units','pixels','position',[890 200 150 20],'style','popupmenu',...
        'String',{'parula','jet','gray','bone','i-parula','i-jet','i-gray','i-bone'},'Value',1,'Callback',{@change_colormap,'obj_phase'});
    
    function change_colormap(src,eventdata,target)
        cmap = src.String{src.Value};
        if strcmp(cmap(1:2),'i-')
            cmap = eval(sprintf('flipud(%s)',cmap(3:end)));
        else
            cmap = eval(cmap);
        end
        switch target
            case 'obj_real'
                colormap(handles.axes_obj_real,cmap);
            case 'obj_abs'
                colormap(handles.axes_obj_abs,cmap);
            case 'obj_phase'
                colormap(handles.axes_obj_phase,cmap);
        end
    end
    
    
    exp_pos_idx = object_info.exp_pos_idx + probe_info.pos_correct_pixel;
    object_plot_min_row_idx = min(exp_pos_idx(:,1));
    object_plot_max_row_idx = max(exp_pos_idx(:,1));
    object_plot_min_col_idx = min(exp_pos_idx(:,2));
    object_plot_max_col_idx = max(exp_pos_idx(:,2));
    
    obj_real_Ilimit = [min(real(object_info.real_space),[],'all') max(real(object_info.real_space),[],'all')];
    obj_abs_Ilimit = [min(abs(object_info.real_space),[],'all') max(abs(object_info.real_space),[],'all')];
    obj_phase_Ilimit = [min(angle(object_info.real_space),[],'all') max(angle(object_info.real_space),[],'all')];
    
    %[obj_row_size,obj_col_size] = size(object_info.real_space);
    %obj_xaxis = ((1:obj_col_size) - round(obj_col_size/2)) * probe_info.x_res * 1E6; % in [um]
    %obj_yaxis = ((1:obj_row_size) - round(obj_row_size/2)) * probe_info.y_res * 1E6; % in [um]
    obj_xaxis = object_info.real_space_xaxis*1E6; % in [um]
    obj_yaxis = object_info.real_space_yaxis*1E6; % in [um]
    % fixed in 2021/02/17
    % the direction in yaxis is that the postive direction point to the row
    % 1 (meaning the inverse direction of the row direction)
    x_lim_scaned = sort([obj_xaxis(object_plot_min_col_idx) obj_xaxis(object_plot_max_col_idx)]); % in [um]
    y_lim_scaned = sort([obj_yaxis(object_plot_min_row_idx) obj_xaxis(object_plot_max_row_idx)]); % in [um]
    
    handles.plot_real = imagesc(handles.axes_obj_real,obj_xaxis,obj_yaxis,real(object_info.real_space));
    handles.plot_abs = imagesc(handles.axes_obj_abs,obj_xaxis,obj_yaxis,abs(object_info.real_space));
    handles.plot_phase = imagesc(handles.axes_obj_phase,obj_xaxis,obj_yaxis,angle(object_info.real_space));
    handles.axes_obj_real.YAxis.Direction = 'normal';
    handles.axes_obj_abs.YAxis.Direction = 'normal';
    handles.axes_obj_phase.YAxis.Direction = 'normal';
    set(handles.axes_obj_real,{'Xlim','Ylim'},{x_lim_scaned,y_lim_scaned});
    set(handles.axes_obj_abs,{'Xlim','Ylim'},{x_lim_scaned,y_lim_scaned});
    set(handles.axes_obj_phase,{'Xlim','Ylim'},{x_lim_scaned,y_lim_scaned});
    xlabel(handles.axes_obj_real,'\mum'); ylabel(handles.axes_obj_real,'\mum');
    xlabel(handles.axes_obj_abs,'\mum'); ylabel(handles.axes_obj_abs,'\mum');
    xlabel(handles.axes_obj_phase,'\mum'); ylabel(handles.axes_obj_phase,'\mum');
    title(handles.axes_obj_real,'Object (real)');
    title(handles.axes_obj_abs,'Object (abs.)');
    title(handles.axes_obj_phase,'Object (phase)');
    axis(handles.axes_obj_real,'image');
    axis(handles.axes_obj_abs,'image');
    axis(handles.axes_obj_phase,'image');

    linkaxes([handles.axes_obj_real,handles.axes_obj_abs,handles.axes_obj_phase],'xy');
    
    set(handles.slider_obj_real_H,{'min','max','value'},{obj_real_Ilimit(1),obj_real_Ilimit(2),obj_real_Ilimit(2)},'SlidersTep',[0.01 0.1]);
    set(handles.slider_obj_real_L,{'min','max','value'},{obj_real_Ilimit(1),obj_real_Ilimit(2),obj_real_Ilimit(1)},'SlidersTep',[0.01 0.1]);
    set(handles.slider_obj_abs_H,{'min','max','value'},{obj_abs_Ilimit(1),obj_abs_Ilimit(2),obj_abs_Ilimit(2)},'SlidersTep',[0.01 0.1]);
    set(handles.slider_obj_abs_L,{'min','max','value'},{obj_abs_Ilimit(1),obj_abs_Ilimit(2),obj_abs_Ilimit(1)},'SlidersTep',[0.01 0.1]);
    set(handles.slider_obj_phase_H,{'min','max','value'},{obj_phase_Ilimit(1),obj_phase_Ilimit(2),obj_phase_Ilimit(2)},'SlidersTep',[0.01 0.1]);
    set(handles.slider_obj_phase_L,{'min','max','value'},{obj_phase_Ilimit(1),obj_phase_Ilimit(2),obj_phase_Ilimit(1)},'SlidersTep',[0.01 0.1]);
    
    handles.listener_obj_real_H = listener(handles.slider_obj_real_H,'Value','PreSet',@(src,eventdata)change_show_range(src,eventdata,'real_H'));
    handles.listener_obj_real_L = listener(handles.slider_obj_real_L,'Value','PreSet',@(src,eventdata)change_show_range(src,eventdata,'real_L'));
    handles.listener_obj_abs_H = listener(handles.slider_obj_abs_H,'Value','PreSet',@(src,eventdata)change_show_range(src,eventdata,'abs_H'));
    handles.listener_obj_abs_L = listener(handles.slider_obj_abs_L,'Value','PreSet',@(src,eventdata)change_show_range(src,eventdata,'abs_L'));
    handles.listener_obj_phase_H = listener(handles.slider_obj_phase_H,'Value','PreSet',@(src,eventdata)change_show_range(src,eventdata,'phase_H'));
    handles.listener_obj_phase_L = listener(handles.slider_obj_phase_L,'Value','PreSet',@(src,eventdata)change_show_range(src,eventdata,'phase_L'));
    
    function change_show_range(src,eventdata,target)
        switch target
            case 'real_H'
                CLim = get(handles.axes_obj_real,'CLim');
                CLim(2) = get(handles.slider_obj_real_H,'Value');
                set(handles.axes_obj_real,'CLim',CLim);
            case 'real_L'
                CLim = get(handles.axes_obj_real,'CLim');
                CLim(1) = get(handles.slider_obj_real_L,'Value');
                set(handles.axes_obj_real,'CLim',CLim);
            case 'abs_H'
                CLim = get(handles.axes_obj_abs,'CLim');
                CLim(2) = get(handles.slider_obj_abs_H,'Value');
                set(handles.axes_obj_abs,'CLim',CLim);
            case 'abs_L'
                CLim = get(handles.axes_obj_abs,'CLim');
                CLim(1) = get(handles.slider_obj_abs_L,'Value');
                set(handles.axes_obj_abs,'CLim',CLim);
            case 'phase_H'
                CLim = get(handles.axes_obj_phase,'CLim');
                CLim(2) = get(handles.slider_obj_phase_H,'Value');
                set(handles.axes_obj_phase,'CLim',CLim);
            case 'phase_L'
                CLim = get(handles.axes_obj_phase,'CLim');
                CLim(1) = get(handles.slider_obj_phase_L,'Value');
                set(handles.axes_obj_phase,'CLim',CLim);
        end
    end

    % probe amplitude part ratio
    NProbe = probe_info.Mp;
    probe_mode_amp = zeros(NProbe,1);
    for probe_sn = 1:NProbe
        probe_mode_amp(probe_sn,1) = sum(abs(probe_info.real_space(:,:,probe_sn)),'all');
    end
    probe_mode_sum_up = sum(probe_mode_amp);
    probe_mode_amp_ratio = probe_mode_amp / probe_mode_sum_up;
    
    % probe intensity part ratio
    probe_mode_int = zeros(NProbe,1);
    for probe_sn = 1:NProbe
        probe_mode_int(probe_sn,1) = sum(abs(probe_info.real_space(:,:,probe_sn)).^2,'all');
    end
    probe_mode_sum_up = sum(probe_mode_int,'all');
    probe_mode_int_ratio = probe_mode_int / probe_mode_sum_up;
    
    TGrapPosition{1} = [50 360 250 250];
    TGrapPosition{2} = [470 360 250 250];
    TGrapPosition{3} = [890 360 250 250];
    BGrapPosition{1} = [50 50 250 250];
    BGrapPosition{2} = [470 50 250 250];
    BGrapPosition{3} = [890 50 250 250];
 
    x_axis = probe_info.real_space_xaxis*1E6; % in [um]
    y_axis = probe_info.real_space_yaxis*1E6; % in [um]
    x_axis_lim = [min(x_axis) max(x_axis)];
    y_axis_lim = [min(y_axis) max(y_axis)];
    
    % plot probes
    for ProbeSN = 1:min(NProbe,3)
        % create axes for plot
        handles.axes_probe_amp{ProbeSN} = axes(handles.tab_probe_amp,'Units','Pixels','Position',TGrapPosition{ProbeSN},'DataAspectRatio',[1 1 1]);
        handles.axes_probe_phase{ProbeSN} = axes(handles.tab_probe_amp,'Units','Pixels','Position',BGrapPosition{ProbeSN},'DataAspectRatio',[1 1 1]);
        handles.axes_probe_int{ProbeSN} = axes(handles.tab_probe_int,'Units','Pixels','Position',TGrapPosition{ProbeSN},'DataAspectRatio',[1 1 1]);
        handles.axes_probe_int_phase{ProbeSN} = axes(handles.tab_probe_int,'Units','Pixels','Position',BGrapPosition{ProbeSN},'DataAspectRatio',[1 1 1]);
        
        % for phase sheet
        imagesc(handles.axes_probe_amp{ProbeSN},x_axis,y_axis,abs(probe_info.real_space(:,:,ProbeSN)));
        handles.axes_probe_amp{ProbeSN}.YAxis.Direction = 'normal';
        handles.axes_probe_amp{ProbeSN}.XLabel.String = '\mum';
        handles.axes_probe_amp{ProbeSN}.YLabel.String = '\mum';
        handles.axes_probe_amp{ProbeSN}.XLim = x_axis_lim;
        handles.axes_probe_amp{ProbeSN}.YLim = y_axis_lim;
        
        imagesc(handles.axes_probe_phase{ProbeSN},x_axis,y_axis,angle(probe_info.real_space(:,:,ProbeSN)));
        handles.axes_probe_phase{ProbeSN}.YAxis.Direction = 'normal';
        handles.axes_probe_phase{ProbeSN}.XLabel.String = '\mum';
        handles.axes_probe_phase{ProbeSN}.YLabel.String = '\mum';
        handles.axes_probe_phase{ProbeSN}.XLim = x_axis_lim;
        handles.axes_probe_phase{ProbeSN}.YLim = y_axis_lim;
        
        title(handles.axes_probe_amp{ProbeSN},sprintf('Probe %d %1.0f%%',ProbeSN,probe_mode_amp_ratio(ProbeSN)*100));
        
        % for intensity sheet
        imagesc(handles.axes_probe_int{ProbeSN},x_axis,y_axis,abs(probe_info.real_space(:,:,ProbeSN)).^2);
        handles.axes_probe_int{ProbeSN}.YAxis.Direction = 'normal';
        handles.axes_probe_int{ProbeSN}.XLabel.String = '\mum';
        handles.axes_probe_int{ProbeSN}.YLabel.String = '\mum';
        handles.axes_probe_int{ProbeSN}.XLim = x_axis_lim;
        handles.axes_probe_int{ProbeSN}.YLim = y_axis_lim;
        
        imagesc(handles.axes_probe_int_phase{ProbeSN},x_axis,y_axis,angle(probe_info.real_space(:,:,ProbeSN)));
        handles.axes_probe_int_phase{ProbeSN}.YAxis.Direction = 'normal';
        handles.axes_probe_int_phase{ProbeSN}.XLabel.String = '\mum';
        handles.axes_probe_int_phase{ProbeSN}.YLabel.String = '\mum';
        handles.axes_probe_int_phase{ProbeSN}.XLim = x_axis_lim;
        handles.axes_probe_int_phase{ProbeSN}.YLim = y_axis_lim;
        
        title(handles.axes_probe_int{ProbeSN},sprintf('Probe %d %1.0f%%',ProbeSN,probe_mode_int_ratio(ProbeSN)*100));
        
    end
        
    
    % plot position correction results
    handles.axes_PC = axes(handles.tab_PC,'Units','Pixels','Position',[100 100 400 400],'DataAspectRatio',[1 1 1]);
    imagesc( handles.axes_PC,angle(object_info.real_space));
    hold on
    q = quiver(handles.axes_PC,object_info.exp_pos_idx(:,2),object_info.exp_pos_idx(:,1),probe_info.pos_correct_pixel(:,2),probe_info.pos_correct_pixel(:,1));
    q.Color = 'Y';
    hold off
    ax = gca;
    ax.DataAspectRatio = [1 1 1];
    ax.XLim = [min(object_info.exp_pos_idx(:,2)) max(object_info.exp_pos_idx(:,2))];
    ax.YLim = [min(object_info.exp_pos_idx(:,1)) max(object_info.exp_pos_idx(:,1))];
    
    handles.PC_info = uicontrol(handles.tab_PC,'Style','Text','Units','Pixels','Position', [550 100 200 200]);
    message = {sprintf('Max. x shift: %.1f nm',(max(probe_info.pos_correct_pixel(:,2))*probe_info.x_res*1E9)),sprintf('Max. y shift: %.1f nm',(max(probe_info.pos_correct_pixel(:,1)*probe_info.y_res)*1E9))};
    set(handles.PC_info,'String',message)
    
    % chi^2
    handles.axes_chi2 = axes(handles.tab_chi2,'Units','Pixels','Position',[100 100 400 400],'DataAspectRatio',[1 1 1]);
    if ~isfield(iteration_para,'FinishedRun')
        iteration_para.FinishedRun = iteration_para.run;
    end
    plot(handles.axes_chi2,iteration_para.chi2(1:iteration_para.FinishedRun));
    handles.axes_chi2.YScale = 'log';
    handles.axes_chi2.XScale = 'log';
   
end