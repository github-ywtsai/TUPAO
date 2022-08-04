function object_info = gen_object_info(init_cond)
    rng(init_cond.rand_seed + 1)
    effective_clip_size = init_cond.effective_clip_size;
    x_res = init_cond.wavelength * init_cond.detector_distance / effective_clip_size / init_cond. effective_x_pixel_size;
    y_res = init_cond.wavelength * init_cond.detector_distance / effective_clip_size / init_cond. effective_y_pixel_size;
    object_info.x_res = x_res;
    object_info.y_res = y_res;
    
    [temp,~] = max(init_cond.exp_pos);
    exp_pos_ymax = temp(1);
    exp_pos_xmax = temp(2);
    [temp,~] = min(init_cond.exp_pos);
    exp_pos_ymin = temp(1);
    exp_pos_xmin = temp(2);
    % exp_pos_ycen = (exp_pos_ymax + exp_pos_ymin)/2;
    % exp_pos_xcen = (exp_pos_xmax + exp_pos_xmin)/2;

    obj_extend_factor = 0.1; % ex: 0.1 meaning 10% probe_area extented
    obj_left_boundary   = exp_pos_xmin - effective_clip_size/2 * x_res - effective_clip_size * x_res * obj_extend_factor;
    obj_right_boundary  = exp_pos_xmax + effective_clip_size/2 * x_res + effective_clip_size * x_res * obj_extend_factor;
    obj_bottom_boundary = exp_pos_ymin - effective_clip_size/2 * y_res - effective_clip_size * y_res * obj_extend_factor;
    obj_top_boundary    = exp_pos_ymax + effective_clip_size/2 * y_res + effective_clip_size * y_res * obj_extend_factor;

    object_row_size = round (( obj_top_boundary - obj_bottom_boundary ) / y_res);
    object_col_size = round (( obj_right_boundary - obj_left_boundary ) / x_res);

    if mod(object_row_size,2) == 0
        object_row_size = object_row_size + 1;
    end
    if mod(object_col_size,2) == 0
        object_col_size = object_col_size + 1;
    end

    object_cen_row_idx = round(object_row_size/2); % row and col index
    object_cen_col_idx = round(object_col_size/2); % row and col index

    object_info.real_space = zeros(object_row_size,object_col_size);

    [object_row_idx , object_col_idx , ~] = find(object_info.real_space == 0);
    object_row_idx = reshape(object_row_idx,object_row_size,object_col_size);
    object_col_idx = reshape(object_col_idx,object_row_size,object_col_size);

    object_info.real_space_xaxis = (object_col_idx(1,:) - object_cen_col_idx) * x_res;
    object_info.real_space_yaxis = -(object_row_idx(:,1) - object_cen_row_idx) * y_res;
    % 20180115 fix the y axis from
    % guessed_object.real_space_yaxis = (object_row_idx(:,1) - object_cen_row_idx) * y_res;

    object_info = pre_assign_object_intensity_flat(init_cond, object_info,0.5);
end


function object_info = pre_assign_object_intensity_flat(init_cond, object_info,value)
    n_of_data = init_cond.n_of_data;
    exp_pos = init_cond.exp_pos;
    real_space_xaxis = object_info.real_space_xaxis;
    real_space_yaxis = object_info.real_space_yaxis;
    real_space = object_info.real_space;
    [real_space_row_size,real_space_col_size] = size(real_space);

    exp_pos_idx = zeros(n_of_data,2);

    for i = 1 : n_of_data

        exp_pos_x = exp_pos(i,2);
        exp_pos_y = exp_pos(i,1);
        [~,target_col_idx] = min(abs(real_space_xaxis - exp_pos_x));
        [~,target_row_idx] = min(abs(real_space_yaxis - exp_pos_y));

        exp_pos_idx(i,:) = [target_row_idx ,target_col_idx];

    end
    
    rng(init_cond.rand_seed + 1)
    BaseObj = ones(real_space_row_size,real_space_col_size)*value;
    RandObj = (rand(real_space_row_size,real_space_col_size)-0.5)*value*0.1;
    RandPhase = exp(1i*rand(real_space_row_size,real_space_col_size)*pi); % phase 0 ~ pi, make image part alway > 0;
    
    object_info.real_space = abs(single(BaseObj + RandObj)).*RandPhase;
    object_info.exp_pos_idx = single(exp_pos_idx);


end

