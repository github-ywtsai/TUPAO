function output_package = flat_background(ptycho_package,reference_pos,mode)
% reference_pos: the empty position on image noted by [x,y] in meter
% example:
% [x1,y1;
%  x2,y2;
%  ....
% ]
% mode are 'phase', 'amplitude' or 'both'
num_reference_pos = size(reference_pos,1);
img = ptycho_package.object_info.real_space;
x_axis = ptycho_package.object_info.real_space_xaxis;
y_axis = ptycho_package.object_info.real_space_yaxis;
[X,Y] = meshgrid(x_axis,y_axis);

x = zeros(num_reference_pos,1);
y = zeros(num_reference_pos,1);
phase_z = zeros(num_reference_pos,1);
amplitude_z = zeros(num_reference_pos,1);
for ii = 1:num_reference_pos
    [closest_row_idx,closest_col_idx] = find_closest_idx(x_axis,y_axis,reference_pos(ii,1),reference_pos(ii,2));
    x(ii) = x_axis(closest_col_idx);
    y(ii) = y_axis(closest_row_idx);
    phase_z(ii) = angle(img(closest_row_idx,closest_col_idx));
    amplitude_z(ii) = abs(img(closest_row_idx,closest_col_idx));
end

% fitting the phase surface
A = [x,y,ones(length(x),1)];
phase_params = A\phase_z;
amplitude_params = A\amplitude_z;
phase_backgrpund = phase_params(1)*X + phase_params(2)*Y + phase_params(3);
amplitude_backgrpund = amplitude_params(1)*X + amplitude_params(2)*Y + amplitude_params(3);
if strcmp(mode , 'both')
        back_ground = amplitude_backgrpund .* exp(1i*phase_backgrpund);
elseif strcmp(mode , 'phase')
        back_ground = exp(1i*phase_backgrpund);
elseif strcmp(mode , 'amplitude')
        back_ground = amplitude_backgrpund;
end

output_package = ptycho_package;
output_package.object_info.real_space = output_package.object_info.real_space./back_ground;



function [closest_row_idx,closest_col_idx] = find_closest_idx(x_axis,y_axis,x,y)
% x_axis and y_axis are the axes for the image
% x,y are the target position
[closest_y_value, closest_row_idx] = min(abs(y-y_axis));
[closest_x_value, closest_col_idx] = min(abs(x-x_axis));

