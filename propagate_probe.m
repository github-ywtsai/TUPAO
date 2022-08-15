function [propagated_probe,propagated_x_axis,propagated_y_axis]  = propagate_probe(z,probe,wavelength,xi_axis,eta_axis)
    % xi_axis and et_axis are the probe_info.real_space_xaxis and probe_info.real_space_yaxis
    % wavelength in meter
    % z and propagating direction:
    % in formula, +z means propagating to "downstream" and -z to
    % "upstream". So, when input, +z means "downstream" and to increase the distance from the
    % source coming from the undulator. -z means "upstream" and to decrease

    % preset parameters
    lambda = wavelength; %[m]
    k = 2*pi/lambda;
    
    xi_axis_res = abs(xi_axis(2)-xi_axis(1));
    eta_axis_res = abs(eta_axis(2)-eta_axis(1));
    [xi,eta] = meshgrid(xi_axis,eta_axis); %[m]
    [U_measured_eta_size,U_measured_xi_size] = size(probe(:,:,1));
    xp_res = abs(lambda*z/U_measured_xi_size/xi_axis_res);
    yp_res = abs(lambda*z/U_measured_eta_size/eta_axis_res);
    [U_propagated_yp_size,U_propagated_xp_size] = size(probe(:,:,1));
    xp_axis = ((1:U_propagated_xp_size)-round(U_propagated_xp_size/2))*xp_res;
    yp_axis = ((1:U_propagated_yp_size)-round(U_propagated_yp_size/2))*yp_res;
    [xp,yp] = meshgrid(xp_axis,yp_axis);
    
    % data prepare for GPU case
    if strcmpi(class(probe),'gpuArray')
        xi = gpuArray(xi);
        eta = gpuArray(eta);
        xp = gpuArray(xp);
        yp = gpuArray(yp);
    end

    if z<0
        U_measured = rot90(probe,2);
    else
        U_measured = probe;
    end
    % calculate propagating

    temp = U_measured.*exp(1i*k/(2*z)*(xi.^2 +eta.^2));
    fft_temp = fftshift(fft2(ifftshift(temp)))*xi_axis_res*eta_axis_res;
    U_propagated = exp(1i*k*z)/(1i*lambda*z)*exp(1i*k/(2*z)*(xp.^2+yp.^2)).*fft_temp;

    propagated_probe = U_propagated;
    propagated_x_axis = xp_axis;
    propagated_y_axis = yp_axis;