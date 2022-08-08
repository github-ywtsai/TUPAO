function output = propagate_probe(z,probe_info)
    % 20220727
    % fix wrong useing of fft2 cause in incorrect phase information
    % worng usage: fft2(A)
    % corrected usage: fft2(ifftshift(A))
    % try in 20190909
    % z and propagating direction:
    % in formula, +z means propagating to "downstream" and -z to
    % "upstream". So, when input, +z means "downstream" and to increase the distance from the
    % source coming from the undulator. -z means "upstream" and to decrease
    % the distace from the source coming from the undulator.
    % Different mode between fft mode and ifft mode.

    % preset parameters
    lambda = probe_info.wavelength; %[m]
    k = 2*pi/lambda;
    ModeNum = size(probe_info.real_space,3); 
    
    xi_axis =  probe_info.real_space_xaxis;
    eta_axis =  probe_info.real_space_yaxis';
    xi_axis_res = abs(xi_axis(2)-xi_axis(1));
    eta_axis_res = abs(eta_axis(2)-eta_axis(1));
    [xi,eta] = meshgrid(xi_axis,eta_axis); %[m]
    
    [U_measured_eta_size,U_measured_xi_size] = size(probe_info.real_space(:,:,1));
    xp_res = abs(lambda*z/U_measured_xi_size/xi_axis_res);
    yp_res = abs(lambda*z/U_measured_eta_size/eta_axis_res);
    [U_propagated_yp_size,U_propagated_xp_size] = size(probe_info.real_space(:,:,1));
    xp_axis = ((1:U_propagated_xp_size)-round(U_propagated_xp_size/2))*xp_res;
    yp_axis = ((1:U_propagated_yp_size)-round(U_propagated_yp_size/2))*yp_res;
    [xp,yp] = meshgrid(xp_axis,yp_axis);
    
    U_propagated = zeros(U_propagated_yp_size,U_propagated_xp_size,ModeNum);
    
    % calculate propagating
    for ModeSN = 1:ModeNum
        U_measured = probe_info.real_space(:,:,ModeSN);
        if z<0
            U_measured = rot90(U_measured,2);
        end
        
        temp = U_measured.*exp(1i*k/(2*z)*(xi.^2 +eta.^2));
        %fft_temp = fftshift(fft2(temp))*(1/xiSize)*(1/etaSize);
        % Above formula is formula without the area normazilation factor.
        % Used before 20201102
        fft_temp = fftshift(fft2(ifftshift(temp)))*xi_axis_res*eta_axis_res;
        % Above formula is formula with the area normazilation factor.
        % Used after 20201102
        U_propagated(:,:,ModeSN) = exp(1i*k*z)/(1i*lambda*z)*exp(1i*k/(2*z)*(lambda*z)^2*(xp.^2+yp.^2)).*fft_temp;
    end
    
    output.propagated_probe = U_propagated;
    output.x_axis = xp_axis;
    output.y_axis = yp_axis;
    output.x_res = xp_res;
    output.y_res = yp_res;
    output.propagating_distance = z; 

end