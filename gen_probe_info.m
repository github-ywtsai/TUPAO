function probe_info = gen_probe_info(init_cond)
    ProbeConf = get_ProbeConf(init_cond);

    % generate empty probe   
    real_space = zeros(ProbeConf.clip_size,ProbeConf.clip_size,ProbeConf.Mp);
    
    switch ProbeConf.ProbGenMode
        case 'Gaussian'
            disp('Generating probe using Gaussian profile...')
            probe_temp = generate_probe_from_Gaussian(ProbeConf);
        case 'Zone Plate'
            disp('Generating probe using simulation of zone plate...')
            probe_temp = generate_probe_from_ZonePlate(ProbeConf);
        case 'Adapt from another result'
            disp('Generating probe from ptychography results...')
            probe_temp = generate_probe_from_sectionfile(init_cond,ProbeConf);
    end
    fprintf('\tDone.\n')
    %normalize
    Isum = sum(abs(probe_temp).^2,'all');
    ratio = sqrt(ProbeConf.photon_flux*init_cond.count_time/Isum);
    probe_temp = probe_temp*ratio;

    real_space(:,:,1) = probe_temp;

    if ProbeConf.Mp > 1
        main_probe_weight = 0.8;
        main_probe_amp_sum = sum(sum(real_space(:,:,1)));

        real_space(:,:,2:ProbeConf.Mp) = abs(rand(ProbeConf.clip_size,ProbeConf.clip_size,ProbeConf.Mp-1));
        other_probe_amp_sum = sum(sum(sum(real_space(:,:,2:ProbeConf.Mp),3)));
        factor = (1-main_probe_weight)/main_probe_weight*main_probe_amp_sum/other_probe_amp_sum;
        real_space(:,:,2:ProbeConf.Mp) = real_space(:,:,2:ProbeConf.Mp) * factor;        
    end
    
    probe_info.ProbeConf = ProbeConf;
    probe_info.real_space = single(real_space);
    probe_info.real_space_xaxis = single(ProbeConf.real_space_xaxis);
    probe_info.real_space_yaxis = single(ProbeConf.real_space_yaxis);
    
    if ProbeConf.AdapPosCorr == 0
        pos_correct_pixel = single(zeros(ProbeConf.n_of_data,2)); % in col and row, shift unit is pixel.
    else
        pos_correct_pixel = adapt_PosCorr_info(ProbeConf);
    end
    
    probe_info.pos_correct_pixel = pos_correct_pixel;
    probe_info.x_res = ProbeConf.x_res;
    probe_info.y_res = ProbeConf.y_res;
    probe_info.Mp = ProbeConf.Mp;
    probe_info.wavelength = init_cond.wavelength;

end

function RefPosCorrPixel = adapt_PosCorr_info(ProbeConf)
    if ~exist(ProbeConf.RefProbeFP,'File')
        return
    end
    Temp = load(ProbeConf.RefProbeFP,'probe_info');
    RefProbeInfo = Temp.probe_info;

    % pos_correct_pixel = [row shift 1, col shift 1; row shift 2, col shift 2; ...]

    ScanPts = ProbeConf.n_of_data;
    [RefScanPts,~] = size(RefProbeInfo.pos_correct_pixel);
    
    if ScanPts ~= RefScanPts
        return
    end
    
    x_res = ProbeConf.x_res;
    y_res = ProbeConf.y_res;

    Ref_x_res = RefProbeInfo.x_res;
    Ref_y_res = RefProbeInfo.y_res;

    RefPosCorrPixel = round([RefProbeInfo.pos_correct_pixel(:,1)*Ref_y_res/y_res, RefProbeInfo.pos_correct_pixel(:,2)*Ref_x_res/x_res]);
    
end

function probe_temp = generate_probe_from_Gaussian(ProbeConf)    
    amplitude = ProbeConf.photon_flux;
    HBSize = ProbeConf.GaussainBeamHSize;
    VBSize = ProbeConf.GaussainBeamVSize;
    probe_col_idx = ProbeConf.probe_col_idx;
    probe_row_idx = ProbeConf.probe_row_idx;
    real_space_col_cen_idx = ProbeConf.real_space_col_cen_idx;
    real_space_row_cen_idx = ProbeConf.real_space_row_cen_idx;
    x_res = ProbeConf.x_res;
    y_res = ProbeConf.y_res;
    clip_size = ProbeConf.clip_size;
    real_space_xaxis_matrix = (probe_col_idx - real_space_col_cen_idx) * x_res;
    real_space_yaxis_matrix = (probe_row_idx - real_space_row_cen_idx) * y_res;
    probe_temp = amplitude .* exp(-power(real_space_xaxis_matrix,2)/2/HBSize^2) .*  exp(-power(real_space_yaxis_matrix,2)/2/VBSize^2);
    
    if ProbeConf.GaussainBroken == 1
        SpeckleVSize = 0.1*VBSize; % [m]
        SpeckleHSize = 0.1*HBSize; % [m]
        Speckle = exp(-power(real_space_xaxis_matrix,2)/2/SpeckleHSize^2) .*  exp(-power(real_space_yaxis_matrix,2)/2/SpeckleVSize^2);
        SeedNum = round((x_res*clip_size)^2/SpeckleVSize/SpeckleHSize);
        Seed = randi([1 ProbeConf.clip_size^2],SeedNum,1);
        Temp = zeros(ProbeConf.clip_size);
        Temp(Seed) = 1;
        Temp = conv2(Temp,Speckle,'same');
        probe_temp = probe_temp.*Temp.*exp(1i*Temp); % break amp and
        % phase
        %probe_temp = probe_temp.*exp(1i*Temp); % only break phase
    end
end

function main_probe = generate_probe_from_ZonePlate(ProbeConf)
    % the min. pitch dr and the period N for the zone plate design
    % using formula from Attwood Soft X-rays and Extereme Ultraviolet Radiation
    %dr = 30e-9; N = 800; %  96 um zoneplate in TPS 25A
    %dr = 50e-9; N = 300; %  60 um zoneplate in TPS 25A
    dr = 70e-9; N = 286; %  60 um zoneplate in TPS 25A
    % diameter of the zone plate
    D = 4*N*dr;
    lambda = ProbeConf.wavelength;
    f = (4*N*dr^2)/lambda;
    
    D_cs = 35e-6; % diameter of the central stop [m]
    
    
    %% define material characterastic
    % auto calculating: not yet
    material = 'Au';
    thickness = 1500e-9; % [m]
    
    %% calculate refractivity
    % n(lambda) = 1 - delta(lambda) + i beta(lambda)
    f1 = 74.9382; % [e/atom] from NIST
    f2 = 6.0726; % [e/atom] % from NIST
    re = 2.8179403E-15; % classical electron rauius [m]
    Na = 6.02214129E23; % Avogadro constant [1/mol]
    Ma = 196.966569; % molar mass [g/mol]
    rho = 19.32; % density [g/cm^3]
    rho = rho / (1E-2)^3; % convert unit from [g cm-3] to [g m-3]
    na = rho*Na/Ma; % number density
    delta = na*re*lambda^2/2/pi*f1;
    beta =  na*re*lambda^2/2/pi*f2;
    n_refractivity = 1 - delta + 1i * beta;
    
    % E(z) = E0 * exp(i * n k z)
    % modulator = E(z)/E0 = exp(i * n k z)
    zp_modulation_factor = exp(1i * n_refractivity * 2*pi/lambda * thickness)/exp(1i * 1 * 2*pi/lambda * thickness);
    
    % central stop modulation facotr
    cs_modulation_factor = 0;
    
    %% create zone plate with central stop
    n = 0:N;
    rn = sqrt(n*lambda*f+(n*lambda).^2/4);
    
    % create zp
    pix_res = dr/2;
    range = D*1.5;
    pix_num = round(range/pix_res);
    if mod(pix_num,2)==0
        pix_num = pix_num + 1;
    end
    cen_idx = (pix_num+1)/2;
    zp = zeros(pix_num);
    x_axis = ((1:pix_num)-cen_idx)*pix_res;
    y_axis = x_axis*-1; % the direction of y is inverse of the row axis
    [x_matrix,y_matrix] = meshgrid(x_axis,y_axis);
    distance_map = sqrt(x_matrix.^2 + y_matrix.^2);
    
    for n_sn = flip(0:N)
        if mod(n_sn,2) == 0
            zp(distance_map<rn(n_sn+1)) = zp_modulation_factor;
        else
            zp(distance_map<rn(n_sn+1)) = 0;
        end
    end

    % create cs
    cs = ones(pix_num);
    cs(distance_map<D_cs/2) = cs_modulation_factor;
    
    exitwave = zp.*cs;
    
    %{
    %% gen probe old version, before 2025/07/10
    % gen on traget probe
    propagating_distance = f + ProbeConf.ZoneplateOffFocal;
    [probe_temp,probe_temp_x_axis,probe_temp_y_axis] = propagate_probe(propagating_distance,exitwave,lambda,x_axis,y_axis);
    
    % rescale on target probe
    probe_temp_res = abs(probe_temp_x_axis(2) - probe_temp_x_axis(1));
    rescale_factor = probe_temp_res/ProbeConf.x_res;
    rescale_pix_num = round(pix_num*rescale_factor);
    if mod(rescale_pix_num,2) == 0
        rescale_pix_num = rescale_pix_num + 1;
    end
    probe_temp_rescale = imresize(probe_temp,[rescale_pix_num,rescale_pix_num]);
    
    rescale_cen = (rescale_pix_num + 1)/2;
    extend_range = (ProbeConf.clip_size-1)/2;
    clip_range = rescale_cen-extend_range:rescale_cen+extend_range;
    main_probe = probe_temp_rescale(clip_range,clip_range);
    
    % blur probe
    blur_width = round(ProbeConf.clip_size/200);
    [row_ind,col_ind] = find(ones(ProbeConf.clip_size,ProbeConf.clip_size));
    cen_ind = (ProbeConf.clip_size+1)/2;
    row_ind = reshape(row_ind,ProbeConf.clip_size,ProbeConf.clip_size)-cen_ind;
    col_ind = reshape(col_ind,ProbeConf.clip_size,ProbeConf.clip_size)-cen_ind;
    blur_matrix = exp(-row_ind.^2/2/pi/blur_width^2) .* exp(-col_ind.^2/2/pi/blur_width^2);
    main_probe = conv2(main_probe,blur_matrix,'same');
    %}
    
    %% gen probe new version, after 2025/07/10
    propagating_distance = f + ProbeConf.ZoneplateOffFocal;
    main_probe = propagating_rescaling_wavefield(propagating_distance,exitwave,lambda,x_axis,y_axis,ProbeConf.x_res,ProbeConf.clip_size);
end

function propagted_and_matched_probe = propagating_rescaling_wavefield(propagating_distance,wavefield,lambda,x_axis,y_axis,target_res,target_size)
% This function is designed for the probe propagating during the generation
% of probes. The input is the the propagating distance and the wavefield.
% After the wavefield propageted, the resolution and field of view is
% changed. Then, using interp method to math the pixel size of the
% wavefiled to the target size and clip it to the purper size.
% The propagating_distance is the propagating distance of the wavefield in
% meter. The wavefield is the candidate for the propagating. The
% lambda,x_axis, y_axis are the related information of the wavefield. the
% target resolution and the target size are the pixel size and clip size of
% the propagated wavefield.
    [pix_num,~] = size(wavefield);
    % gen on traget probe
    [probe_temp,probe_temp_x_axis,probe_temp_y_axis] = propagate_probe(propagating_distance,wavefield,lambda,x_axis,y_axis);
    % rescale on target probe
    probe_temp_res = abs(probe_temp_x_axis(2) - probe_temp_x_axis(1));
    rescale_factor = probe_temp_res/target_res;
    rescale_pix_num = round(pix_num*rescale_factor);
    if mod(rescale_pix_num,2) == 0
        rescale_pix_num = rescale_pix_num + 1;
    end
    probe_temp_rescale = imresize(probe_temp,[rescale_pix_num,rescale_pix_num]);
    
    rescale_cen = (rescale_pix_num + 1)/2;
    extend_range = (target_size-1)/2;
    clip_range = rescale_cen-extend_range:rescale_cen+extend_range;
    main_probe = probe_temp_rescale(clip_range,clip_range);
    
    % blur probe
    blur_width = round(target_size/200);
    [row_ind,col_ind] = find(ones(target_size,target_size));
    cen_ind = (target_size+1)/2;
    row_ind = reshape(row_ind,target_size,target_size)-cen_ind;
    col_ind = reshape(col_ind,target_size,target_size)-cen_ind;
    blur_matrix = exp(-row_ind.^2/2/pi/blur_width^2) .* exp(-col_ind.^2/2/pi/blur_width^2);
    propagted_and_matched_probe = conv2(main_probe,blur_matrix,'same');
end

    

function probe_temp = generate_probe_from_ZonePlate_old(ProbeConf)    
    % using formula A = (1+cos(kr^2) )/2
    Beamsize = ProbeConf.GaussainBeamHSize; % [um]
    Broken = ProbeConf.GaussainBroken;
    r0 = Beamsize/2; % [um]
    r0 = r0*1E-6; % in [m]
    k = pi/2/r0^2;
    probe_col_idx = ProbeConf.probe_col_idx;
    probe_row_idx = ProbeConf.probe_row_idx;
    real_space_col_cen_idx = ProbeConf.real_space_col_cen_idx;
    real_space_row_cen_idx = ProbeConf.real_space_row_cen_idx;
    x_res = ProbeConf.x_res;
    y_res = ProbeConf.y_res;
    clip_size = ProbeConf.clip_size;
    real_space_xaxis_matrix = (probe_col_idx - real_space_col_cen_idx) * x_res;
    real_space_yaxis_matrix = (probe_row_idx - real_space_row_cen_idx) * y_res;
    real_space_r_matrix = sqrt(real_space_xaxis_matrix.^2+real_space_yaxis_matrix.^2);
    probe_amp = (1+cos(k*real_space_r_matrix.^2))/2;
    probe_phase = (1-probe_amp)*pi;

    probe = probe_amp .* exp(1i*probe_phase);
    probe = probe.*exp(-real_space_r_matrix.^2/(2*r0^2));
    
    if Broken== 1
        SpeckleSize = 0.1*r0; % [m]
        Speckle = exp(-power(real_space_xaxis_matrix,2)/2/SpeckleSize^2) .*  exp(-power(real_space_yaxis_matrix,2)/2/SpeckleSize^2);
        SeedNum = round((x_res*clip_size)^2/SpeckleSize/SpeckleSize);
        Seed = randi([1 ProbeConf.clip_size^2],SeedNum,1);
        Temp = zeros(ProbeConf.clip_size);
        Temp(Seed) = 1;
        Temp = conv2(Temp,Speckle,'same');
        probe = probe.*Temp.*exp(1i*Temp);
    end

    probe_temp = probe;
end

function probe_temp = generate_probe_from_sectionfile(init_cond,ProbeConf)
    clip_size = init_cond.effective_clip_size;
    detector_distance = init_cond.detector_distance;
    x_pixel_size = init_cond. effective_x_pixel_size;
    y_pixel_size = init_cond. effective_y_pixel_size;
    n_of_data = init_cond.n_of_data;
    amplitude = ProbeConf.photon_flux;
    wavelength = init_cond.wavelength;
    
    x_res = wavelength * detector_distance / clip_size / x_pixel_size;
    y_res = wavelength * detector_distance / clip_size / y_pixel_size;
    
    real_space_row_cen_idx = (clip_size+1)/2;
    real_space_col_cen_idx = (clip_size+1)/2;
    
    [probe_row_idx,probe_col_idx,~] = find(ones(clip_size,clip_size));
    probe_row_idx = reshape(probe_row_idx,clip_size,clip_size);
    probe_col_idx = reshape(probe_col_idx,clip_size,clip_size);
    real_space_xaxis = ProbeConf.real_space_xaxis;
    real_space_yaxis = ProbeConf.real_space_xaxis;
    
    probe_temp = zeros(clip_size,clip_size);
    
    if ~exist(ProbeConf.RefProbeFP,'File')
        error('!!!-----Section file for probe adapting doesn''t exist.----!!!')
        return
    end
    Temp = load(ProbeConf.RefProbeFP,'probe_info');
    AdaptProbeInfo = Temp.probe_info;
    if strcmpi(ProbeConf.AdaptProbeMode,'All')
        AdaptProbe= sum(AdaptProbeInfo.real_space,3);
    else
        AdaptProbe= AdaptProbeInfo.real_space(:,:,ProbeConf.AdaptProbeMode);
    end
    AdaptProbeXAxis = AdaptProbeInfo.real_space_xaxis;
    AdaptProbeYAxis = AdaptProbeInfo.real_space_yaxis;
    AdaptProbeXRes = AdaptProbeInfo.x_res;
    AdaptProbeYRes = AdaptProbeInfo.y_res;

    if abs(ProbeConf.AdaptPropagating) > 50E-6
        AdaptProbe = propagating_rescaling_wavefield(ProbeConf.AdaptPropagating,AdaptProbe,wavelength,AdaptProbeXAxis,AdaptProbeYAxis,AdaptProbeXRes,clip_size);
    else
        fprintf('!!!-----Ignore propagating, the propagating distance is small than 100 um----!!!')
    end
    
    ScalingRatio = x_res/AdaptProbeXRes;
    
    ScaledAdpatProbe = imresize(AdaptProbe,1/ScalingRatio);
    % after scaling, the x_res and y_res are the same between probe and
    % adaptprobe
    [ScaledAdpatProbeRowSize, ScaledAdpatProbeColSize] = size(ScaledAdpatProbe);
    ScaledAdpatProbeRowCen = round((ScaledAdpatProbeRowSize+1)/2);
    ScaledAdpatProbeColCen = round((ScaledAdpatProbeColSize+1)/2);
    
    if ScaledAdpatProbeRowSize>clip_size
        RowStart = ScaledAdpatProbeRowCen - (clip_size-1)/2;
        RowEnd = ScaledAdpatProbeRowCen + (clip_size-1)/2;
        ColStart = ScaledAdpatProbeColCen - (clip_size-1)/2;
        ColEnd = ScaledAdpatProbeColCen + (clip_size-1)/2;
        Probe = ScaledAdpatProbe(RowStart:RowEnd,ColStart:ColEnd);
        probe_temp = Probe;
    else
        if mod(ScaledAdpatProbeRowSize,2) == 0
            ScaledAdpatProbe(end,:) = [];
        end
        if mod(ScaledAdpatProbeColSize,2) == 0
            ScaledAdpatProbe(:,end) = [];
        end
        [ScaledAdpatProbeRowSize, ScaledAdpatProbeColSize] = size(ScaledAdpatProbe);
        RowStart = real_space_row_cen_idx - (ScaledAdpatProbeRowSize-1)/2;
        RowEnd = real_space_row_cen_idx + (ScaledAdpatProbeRowSize-1)/2;
        ColStart = real_space_col_cen_idx - (ScaledAdpatProbeColSize-1)/2;
        ColEnd = real_space_col_cen_idx + (ScaledAdpatProbeColSize-1)/2;
        probe_temp(RowStart:RowEnd,ColStart:ColEnd) = ScaledAdpatProbe;
    end
    
    % clear area around center
    
    boundary = (clip_size-1)/2;
    [X,Y] = meshgrid(-boundary:boundary,-boundary:boundary);
    w_pixel = boundary/x_res; % width for gaussian in pixel
    Gaussian = exp(-0.5*(X/w_pixel).^2) .* exp(-0.5*(Y/w_pixel).^2);
    probe_temp = probe_temp.*Gaussian;
    
end

function ProbeConf = get_ProbeConf(init_cond)
    config_probe_table = init_cond.config_tables.config_probe_table;

    ProbeConf.photon_flux = config_probe_table{'photon_flux','Value'}{1};
    ProbeConf.Mp = config_probe_table{'mixture_statue','Value'}{1};
    Method = config_probe_table{'probe_generate_method','Value'}{1};
    if strcmpi(Method,'G')
        Method = 'Gaussian';
    elseif strcmpi(Method,'Z')
        Method = 'Zone Plate';
    elseif strcmpi(Method,'A')
        Method = 'Adapt from another result';        
    end
    ProbeConf.ProbGenMode = Method;
    
    %% Define by parameters (gaussian)
    ProbeConf.GaussainBeamVSize = config_probe_table{'gaussian_ver_beamsize','Value'}{1}*1E-6; % from um to m
    ProbeConf.GaussainBeamHSize = config_probe_table{'gaussian_hor_beamsize','Value'}{1}*1E-6; % from um to m
    ProbeConf.GaussainBroken = logical(config_probe_table{'gaussian_broken_profile','Value'}{1}); 
    
    %% Define by parameters (zone plate)
    ProbeConf.ZoneplateOffFocal = config_probe_table{'zoneplate_off_focal_um','Value'}{1}*1E-6 ; % from um to m
    
    %% Adapt from another result
    ProbeConf.RefProbeFP = config_probe_table{'adapt_section_file','Value'}{1};
    ProbeConf.AdaptProbeMode = config_probe_table{'adapt_mode_index','Value'}{1};
    ProbeConf.AdapPosCorr = logical(config_probe_table{'adapt_pos_corr_from_file','Value'}{1});
    ProbeConf.AdaptPropagating = config_probe_table{'adapt_probe_propagating_um','Value'}{1}*1E-6; % from um to meter
    
    %% probe upstream constrain
    ProbeConf.UpStreamConstrain = logical(config_probe_table{'probe_upstream_constrain','Value'}{1});
    ProbeConf.ApertureDist = logical(config_probe_table{'aperture_distance_m','Value'}{1}); % in meter
    ProbeConf.ApertureSize = config_probe_table{'aperture_size_um','Value'}{1}*1E-6; % from um to meter
    
    clip_size = init_cond.effective_clip_size;
    n_of_data = init_cond.n_of_data;
    
    wavelength = init_cond.wavelength;
    detector_distance = init_cond.detector_distance;
    x_pixel_size = init_cond. effective_x_pixel_size;
    y_pixel_size = init_cond. effective_y_pixel_size;
    x_res = wavelength * detector_distance / clip_size / x_pixel_size;
    y_res = wavelength * detector_distance / clip_size / y_pixel_size;

    real_space_row_cen_idx = (clip_size+1)/2;
    real_space_col_cen_idx = (clip_size+1)/2;

    [probe_row_idx,probe_col_idx,~] = find(ones(clip_size,clip_size));
    probe_row_idx = reshape(probe_row_idx,clip_size,clip_size);
    probe_col_idx = reshape(probe_col_idx,clip_size,clip_size);
    real_space_xaxis = (probe_col_idx(1,:) - real_space_col_cen_idx) * x_res;
    real_space_yaxis = -(probe_row_idx(:,1) - real_space_row_cen_idx) * y_res;
    
    ProbeConf.wavelength = wavelength;
    ProbeConf.clip_size = clip_size;
    ProbeConf.n_of_data = n_of_data;
    ProbeConf.x_res = x_res;
    ProbeConf.y_res = y_res;
    ProbeConf.real_space_row_cen_idx = real_space_row_cen_idx;
    ProbeConf.real_space_col_cen_idx = real_space_col_cen_idx;
    ProbeConf.probe_row_idx = probe_row_idx;
    ProbeConf.probe_col_idx = probe_col_idx;
    ProbeConf.real_space_xaxis = real_space_xaxis;
    ProbeConf.real_space_yaxis = real_space_yaxis;
    
    ProbeConf.upstream_ROI = genUpstreamROI(ProbeConf);
end

function Upstream_ROI = genUpstreamROI(ProbeConf)
    probe = zeros(ProbeConf.clip_size);
    z = -ProbeConf.ApertureDist; % in meter
    ApertureSize = ProbeConf.ApertureSize; % in meter
    wavelength = ProbeConf.wavelength;
    xi_axis =  ProbeConf.real_space_xaxis;
    eta_axis = ProbeConf.real_space_yaxis;
    
    [propagated_probe,propagated_x_axis,propagated_y_axis]  = propagate_probe(z,probe,wavelength,xi_axis,eta_axis);
    [propagated_x, propagated_y] = meshgrid(propagated_x_axis,propagated_y_axis);
    propagated_r_matrix = sqrt(propagated_x.^2 + propagated_y.^2);
    Upstream_ROI = ones(ProbeConf.clip_size);
    Upstream_ROI(propagated_r_matrix >  (ApertureSize/2)) = 0.5;
end