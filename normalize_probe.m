function probe_info = normalize_probe(probe_info,measurement_info)

% assume the object is in ones, the sum of the probe should be in consist
% on the target and detector.

[~,~,num_diffrcation_pattern] = size(measurement_info.measured_amp);
individual_roi = ones(size(measurement_info.individual_mask));
individual_roi(measurement_info.individual_mask==true) = nan;
overall_roi = sum(individual_roi,3);
flat_diffraction_pattern_amp = sum(measurement_info.measured_amp.*individual_roi,3)/num_diffrcation_pattern;
flat_diffraction_pattern_amp_total_count = sum(flat_diffraction_pattern_amp,'all','omitnan');

main_probe = probe_info.real_space(:,:,1);
main_probe_on_detector = ifftshift(fft2(ifftshift(main_probe)));

main_probe_on_detector_total_count = sum(abs(main_probe_on_detector).*overall_roi,'all','omitnan');

ratio = flat_diffraction_pattern_amp_total_count/main_probe_on_detector_total_count;


probe_info.real_space = probe_info.real_space*ratio;