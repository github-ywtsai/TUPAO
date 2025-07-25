function create_config_tables_from_files(projectFF)

% =========================================================================
% CREATE_CONFIG_TABLES.M
% -------------------------------------------------------------------------
% This script reads three configuration files (config_init_cond.txt, 
% config_probe.txt, config_iteration.txt) and converts them into 
% MATLAB table data types for structured access in subsequent programs.
%
% Each table will contain 'Value' and 'Description' columns, and uses
% manually defined parameter names as RowNames for improved readability.
% =========================================================================

fprintf('Starting to convert config files to MATLAB tables...\n\n');

%% ==================== 1. Process config_init_cond.txt =====================
fprintf('1. Processing config_init_cond.txt...\n');

% --- Manually define parameter names (RowNames) ---
names_init_cond = {
    'folder_path'
    'scan_file_name'
    'mask_file_1'
    'mask_file_2'
    'mask_file_3'
    'mask_file_4'
    'beam_center_X'
    'beam_center_Y'
    'wavelength_A'
    'detector_distance_m'
    'pixel_size_X_um'
    'pixel_size_Y_um'
    'data_clipping_size'
    'probe_extend_factor'
    'random_seed'
    'core_algorithm'
    };

% --- Manually define English descriptions ---
descriptions_init_cond = {
    'The folder of the files'
    'The file name of the scan position file'
    'The file name of the csv mask file (None for skipping)'
    'The file name of the csv mask file (None for skipping)'
    'The file name of the csv mask file (None for skipping)'
    'The file name of the csv mask file (None for skipping)'
    'Beam center X [um] (Auto for loading from header)'
    'Beam center Y [um] (Auto for loading from header)'
    'Wavelength [A] (Auto for loading from header)'
    'Sample to detector distance [m] (Auto for loading from header)'
    'Detector X pixel size [um] (Auto for loading from header)'
    'Detector Y pixel size [um] (Auto for loading from header)'
    'Data clipping size'
    'Probe extending factor'
    'Seed for random (None for random seed)'
    'Core algorithm: ePIE, rPIE, mPIE, ML, DM'
    };


% Call the parsing function and create the table
file_init_cond = fullfile(projectFF,'config_init_cond.txt');
config_init_cond_table = parse_config_file_to_table(file_init_cond, names_init_cond, descriptions_init_cond);

% record projectFF information
config_init_cond_table{'project_folder_path',:} = {tools.get_absolute_path(projectFF),'The folder of the configurations and results'};

% Display the result
disp('Created config_init_cond_table:');
disp(config_init_cond_table);
fprintf('\n');


%% ===================== 2. Process config_probe.txt ========================
fprintf('2. Processing config_probe.txt...\n');

% --- Manually define parameter names (RowNames) ---
names_probe = {
    'photon_flux'
    'mixture_statue'
    'probe_generate_method'
    'gaussian_ver_beamsize'
    'gaussian_hor_beamsize'
    'gaussian_broken_profile'
    'zoneplate_off_focal_um'
    'adapt_section_file'
    'adapt_mode_index'
    'adapt_pos_corr_from_file'
    'adapt_probe_propagating_um'
    'probe_upstream_constrain'
    'aperture_distance_m'
    'aperture_size_um'
    };

% --- Manually define English descriptions ---
descriptions_probe = {
    'Photon flux / sec'
    'Mixture statue of probe'
    'Probe generate method (G: Gaussian, Z: ZonePlate, A:Adapt)'
    '(Gaussian) Ver. beamsize [um]'
    '(Gaussian) Hor. beamsize [um]'
    '(Gaussian) Broken profile (1: yes, 0: no)'
    '(Zoneplate) off focal distance [um]'
    '(Adapt) Section File Path'
    '(Adapt) mixture-state mode index for adapting (All for sum up)'
    '(Adapt) Adapt position corr. from section file (1: yes, 0: no)'
    '(Adapt) probe propagating [um]'
    'Probe upstream constrain (1: yes, 0: no)'
    'Aperture distance [m]'
    'Aperture size (the diameter in [um])'
    };

% Call the parsing function and create the table
file_probe = fullfile(projectFF,'config_probe.txt');
config_probe_table = parse_config_file_to_table(file_probe, names_probe, descriptions_probe);

% Display the result
disp('Created config_probe_table:');
disp(config_probe_table);
fprintf('\n');


%% =================== 3. Process config_iteration.txt ======================
fprintf('3. Processing config_iteration.txt...\n');

% --- Manually define parameter names (RowNames) ---
names_iteration = {
    'max_iteration'
    'save_results_period'
    'alpha_object_update'
    'probe_update_start'
    'beta_probe_update'
    'real_space_constraint_start'
    'real_space_constraint_factor'
    'pos_corr_start'
    'pos_corr_period'
    'pos_corr_points'
    'pos_corr_range_nm'
    'probe_deny_area_factor'
    'probe_deny_reduce_ratio'
    'draw_results_period'
    };
    
% --- Manually define English descriptions ---
descriptions_iteration = {
    'Max Iteration Number'
    'Save Results period'
    'Alpha for object update'
    'Probe update start run'
    'Beta for probe update'
    'Real space constraint start run'
    'Real space constraint factor'
    'Position correction start run'
    'Position correction period'
    'Position correction points'
    'Position correction range (nm)'
    'Probe deny area factor'
    'Probe deny reducing ratio'
    'Draw results period'
    };

% Call the parsing function and create the table
file_iteration = fullfile(projectFF,'config_iteration.txt');
config_iteration_table = parse_config_file_to_table(file_iteration, names_iteration, descriptions_iteration);

% Display the result
disp('Created config_iteration_table:');
disp(config_iteration_table);
fprintf('\n');


%% ======================= 4. Save Tables        ==========================
save(fullfile(projectFF,'config_tables.mat'),'config_init_cond_table','config_probe_table','config_iteration_table');

%% ======================= 5. Clean Up Workspace ==========================
% This section cleans up all temporary variables used for setup,
% leaving only the final config tables in the workspace.

%fprintf('Cleaning up temporary variables...\n');
%clearvars -except config_init_cond_table config_probe_table config_iteration_table

%fprintf('All tasks completed. Workspace is clean.\n');
end

%% ======================== Core Parsing Function =================================
function config_table = parse_config_file_to_table(filepath, row_names, descriptions)
    % This function reads the specified config file and converts its content
    % into a MATLAB table.
    
    values = {};

    % Open the file for reading
    fid = fopen(filepath, 'r');
    if fid == -1
        error('Cannot open file: %s', filepath);
    end

    % Skip the header line
    fgetl(fid); 

    % Read and parse line by line
    while ~feof(fid)
        line = fgetl(fid);
        
        comment_pos = strfind(line, '%');

        if isempty(comment_pos)
            val_str = strtrim(line);
        else
            val_str = strtrim(line(1:comment_pos(1)-1));
        end

        % Only add a row if a 'Value' is present
        if ~isempty(val_str)
            % convert to number if possible
            [num, status] = str2num(val_str); % Use str2num for scientific notation like 1E3
            if status
                values{end+1, 1} = num;
            else
                values{end+1, 1} = val_str; % Keep as string if conversion fails
            end
            
            
        end
    end
    fclose(fid);

    % Check if the number of parsed values matches the number of defined names
    if length(values) ~= length(row_names)
        warning('In %s, number of parsed values (%d) does not match number of defined names (%d). Please check file content or name definitions.', ...
                filepath, length(values), length(row_names));
        % Adjust arrays to the minimum length to prevent erroring
        min_len = min(length(values), min(length(row_names), length(descriptions)));
        values = values(1:min_len);
        row_names = row_names(1:min_len);
        descriptions = descriptions(1:min_len);
    end

    % Create the final table
    config_table = table(values, descriptions, 'RowNames', row_names, 'VariableNames', {'Value', 'Description'});
end