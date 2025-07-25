function absolute_path = get_absolute_path(input_path)

if ~isstring(input_path) && ~ischar(input_path)
    error('Input path must be a string.')
end

file_obj = java.io.File(input_path);

if file_obj.isAbsolute()
    path_to_normalize = input_path;
else
    current_matlab_dir = pwd;
    path_to_normalize = fullfile(current_matlab_dir,input_path);
end


final_file_obj = java.io.File(path_to_normalize);
absolute_path = char(final_file_obj.getCanonicalPath());
