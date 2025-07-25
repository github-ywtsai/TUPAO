classdef tools
    methods (Static = true)
        create_config_tables_from_files(projectFF)
        idle_GPU_index = find_idle_GPU()
        absolute_path = get_absolute_path(input_path)
        ptycho_pachage = initialize(projectFF,options)
        save_section_file(fn,ptycho_package)
    end
end