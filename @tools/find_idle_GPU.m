    function idle_GPU_index = find_idle_GPU(options)
    arguments
        options.release_GPU = true;
    end
    if options.release_GPU
        fprintf('Releasing occupied GPU...\n')
        gpuDevice([]);
    end
        idle_GPU_index = [];
        command = 'nvidia-smi pmon -c 1';
        [status, cmdout] = system(command);
        lines = strsplit(cmdout,'\n');
        for ii = 1:length(lines)
            line = strtrim(lines{ii});
            if isempty(line) || startsWith(line,'#')
                continue;
            end
            temp = strsplit(line,' ');
            GPU_index = str2double(temp{1})+1;
            name = temp{end};
            if strcmp(name,'-')
                idle_GPU_index = GPU_index;
                fprintf('Find idle GPU %d.\n',idle_GPU_index);
                break
            end
        end
        if isempty(idle_GPU_index)
            error('No idle GPU found.\n')
        end
    end