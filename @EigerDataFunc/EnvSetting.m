function EnvSettingCheck = EnvSetting()
% EnvSettingCheck = true;

if ispc
    OS = 'PC';
elseif isunix
    OS = 'UNIX';
elseif ismac
    OS = 'MAC';
end

switch OS
    case 'PC'
        DLLCheck = CheckingDLL();
        if DLLCheck == 1
            % disp('Environment configuration checking passed.')
            EnvSettingCheck = true;
        else
            DLLPATH = fullfile(pwd,'./@EigerDataFunc');
            cmd = sprintf('setx HDF5_PLUGIN_PATH "%s', DLLPATH);
            system(cmd);
            disp('Environment variables are changed.')
            disp('The changing is NOT available unitl restarting Matlab.')
            EnvSettingCheck = false;

        end
    case 'UNIX'
        disp('The Environment checking is ignored.')
        disp('Environment variables must be configured manually in UNIX.')
        setenv('HDF5_PLUGIN_PATH','/blsw/opt/areaDetector/root/usr/lib/h5plugin');
        EnvSettingCheck = true;
    case 'MAC'
        disp('The Environment checking is ignored.')
        disp('Environment variables must be configured manually in MAC.')
        EnvSettingCheck = true;
end

function DLLCheck = CheckingDLL()
DLLMemberList = {'libh5blosc.dll','libh5bz2.dll','libh5lz4.dll','libh5lzf.dll','libh5mafisc.dll','libh5zfp.dll'};
CheckingList = zeros(length(DLLMemberList),1);
Default_HDF5_PLUGIN_PATH = getenv('HDF5_PLUGIN_PATH');
for DLLSN = 1:length(CheckingList)
    DLLFP = fullfile(Default_HDF5_PLUGIN_PATH,DLLMemberList{DLLSN});
    CheckingList(DLLSN) = logical(exist(DLLFP,'file'));
end
DLLCheck = all(CheckingList);
