classdef EigerDataFunc
    % 2021/07/15
    methods (Static = true)
        MasterInfo =ReadEigerHDF5Master(MasterFP);
        DataOutput = ReadEigerHDF5Data(MasterInfo,RequestSN,XRange,YRange);
         EnvSettingCheck = EnvSetting();
    end
end