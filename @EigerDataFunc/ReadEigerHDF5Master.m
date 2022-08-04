function MasterInfo =ReadEigerHDF5Master(MasterFP)

[MasterFF,MasterFN] = AnalyzeMasterFP(MasterFP);

MasterInfo.MasterFF = MasterFF;
MasterInfo.MasterFN = MasterFN;
MasterInfo.MasterFP = MasterFP;
MasterInfo.BitDepthImage = double(h5read(MasterFP,'/entry/instrument/detector/bit_depth_image'));
MasterInfo.BadPointValue = power(2,MasterInfo.BitDepthImage)-1;
MasterInfo.XPixelsInDetector = double(h5read(MasterFP,'/entry/instrument/detector/detectorSpecific/x_pixels_in_detector'));
MasterInfo.YPixelsInDetector = double(h5read(MasterFP,'/entry/instrument/detector/detectorSpecific/y_pixels_in_detector'));
MasterInfo.CountTime = double(h5read(MasterFP,'/entry/instrument/detector/count_time'));
MasterInfo.DetectorDistance = double(h5read(MasterFP,'/entry/instrument/detector/detector_distance')); % [m]
MasterInfo.XPixelSize = double(h5read(MasterFP,'/entry/instrument/detector/x_pixel_size')); % [m]
MasterInfo.YPixelSize = double(h5read(MasterFP,'/entry/instrument/detector/y_pixel_size')); % [m]
MasterInfo.Wavelength = double(h5read(MasterFP,'/entry/instrument/beam/incident_wavelength'))*1E-10; % read as [A], save as [m]
MasterInfo.BeamCenterX= round(double(h5read(MasterFP,'/entry/instrument/detector/beam_center_x')));
MasterInfo.BeamCenterY= round(double(h5read(MasterFP,'/entry/instrument/detector/beam_center_y')));
MasterInfo.PixelMask = logical(transpose(h5read(MasterFP,'/entry/instrument/detector/detectorSpecific/pixel_mask')));
MasterInfo.AveragedDataSheetNum = 1; % Record the averaged data sheets 

temp = h5info(MasterFP,'/entry/data');
Links = temp.Links;
NLinks = length(Links); 
for LinkIdx = 1:NLinks
    LinkedFN = Links(LinkIdx).Value{1};
    LinkedFP = fullfile(MasterInfo.MasterFF, LinkedFN);
    if exist(LinkedFP,'file')
        MasterInfo.Links(LinkIdx).FN = LinkedFN;
        MasterInfo.Links(LinkIdx).FF = MasterInfo.MasterFF;
        MasterInfo.Links(LinkIdx).FP = LinkedFP;
        MasterInfo.Links(LinkIdx).Location = Links(LinkIdx).Value{2};
        MasterInfo.Links(LinkIdx).ImageNrLow = h5readatt(MasterInfo.Links(LinkIdx).FP,MasterInfo.Links(LinkIdx).Location,'image_nr_low');
        MasterInfo.Links(LinkIdx).ImageNrHigh = h5readatt(MasterInfo.Links(LinkIdx).FP,MasterInfo.Links(LinkIdx).Location,'image_nr_high');
    else
        break
    end
end
MasterInfo.DataSheetNum = MasterInfo.Links(end).ImageNrHigh;




function [MasterFF,MasterFN] = AnalyzeMasterFP(MasterFP)

[MasterFF,MasterFN,MasterEXT] = fileparts(MasterFP);
MasterFN = [MasterFN MasterEXT];

if ~strcmpi(MasterFP, fullfile(MasterFF,MasterFN))
    error('Master file name analyze fault.')
end