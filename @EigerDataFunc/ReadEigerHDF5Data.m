function DataOutput = ReadEigerHDF5Data(MasterInfo,RequestSN,XRange,YRange)
% ***** Output data format is UINT32 *****

if isempty(XRange)
    XRange = [1 MasterInfo.XPixelsInDetector];
end
if isempty(YRange)
    YRange = [1 MasterInfo.YPixelsInDetector];
end
% find linked data file
NLinkFile = length(MasterInfo.Links);
LinkFileSN = 0;
for LinkFileIdx = 1:NLinkFile
    if and(RequestSN >= MasterInfo.Links(LinkFileIdx).ImageNrLow, RequestSN <= MasterInfo.Links(LinkFileIdx).ImageNrHigh)
        LinkFileSN = LinkFileIdx;
        RequestSNinLinkFile = double(RequestSN - MasterInfo.Links(LinkFileIdx).ImageNrLow + 1);
        break
    end
end

% Request SN out of range
if LinkFileSN == 0
    DataOutput = [];
    return
end

DataOutput = h5read(MasterInfo.Links(LinkFileSN).FP,MasterInfo.Links(LinkFileSN).Location,[XRange(1),YRange(1),RequestSNinLinkFile],[XRange(2)-XRange(1)+1,YRange(2)-YRange(1)+1,1]);
DataOutput = transpose(DataOutput);