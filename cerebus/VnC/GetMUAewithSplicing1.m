function MUAe = GetMUAewithSplicing(RawData,RawDataSampleRate,MUAparameters,interpolationTimepoints)
% RawData must be in 2D matrix, time x trial format

% parameters


MUAesamplerate = MUAparameters.MUAeSamplingrate;
MUAebpFreq = MUAparameters.MUAeBandpassFreq;
MUAelpFreq = MUAparameters.MUAeLowpassFreq;


% MUAesupersamplingRatio = lcm(RawDataSampleRate, MUAesamplerate) / RawDataSampleRate;
MUAedownsamplingRatio = RawDataSampleRate / MUAesamplerate;

% interpolationLength = interpolationLength / 1000;

% build filter

% band pass 500-5000Hz
Fn = RawDataSampleRate/2;
Fbp=[MUAebpFreq(1),MUAebpFreq(2)];
N  = 2;    % filter order
[Bbp, Abp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]); % BandPass

% low pass 200Hz
Fl= MUAelpFreq;
Fn = RawDataSampleRate/2;
N = 2;
[BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass

% spilice data
idx0 = isnan(RawData);
idx0 = sum(idx0,2) > 0;
xx0 = 1:size(RawData,1);
x0 = xx0;
x0(idx0) = [];
% idx1 = repmat(idx0,1,size(RawData,2));
tmpdata = RawData;
RawData = [];
for thisTrial = 1:size(tmpdata,2)
    y = tmpdata(:,thisTrial);
    y(idx0) = [];
    RawData(:,thisTrial) = y;
end


% band pass filtering 500-5000Hz
BPdata = filtfilt(Bbp, Abp, RawData);
clearvars RawData

% fullwave rectification
FRdata = abs(BPdata);
clearvars BPdata

% linear interpolation
% if interpolationLength > 0
%     % calculate interpolation parameters
%     NumTrials = size(FRdata,2);
%     xx = 1:size(FRdata,1);
%     NumGaps = size(interpolationTimepoints,1);
%     for thisTrial = 1:NumTrials
%         if sum(interpolationTimepoints(:,thisTrial) == 0) > 0
%             continue
%         end
%         interpolationStartTimepoints = interpolationTimepoints(:,thisTrial) - round(interpolationLength * RawDataSampleRate / 2);
%         idx = 0:(round(interpolationLength * RawDataSampleRate)-1);
%         idx = repmat(interpolationStartTimepoints,1,round(interpolationLength * RawDataSampleRate)) + ...
%             repmat(idx,NumGaps,1);
%         idx = reshape(idx,[],1);
%         x = xx;
%         x(idx) = [];
%         
%         tmpdata = FRdata(:,thisTrial);
%         tmpdata(idx) = [];
%         FRdata(:,thisTrial) = interp1(x,tmpdata,xx,'linear');
%     end
% end

% MUAe supersampling
% if MUAesupersamplingRatio > 1
%     datalength = size(MUAeData,1);
%     x = 1:datalength;
%     xx = 1:(1/MUAesupersamplingRatio):datalength;
%     MUAeData = interp1(x,MUAeData,xx,'spline');
% %     MUAeData = spline(x, MUAeData', xx)';
% end
% if MUAesupersamplingRatio > 1
%     datalength = size(FRdata,1);
%     x = 1:datalength;
%     xx = 1:(1/MUAesupersamplingRatio):datalength;
%     FRdata = interp1(x,FRdata,xx,'linear');
% %     MUAeData = spline(x, MUAeData', xx)';
% end

% MUAe downsampling % moved to in front of low pass filtering. 
% FRdata = downsample(FRdata,MUAedownsamplingRatio);

% low pass filtering 200Hz
MUAeData = filtfilt(BMUAe, AMUAe, FRdata);
tmpdata = MUAeData;
MUAeData = [];

% linear interpolate gaps

for thisTrial = 1:size(tmpdata,2)
    
    MUAeData(:,thisTrial) = interp1(x0,tmpdata(:,thisTrial),xx0,'linear');
end

% MUAe downsampling
% MUAeData = downsample(MUAeData,MUAedownsamplingRatio);
downsampleIdx = round(1:MUAedownsamplingRatio:size(MUAeData,1));
tmpdata = MUAeData;
MUAeData = [];
for thisTrial = 1:size(tmpdata,2)
    MUAeData(:,thisTrial) = tmpdata(downsampleIdx,thisTrial);
end

MUAe.data = MUAeData;
MUAe.MUAesamplerate = MUAesamplerate;
MUAe.MUAebpFreq = MUAebpFreq;
MUAe.MUAelpFreq = MUAelpFreq;