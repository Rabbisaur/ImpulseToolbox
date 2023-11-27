function MUAe = GetMUAeLFPwithInterpolationSubroutine4Shapley(RawData,RawDataSampleRate,MUAparameters,AllInterpolationPoints,trialIdx,LFPparameters)
% RawData must be in 2D matrix, time x trial format

% parameters
LFPsamplerate = LFPparameters.LFPsamplingrate;
% LFPlowpassFreq = LFPparameters.LFPlowpassFreq;

MUAesamplerate = MUAparameters.MUAeSamplingrate;
MUAebpFreq = MUAparameters.MUAeBandpassFreq;
% MUAelpFreq = MUAparameters.MUAeLowpassFreq;

% LFPsupersamplingRatio = lcm(RawDataSampleRate, LFPsamplerate) / RawDataSampleRate;
MUAesupersamplingRatio = lcm(RawDataSampleRate, MUAesamplerate) / RawDataSampleRate;

% LFPdownsamplingRatio = RawDataSampleRate * LFPsupersamplingRatio / LFPsamplerate;
MUAedownsamplingRatio = RawDataSampleRate * MUAesupersamplingRatio / MUAesamplerate;

% interpolationLength = interpolationLength / 1000;

% build filter

% band pass 1000-5000Hz
Fn = RawDataSampleRate/2;
Fbp=[MUAebpFreq(1),MUAebpFreq(2)];
N  = 2;    % filter order
[Bbp, Abp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]); % BandPass

% low pass 200Hz
% Fl= MUAelpFreq;
% Fn = MUAesamplerate/2;
% N = 2;
% [BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass

% Fl= MUAelpFreq;
% Fn = MUAesamplerate/2;
% N = 2;
% [BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass

% low pass 150Hz
% Fl = LFPlowpassFreq;
% Fn = RawDataSampleRate/2;
% N = 2;
% [BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass

% band stop 50Hz
Fn = RawDataSampleRate/2;
Fbp=[49,51];
N  = 2;    % filter order
[Bbs, Abs] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn],'stop'); % Bandstop

% band stop filtering 50Hz
BSdata = filtfilt(Bbs, Abs, RawData);

% LFP filtering
% LFPdata = filtfilt(BLFP, ALFP, BSdata);

% LFP supersampling
% if LFPsupersamplingRatio > 1
%     datalength = size(LFPdata,1);
%     x = 1:datalength;
%     xx = 1:(1/LFPsupersamplingRatio):datalength;
%     LFPdata = interp1(x,LFPdata,xx,'linear');
% %     LFPdata = spline(x, LFPdata, xx);
% end

% LFP downsampling
% LFPdata = downsample(LFPdata,LFPdownsamplingRatio);

% band pass filtering 500-5000Hz
BPdata = filtfilt(Bbp, Abp, BSdata);
clearvars BSdata

% half wave rectification
BPdata(BPdata<0) = 0;
FRdata = BPdata;
clearvars BPdata

% linear interpolation
if isempty(AllInterpolationPoints)
    % skip this step
else
    % linear interpolate just before and after each AR interpolation
    % points
    NumTrials = sum(trialIdx);
    InterpolationData = FRdata(:,trialIdx);
    for thisTrial = 1:NumTrials
        trialARinterPoints = AllInterpolationPoints{thisTrial};
        trialARinterPoints(1:2) = false; % set the first two and last two to be false
        trialARinterPoints(end-1:end) = false;
        tmpDiff = diff([trialARinterPoints;0]);
        trialMUAinterPoints = trialARinterPoints;
        trialMUAinterPoints(tmpDiff == 1) = true;
        trialMUAinterPoints(find(tmpDiff == -1)+1) = true;
        y = InterpolationData(:,thisTrial);
        x = 1:numel(y);
        xx = 1:numel(y);
        x(trialMUAinterPoints) = [];
        y(trialMUAinterPoints) = [];
        yy = interp1(x,y,xx,'linear');
        InterpolationData(:,thisTrial) = yy;
    end
    FRdata(:,trialIdx) = InterpolationData;
    % find Interpolation start point
    
end

baseline = std(reshape(FRdata,[],1)) * 2;
FRdata = FRdata > baseline;
FRdata = diff(FRdata,1,1);
FRdata = FRdata == 1;
windowSize = round(RawDataSampleRate/MUAesamplerate);

MUAeData = movsum(FRdata,windowSize,1);
% % band stop filtering 50Hz
% MUAeData = filtfilt(Bbs, Abs, FRdata);

% MUAe supersampling
if MUAesupersamplingRatio > 1
    datalength = size(MUAeData,1);
    x = 1:datalength;
    xx = 1:(1/MUAesupersamplingRatio):datalength;
    MUAeData = interp1(x,MUAeData,xx,'spline');
%     MUAeData = spline(x, MUAeData', xx)';
end
% if MUAesupersamplingRatio > 1
%     datalength = size(FRdata,1);
%     x = 1:datalength;
%     xx = 1:(1/MUAesupersamplingRatio):datalength;
%     FRdata = interp1(x,FRdata,xx,'linear');
% %     MUAeData = spline(x, MUAeData', xx)';
% end

% LFPsamplerate

% % MUAe downsampling % moved to in front of low pass filtering. 
% FRdata = downsample(FRdata,MUAedownsamplingRatio);

% low pass filtering 200Hz
% MUAeData = filtfilt(BMUAe, AMUAe, FRdata);

% MUAe downsampling
MUAeData = downsample(MUAeData,MUAedownsamplingRatio);


MUAe.data = MUAeData;
MUAe.MUAesamplerate = MUAesamplerate;
% MUAe.MUAebpFreq = MUAebpFreq;
% MUAe.MUAelpFreq = MUAelpFreq;

% LFP.data = LFPdata;
% LFP.LFPsamplerate = LFPsamplerate;
% LFP.LFPlowpassFreq = LFPlowpassFreq;