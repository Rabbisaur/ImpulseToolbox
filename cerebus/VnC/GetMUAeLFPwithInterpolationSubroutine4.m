function [MUAe, LFP] = GetMUAeLFPwithInterpolationSubroutine4(RawData,RawDataSampleRate,MUAparameters,AllInterpolationPoints,trialIdx,LFPparameters)
% RawData must be in 2D matrix, time x trial format

% parameters
LFPsamplerate = LFPparameters.LFPsamplingrate;
LFPlowpassFreq = LFPparameters.LFPlowpassFreq;

MUAesamplerate = MUAparameters.MUAeSamplingrate;
MUAebpFreq = MUAparameters.MUAeBandpassFreq;
MUAelpFreq = MUAparameters.MUAeLowpassFreq;

LFPsupersamplingRatio = lcm(RawDataSampleRate, LFPsamplerate) / RawDataSampleRate;
MUAesupersamplingRatio = lcm(RawDataSampleRate, MUAesamplerate) / RawDataSampleRate;

LFPdownsamplingRatio = RawDataSampleRate * LFPsupersamplingRatio / LFPsamplerate;
MUAedownsamplingRatio = RawDataSampleRate * MUAesupersamplingRatio / MUAesamplerate;

% interpolationLength = interpolationLength / 1000;

% build filter

% band pass 1000-5000Hz
Fn = RawDataSampleRate/2;
Fbp=[MUAebpFreq(1),MUAebpFreq(2)];
N  = 2;    % filter order
[Bbp, Abp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]); % BandPass

% low pass 200Hz
Fl= MUAelpFreq;
Fn = RawDataSampleRate/2;
N = 2;
[BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass

% Fl= MUAelpFreq;
% Fn = MUAesamplerate/2;
% N = 2;
% [BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass

% low pass 150Hz
Fl = LFPlowpassFreq;
Fn = RawDataSampleRate/2;
N = 2;
[BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass

% band stop 50Hz
Fn = RawDataSampleRate/2;
Fbp=[49,51];
N  = 2;    % filter order
[Bbs, Abs] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn],'stop'); % Bandstop

% band stop filtering 50Hz
BSdata = filtfilt(Bbs, Abs, RawData);

% LFP filtering
LFPdata = filtfilt(BLFP, ALFP, BSdata);

% LFP supersampling
if LFPsupersamplingRatio > 1
    datalength = size(LFPdata,1);
    x = 1:datalength;
    xx = 1:(1/LFPsupersamplingRatio):datalength;
    LFPdata = interp1(x,LFPdata,xx,'linear');
%     LFPdata = spline(x, LFPdata, xx);
end

% LFP downsampling
NumElec = size(LFPdata,3);
tmpdata = LFPdata;
LFPdata = zeros(floor(size(tmpdata,1)/LFPdownsamplingRatio)+1,size(tmpdata,2),size(tmpdata,3));
% for thisElec = 1:NumElec
%     LFPdata(:,:,thisElec) = downsample(tmpdata(:,:,thisElec),LFPdownsamplingRatio);
% end

LFPdata = downsample(tmpdata,LFPdownsamplingRatio);

% band pass filtering 500-5000Hz
BPdata = filtfilt(Bbp, Abp, BSdata);
clearvars BSdata

% fullwave rectification
FRdata = abs(BPdata);
clearvars BPdata

% linear interpolation
if isempty(AllInterpolationPoints)
    % skip this step
else
    % linear interpolate just before and after each AR interpolation
    % points
    NumTrials = size(FRdata,2);
    FRdata2 = FRdata;
%     parfor thisElec = 1:size(FRdata,3)
%         InterpolationData = FRdata(:,:,thisElec);
%         for thisTrial = 1:NumTrials
%             trialARinterPoints = AllInterpolationPoints{thisTrial};
%             %         trialARinterPoints(1:2) = false; % set the first two and last two to be false
%             %         trialARinterPoints(end-1:end) = false;
%             %         tmpDiff = diff([trialARinterPoints;0]);
%             %         trialMUAinterPoints = trialARinterPoints;
%             %         trialMUAinterPoints(tmpDiff == 1) = true;
%             %         trialMUAinterPoints(find(tmpDiff == -1)+1) = true;
%             y = InterpolationData(:,thisTrial);
%             x = 1:numel(y);
%             xx = 1:numel(y);
%             x(trialARinterPoints) = [];
%             y(trialARinterPoints) = [];
%             yy = interp1(x,y,xx,'linear');
%             InterpolationData(:,thisTrial) = yy;
%         end
%         FRdata(:,:,thisElec) = InterpolationData;
%     end
    
% improved algorithm
    for thisTrial = 1:NumTrials
        trialARinterPoints = AllInterpolationPoints{thisTrial};
        InterpolationData = FRdata2(:,thisTrial,:);
        x = 1:size(InterpolationData,1);
        xx = 1:size(InterpolationData,1);
        x(trialARinterPoints) = [];
        InterpolationData(trialARinterPoints,:,:) = [];
        InterpolationData = interp1(x,InterpolationData,xx,'linear');
        
        FRdata2(:,thisTrial,:) = InterpolationData;
    end
    
    
    % find Interpolation start point
    
end

% low pass filtering 200Hz
idx = ~isfinite(FRdata);
FRdata(idx) = 0;
FRdata = filtfilt(BMUAe, AMUAe, FRdata);

% MUAe supersampling
% if MUAesupersamplingRatio > 1
%     datalength = size(MUAeData,1);
%     x = 1:datalength;
%     xx = 1:(1/MUAesupersamplingRatio):datalength;
%     MUAeData = interp1(x,MUAeData,xx,'spline');
% %     MUAeData = spline(x, MUAeData', xx)';
% end
if MUAesupersamplingRatio > 1
    datalength = size(FRdata,1);
    x = 1:datalength;
    xx = 1:(1/MUAesupersamplingRatio):datalength;
    FRdata = interp1(x,FRdata,xx,'linear');
%     MUAeData = spline(x, MUAeData', xx)';
end

% MUAe downsampling % moved to in front of low pass filtering. 
% FRdata = downsample(FRdata,MUAedownsamplingRatio);

NumElec = size(FRdata,3);
tmpdata = FRdata;
FRdata = zeros(floor(size(tmpdata,1)/MUAedownsamplingRatio)+1,size(FRdata,2),size(FRdata,3));
% for thisElec = 1:NumElec
%     FRdata(:,:,thisElec) = downsample(tmpdata(:,:,thisElec),MUAedownsamplingRatio);
% end
FRdata = downsample(tmpdata,MUAedownsamplingRatio);

% MUAe downsampling
% MUAeData = downsample(MUAeData,MUAedownsamplingRatio);


MUAe.data = FRdata;
MUAe.MUAesamplerate = MUAesamplerate;
MUAe.MUAebpFreq = MUAebpFreq;
MUAe.MUAelpFreq = MUAelpFreq;

LFP.data = LFPdata;
LFP.LFPsamplerate = LFPsamplerate;
LFP.LFPlowpassFreq = LFPlowpassFreq;