function [MUAe, LFP] = GetMUAeLFPwithInterpolationSubroutine(RawData,RawDataSampleRate,MUAparameters,interpolationTimepoints,interpolationLength,LFPparameters)
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

interpolationLength = interpolationLength / 1000;

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
    LFPdata = spline(x, LFPdata, xx);
end

% LFP downsampling
LFPdata = downsample(LFPdata,LFPdownsamplingRatio);

% band pass filtering 500-5000Hz
BPdata = filtfilt(Bbp, Abp, BSdata);

% fullwave rectification
FRdata = abs(BPdata);
clearvars BPdata

% linear interpolation
if interpolationLength > 0
    % calculate interpolation parameters
    xx = 1:size(FRdata,1);
    if size(interpolationTimepoints,1)<size(interpolationTimepoints,2)
        interpolationTimepoints = interpolationTimepoints';
    end
    NumGaps = numel(interpolationTimepoints);
    interpolationStartTimepoints = interpolationTimepoints - round(interpolationLength * RawDataSampleRate / 2);
    idx = 0:(round(interpolationLength * RawDataSampleRate)-1);
    interpolationTimepoints = repmat(interpolationStartTimepoints,1,round(interpolationLength * RawDataSampleRate)) + ...
        repmat(idx,NumGaps,1);
    interpolationTimepoints = reshape(interpolationTimepoints,[],1);
    x = xx;
    x(interpolationTimepoints) = [];
    
    NumTrials = size(FRdata,2);
    for thisTrial = 1:NumTrials
        tmpdata = FRdata(:,thisTrial);
        tmpdata(interpolationTimepoints) = [];
        FRdata(:,thisTrial) = interp1(x,tmpdata,xx,'linear');
    end
end

% low pass filtering 200Hz
MUAeData = filtfilt(BMUAe, AMUAe, FRdata);

% MUAe supersampling
if MUAesupersamplingRatio > 1
    datalength = size(MUAeData,1);
    x = 1:datalength;
    xx = 1:(1/MUAesupersamplingRatio):datalength;
    MUAeData = spline(x, MUAeData', xx)';
end

% MUAe downsampling
MUAeData = downsample(MUAeData,MUAedownsamplingRatio);

MUAe.data = MUAeData;
MUAe.MUAesamplerate = MUAesamplerate;
MUAe.MUAebpFreq = MUAebpFreq;
MUAe.MUAelpFreq = MUAelpFreq;

LFP.data = LFPdata;
LFP.LFPsamplerate = LFPsamplerate;
LFP.LFPlowpassFreq = LFPlowpassFreq;