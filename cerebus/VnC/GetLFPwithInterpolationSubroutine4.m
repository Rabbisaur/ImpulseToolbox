function [LFP] = GetLFPwithInterpolationSubroutine4(RawData,RawDataSampleRate,LFPparameters)
% RawData must be in 2D matrix, time x trial format

% parameters
LFPsamplerate = LFPparameters.LFPsamplingrate;
LFPlowpassFreq = LFPparameters.LFPlowpassFreq;


LFPsupersamplingRatio = lcm(RawDataSampleRate, LFPsamplerate) / RawDataSampleRate;
% MUAesupersamplingRatio = lcm(RawDataSampleRate, MUAesamplerate) / RawDataSampleRate;

LFPdownsamplingRatio = RawDataSampleRate * LFPsupersamplingRatio / LFPsamplerate;
% MUAedownsamplingRatio = RawDataSampleRate * MUAesupersamplingRatio / MUAesamplerate;

% interpolationLength = interpolationLength / 1000;

% build filter

% band pass 1000-5000Hz
% Fn = RawDataSampleRate/2;
% Fbp=[MUAebpFreq(1),MUAebpFreq(2)];
% N  = 2;    % filter order
% [Bbp, Abp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]); % BandPass

% low pass 200Hz
% Fl= MUAelpFreq;
% Fn = RawDataSampleRate/2;
% N = 2;
% [BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass

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
% LFPdata = zeros(floor(size(tmpdata,1)/LFPdownsamplingRatio)+1,size(tmpdata,2),size(tmpdata,3));
% for thisElec = 1:NumElec
%     LFPdata(:,:,thisElec) = downsample(tmpdata(:,:,thisElec),LFPdownsamplingRatio);
% end

LFPdata = downsample(tmpdata,LFPdownsamplingRatio);

LFP.data = LFPdata;
LFP.LFPsamplerate = LFPsamplerate;
LFP.LFPlowpassFreq = LFPlowpassFreq;