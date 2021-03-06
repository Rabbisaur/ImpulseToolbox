function MUAe = GetMUAewithInterpolationSubroutine(RawData,RawDataSampleRate,MUAparameters,interpolationTimepoints,interpolationLength)
% RawData must be in 2D matrix, time x trial format

% parameters


MUAesamplerate = MUAparameters.MUAeSamplingrate;
MUAebpFreq = MUAparameters.MUAeBandpassFreq;
MUAelpFreq = MUAparameters.MUAeLowpassFreq;


MUAesupersamplingRatio = lcm(RawDataSampleRate, MUAesamplerate) / RawDataSampleRate;
MUAedownsamplingRatio = RawDataSampleRate * MUAesupersamplingRatio / MUAesamplerate;

% interpolationLength = interpolationLength / 1000;

% build filter

% band pass 500-5000Hz
Fn = RawDataSampleRate/2;
Fbp=[MUAebpFreq(1),MUAebpFreq(2)];
N  = 2;    % filter order
[Bbp, Abp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]); % BandPass

% low pass 200Hz
Fl= MUAelpFreq;
% Fn = MUAesamplerate/2;
Fn = RawDataSampleRate/2;
N = 2;
[BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass


% band pass filtering 500-5000Hz
BPdata = filtfilt(Bbp, Abp, RawData);
clearvars RawData

% fullwave rectification
FRdata = abs(BPdata);
clearvars BPdata

% linear interpolation
if false
    figure
    hold on
    plot(FRdata,'blue')
    MUAeData2 = filtfilt(BMUAe, AMUAe, FRdata);
end
if interpolationLength > 0
    % calculate interpolation parameters
    NumTrials = size(FRdata,2);
    xx = 1:size(FRdata,1);
    NumGaps = size(interpolationTimepoints,1);
    for thisTrial = 1:NumTrials
        if sum(interpolationTimepoints(:,thisTrial) == 0) > 0
            continue
        end
        interpolationStartTimepoints = interpolationTimepoints(:,thisTrial) - floor(interpolationLength * RawDataSampleRate / 2);
        idx = 0:(round(interpolationLength * RawDataSampleRate)-1);
        idx = repmat(interpolationStartTimepoints,1,round(interpolationLength * RawDataSampleRate)) + ...
            repmat(idx,NumGaps,1);
        idx = reshape(idx,[],1);
        x = xx;
        x(idx) = [];
        
        tmpdata = FRdata(:,thisTrial);
        tmpdata(idx) = [];
        FRdata(:,thisTrial) = interp1(x,tmpdata,xx,'linear');
    end
end

if false
    plot(FRdata,'red')
    xlim([1000 size(FRdata,1)-1000])
end
% low pass filtering 200Hz
MUAeData = filtfilt(BMUAe, AMUAe, FRdata);

if false
    figure
    hold on
    plot(mean(MUAeData2,2),'blue')
    plot(mean(MUAeData,2),'red')
    xlim([1000 size(FRdata,1)-1000])
end

if MUAesupersamplingRatio > 1
    datalength = size(MUAeData,1);
    x = 1:datalength;
    xx = 1:(1/MUAesupersamplingRatio):datalength;
    MUAeData = interp1(x,MUAeData,xx,'linear');
%     MUAeData = spline(x, MUAeData', xx)';
end

% MUAe downsampling
MUAeData = downsample(MUAeData,MUAedownsamplingRatio);


MUAe.data = MUAeData;
MUAe.MUAesamplerate = MUAesamplerate;
MUAe.MUAebpFreq = MUAebpFreq;
MUAe.MUAelpFreq = MUAelpFreq;
MUAe.interpolationLength = interpolationLength;