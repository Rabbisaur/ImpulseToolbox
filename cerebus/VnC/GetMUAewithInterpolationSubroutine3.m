function [MUAe, BPdata]= GetMUAewithInterpolationSubroutine3(RawData,RawDataSampleRate,MUAparameters,interpolationTimepoints,interpolationLength,BPdata)
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

if isempty(BPdata)
    % band pass filtering 500-5000Hz
    BPdata = filtfilt(Bbp, Abp, RawData);
    clearvars RawData
else
    % do nothing
end
% fullwave rectification
FRdata = abs(BPdata);
% clearvars BPdata




% linear interpolation
if false
    figure
    hold on
    plot(FRdata,'blue')
    MUAeData2 = filtfilt(BMUAe, AMUAe, FRdata);
end
if interpolationLength > 0
    paddingSize = 100;
    padding = zeros(paddingSize,1);
    % calculate interpolation parameters
    NumTrials = size(FRdata,2);
    
    NumGaps = size(interpolationTimepoints,1);
    for thisTrial = 1:NumTrials
        if sum(interpolationTimepoints(:,thisTrial) == 0) > 0
            continue
        end
        interpolationStartTimepoints = interpolationTimepoints(:,thisTrial) - floor(interpolationLength * RawDataSampleRate / 2);
        idx = 0:(round(interpolationLength * RawDataSampleRate)-1);
        idx = repmat(interpolationStartTimepoints,1,round(interpolationLength * RawDataSampleRate)) + ...
            repmat(idx,NumGaps,1);
        idx = reshape(idx',[],1);
        
        
        tmpdata = FRdata(:,thisTrial);
        
        tmpdata = [padding;tmpdata;padding];
        idx = idx + paddingSize;
        
        xx = 1:size(tmpdata);
        x = xx;
        x(idx) = [];
        
        tmpdata(idx) = [];
        tmpdata = interp1(x,tmpdata,xx,'linear');
        tmpdata = tmpdata((paddingSize+1):(end-paddingSize));
        FRdata(:,thisTrial) = tmpdata;
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

% datafilt = MUAeData;
% FsD = MUAesamplerate;
% Fn = FsD/2; % Downsampled Nyquist frequency
% for v = [linenoiseFreq linenoiseFreq*2 linenoiseFreq*3]
%     Fbp = [v-2,v+2];
%     [Blp, Alp] = butter(4, [min(Fbp)/Fn max(Fbp)/Fn],'stop'); % compute filter coefficients
%     datafilt = filtfilt(Blp, Alp, datafilt);
% end
% MUAeData = datafilt;

MUAe.data = MUAeData;
MUAe.MUAesamplerate = MUAesamplerate;
MUAe.MUAebpFreq = MUAebpFreq;
MUAe.MUAelpFreq = MUAelpFreq;
MUAe.interpolationLength = interpolationLength;