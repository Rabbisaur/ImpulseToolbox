function [ARdata,ArtifactLocation] = VnCRemoveArtifactsSubroutine3(RawTrace,StimulationParam,StimStartEndTime,SampleRate,SupersamplingRatio,debug)

% inputs

% RawTrace = AMF';
% SampleRate = 30000;
% SupersamplingRatio = 16;

% StimulationParam.numberofStimPulses = 50;
% StimulationParam.PulsesFrequency = 300; %Hz
% StimulationParam.stimulationwaveformwidth = 0.4; % ms
% StimulationParam.WaveDPbefore = 12; %number of datapoints, relative to peak of artifact waveform
% StimulationParam.WaveDPafter = 12; %number of datapoints, relative to peak of artifact waveform
dataAmpResolution = 256;
AvgWindowSize = 3;

% debug = 1;
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
% parameters
% Gain = 100;
% AverageTrace = mean(RawTrace,2);
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
StimulationFreq = 1/(StimulationParam.stimulationwaveformwidth/1000);
numberofStimPulses = StimulationParam.numberofStimPulses;
stimulationwaveformwidth = StimulationParam.stimulationwaveformwidth * SampleRate / 1000 * SupersamplingRatio;

WaveDPbefore = StimulationParam.WaveDPbefore;
WaveDPafter = StimulationParam.WaveDPafter;

SupersamplingRatio = round(SupersamplingRatio);
NumTrials = size(RawTrace,2);

% trialAlignPoint = round(abs(trialParam.startT)*SampleRate*SupersamplingRatio);

% intervalthreshold = 0.03; % 3 percent
% AlignErrorThreshold = 1; % data points

% construct filter
Fl = round(StimulationFreq*0.8);
Fn = SampleRate*SupersamplingRatio/2;
N = 2;
[Bartifact, Aartifact] = butter(N,Fl/Fn,'high'); % Artifact highpass filter
% 
Fl = StimulationPulsesFrequency * 1.5;
Fn = SampleRate*SupersamplingRatio/2;
N = 2;
[Bpulse, Apulse] = butter(N,Fl/Fn,'low'); % pulse low pass filter

% supersampling
if SupersamplingRatio > 1
    WaveDPbefore = round(WaveDPbefore * SupersamplingRatio);
    WaveDPafter = round(WaveDPafter * SupersamplingRatio);
    
    xdata = 1:size(RawTrace,1);
    xdata = xdata';
    xdata2 = 1:(1/SupersamplingRatio):size(RawTrace,1);
    xdata2 = xdata2';
    
    RawTrace = interp1(xdata,RawTrace,xdata2,'spline');
    StimStartEndTime = StimStartEndTime * SupersamplingRatio;
    
    %     AverageTrace = spline(xdata,AverageTrace,xdata2);
    
elseif SupersamplingRatio < 1
    error('SupersamplingRatio must be greater or equal than 1.')
end

% filter the trace
% HighPassAverage = filtfilt(Bartifact,Aartifact,AverageTrace);
HighPassdata = filtfilt(Bartifact,Aartifact,RawTrace);
% HighPassdata = HighPassdata * Gain;
HighPassdata = filtfilt(Bpulse,Apulse,HighPassdata);
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
if size(HighPassdata,2) > size(HighPassdata,1)
    HighPassdata = HighPassdata';
end


% Detection
% determin a threshold
% ThresholdData = zeros(size(HighPassdata,1),size(HighPassdata,2),2);
ValidNumberPulses = zeros(NumTrials,2);
for thisDirection = 1:2
    HighPassdata = -HighPassdata;
    for thisTrial = 1:NumTrials
        trialdata = HighPassdata(:,thisTrial);
        trialdata(1:(StimStartEndTime(1,thisTrial)-WaveDPbefore)) = 0;
        trialdata((StimStartEndTime(2,thisTrial)+WaveDPafter):end) = 0;
        trialdata = trialdata/max(trialdata)*dataAmpResolution;
        thresholdhi = max(trialdata);
        thresholdlow = thresholdhi/2;
        % search for artifact
        PulsesEqualFlag = 0;
%         for thisThreshold = threshold:-2:round(dataAmpResolution*0.5)
        thisThreshold = (thresholdhi + thresholdlow)/2;
%         movedirectionflag = 0;
        counter = 0;
        while 1
            counter = counter + 1;
            tmp = trialdata > thisThreshold;
            tmp = diff(tmp);
            tmp = tmp == -1;
            NumTrans = sum(tmp);
            % break search when then number of transgression equals number of
            % stimuluation pulses
            if NumTrans == numberofStimPulses
                PulsesEqualFlag = 1;
                break;
            elseif NumTrans < numberofStimPulses % move down
                    thresholdhi = thisThreshold;
            elseif NumTrans > numberofStimPulses % move up
                    thresholdlow = thisThreshold;
            end
            thisThreshold = (thresholdhi + thresholdlow)/2;
            if thresholdhi <= thresholdlow+1 % cannot find
                PulsesEqualFlag = 0;
                break
            end
        end
        if PulsesEqualFlag
            ValidNumberPulses(thisTrial,3-thisDirection) = 1;
        else
            
        end
    end
end

% Find stimulation exact locations, get waveform
ArtifactLocation = zeros(numberofStimPulses,NumTrials);
ARdata = RawTrace;
for thisTrial = 1:NumTrials
    if sum(ValidNumberPulses(thisTrial,:)) == 2
        data = HighPassdata(:,thisTrial);
        data(1:((StimStartEndTime(1,thisTrial) - WaveDPbefore))) = 0;
        data((StimStartEndTime(2,thisTrial) + WaveDPafter):end) = 0;
        loc = getArtifactTime(data,SampleRate,SupersamplingRatio,StimulationPulsesFrequency,stimulationwaveformwidth,numberofStimPulses);
    else
        loc = [];
    end
    if isempty(loc)
        ValidNumberPulses(thisTrial,:) = 0;
        predictedLoc = (0:(numberofStimPulses-1)) * 1/StimulationPulsesFrequency * SampleRate * SupersamplingRatio;
        predictedLoc = predictedLoc + StimStartEndTime(1,thisTrial);
        predictedLoc = predictedLoc + StimulationParam.stimulationwaveformwidth/1000 * SampleRate * SupersamplingRatio;
        predictedLoc = round(predictedLoc);
        loc = predictedLoc;
    end
    ArtifactLocation(:,thisTrial) = loc;
    wave = getArtifactWaveform(RawTrace(:,thisTrial),loc,numberofStimPulses,WaveDPbefore,WaveDPafter);
    
    % calculate rolling average and substract from the current pulse
    for thisLoc = 1:numberofStimPulses
        if thisLoc < AvgWindowSize/2
            pulseIdxStart = 1;
            pulseIdxEnd = 1 + AvgWindowSize - 1;
        elseif numberofStimPulses - thisLoc < AvgWindowSize/2
            pulseIdxStart = numberofStimPulses - AvgWindowSize + 1;
            pulseIdxEnd = numberofStimPulses;
        else
            pulseIdxStart = thisLoc - ceil(AvgWindowSize/2) + 1;
            pulseIdxEnd = thisLoc + floor(AvgWindowSize/2);
        end
        meanWave = mean(wave(:,pulseIdxStart:pulseIdxEnd),2);
        idx = (loc(thisLoc)-WaveDPbefore):(loc(thisLoc)+WaveDPafter-1);
        ARdata(idx,thisTrial) = ARdata(idx,thisTrial) - meanWave;
    end
end

if debug
    colors = lines(NumTrials);
    figure
    subplot(2,3,[1 2]),hold on
    plot(RawTrace)
    for thisTrial = 1:NumTrials
        idx = ArtifactLocation(:,thisTrial);
        plot(idx,RawTrace(idx,thisTrial),'ro')
    end
    subplot(2,3,3),hold on
    for thisTrial = 1:NumTrials
        wave = getArtifactWaveform(RawTrace(:,thisTrial),ArtifactLocation(:,thisTrial),numberofStimPulses,WaveDPbefore,WaveDPafter);
        plot(wave,'color',colors(thisTrial,:))
    end
    
    subplot(2,3,[4 5]),hold on
    plot(ARdata)
    for thisTrial = 1:NumTrials
        idx = ArtifactLocation(:,thisTrial);
        plot(idx,ARdata(idx,thisTrial),'ro')
    end
    
    subplot(2,3,6),hold on
    for thisTrial = 1:NumTrials
        wave = getArtifactWaveform(ARdata(:,thisTrial),ArtifactLocation(:,thisTrial),numberofStimPulses,WaveDPbefore,WaveDPafter);
        plot(wave,'color',colors(thisTrial,:))
    end
    input('continue?','s')
    close all
end

% down sample the raw trace
ARdata = downsample(ARdata,SupersamplingRatio);
ArtifactLocation = round(ArtifactLocation/SupersamplingRatio);
end

% supporting subroutines
function loc = getArtifactTime(trace,SampleRate,SupersamplingRatio,StimulationPulsesFrequency,stimulationwaveformwidth,numberofStimPulses)
% [~,tmp,~,~] = findpeaks(trace,'WidthReference','halfheight','MinPeakDistance',1/StimulationPulsesFrequency*SampleRate*SupersamplingRatio*0.95,'SortStr','descend','Annotate','extents');
locvalidFlag = ones(2,1);
sortflag = 1;
loc1 = peakseek(trace,1/StimulationPulsesFrequency*SampleRate*SupersamplingRatio*0.95,stimulationwaveformwidth,sortflag);
if numel(loc1) < numberofStimPulses
    locvalidFlag(1) = 0;
else
    loc1 = loc1(1:numberofStimPulses);
    loc1 = sort(loc1,'ascend'); % positive peaks
end
% [~,loc2,~,~] = findpeaks(-trace,'WidthReference','halfheight','MinPeakDistance',1/StimulationPulsesFrequency*SampleRate*SupersamplingRatio*0.95,'SortStr','descend','Annotate','extents');
loc2 = peakseek(-trace,1/StimulationPulsesFrequency*SampleRate*SupersamplingRatio*0.95,stimulationwaveformwidth,sortflag);
if numel(loc2) < numberofStimPulses
    locvalidFlag(2) = 0;
else
    loc2 = loc2(1:numberofStimPulses);
    loc2 = sort(loc2,'ascend'); % negative peaks
end

if sum(locvalidFlag) == 0
    loc = [];
    return
elseif sum(locvalidFlag) == 1
    if locvalidFlag(1)
        loc = loc1;
    else
        loc = loc2;
    end
    return
end

loc = loc1;
for thisLoc = 1:numberofStimPulses
    if loc1(thisLoc) < loc2(thisLoc)
        idx = loc1(thisLoc) : loc2(thisLoc);
        tmp = trace(idx);
        [~,idx] = min(abs(tmp));
        loc(thisLoc) = loc1(thisLoc) + idx;
    else
        idx = loc2(thisLoc) : loc1(thisLoc);
        tmp = trace(idx);
        [~,idx] = min(abs(tmp));
        loc(thisLoc) = loc2(thisLoc) + idx;
    end
end
end

function wave = getArtifactWaveform(trace,loc,numberofStimPulses,WaveDPbefore,WaveDPafter)
wave = zeros(WaveDPbefore + WaveDPafter,numberofStimPulses);
for thisLoc = 1:numberofStimPulses
    idx = (loc(thisLoc)-WaveDPbefore):(loc(thisLoc)+WaveDPafter-1);
    tmp = trace(idx);
    tmp = tmp - mean([tmp(1);tmp(end)]);
%     tmp = tmp - mean([tmp(1:round(WaveDPbefore/5));tmp(round(end-WaveDPafter/5):end)]);
%     tmp = tmp - tmp(WaveDPbefore+1);
    wave(:,thisLoc) = tmp;
end
end