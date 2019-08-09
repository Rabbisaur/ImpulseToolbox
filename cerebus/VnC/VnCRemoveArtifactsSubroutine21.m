function [ARdata,ArtifactLocation] = VnCRemoveArtifactsSubroutine2(RawTrace,StimulationParam,trialParam,StimStartEndTime,SampleRate,SupersamplingRatio,debug)

% inputs

% RawTrace = AMF';
% SampleRate = 30000;
% SupersamplingRatio = 16;

% StimulationParam.numberofStimPulses = 50;
% StimulationParam.PulsesFrequency = 300; %Hz
% StimulationParam.stimulationwaveformwidth = 0.4; % ms
% StimulationParam.WaveDPbefore = 12; %number of datapoints, relative to peak of artifact waveform
% StimulationParam.WaveDPafter = 12; %number of datapoints, relative to peak of artifact waveform

% debug = 1;
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
% parameters
Gain = 100;
% AverageTrace = mean(RawTrace,2);
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
StimulationFreq = 1/(StimulationParam.stimulationwaveformwidth/1000);
numberofStimPulses = StimulationParam.numberofStimPulses;
stimulationwaveformwidth = StimulationParam.stimulationwaveformwidth * SampleRate / 1000 * SupersamplingRatio;

WaveDPbefore = StimulationParam.WaveDPbefore;
WaveDPafter = StimulationParam.WaveDPafter;

SupersamplingRatio = round(SupersamplingRatio);
NumTrials = size(RawTrace,2);

trialAlignPoint = round(abs(trialParam.startT)*SampleRate*SupersamplingRatio);

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
HighPassdata = HighPassdata * Gain;
HighPassdata = filtfilt(Bpulse,Apulse,HighPassdata);
% HighPassdata = RawTrace;
% HighPassdata = filtfilt(Bpulse,Apulse,RawTrace);
% HighPassdata = HighPassdata * Gain;
% HighPassdata = filtfilt(Bartifact,Aartifact,HighPassdata);
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
if size(HighPassdata,2) > size(HighPassdata,1)
    HighPassdata = HighPassdata';
end

for thisTrial = 1:NumTrials
    HighPassdata(:,thisTrial) = HighPassdata(:,thisTrial) - mean(HighPassdata(:,thisTrial));
end


% find the peak from the AverageTrace
% loc = getArtifactTime(HighPassAverage,SampleRate,SupersamplingRatio,StimulationPulsesFrequency,numberofStimPulses);
% trialAlignPoint = mean(loc);
% AvgLoc = loc;
%
% % determine if the trials are aligned good enough
% % check for variation of the peaks
% artifactintervals = diff(loc) / SampleRate / SupersamplingRatio;
% if max(artifactintervals) - min(artifactintervals) > 1/StimulationPulsesFrequency * intervalthreshold
%     disp('The artifact pulses on the average trace has too much variation.')
%     % get artifact waveforms from the average
%     wave = getArtifactWaveform(HighPassAverage,loc,numberofStimPulses,WaveDPbefore,WaveDPafter);
%     figure
%     subplot(1,3,[1 2])
%     hold on
%     plot(HighPassAverage)
%     plot(loc,HighPassAverage(loc),'ro')
%     subplot(1,3,3)
%     plot(wave)
%     s = input('continue? (Y/N)','s');
%     if s == 'Y' || 'y'
%         % do nothing
%     else
%         return
%     end
% end
% cut off 1/3 of the artifact pulse flanking the align points
% idx = zeros(numberofStimPulses,NumTrials);
% for thisTrial = 1:NumTrials
%     trace = HighPassdata(:,thisTrial);
%     wave = getArtifactWaveform(trace,loc,numberofStimPulses,round(stimulationwaveformwidth/6) * SupersamplingRatio,round(stimulationwaveformwidth/6)*SupersamplingRatio);
%     [~,idx(:,thisTrial)] = min(abs(wave),[],1);
% end
% realignFlag = 0;
% AlignError = (max(max(idx)) - min(min(idx)))/SupersamplingRatio;
% if  AlignError > AlignErrorThreshold
% %     disp(['Trials are not very well aligned (error = ',num2str(AlignError),'DPs), performing realign'])
%     realignFlag = 1;
% end

% realign the trials using RawTrace if needed
% if realignFlag
% Realign
meanPeakPoint = zeros(NumTrials,1);
if debug
    figure
    colors = lines(NumTrials);
end

for thisTrial = 1:NumTrials
    data = HighPassdata(:,thisTrial);
    data(1:(StimStartEndTime(1,thisTrial)-WaveDPbefore)) = 0;
    data((StimStartEndTime(2,thisTrial)+WaveDPafter):end) = 0;
    predictedLoc = StimStartEndTime(1,thisTrial) : 1/StimulationPulsesFrequency * SampleRate * SupersamplingRatio: StimStartEndTime(1,thisTrial) + 1/StimulationPulsesFrequency * SampleRate * SupersamplingRatio *numberofStimPulses;
    predictedLoc = predictedLoc + StimulationParam.stimulationwaveformwidth/1000 * SampleRate * SupersamplingRatio;
    predictedLoc = round(predictedLoc);
    loc = getArtifactTime(data,SampleRate,SupersamplingRatio,StimulationPulsesFrequency,stimulationwaveformwidth,numberofStimPulses);
    if isempty(loc)
        loc = predictedLoc;
    end
    meanPeakPoint(thisTrial) = loc(1);
end
AlignmentOffset = round(meanPeakPoint - trialAlignPoint);


for thisTrial = 1:NumTrials
    HighPassdata(:,thisTrial) = circshift(HighPassdata(:,thisTrial),-AlignmentOffset(thisTrial));
    HighPassdata(1:abs(AlignmentOffset(thisTrial)),thisTrial) = HighPassdata(abs(AlignmentOffset(thisTrial))+1,thisTrial);
    HighPassdata((end-abs(AlignmentOffset(thisTrial))) : end,thisTrial) = HighPassdata(end-abs(AlignmentOffset(thisTrial))-1,thisTrial);
    
    RawTrace(:,thisTrial) = circshift(RawTrace(:,thisTrial),-AlignmentOffset(thisTrial));
    RawTrace(1:abs(AlignmentOffset(thisTrial)),thisTrial) = RawTrace(abs(AlignmentOffset(thisTrial))+1,thisTrial);
    RawTrace((end-abs(AlignmentOffset(thisTrial))) : end,thisTrial) = RawTrace(end-abs(AlignmentOffset(thisTrial))-1,thisTrial);
    
%     if debug
%         subplot(1,3,2),hold on
%         wave = getArtifactWaveform(RawTrace(:,thisTrial),AvgLoc,numberofStimPulses,WaveDPbefore,WaveDPafter);
%         plot(wave,'color',colors(thisTrial,:))
%     end
    StimStartEndTime(:,thisTrial) = StimStartEndTime(:,thisTrial) - AlignmentOffset(thisTrial);
end
% end

% get average artifact % substract artifact from the raw trace
ArtifactLocation = zeros(numberofStimPulses,NumTrials);
ARdata = RawTrace;
for thisTrial = 1:NumTrials
    data = HighPassdata(:,thisTrial);
    data(1:(StimStartEndTime(1,thisTrial)-WaveDPbefore)) = 0;
    data((StimStartEndTime(2,thisTrial)+WaveDPafter):end) = 0;
    
    
    
    loc = getArtifactTime(data,SampleRate,SupersamplingRatio,StimulationPulsesFrequency,stimulationwaveformwidth,numberofStimPulses);
    if isempty(loc)
        predictedLoc = StimStartEndTime(1,thisTrial) : 1/StimulationPulsesFrequency * SampleRate * SupersamplingRatio: StimStartEndTime(1,thisTrial) + 1/StimulationPulsesFrequency * SampleRate * SupersamplingRatio *(numberofStimPulses-1);
        predictedLoc = predictedLoc + StimulationParam.stimulationwaveformwidth/1000 * SampleRate * SupersamplingRatio;
        predictedLoc = round(predictedLoc);
        loc = predictedLoc;
    end
    ArtifactLocation(:,thisTrial) = loc;
    wave = getArtifactWaveform(RawTrace(:,thisTrial),loc,numberofStimPulses,WaveDPbefore,WaveDPafter);
    if debug
        subplot(1,3,3),hold on
        plot(wave,'color',colors(thisTrial,:))
    end
    meanWave = mean(wave,2);
    for thisLoc = 1:numberofStimPulses
        idx = (loc(thisLoc)-WaveDPbefore):(loc(thisLoc)+WaveDPafter-1);
        ARdata(idx,thisTrial) = ARdata(idx,thisTrial) - meanWave;
    end
end

if debug
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
        %     [~,idx] = min(abs(trace(idx)-mean([trace(idx(1)),trace(idx(end))])));
        tmp = trace(idx);
%         tmp = tmp - mean(tmp);
        [~,idx] = min(abs(tmp));
        loc(thisLoc) = loc1(thisLoc) + idx;
    else
        idx = loc2(thisLoc) : loc1(thisLoc);
        %     [~,idx] = min(abs(trace(idx)-mean([trace(idx(1)),trace(idx(end))])));
        tmp = trace(idx);
%         tmp = tmp - mean(tmp);
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
    tmp = tmp - mean([tmp(1),tmp(end)]);
    wave(:,thisLoc) = tmp;
end
end