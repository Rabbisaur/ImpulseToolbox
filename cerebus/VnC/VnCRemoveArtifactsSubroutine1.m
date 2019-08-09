function [ARtrace, loc] = VnCRemoveArtifactsSubroutine(AlignedMeanTrace, RawTrace,StimulationParam,SampleRate,SupersamplingRatio,debug)

% RawTrace should be in the format of time by trials

NumTrials = size(RawTrace,2);

WaveDPbefore = StimulationParam.WaveDPbefore;
WaveDPafter = StimulationParam.WaveDPafter;
numberofStimPulses = StimulationParam.numberofStimPulses;

% supersampling
if SupersamplingRatio ~= 1
    WaveDPbefore = round(WaveDPbefore * SupersamplingRatio);
    WaveDPafter = round(WaveDPafter * SupersamplingRatio);
    
    xdata = 1:size(RawTrace,1);
    xdata2 = 1:(1/SupersamplingRatio):size(RawTrace,1);
    
    RawTrace = spline(xdata,RawTrace',xdata2);
    RawTrace = RawTrace';
    
    
    AlignedMeanTrace = spline(xdata,AlignedMeanTrace,xdata2);
    
end

% AlignedMeanTrace = mean(RawTrace,2); % calculate the mean trace as trial average;

StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
StimulationFreq = 1/(StimulationParam.stimulationwaveformwidth/1000);



Fl = StimulationFreq/2;
Fn = SampleRate*SupersamplingRatio/2;
N = 2;
[BMUAe, AMUAe] = butter(N,Fl/Fn,'high'); % MUAe low pass



tmpdata2 = filtfilt(BMUAe,AMUAe,AlignedMeanTrace);
tmpdata2 = tmpdata2 - mean(tmpdata2);

% readjust peak location
[~,loc1,~,~] = findpeaks(tmpdata2,'WidthReference','halfheight','MinPeakDistance',1/StimulationPulsesFrequency*SampleRate*SupersamplingRatio*0.95,'SortStr','descend','Annotate','extents');
loc1 = loc1(1:numberofStimPulses);
loc = sort(loc1,'ascend'); % positive peaks
% [~,loc2,~,~] = findpeaks(-tmpdata2,'WidthReference','halfheight','MinPeakDistance',1/StimulationPulsesFrequency*SampleRate*SupersamplingRatio*0.95,'SortStr','descend','Annotate','extents');
% loc2 = loc2(1:numberofStimPulses);
% loc2 = sort(loc2,'ascend'); % negative peaks
% 
% loc = loc1;
% for thisLoc = 1:numberofStimPulses
%     if loc1(thisLoc) < loc2(thisLoc)
%         idx = loc1(thisLoc) : loc2(thisLoc);
%         %     [~,idx] = min(abs(trace(idx)-mean([trace(idx(1)),trace(idx(end))])));
%         [~,idx] = min(abs(tmpdata2(idx)));
%         loc(thisLoc) = loc1(thisLoc) + idx;
%     else
%         idx = loc2(thisLoc) : loc1(thisLoc);
%         %     [~,idx] = min(abs(trace(idx)-mean([trace(idx(1)),trace(idx(end))])));
%         [~,idx] = min(abs(tmpdata2(idx)));
%         loc(thisLoc) = loc2(thisLoc) + idx;
%     end
% end


% get waveforms
numLoc = numel(loc);

if debug
    wave = zeros(WaveDPbefore+WaveDPafter,numLoc);
    
    for thisLoc = 1:numLoc
        idx = (loc(thisLoc)-WaveDPbefore):(loc(thisLoc)+WaveDPafter-1);
        tmp = tmpdata2(idx);
        tmp = tmp - mean([tmp(1),tmp(end)]);
        wave(:,thisLoc) = tmp;
    end
end

% get waveforms and calculate mean from the rawtrace
waveOri = zeros(WaveDPbefore+WaveDPafter,numLoc,NumTrials);
meanWaveOri = zeros(WaveDPbefore+WaveDPafter,NumTrials);
for thisTrial = 1:NumTrials
    for thisLoc = 1:numLoc
        idx = (loc(thisLoc)-WaveDPbefore):(loc(thisLoc)+WaveDPafter-1);
        tmp = RawTrace(idx,thisTrial);
        tmp = tmp - mean([tmp(1),tmp(end)]);
        waveOri(:,thisLoc,thisTrial) = tmp;
    end
    meanWaveOri(:,thisTrial) = mean(waveOri(:,:,thisTrial),2);
end
clearvars waveOri

% artifact removal by substracting average at each artifact pulse
ARtrace = RawTrace;
if ~debug
    clearvars RawTrace
end
for thisTrial = 1:NumTrials
    for thisLoc = 1:numLoc
        idx = (loc(thisLoc)-WaveDPbefore):(loc(thisLoc)+WaveDPafter-1);
        ARtrace(idx,thisTrial) = ARtrace(idx,thisTrial) - meanWaveOri(:,thisTrial);
    end
end

if debug
    Residuals = zeros(WaveDPbefore+WaveDPafter,numLoc,NumTrials);
    % meanResiduals = zeros(WaveDPbefore+WaveDPafter,NumTrials);
    for thisTrial = 1:NumTrials
        for thisLoc = 1:numLoc
            idx = (loc(thisLoc)-WaveDPbefore):(loc(thisLoc)+WaveDPafter-1);
            tmp = ARtrace(idx,thisTrial);
            tmp = tmp - mean([tmp(1),tmp(end)]);
            Residuals(:,thisLoc,thisTrial) = tmp;
        end
    end
    meanResiduals = mean(Residuals,3);
end

if debug
    figure
    % high passed average trace
    subplot(3,4,[1:3]),hold on
    plot(tmpdata2)
    plot(loc, tmpdata2(loc), 'ro')
    
    % waveforms on high passed average trace
    subplot(3,4,4),hold on
    plot(wave)
    waveMax = max(max(wave));
    waveMin = min(min(wave));
    waveylimGain = 1.5;
    ylim([waveMin*waveylimGain waveMax*waveylimGain])
    
    % raw trace
    subplot(3,4,[5:7]),hold on
    plot(RawTrace)
    plot(loc, RawTrace(loc,:), 'ro')
    subplot(3,4,8),hold on
    plot(meanWaveOri)
    ylim([waveMin*waveylimGain waveMax*waveylimGain])
    
    % artifact removed trace
    subplot(3,4,[9:11]),hold on
    plot(ARtrace)
    plot(loc, ARtrace(loc,:), 'ro')
    subplot(3,4,12),hold on
    plot(meanResiduals)
    ylim([waveMin*waveylimGain waveMax*waveylimGain])
    input('Continue?')
end

% down sampling
ARtrace = downsample(ARtrace,SupersamplingRatio);
loc = round(loc/SupersamplingRatio);
end