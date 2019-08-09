function [ARdata,ArtefactLocation] = VnCRemoveArtifactsSubroutine17(RawTrace,StimulationParam,StimStartTime,SampleRate)
% supersampling, fine alignment and SALPA on segement

% inputs
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
NumTrials = size(RawTrace,2);

tau = 0.0015; % sec
interpolationLengthFactor = 2;


% parameters
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
stimulationwaveformwidth = StimulationParam.stimulationwaveformwidth*SampleRate;
StimulationPulseInterval = 1/StimulationPulsesFrequency;

% padding
NumDPpadding = 100;
RawTrace = [zeros(NumDPpadding,NumTrials);RawTrace;zeros(NumDPpadding,NumTrials)];
NumDP = size(RawTrace,1);

% for each trial, remove the piece of data at the moment of stimulation and
% perform salpa at inter stimulation interval, then perform linear
% interpolation to bridge the gap


% get predicted stimulation locations
StimStartTime = StimStartTime + NumDPpadding + 1;
NumPulsesBefore = floor(StimStartTime/SampleRate / StimulationPulseInterval);
NumPulsesAfter = floor((NumDP - StimStartTime)/SampleRate / StimulationPulseInterval);
numberofStimPulses = NumPulsesBefore + NumPulsesAfter + 1;

FirstPulseTimepoint = round(StimStartTime - NumPulsesBefore * StimulationPulseInterval* SampleRate);

predictedLoc = (0:(numberofStimPulses-1)) * StimulationPulseInterval * SampleRate;
predictedLoc = predictedLoc + FirstPulseTimepoint;
loc0 = round(predictedLoc);

WaveDPbefore = 0;
WaveDPafter = stimulationwaveformwidth*interpolationLengthFactor;

% determine AR interpolation length as 1.5 times pulse width
% disp(['AR interpolation length = ',num2str((WaveDPbefore+WaveDPafter)/SampleRate), 'ms'])

Interpolateddata = RawTrace;
% clearvars RawTrace

% calculate artefect locations
ArtefactLocation = zeros(numberofStimPulses,NumTrials);
idx = round(-WaveDPbefore) : round(WaveDPafter);

numinterpDP = numel(idx);
if size(idx,2)>size(idx,1) % column
    idx = idx';
end
idx = repmat(idx,1,numel(loc0));
if size(loc0,1)>size(loc0,2) % row
    loc0 = loc0';
end
loc = repmat(loc0,numinterpDP,1);
idx = idx + loc;

for thisTrial = 1:NumTrials
    ArtefactLocation(:,thisTrial) = loc0;
    
    y = Interpolateddata(:,thisTrial);
    
    % perform Salpa here
%     y = salpa(y,'tau',round(tau*SampleRate));
%     y(isnan(y)) = 0;
    
%     yy = y(1:idx(1,1));
%     yy = salpa(yy,'tau',round(tau*SampleRate));
%     y(1:idx(1,1)) = yy;
%     figure
    for thisArtefect = 2:(numberofStimPulses-2)
        % perform salpa for each piece of data
        yidx = (idx(end,thisArtefect)+1):(idx(1,thisArtefect+1)-1);
        yy = y(yidx);
%         subplot(1,3,1)
%         ylim([-2500 2500])
%         hold on
%         plot(yy)

        % extrapolation on both end
        %         Csize = numel(yy);
        %         yy = [flipud(yy);yy;flipud(yy)];
        padlength = round(tau*1.2*SampleRate);
        yy = [yy(1:padlength)-yy(padlength) + yy(1)+yy(1)-yy(2); yy; yy(end-padlength+1:end)-yy(end-padlength+1)+yy(end)+yy(end)-yy(end-1)];
        
%         subplot(1,3,2)
%         ylim([-2500 2500])
%         hold on
%         plot(yy)
        
        yy = salpa(yy,'tau',round(tau*SampleRate));
        yy = yy(padlength+1:padlength+numel(yidx));
%         subplot(1,3,3)
%         ylim([-2500 2500])
%         hold on
%         plot(yy)
%         xx = (1:numel(yy))';
%         f2 = fit(xx,yy,'exp1');
%         yy = yy - f2(xx);
        
        %         yy = yy((Csize+1):2*Csize);
        y(yidx) = yy;
    end
    y(isnan(y)) = 0;
    % perform linear interpolation
    x = 1:numel(y);
    xx = x;
    x(idx) = [];
    y(idx) = [];
    Interpolateddata(:,thisTrial) = interp1(x,y,xx,'linear');
end
ARdata = Interpolateddata;
clearvars Interpolateddata;

% remove padding
if NumDPpadding > 0
    ARdata = ARdata((NumDPpadding+1):(end-NumDPpadding),:);
    ArtefactLocation = ArtefactLocation - NumDPpadding;
end
end