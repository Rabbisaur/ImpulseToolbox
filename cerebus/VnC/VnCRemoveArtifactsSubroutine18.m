function [ARdata,ArtefactLocation,Interpolateddata] = VnCRemoveArtifactsSubroutine18(RawTrace,StimulationParam,StimStartTime,SampleRate)
% SALPA on full trial, 1ms tau

% inputs
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
NumTrials = size(RawTrace,2);


interpolationLengthFactor = [2.2 2.3];


% parameters
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
stimulationwaveformwidth = StimulationParam.stimulationwaveformwidth*SampleRate;
StimulationPulseInterval = 1/StimulationPulsesFrequency;
tau = 0.001; % sec
% tau = 1/StimulationPulsesFrequency; % sec

% padding
NumDPpadding = 0;
if NumDPpadding > 0
    RawTrace = [zeros(NumDPpadding,NumTrials);RawTrace;zeros(NumDPpadding,NumTrials)];
end
NumDP = size(RawTrace,1);


% get predicted stimulation locations
StimStartTime = StimStartTime + NumDPpadding;
NumPulsesBefore = floor(StimStartTime/SampleRate / StimulationPulseInterval);
NumPulsesAfter = floor((NumDP - StimStartTime)/SampleRate / StimulationPulseInterval);
numberofStimPulses = NumPulsesBefore + NumPulsesAfter + 1;

FirstPulseTimepoint = floor(StimStartTime - NumPulsesBefore * StimulationPulseInterval* SampleRate);

predictedLoc = (0:(numberofStimPulses-1)) * StimulationPulseInterval * SampleRate;
predictedLoc = predictedLoc + FirstPulseTimepoint;
loc0 = round(predictedLoc);

WaveDPbefore = 0;
WaveDPafter1 = stimulationwaveformwidth*interpolationLengthFactor(1);
WaveDPafter2 = stimulationwaveformwidth*interpolationLengthFactor(2);

% determine AR interpolation length as 1.5 times pulse width
% disp(['AR interpolation length = ',num2str((WaveDPbefore+WaveDPafter)/SampleRate), 'ms'])

Interpolateddata = RawTrace;
ARdata = Interpolateddata;

% calculate artefect locations
ArtefactLocation = zeros(numberofStimPulses,NumTrials);
idx1 = round(-WaveDPbefore) : round(WaveDPafter1);

% numinterpDP = numel(idx1);
if size(idx1,2)>size(idx1,1) % column
    idx1 = idx1';
end
if size(loc0,1)>size(loc0,2) % row
    loc0 = loc0';
end
idx1 = idx1 + loc0;

idx2 = round(-WaveDPbefore) : round(WaveDPafter2);

% numinterpDP = numel(idx2);
if size(idx2,2)>size(idx2,1) % column
    idx2 = idx2';
end
if size(loc0,1)>size(loc0,2) % row
    loc0 = loc0';
end
idx2 = idx2 + loc0;

for thisTrial = 1:NumTrials
    ArtefactLocation(:,thisTrial) = loc0;
    
    y = RawTrace(:,thisTrial);
    
    % perform linear interpolation 1
    x = (1:numel(y))';
    xx = x;
    x(idx1) = [];
    y(idx1) = [];
    y = interp1(x,y,xx,'linear');
    
    Interpolateddata(:,thisTrial)=y;
    % perform Salpa here
    y = salpa(y,'tau',round(tau*SampleRate));
    y(isnan(y)) = 0;
    
%     yy = y(1:idx(1,1));
%     yy = salpa(yy,'tau',round(tau*SampleRate));
%     y(1:idx(1,1)) = yy;
%     figure
%     for thisArtefect = 2:(numberofStimPulses-2)
%         % perform salpa for each piece of data
%         yidx = (idx(end,thisArtefect)+1):(idx(1,thisArtefect+1)-1);
%         yy = y(yidx);
% %         subplot(1,3,1)
% %         ylim([-2500 2500])
% %         hold on
% %         plot(yy)
% 
%         % extrapolation on both end
%         %         Csize = numel(yy);
%         %         yy = [flipud(yy);yy;flipud(yy)];
%         padlength = round(tau*1.2*SampleRate);
%         yy = [yy(1:padlength)-yy(padlength) + yy(1)+yy(1)-yy(2); yy; yy(end-padlength+1:end)-yy(end-padlength+1)+yy(end)+yy(end)-yy(end-1)];
%         
% %         subplot(1,3,2)
% %         ylim([-2500 2500])
% %         hold on
% %         plot(yy)
%         
%         yy = salpa(yy,'tau',round(tau*SampleRate));
%         yy = yy(padlength+1:padlength+numel(yidx));
% %         subplot(1,3,3)
% %         ylim([-2500 2500])
% %         hold on
% %         plot(yy)
% %         xx = (1:numel(yy))';
% %         f2 = fit(xx,yy,'exp1');
% %         yy = yy - f2(xx);
%         
%         %         yy = yy((Csize+1):2*Csize);
%         y(yidx) = yy;
%     end
%     y(isnan(y)) = 0;
    % perform linear interpolation 2
    x = (1:numel(y))';
    xx = x;
    x(idx2) = [];
    y(idx2) = [];
    ARdata(:,thisTrial) = interp1(x,y,xx,'linear');
end


% remove padding
if NumDPpadding > 0
    ARdata = ARdata((NumDPpadding+1):(end-NumDPpadding),:);
    ArtefactLocation = ArtefactLocation - NumDPpadding;
end
end