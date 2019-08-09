function [ARdata,ArtifactLocation] = VnCRemoveArtifactsSubroutine5(RawTrace,StimulationParam,StimStartTime,SampleRate)

tau = 0.001; % sec

% inputs
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
% parameters
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
numberofStimPulses = StimulationParam.numberofStimPulses;
WaveDPbefore = StimulationParam.WaveDPbefore;
WaveDPafter = StimulationParam.WaveDPafter;

% build a 10Hz high pass filter
% Fn = SampleRate/2;
% Fhp= 10;
% N  = 4;    % filter order
% [Bhp, Ahp] = butter(N, Fhp/Fn,'high'); % high pass
% RawTrace = filtfilt(Bhp, Ahp,RawTrace);

NumTrials = size(RawTrace,2);

predictedLoc = (0:(numberofStimPulses-1)) * 1/StimulationPulsesFrequency * SampleRate;
predictedLoc = predictedLoc + StimStartTime;
predictedLoc = predictedLoc + StimulationParam.stimulationwaveformwidth * SampleRate;
predictedLoc = floor(predictedLoc);
loc0 = predictedLoc;
% construct interpolation idx matrix
idx = -WaveDPbefore : WaveDPafter;
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

% Find stimulation exact locations, get waveform
ArtifactLocation = zeros(numberofStimPulses,NumTrials);
ARdata = RawTrace;
Interpolateddata = RawTrace;

for thisTrial = 1:NumTrials
    ArtifactLocation(:,thisTrial) = loc0;
    % linear interpolate
    y = RawTrace(:,thisTrial);
    x = 1:numel(y);
    xx = x;
    x(idx) = [];
    y(idx) = [];
    Interpolateddata(:,thisTrial) = interp1(x,y,xx,'linear');
end

% inspect interpolation result
if true
    figure
    hold on
    plot(Interpolateddata)
    plot(loc0, max(max(Interpolateddata))*ones(size(loc0)),'rx')
end
    
for thisTrial = 1:NumTrials
    % use salpa algorithm to remove exponential decay.
    y = salpa(yy,'tau',round(tau*SampleRate));
    
    nanidx = isnan(y);
    y(nanidx) = RawTrace(nanidx,thisTrial);
    ARdata(:,thisTrial) = y;
    
%     if thisTrial == 3
%         figure
%         hold on
%         plot(RawTrace(:,3))
%         plot(yy)
%         plot(y)
%     end
end
ArtifactLocation = round(ArtifactLocation) + WaveDPafter - WaveDPbefore;
end