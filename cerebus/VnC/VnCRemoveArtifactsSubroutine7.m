function [ARdata,ArtifactLocation,Interpolateddata,WaveDPbefore,WaveDPafter] = VnCRemoveArtifactsSubroutine7(RawTrace,StimulationParam,StimStartTime,SampleRate)

tau = 0.001; % sec
SSrate = 1; % super sampling rate
cycleLimit = 100;
% ThresholdGain = 3;

alpha = 0.1;

% inputs
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
% parameters
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
numberofStimPulses = StimulationParam.numberofStimPulses;
WaveDPbefore = StimulationParam.WaveDPbefore;
WaveDPafter = StimulationParam.WaveDPafter;
% WaveformFreq = 1/StimulationParam.stimulationwaveformwidth;

% Fl=WaveformFreq/10;
Fl = 200;
Fn = SampleRate/2;
N = 4;
[BWave, AWave] = butter(N,Fl/Fn,'high'); % Waveform hi pass
HipassTrace = filtfilt(BWave, AWave,RawTrace);

% Fl=StimulationPulsesFrequency*2;
% Fn = SampleRate/2;
% N = 4;
% [BPulse, APulse] = butter(N,Fl/Fn,'low'); % Waveform hi pass
% HipassTrace = filtfilt(BWave, AWave,RawTrace);
HipassTrace = abs(HipassTrace);
% HipassTrace = filtfilt(BPulse, APulse,HipassTrace);


NumTrials = size(RawTrace,2);

predictedLoc = (0:(numberofStimPulses-1)) * 1/StimulationPulsesFrequency * SampleRate;
predictedLoc = predictedLoc + StimStartTime;
predictedLoc = predictedLoc + StimulationParam.stimulationwaveformwidth*0.95 * SampleRate;
% predictedLoc = floor(predictedLoc);

loc0 = round(predictedLoc * SSrate);

% Find stimulation exact locations, get waveform
numSamplePoints = size(HipassTrace,1);
x = 1:numSamplePoints;
xx = 1:1/SSrate:numSamplePoints;
HipassTrace = spline(x,HipassTrace',xx)';

RawTrace = spline(x,RawTrace',xx)';

Threshold = median(abs(reshape(HipassTrace,[],1))/0.6745,1);

% search for interpolation parameters
cycle = 1;
terminateFlag = 0;
% GoFlag = zeros(numberofStimPulses,1);
while cycle < cycleLimit && terminateFlag == 0
% construct interpolation idx matrix
idx = round(-WaveDPbefore* SSrate : WaveDPafter* SSrate);

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


% ARdata = HipassTrace;
Interpolateddata = HipassTrace;

for thisTrial = 1:NumTrials
%     ArtifactLocation(:,thisTrial) = loc0;
    % linear interpolate
    y = HipassTrace(:,thisTrial);
    x = 1:numel(y);
    xx = x;
    x(idx) = [];
    y(idx) = [];
    Interpolateddata(:,thisTrial) = interp1(x,y,xx,'linear');
end

idx2 = Interpolateddata(idx,:)>Threshold* 3;
if sum(sum(idx2 == 0)) / numel(idx2) >= 0.999
    terminateFlag = 1;
else
    WaveDPbefore = WaveDPbefore * (1+alpha);
    WaveDPafter = WaveDPafter * (1+alpha);
end
if (WaveDPbefore+WaveDPafter) /SampleRate > (1/StimulationPulsesFrequency)/2
    break
end
cycle = cycle + 1;
disp(['cycle=',num2str(cycle),'; Time =',num2str((WaveDPbefore+WaveDPafter)/30),'ms;'])
end % end of interpolation parameter search

Interpolateddata = RawTrace;
ArtifactLocation = zeros(numberofStimPulses,NumTrials);
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
% if true
%     figure
%     hold on
%     plot(mean(Interpolateddata,2))
%     plot(loc0, max(mean(Interpolateddata,2))*ones(size(loc0)),'rx')
%     plot(idx,max(mean(Interpolateddata,2))*ones(size(idx)),'r.')
% end
    
for thisTrial = 1:NumTrials
    % use salpa algorithm to remove exponential decay.
    yy = Interpolateddata(:,thisTrial);
    y = salpa(yy,'tau',round(tau*SampleRate*SSrate));
    
    nanidx = isnan(y);
    y(nanidx) = RawTrace(nanidx,thisTrial);
    ARdata(:,thisTrial) = y;
end
% if true
%     figure
%     hold on
%     plot(mean(ARdata,2))
%     plot(loc0, max(mean(ARdata,2))*ones(size(loc0)),'rx')
%     plot(idx,80*ones(size(idx)),'r.')
% end

% downsample
ARdata = downsample(ARdata,SSrate);
Interpolateddata = downsample(Interpolateddata,SSrate);

ArtifactLocation = round(ArtifactLocation/SSrate);
end