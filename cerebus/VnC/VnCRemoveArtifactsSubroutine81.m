function [ARdata,ArtifactLocation,Interpolateddata,WaveDPbefore,WaveDPafter] = VnCRemoveArtifactsSubroutine8(RawTrace,StimulationParam,StimStartTime,SampleRate)

tau = 0.001; % sec
SSrate = 1; % super sampling rate
% cycleLimit = 100;
% ThresholdGain = 3;

% alpha = 0.1;

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
Fl = 100;
Fn = SampleRate/2;
N = 4;
[BWave, AWave] = butter(N,Fl/Fn,'high'); % Waveform hi pass
HipassTrace = filtfilt(BWave, AWave,RawTrace);

% Fl=StimulationPulsesFrequency*2;
% Fn = SampleRate/2;
% N = 4;
% [BPulse, APulse] = butter(N,Fl/Fn,'low'); % Waveform hi pass
% HipassTrace = filtfilt(BWave, AWave,RawTrace);
% HipassTrace = abs(HipassTrace);
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

StimWidth = 1/StimulationPulsesFrequency/2;
idx = (0:round(StimWidth*SampleRate*SSrate-1))';
WaveLength = numel(idx);
idx = idx - floor(max(idx)/2);
idx1 = repmat(idx,1,numberofStimPulses);
idx1 = idx1 + repmat(loc0,numel(idx),1);
d = zeros(WaveLength,NumTrials);
for thisTrial = 1:NumTrials
    Waveform = HipassTrace(idx1,thisTrial);
    Waveform = reshape(Waveform,numel(idx),numberofStimPulses);
%     figure
%     hold on
    %     plot(Waveform,'black')
%     plot(mean(Waveform,2),'red','linewidth',2)
    meanWaveform = mean(Waveform,2);
    %     K = (meanWaveform(2:end) - meanWaveform(1))./(1:(numel(meanWaveform)-1))'*;
    %     plot(K,'blue','linewidth',2)
%     counter = 0;    
    for InterpWidth = 2:2:(WaveLength-2)
%         counter = counter + 1;
        interpIdx = floor(WaveLength / 2 + (0:(InterpWidth-1)) - (InterpWidth/2));
        x = 1:WaveLength;
        x(interpIdx) = [];
        tmp = meanWaveform;
        tmp(interpIdx) = [];
        xx = 1:WaveLength;
        tmp = interp1(x,tmp,xx);
        d(InterpWidth,thisTrial) = sum((tmp' - mean(meanWaveform)).^2);
    end
    if mod(size(d,1),2) == 0
        d(1:2:end-1,thisTrial) = d(2:2:end,thisTrial);
    else
        d(1:2:end-1,thisTrial) = d(2:2:end,thisTrial);
    end
end


% figure
% hold on
% plot(d)
% plot(mean(d,2),'blue','linewidth',2)
% plot(max(d,[],2),'red','linewidth',2)
d(isnan(d)) = 0;
Dmax = max(d,[],2);
Dmax = Dmax - mean(Dmax(round(end/3):end));
threshold = std(Dmax(round(end/3):end));
idx = find(Dmax<threshold/3,1,'first');
% plot([idx idx],[min(Dmax) max(Dmax)])

WaveDPbefore = idx*3/5/SSrate;
WaveDPafter = WaveDPbefore;

disp(['Interpolation length = ',num2str(WaveDPbefore*2/30), 'ms'])

Interpolateddata = RawTrace;
ArtifactLocation = zeros(numberofStimPulses,NumTrials);

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

ARdata = Interpolateddata;

% for thisTrial = 1:NumTrials
%     % use salpa algorithm to remove exponential decay.
%     yy = Interpolateddata(:,thisTrial);
%     y = salpa(yy,'tau',round(tau*SampleRate*SSrate));
%     
%     nanidx = isnan(y);
%     y(nanidx) = Interpolateddata(nanidx,thisTrial);
%     ARdata(:,thisTrial) = y;
% end
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
WaveDPbefore = round(WaveDPbefore/SSrate);
WaveDPafter = WaveDPbefore;
ArtifactLocation = round(ArtifactLocation/SSrate);
end