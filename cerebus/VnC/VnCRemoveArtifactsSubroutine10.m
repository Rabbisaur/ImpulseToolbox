function [ARdata,ArtifactLocation,Interpolateddata,WaveDPbefore,WaveDPafter] = VnCRemoveArtifactsSubroutine10(RawTrace,StimulationParam,StimStartTime,SampleRate,EID,thisInstance,instanceinfo,GoodMStrialIdx,FineAlignFlag)

if ispc
    slash = '\';
else
    slash = '/';
end

NumTrials = size(RawTrace,2);
NumDP = size(RawTrace,1);
tau = 0.0008; % sec
SSrate = 1; % super sampling rate

% FineAlignFlag = 1;

if FineAlignFlag
%     shiftvalue0 = instanceinfo(thisInstance).trialInfo.fineAlign.shiftvalue0;
%     SSrate = instanceinfo(thisInstance).trialInfo.fineAlign.SupersamplingRatio;
%     GoodMStrialIdx = instanceinfo(thisInstance).trialInfo.fineAlign.GoodMStrialIdx;
    SSrate = 8;
end

% inputs
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
% parameters
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
numberofStimPulses = StimulationParam.numberofStimPulses;

Fl = 100;
Fn = SampleRate/2;
N = 4;
[BWave, AWave] = butter(N,Fl/Fn,'high'); % Waveform hi pass
HipassTrace = filtfilt(BWave, AWave,RawTrace);

% super sampling
numSamplePoints = size(HipassTrace,1);
x = 1:numSamplePoints;
xx = 1:1/SSrate:numSamplePoints;
HipassTrace = spline(x,HipassTrace',xx)';
RawTrace = spline(x,RawTrace',xx)';

% padding
NumDPpadding = 100*SSrate;
RawTrace = [zeros(NumDPpadding,NumTrials);RawTrace;zeros(NumDPpadding,NumTrials)];
HipassTrace = [zeros(NumDPpadding,NumTrials);HipassTrace;zeros(NumDPpadding,NumTrials)];


StimulationPulseInterval = 1/StimulationPulsesFrequency;

StimStartTime = StimStartTime + NumDPpadding/SSrate;
predictedLoc = (0:(numberofStimPulses-1)) * 1/StimulationPulsesFrequency * SampleRate;
predictedLoc = predictedLoc + StimStartTime;
predictedLoc = predictedLoc + StimulationParam.stimulationwaveformwidth*0.95 * SampleRate;

loc0 = round(predictedLoc * SSrate);


% Find stimulation exact locations, get waveform


StimWidth = 1/StimulationPulsesFrequency/2;
idx = (0:round(StimWidth*SampleRate*SSrate-1))';
WaveLength = numel(idx);
idx = idx - floor(max(idx)/2);
idx1 = repmat(idx,1,numberofStimPulses);
idx1 = idx1 + repmat(loc0,numel(idx),1);
d = zeros(WaveLength,NumTrials);

% fine alignment
NumTrials2 = sum(GoodMStrialIdx);
HipassTrace2 = HipassTrace(:,GoodMStrialIdx);
RawTrace2 = RawTrace(:,GoodMStrialIdx);
meanWaveform = zeros(WaveLength,NumTrials2);
h = figure('position',[100 100 900 700],'visible','off');
for thisTrial = 1:NumTrials2
    Waveform = HipassTrace2(idx1,thisTrial);
    Waveform = reshape(Waveform,numel(idx),numberofStimPulses);
    meanWaveform(:,thisTrial) = mean(Waveform,2);
end
for thisTrial = 1:NumTrials2
    meanWaveform(:,thisTrial) = meanWaveform(:,thisTrial) - mean([meanWaveform(1,thisTrial),meanWaveform(end,thisTrial)]);
end

subplot(2,1,1)
hold on
ydata = meanWaveform;
ydata = downsample(ydata,SSrate);
plot(ydata)
xlabel('Sample points')
ylabel('Voltage')
[~,maxIdx] = max(meanWaveform,[],1);
[~,minIdx] = min(meanWaveform,[],1);
zerocrossingpoint = minIdx;
for thisTrial = 1:NumTrials2
    [~,zerocrossingpoint(thisTrial)] = min((meanWaveform(min([maxIdx(thisTrial),minIdx(thisTrial)]):max([maxIdx(thisTrial),minIdx(thisTrial)]),thisTrial)-0).^2);
    zerocrossingpoint(thisTrial) = zerocrossingpoint(thisTrial) + min([maxIdx(thisTrial),minIdx(thisTrial)])-1;
    plot(zerocrossingpoint(thisTrial)/SSrate,meanWaveform(zerocrossingpoint(thisTrial),thisTrial),'rx')
end

middleIdx = median(zerocrossingpoint);
plot([middleIdx middleIdx]/SSrate,[min(min(meanWaveform)),max(max(meanWaveform))],'black','linewidth',2)
shiftvalue0 = zerocrossingpoint - middleIdx;
meanWaveformShifted = meanWaveform;
for thisTrial = 1:NumTrials2
    meanWaveformShifted(:,thisTrial) = circshift(meanWaveform(:,thisTrial),-round(shiftvalue0(thisTrial)));
end

subplot(2,1,2)
hold on
xlabel('Sample points')
ylabel('Voltage')
ydata = meanWaveformShifted;
ydata = downsample(ydata,SSrate);
plot(ydata)
[~,maxIdx] = max(meanWaveformShifted,[],1);
[~,minIdx] = min(meanWaveformShifted,[],1);
for thisTrial = 1:NumTrials2
    [~,zerocrossingpoint(thisTrial)] = min((meanWaveformShifted(min([maxIdx(thisTrial),minIdx(thisTrial)]):max([maxIdx(thisTrial),minIdx(thisTrial)]),thisTrial)-0).^2);
    zerocrossingpoint(thisTrial) = zerocrossingpoint(thisTrial) + min([maxIdx(thisTrial),minIdx(thisTrial)])-1;
    plot(zerocrossingpoint(thisTrial)/SSrate,meanWaveformShifted(zerocrossingpoint(thisTrial),thisTrial),'go')
end
middleIdx2 = median(zerocrossingpoint);
plot([middleIdx2 middleIdx2]/SSrate,[min(min(meanWaveform)),max(max(meanWaveform))],'black','linewidth',2)

for thisTrial = 1:NumTrials2
    HipassTrace2(:,thisTrial) = circshift(HipassTrace2(:,thisTrial),-round(shiftvalue0(thisTrial)));
    RawTrace2(:,thisTrial) = circshift(RawTrace2(:,thisTrial),-round(shiftvalue0(thisTrial)));
end


HipassTrace(:,GoodMStrialIdx)=HipassTrace2;
RawTrace(:,GoodMStrialIdx)=RawTrace2;
idx2 = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,3,'last');
idx2 = idx2(1);
dirpath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx2),'debugimg'];
if exist(dirpath,'dir')==7
    
else
    mkdir(dirpath)
end
savepath = [dirpath,'/',num2str(EID),'FineAlignment.png'];
saveas(h,savepath)
% savepath = [dirpath,'/',num2str(EID),'FineAlignment.fig'];
% saveas(h,savepath)
close(h)

% perform measurment calculate interpolation length
h = figure('position',[100 100 900 700],'visible','off');

for thisTrial = 1:NumTrials2
    Waveform = HipassTrace2(idx1,thisTrial);
    Waveform = reshape(Waveform,numel(idx),numberofStimPulses);

    subplot(2,2,1)
    hold on
    xlabel('Sample points')
    ylabel('Voltage')
    ydata = mean(Waveform,2);
    ydata = downsample(ydata,SSrate);
    plot(ydata)
    meanWaveform = mean(Waveform,2);

    for InterpWidth = 1:2:WaveLength-2
        interpIdx = floor(WaveLength / 2 + (0:(InterpWidth-1)) - (InterpWidth/2)+1);
        x = 1:WaveLength;
        x(interpIdx) = [];
        tmp = meanWaveform;
        tmp(interpIdx) = [];
        xx = 1:WaveLength;
        tmp = interp1(x,tmp,xx);
        d(InterpWidth,thisTrial) = sum((tmp' - mean([meanWaveform(1),meanWaveform(end)])).^2);
    end
    if mod(size(d,1),2) == 0
        d(2:2:end,thisTrial) = d(1:2:end,thisTrial);
    else
        d(2:2:end,thisTrial) = d(1:2:end-1,thisTrial);
    end
    
end

d(isnan(d)) = 0;
Dmax = max(d,[],2);
subplot(2,2,2)
xlabel('Sample points')
ylabel('SS error')
hold on
ydata = d;
ydata = downsample(ydata,SSrate);
plot(ydata)
ydata = Dmax;
ydata = downsample(ydata,SSrate);
plot(ydata,'red','linewidth',2)
Dmax = smooth(Dmax,5*SSrate);
Dmax = Dmax(1:end-3);
% K = (Dmax-Dmax(end))./(numel(Dmax):-1:1)';
x = round(numel(Dmax)*2/3):numel(Dmax);
g = fit(x',Dmax(x),'poly1');
base = g(1:numel(Dmax));
idx2 = find((Dmax-base)>0);
if isempty(idx2)
    idx2 = numel(Dmax);
else
    tmp = find(diff(idx2)~=1,1,'first');
    if isempty(tmp)
        idx2 = idx2(end);
    else
        idx2 = idx2(tmp)+1;
    end
end
subplot(2,2,3)
xlabel('Sample points')
ylabel('Max SS error')
hold on
plot(ydata)
ydata = base;
ydata = downsample(ydata,SSrate);
plot(ydata)
plot([idx2 idx2]/SSrate,[min(Dmax) max(Dmax)])

WaveDPbefore = (idx2/2+1)/SSrate;
WaveDPafter = WaveDPbefore;

disp(['AR interpolation length = ',num2str(WaveDPbefore*2/SampleRate*1000), 'ms'])

NumPulsesBefore = floor(StimStartTime/SampleRate / StimulationPulseInterval);
NumPulsesAfter = floor((NumDP - StimStartTime)/SampleRate / StimulationPulseInterval);
FirstPulseDP = round(StimStartTime - NumPulsesBefore * StimulationPulseInterval * SampleRate);

numberofStimPulses = NumPulsesBefore + NumPulsesAfter + 1;
StimStartTime = FirstPulseDP;



StimStartTime = StimStartTime + NumDPpadding/SSrate;
predictedLoc = (0:(numberofStimPulses-1)) * 1/StimulationPulsesFrequency * SampleRate;
predictedLoc = predictedLoc + StimStartTime;
predictedLoc = predictedLoc + StimulationParam.stimulationwaveformwidth*0.95 * SampleRate;

loc0 = round(predictedLoc * SSrate);


Interpolateddata = RawTrace;
ArtifactLocation = zeros(numberofStimPulses,NumTrials);

idx = round(-WaveDPbefore* SSrate) : round(WaveDPafter* SSrate);

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

subplot(2,2,4)
xlabel('Sample points')
ylabel('Voltage')
hold on
for thisTrial = 1:NumTrials
    ArtifactLocation(:,thisTrial) = loc0;
    % linear interpolate
    y = RawTrace(:,thisTrial);
    x = 1:numel(y);
    xx = x;
    x(idx) = [];
    y(idx) = [];
    waveform = Interpolateddata(idx,thisTrial);
    waveform = reshape(waveform',[],numberofStimPulses);
    maxvalue = max(waveform,[],1);
    minvalue = min(waveform,[],1);
    [~,idx2]=max(maxvalue - minvalue);
    ydata = waveform(:,idx2)-mean([waveform(1,idx2),waveform(end,idx2)]);
    ydata = downsample(ydata,SSrate);
    plot(ydata,'black');
    Interpolateddata(:,thisTrial) = interp1(x,y,xx,'linear');
    waveform = Interpolateddata(idx,thisTrial);
    waveform = reshape(waveform',[],numberofStimPulses);
    ydata = waveform(:,idx2)-mean([waveform(1,idx2),waveform(end,idx2)]);
    ydata = downsample(ydata,SSrate);
    plot(ydata,'red');
end

% save figure
idx2 = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,3,'last');
idx2 = idx2(1);
dirpath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx2),'debugimg'];
if exist(dirpath,'dir')==7
else
    mkdir(dirpath)
end
savepath = [dirpath,'/',num2str(EID),'AR.png'];
saveas(h,savepath)
% savepath = [dirpath,'/',num2str(EID),'AR.fig'];
% saveas(h,savepath)
close(h)

ARdata = Interpolateddata;

for thisTrial = 1:NumTrials
    % use salpa algorithm to remove exponential decay.
    yy = Interpolateddata(:,thisTrial);
    y = salpa(yy,'tau',round(tau*SampleRate*SSrate));
    
    nanidx = isnan(y);
    y(nanidx) = Interpolateddata(nanidx,thisTrial);
    ARdata(:,thisTrial) = y;
end
if true
    [~,idx3] = max(max(HipassTrace,[],1)-min(HipassTrace,[],1));
    h=figure('position',[100 100 900 700],'visible','off');
    hold on
    xlabel('Samples')
    ylabel('Voltage')
    ydata = RawTrace(:,idx3);
    ydata = downsample(ydata,SSrate);
    plot(ydata)
    ydata = Interpolateddata(:,idx3);
    ydata = downsample(ydata,SSrate);
    plot(ydata)
    ydata = ARdata(:,idx3);
    ydata = downsample(ydata,SSrate);
    plot(ydata)
    
    plot(loc0/SSrate, max(mean(ARdata,2))*ones(size(loc0))*10,'rx')
    plot(idx/SSrate,80*ones(size(idx))*10,'r.')
    % save figure
    idx2 = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,3,'last');
    idx2 = idx2(1);
    dirpath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx2),'debugimg'];
    if exist(dirpath,'dir')==7
    else
        mkdir(dirpath)
    end
    savepath = [dirpath,'/',num2str(EID),'biggestTrialCompare.png'];
    saveas(h,savepath)
    savepath = [dirpath,'/',num2str(EID),'biggestTrialCompare.fig'];
    saveas(h,savepath)
    close(h)
end

% remove padding
if NumDPpadding > 0
    ARdata = ARdata((NumDPpadding+1):(end-NumDPpadding),:);
    Interpolateddata = Interpolateddata((NumDPpadding+1):(end-NumDPpadding),:);
    ArtifactLocation = ArtifactLocation - NumDPpadding;
end

% downsample
ARdata = downsample(ARdata,SSrate);
Interpolateddata = downsample(Interpolateddata,SSrate);
ArtifactLocation = round(ArtifactLocation/SSrate);
WaveDPbefore = round(WaveDPbefore);
WaveDPafter = WaveDPbefore;

end