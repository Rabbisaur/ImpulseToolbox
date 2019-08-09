function electrodeTrialData = VnCGetTrialDataAR(instanceinfo, instanceNumber, ElecNumber,TrialAlignST,trialParam,filterParam)

% if ispc
%     slash = '\';
% else
%     slash = '/';
% end

numSamplesBefore = 6;
numSamplesAfter = 6;

thisInstance = instanceNumber;
% build basic trial index
SampleRate = double(instanceinfo(thisInstance).samplerate);
TrialTimeIndex = trialParam.startT * SampleRate : trialParam.endT * SampleRate;

numValidTrials = numel(TrialAlignST);
numDPinaTrial = numel(TrialTimeIndex);

TrialTimeIndex = repmat(TrialTimeIndex,numValidTrials,1);
TrialAlignST = repmat(TrialAlignST,1,numDPinaTrial);
% build the full trial index
% [~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
TrialTimeIndex = TrialTimeIndex+TrialAlignST;

% build trial matrix

% load electrode data back and build aligned trial data matrix
% numDP = instanceinfo(thisInstance).trialInfo.NumDP(trialIdx);
thisElec = find(instanceinfo(thisInstance).ElecOrder == ElecNumber);

% load data
elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{ElecNumber},'r');
tmpdata = fread(elecfp,inf,'int16=>double');
fclose(elecfp);

% remove 50Hz
tmpdata = filtfilt(filterParam.Bbs, filterParam.Abs, tmpdata);

% high pass filtering
tmpdata = filtfilt(filterParam.Bhigh, filterParam.Ahigh, tmpdata);

threshold = 8* median(abs(tmpdata)/0.6745);

idx = tmpdata >= threshold | tmpdata <= -threshold;
idx1 = find(diff(idx) == 1);
idx2 = find(diff(idx) == -1);

idx1 = idx1 - numSamplesBefore;
idx2 = idx2 + numSamplesAfter;

for i = 1:numel(idx1)
    tmpdata(idx1(i):idx2(i)) = 0;
end

tmpdata1 = salpa(tmpdata,'tau',round(0.001*SampleRate));

trialRawData = tmpdata1(TrialTimeIndex);




% output
electrodeTrialData.trialData = trialRawData;
electrodeTrialData.SampleRate = SampleRate;
electrodeTrialData.trialParam = trialParam;
electrodeTrialData.EID = ElecNumber;
electrodeTrialData.ElecIdx = thisElec;
electrodeTrialData.AlignStamps = TrialAlignST;
end