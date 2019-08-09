function electrodeTrialData = TSLGetTrialData(instanceinfo, ElecID,TrialAlignST,trialParam)

% build basic trial index
SampleRate = double(instanceinfo.samplerate);
TrialTimeIndex = trialParam.startT * SampleRate : trialParam.endT * SampleRate;

numValidTrials = numel(TrialAlignST);
numDPinaTrial = numel(TrialTimeIndex);

TrialTimeIndex = repmat(TrialTimeIndex,numValidTrials,1);
if size(TrialAlignST,2) > size(TrialAlignST,1)
    TrialAlignST = TrialAlignST';
end
electrodeTrialData.AlignStamps = TrialAlignST;
TrialAlignST = repmat(TrialAlignST,1,numDPinaTrial);
% build the full trial index
TrialTimeIndex = TrialTimeIndex+TrialAlignST;

% build trial matrix

% load electrode data back and build aligned trial data matrix
thisElec = find(instanceinfo.ElecOrder == ElecID);

% load data
elecfp = fopen(instanceinfo.electrodeCachePath{ElecID},'r');
tmpdata = fread(elecfp,inf,'int16=>double');
fclose(elecfp);
trialRawData = tmpdata(TrialTimeIndex);
if size(trialRawData,2) > size(trialRawData,1)
    trialRawData = trialRawData';
end
% output
electrodeTrialData.trialData = trialRawData;
electrodeTrialData.SampleRate = SampleRate;
electrodeTrialData.trialParam = trialParam;
electrodeTrialData.EID = ElecID;
electrodeTrialData.ElecIdx = thisElec;
end