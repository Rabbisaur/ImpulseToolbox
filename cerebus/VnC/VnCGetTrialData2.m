function electrodeTrialData = VnCGetTrialData2(instanceinfo, instanceNumber, ElecNumber,TrialAlignST,trialParam)

% if ispc
%     slash = '\';
% else
%     slash = '/';
% end

thisInstance = instanceNumber;
% build basic trial index
SampleRate = double(instanceinfo(thisInstance).samplerate);
TrialTimeIndex = trialParam.startT * SampleRate : trialParam.endT * SampleRate;

numValidTrials = numel(TrialAlignST);
numDPinaTrial = numel(TrialTimeIndex);

TrialTimeIndex = repmat(TrialTimeIndex,numValidTrials,1);
TrialAlignSTmat = repmat(TrialAlignST,1,numDPinaTrial);
% build the full trial index
% [~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
TrialTimeIndex = TrialTimeIndex+TrialAlignSTmat;

% build trial matrix

% load electrode data back and build aligned trial data matrix
% numDP = instanceinfo(thisInstance).trialInfo.NumDP(trialIdx);
thisElec = find(instanceinfo(thisInstance).ElecOrder == ElecNumber);
    
    % load data
    elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{ElecNumber},'r');
    tmpdata = fread(elecfp,inf,'int16=>double');
    fclose(elecfp);
    trialRawData = tmpdata(TrialTimeIndex);
    % output
    electrodeTrialData.trialData = trialRawData;
    electrodeTrialData.SampleRate = SampleRate;
    electrodeTrialData.trialParam = trialParam;
    electrodeTrialData.EID = ElecNumber;
    electrodeTrialData.ElecIdx = thisElec;
    electrodeTrialData.AlignStamps = TrialAlignST;
end