function electrodeTrialData = VnCGetTrialData3(instanceinfo, ElecNumber,TrialAlignST,trialParam,basepath)

% if ispc
%     slash = '\';
% else
%     slash = '/';
% end

% thisInstance = instanceNumber;
% build basic trial index
SampleRate = double(instanceinfo.samplerate);
TrialTimeIndex = trialParam.startT * SampleRate : trialParam.endT * SampleRate;

TrialTimeIndex = bsxfun(@plus,TrialTimeIndex,TrialAlignST');    

% build trial matrix

% load electrode data back and build aligned trial data matrix
% numDP = instanceinfo(thisInstance).trialInfo.NumDP(trialIdx);
thisElec = find(instanceinfo.ElecOrder == ElecNumber);
    
    % load data
    datapath = instanceinfo.electrodeCachePath{ElecNumber};
    idx = strfind(datapath,'instance');
    datapath = [basepath,'/',datapath(idx:end)];
    elecfp = fopen(datapath,'r');
    tmpdata = fread(elecfp,inf,'int16=>double');
    fclose(elecfp);
    tmp = TrialTimeIndex < 1;
    TrialTimeIndex(tmp) = 1;
    tmp = TrialTimeIndex > numel(tmpdata);
    TrialTimeIndex(tmp) = numel(tmpdata);
    trialRawData = tmpdata(TrialTimeIndex);
    % output
    electrodeTrialData.trialData = trialRawData;
    electrodeTrialData.SampleRate = SampleRate;
    electrodeTrialData.trialParam = trialParam;
    electrodeTrialData.EID = ElecNumber;
    electrodeTrialData.ElecIdx = thisElec;
    electrodeTrialData.AlignStamps = TrialAlignST;
end