function electrodeTrialData = VnCGetTrialDataCell(instanceinfo, instanceNumber, ElecNumber)

% if ispc
%     slash = '\';
% else
%     slash = '/';
% end

thisInstance = instanceNumber;
% build basic trial index
% SampleRate = double(instanceinfo(thisInstance).samplerate);
% TrialTimeIndex = trialParam.startT * SampleRate : trialParam.endT * SampleRate;

% numValidTrials = numel(TrialAlignST);
% numDPinaTrial = numel(TrialTimeIndex);

% TrialTimeIndex = repmat(TrialTimeIndex,numValidTrials,1);
% TrialAlignST = repmat(TrialAlignST,1,numDPinaTrial);
% build the full trial index
% [~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
% TrialTimeIndex = TrialTimeIndex+TrialAlignST;

% build trial matrix

% load electrode data back and build aligned trial data matrix
% numDP = instanceinfo(thisInstance).trialInfo.NumDP(trialIdx);
thisElec = find(instanceinfo(thisInstance).ElecOrder == ElecNumber);
    
    % load data
    if ~exist(instanceinfo(thisInstance).electrodeCachePath{ElecNumber},'file')
        error('Electrode cache file not exist!')
    end
    elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{ElecNumber},'r');
    tmpdata = fread(elecfp,inf,'int16=>double');
    fclose(elecfp);
    
    numTrials = numel(instanceinfo(thisInstance).trialInfo.beginST);
    trialRawData = cell(numTrials,1);
    for thisTrial = 1:numTrials
        if instanceinfo(thisInstance).trialInfo.beginST(thisTrial) == 0
            if thisTrial > 1
                instanceinfo(thisInstance).trialInfo.beginST(thisTrial) = instanceinfo(thisInstance).trialInfo.endST(thisTrial-1);
            elseif thisTrial == 1
                instanceinfo(thisInstance).trialInfo.beginST(thisTrial) = 1;
            end
        end
        if instanceinfo(thisInstance).trialInfo.endST(thisTrial) == 0
            if thisTrial == numel(instanceinfo(thisInstance).trialInfo.beginST)
                instanceinfo(thisInstance).trialInfo.endST(thisTrial) = numel(tmpdata);
            else
                instanceinfo(thisInstance).trialInfo.endST(thisTrial) = instanceinfo(thisInstance).trialInfo.beginST(thisTrial+1);
            end
        end
        TrialTimeIndex = instanceinfo(thisInstance).trialInfo.beginST(thisTrial) : instanceinfo(thisInstance).trialInfo.endST(thisTrial);
        trialRawData{thisTrial} = tmpdata(TrialTimeIndex);
    end
    % output
    electrodeTrialData.trialData = trialRawData;
    electrodeTrialData.trialBeginST = instanceinfo(thisInstance).trialInfo.beginST;
    electrodeTrialData.trialEndST = instanceinfo(thisInstance).trialInfo.endST;
    electrodeTrialData.EID = ElecNumber;
    electrodeTrialData.ElecIdx = thisElec;
%     electrodeTrialData.AlignStamps = TrialAlignST;
end