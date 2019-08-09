function instanceinfo = VnCfindTrials4(instanceinfo,numTrialinMat,TrialBit,AlignBit)

% define
MaxPossibleNumberTrials = 20000;
xingFlag = 1; % default guess

if numTrialinMat>0
    MaxPossibleNumberTrials = numTrialinMat;
end

trialCode = 2^TrialBit;
disp('Searching for trials...')
numInstances = numel(instanceinfo);
for thisInstance = 1:numInstances
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances),'...'])
    instanceinfo(thisInstance).nev.digitimestamps = double(instanceinfo(thisInstance).nev.digitimestamps);
    instanceinfo(thisInstance).nev.digidata = double(instanceinfo(thisInstance).nev.digidata);
    
    % find Par.trialB, presumbly a starting signal of a trial

    
    % search for trials, serial port
    trials = 0;
    % perform a serial search
    currentPoz = 1;
    trialNumEncodeST = zeros(MaxPossibleNumberTrials,1);
    trialNumberinSerial = zeros(MaxPossibleNumberTrials,1);
    trialDigiCode = cell(MaxPossibleNumberTrials,1);
    trialDigiCodeST = cell(MaxPossibleNumberTrials,1);
    lastTrialNumber = 0;
    for thisTrial = 1:MaxPossibleNumberTrials
        encode=double(num2str(thisTrial));%serial port encodes. e.g. 0 is encoded as 48, 1 as 49, 10 as [49 48], 12 as [49 50]

        % revise
        tmpdata = instanceinfo(thisInstance).nev.digidata(currentPoz : end);
        dataLength = numel(tmpdata);
        encodeLength = numel(encode);
        idx = zeros(dataLength,1);
        for thisByte = 1:encodeLength
            tmp = tmpdata == encode(thisByte);
            tmp = circshift(tmp,-(thisByte-1));
            tmp = tmp';
            idx = idx + tmp;
        end
        
        PozIdx = find(idx == encodeLength,1,'first');
        if ~isempty(PozIdx)
            % found a match
            trials = trials + 1;
            thisPoz = currentPoz + PozIdx - 1;
            trialNumEncodeST(thisTrial) = instanceinfo(thisInstance).nev.digitimestamps(thisPoz+numel(encode)-1);
            if lastTrialNumber == 0
                idx = instanceinfo(thisInstance).nev.digitimestamps > 0 & instanceinfo(thisInstance).nev.digitimestamps <= trialNumEncodeST(thisTrial);
            else
                idx = instanceinfo(thisInstance).nev.digitimestamps > trialNumEncodeST(lastTrialNumber) & instanceinfo(thisInstance).nev.digitimestamps <= trialNumEncodeST(thisTrial);
            end
            trialDigiCode{thisTrial} = instanceinfo(thisInstance).nev.digidata(idx);
            trialDigiCodeST{thisTrial} = instanceinfo(thisInstance).nev.digitimestamps(idx);
            trialNumberinSerial(thisTrial) = thisTrial;
            lastTrialNumber = thisTrial;
            currentPoz = currentPoz + PozIdx + encodeLength - 1; % move to the next position
        end
    end
    
    if trials == 0 % if the data is not coded in Xing's way
        xingFlag = 0;
        idx = instanceinfo(thisInstance).nev.digidata == trialCode;
        TrialTimeStamps = double(instanceinfo(thisInstance).nev.digitimestamps(idx));
        % find trial alignment time
        numTrials = numel(TrialTimeStamps);
        for trials = 1:numTrials
            if trials == numTrials
                idx = instanceinfo(thisInstance).nev.digitimestamps <= instanceinfo(thisInstance).nev.digitimestamps(end) & instanceinfo(thisInstance).nev.digitimestamps > TrialTimeStamps(trials);
            else
                idx = instanceinfo(thisInstance).nev.digitimestamps > TrialTimeStamps(trials) & instanceinfo(thisInstance).nev.digitimestamps < TrialTimeStamps(trials+1);
            end
            trialDigiCode{trials} = instanceinfo(thisInstance).nev.digidata(idx);
            trialDigiCodeST{trials} = instanceinfo(thisInstance).nev.digitimestamps(idx);
        end
    end
    
    if xingFlag
        numTrials = trials;
    else
        
    end
    
    disp(['Got ', num2str(numTrials), ' trials.'])
    
    
    
    
    instanceinfo(thisInstance).trialInfo.trialNumberinSerial = trialNumberinSerial;
    if trials == 0 % if the data is not coded in Xing's way
        instanceinfo(thisInstance).trialInfo.beginST = TrialTimeStamps';
        instanceinfo(thisInstance).trialInfo.endST = [TrialTimeStamps(2:end)';instanceinfo(thisInstance).nev.digitimestamps(end)];
    else
        instanceinfo(thisInstance).trialInfo.beginST = [1;trialNumEncodeST(1:end-1)];
        instanceinfo(thisInstance).trialInfo.endST = trialNumEncodeST;
    end
    
    instanceinfo(thisInstance).trialInfo.trialDigiCode = trialDigiCode;
    instanceinfo(thisInstance).trialInfo.trialDigiCodeST = trialDigiCodeST;
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%% extra steps %%%%%%%%%%%%%%%%%%%%%

% trim all the trial marks to the max number of trials found
if xingFlag
    lastTrialIdx = zeros(numInstances,1);
    for thisInstance = 1:numInstances
        lastTrialIdx(thisInstance) = find(instanceinfo(thisInstance).trialInfo.trialNumberinSerial>0, 1, 'last');
    end
    MaxTrials= max(lastTrialIdx);
    for thisInstance = 1:numInstances
        instanceinfo(thisInstance).trialInfo.trialNumberinSerial = instanceinfo(thisInstance).trialInfo.trialNumberinSerial(1:MaxTrials);
        instanceinfo(thisInstance).trialInfo.beginST = instanceinfo(thisInstance).trialInfo.beginST(1:MaxTrials);
        instanceinfo(thisInstance).trialInfo.endST = instanceinfo(thisInstance).trialInfo.endST(1:MaxTrials);
        instanceinfo(thisInstance).trialInfo.trialDigiCode = instanceinfo(thisInstance).trialInfo.trialDigiCode(1:MaxTrials);
        instanceinfo(thisInstance).trialInfo.trialDigiCodeST = instanceinfo(thisInstance).trialInfo.trialDigiCodeST(1:MaxTrials);
    end
    numTrials = MaxTrials;
end

% Get align point for each trial
for thisInstance = 1:numInstances
    trialAlign = zeros(numTrials,1);
    for thisTrial = 1: numTrials
        idx = find(instanceinfo(thisInstance).trialInfo.trialDigiCode{thisTrial}==2^AlignBit);
        if isempty(idx)
            trialAlign(thisTrial) = nan;
        else
            if numel(idx) > 1
                warning('There are multiple align point in a trial.')
                idx2 = find(instanceinfo(thisInstance).trialInfo.trialDigiCode{thisTrial}==trialCode);
                if numel(idx2) > 1
                    warning('There are multiple trial code with in a trial, trial seeking was wrong!')
                    idx2 = idx2(end);
                end
                if ~isempty(idx2)
                    idx = idx(idx >= idx2);
                end
                if isempty(idx)
                    trialAlign(thisTrial) = nan;
                else
                    idx = idx(1);
                    trialAlign(thisTrial) = instanceinfo(thisInstance).trialInfo.trialDigiCodeST{thisTrial}(idx);
                end
            else
                trialAlign(thisTrial) = instanceinfo(thisInstance).trialInfo.trialDigiCodeST{thisTrial}(idx);
            end
        end
    end
    
    idx = isnan(trialAlign);
    if sum(idx) > 0
        tmp = trialAlign(~idx);
        if isempty(tmp)
            warning(['No digital event available in instance ', num2str(thisInstance),'!'])
        else
            trialAlign(idx) = tmp(1);
        end
    end
    instanceinfo(thisInstance).trialInfo.trialAlignST = trialAlign;
end

% if it is coded in Xing's way and there is a difference in number of
% trials, try guess the trial information from good instances but lable it.
if xingFlag
    badInstanceIdx = zeros(numInstances);
    for thisInstance = 1:numInstances
        idx = instanceinfo(thisInstance).trialInfo.trialNumberinSerial == 0;
        badInstanceIdx(thisInstance) = sum(idx) > 0;
    end
    goodInstanceIdx = badInstanceIdx==0;
    badInstanceIdx = find(badInstanceIdx);
    NumBadinstance = numel(badInstanceIdx);
    if  NumBadinstance > 0
        disp('We found difference in number of trials, trying to repair it.')
        goodInstanceIdx = find(goodInstanceIdx,1,'first');
        disp(['Instance ', num2str(goodInstanceIdx),' is used as good reference.'])
        for thisBadInstance = 1:NumBadinstance
            InstanceIdx = badInstanceIdx(thisBadInstance);
            disp(['Repairing instance ',num2str(InstanceIdx)]);
            % find an anchor point
            idx1 = instanceinfo.instanceID == InstanceIdx;
            idx2 = instanceinfo.instanceID == goodInstanceIdx;
            trialNumberDiff = instanceinfo(idx1).trialInfo.trialNumberinSerial - instanceinfo(idx2).trialInfo.trialNumberinSerial;
            AnchorTrial = find(trialNumberDiff == 0,1,'last');
            disp(['Trial #', num2str(AnchorTrial), ' is used as an anchor for calcuating the time difference between instances'])
            STdiff = instanceinfo(idx1).trialInfo.endST(AnchorTrial) - instanceinfo(goodInstanceIdx).trialInfo.endST(AnchorTrial);
            disp(['Time difference between the two instances is ', num2str(STdiff),' samples.'])
            idx = instanceinfo(idx1).trialInfo.trialNumberinSerial == 0;
            NumTrialsBefore = sum(idx == 0);
            NumTrialsRepaired = sum(idx);
            
            % repair
            instanceinfo(idx1).trialInfo.repairedTrialIdx = idx;
            instanceinfo(idx1).trialInfo.beginST(idx) = instanceinfo(goodInstanceIdx).trialInfo.beginST(idx) + STdiff;
            instanceinfo(idx1).trialInfo.endST(idx) = instanceinfo(goodInstanceIdx).trialInfo.endST(idx) + STdiff;
            instanceinfo(idx1).trialInfo.trialAlignST(idx) = instanceinfo(goodInstanceIdx).trialInfo.trialAlignST(idx) + STdiff;
            idx = find(idx);
            for thisTrial = 1:numel(idx)
                trialID = idx(idx1);
                instanceinfo(idx1).trialInfo.trialDigiCode{trialID} = instanceinfo(goodInstanceIdx).trialInfo.trialDigiCode{trialID};
                instanceinfo(idx1).trialInfo.trialDigiCodeST{trialID} = instanceinfo(goodInstanceIdx).trialInfo.trialDigiCodeST{trialID} + STdiff;
            end
            disp([num2str(NumTrialsRepaired), ' trial(s) were repaired'])
            disp(['Number of trials changed from ', num2str(NumTrialsBefore), ' -> ', num2str(NumTrialsBefore + NumTrialsRepaired)])
        end
    end
end