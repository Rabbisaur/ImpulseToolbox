function instanceinfo = VnCfindTrials2(instanceinfo,matdatapath,TrialBit,AlignBit)

% define
MaxPossibleNumberTrials = 20000;
xingFlag = 1; % default guess

if isempty(matdatapath)
else
    load(matdatapath);
    if exist(allFixT)
        numTrialinMat = numel(allFixT);
        MaxPossibleNumberTrials = numTrialinMat;
    end
end
trialCode = 2^TrialBit;
disp('Searching for trials...')
numInstances = numel(instanceinfo);
for thisInstance = 1:numInstances
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances),'...'])
    instanceinfo(thisInstance).nev.digitimestamps = double(instanceinfo(thisInstance).nev.digitimestamps);
    instanceinfo(thisInstance).nev.digidata = double(instanceinfo(thisInstance).nev.digidata);
%     samplerate = double(instanceinfo(thisInstance).samplerate);
    
    % find Par.trialB, presumbly a starting signal of a trial
    
    %
    %         TrialLengthsST = diff(TrialTimeStamps);
    %         TrialLengthsT = TrialLengthsST/samplerate;
    
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
        for thisPoz = currentPoz : (numel(instanceinfo(thisInstance).nev.digidata)-numel(encode)+1)
            currentDigi = instanceinfo(thisInstance).nev.digidata(thisPoz:(thisPoz+numel(encode)-1));
            if issame(currentDigi,encode)
                % found a match
                trials = trials + 1;
                trialNumEncodeST(thisTrial) = instanceinfo(thisInstance).nev.digitimestamps(thisPoz+numel(encode)-1);
                if lastTrialNumber == 0
                    idx = instanceinfo(thisInstance).nev.digitimestamps > 0 & instanceinfo(thisInstance).nev.digitimestamps <= trialNumEncodeST(thisTrial);
                else
                    idx = instanceinfo(thisInstance).nev.digitimestamps > trialNumEncodeST(lastTrialNumber) & instanceinfo(thisInstance).nev.digitimestamps <= trialNumEncodeST(thisTrial);
                end
                trialDigiCode{thisTrial} = instanceinfo(thisInstance).nev.digidata(idx);
                trialDigiCodeST{thisTrial} = instanceinfo(thisInstance).nev.digitimestamps(idx);
                currentPoz = thisPoz + numel(encode)-1;
                trialNumberinSerial(thisTrial) = thisTrial;
                lastTrialNumber = thisTrial;
                break;
            else
                % continue to the end
            end
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
        instanceinfo(thisInstance).trialInfo.beginST = [0;trialNumEncodeST(1:end-1)];
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
        idx = find(instanceinfo(thisInstance).trialInfo.trialDigiCode{thisTrial}==2^AlignBit,1,'first');
        if isempty(idx)
            trialAlign(thisTrial) = nan;
        else
            trialAlign(thisTrial) = instanceinfo(thisInstance).trialInfo.trialDigiCodeST{thisTrial}(idx);
        end
    end
    
    idx = isnan(trialAlign);
    if sum(idx) > 0
        tmp = trialAlign(~idx);
        trialAlign(idx) = tmp(1);
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
            trialNumberDiff = instanceinfo(InstanceIdx).trialInfo.trialNumberinSerial - instanceinfo(goodInstanceIdx).trialInfo.trialNumberinSerial;
            AnchorTrial = find(trialNumberDiff == 0,1,'first');
            disp(['Trial #', num2str(AnchorTrial), ' is used as an anchor for calcuating the time difference between instances'])
            STdiff = instanceinfo(InstanceIdx).trialinfo.endST(AnchorTrial) - instanceinfo(goodInstanceIdx).trialinfo.endST(AnchorTrial);
            disp(['Time difference between the two instances is ', num2str(STdiff),' samples.'])
            idx = instanceinfo(InstanceIdx).trialInfo.trialNumberinSerial == 0;
            NumTrialsBefore = sum(idx == 0);
            NumTrialsRepaired = sum(idx);
            
            % repair
            instanceinfo(InstanceIdx).trialInfo.repairedTrialIdx = idx;
            instanceinfo(InstanceIdx).trialInfo.beginST(idx) = instanceinfo(goodInstanceIdx).trialinfo.beginST(idx) + STdiff;
            instanceinfo(InstanceIdx).trialInfo.endST(idx) = instanceinfo(goodInstanceIdx).trialinfo.endST(idx) + STdiff;
            instanceinfo(InstanceIdx).trialInfo.AlignST(idx) = instanceinfo(goodInstanceIdx).trialInfo.AlignST(idx) + STdiff;
            idx = find(idx);
            for thisTrial = 1:numel(idx)
                trialID = idx(thisTrial);
                instanceinfo(InstanceIdx).trialInfo.trialDigiCode{trialID} = instanceinfo(goodInstanceIdx).trialInfo.trialDigiCode{trialID};
                instanceinfo(InstanceIdx).trialInfo.trialDigiCodeST{trialID} = instanceinfo(goodInstanceIdx).trialInfo.trialDigiCodeST{trialID} + STdiff;
            end
            disp([num2str(NumTrialsRepaired), 'trial(s) were repaired'])
            disp(['Number of trials changed from ', num2str(NumTrialsBefore), ' -> ', num2str(NumTrialsBefore + NumTrialsRepaired)])
        end
    end
end