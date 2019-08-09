function instanceinfo = VnCfindTrials(instanceinfo,matdatapath,TrialBit,AlignBit)

if isempty(matdatapath)
else
    load(matdatapath);
end
trialCode = 2^TrialBit;
disp('Searching for trials...')
numInstances = numel(instanceinfo);
for thisInstance = 1:numInstances
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances),'...'])
    instanceinfo(thisInstance).nev.digitimestamps = double(instanceinfo(thisInstance).nev.digitimestamps);
    instanceinfo(thisInstance).nev.digidata = double(instanceinfo(thisInstance).nev.digidata);
    samplerate = double(instanceinfo(thisInstance).samplerate);
    
    numTrials = instanceinfo(thisInstance).numTrials;
%     if numTrials == 1 % conintuous recording, split trials based on digicode
        % find Par.trialB, presumbly a starting signal of a trial
        idx = instanceinfo(thisInstance).nev.digidata == trialCode;
        TrialTimeStamps = double(instanceinfo(thisInstance).nev.digitimestamps(idx));
%         
%         TrialLengthsST = diff(TrialTimeStamps);
%         TrialLengthsT = TrialLengthsST/samplerate;
        
        % search for trials, serial port
        trials = 0;
        finish_flag = false;
        trialNumEncodeST = [];
        trialNumberinSerial = [];
        while ~finish_flag
            encode=double(num2str(trials+1));%serial port encodes. e.g. 0 is encoded as 48, 1 as 49, 10 as [49 48], 12 as [49 50]
            tempInd=strfind(instanceinfo(thisInstance).nev.digidata,encode);
            if isempty(tempInd)
                finish_flag = true;
            else % found one trial
                trials = trials + 1;
                trialNumberinSerial(trials) = trials;
                trialNumEncodeST(trials) = instanceinfo(thisInstance).nev.digitimestamps(tempInd(1));
                if trials == 1
                    idx = instanceinfo(thisInstance).nev.digitimestamps > 0 & instanceinfo(thisInstance).nev.digitimestamps < trialNumEncodeST(trials);
                else
                    % idx = instanceinfo(thisInstance).nev.digitimestamps > TrialTimeStamps(trials) & instanceinfo(thisInstance).nev.digitimestamps < trialNumEncodeTime(trials);
                    idx = instanceinfo(thisInstance).nev.digitimestamps > trialNumEncodeST(trials-1) & instanceinfo(thisInstance).nev.digitimestamps < trialNumEncodeST(trials);
                end
                trialDigiCode{trials} = instanceinfo(thisInstance).nev.digidata(idx);
                trialDigiCodeST{trials} = instanceinfo(thisInstance).nev.digitimestamps(idx);
                
                
            end
        end
        if trials == 0 % if the data is not coded in Xing's way
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
%         if numel(TrialTimeStamps)~=numel(trialNumEncodeST)
%             warning('number of trial begining and end does NOT match')
%         end
        
       
        numNevTrials = numel(trialNumberinSerial);
        if ~isempty(matdatapath)
            if exist('performance','var')
                numMatTrials = numel(performance);
                
                if numNevTrials ~= numMatTrials
                    warning('Number of trials mismatch between Nev and Mat');
                    % truncate the nev trials to match with mat trials
                    numNevTrials = min(numMatTrials,numNevTrials);
                end
                
                fixationKeptTrials = performance ~= 0; % performance == 0 means fixation broken
            end
        end
        
        
        %         microstimtrialflag = allCurrentLevel > 0;
        trialAlign = zeros(numNevTrials,1);
        for thisTrial = 1: numNevTrials
            
            if ~isempty(matdatapath) && exist('fixationKeptTrials','var')
                if fixationKeptTrials(thisTrial) == 1
                    idx = find(trialDigiCode{thisTrial}==2^AlignBit,1,'first');
                    trialAlign(thisTrial) = trialDigiCodeST{thisTrial}(idx);
                else
                    trialAlign(thisTrial) = -1;
                end
                
            else
                idx = find(trialDigiCode{thisTrial}==2^AlignBit,1,'first');
                if isempty(idx)
                    trialAlign(thisTrial) = -1;
                else
                    trialAlign(thisTrial) = trialDigiCodeST{thisTrial}(idx);
                end
            end
        end
        
        if ~isempty(matdatapath)
            if exist('peformance','var') == 1
                validtrialIdx = trialAlign ~=-1 & performance' ~= 0;
            else
                validtrialIdx = trialAlign ~=-1;
            end
        else
            validtrialIdx = trialAlign ~=-1;
        end
        instanceinfo(thisInstance).trialInfo.validtrialIdx = find(validtrialIdx);
        if exist('performance','var')
            instanceinfo(thisInstance).trialInfo.validperformance = performance(validtrialIdx)';
        end
        instanceinfo(thisInstance).trialInfo.validTrialAlignST = trialAlign(validtrialIdx);
        if exist('allCurrentLevel','var')
            instanceinfo(thisInstance).trialInfo.validAllCurrentLevel = allCurrentLevel(validtrialIdx);
        end
                instanceinfo(thisInstance).trialInfo.validbeginST = TrialTimeStamps(validtrialIdx);
        instanceinfo(thisInstance).trialInfo.validtrialNumberinSerial = trialNumberinSerial(validtrialIdx);
        instanceinfo(thisInstance).trialInfo.validEndST = trialNumEncodeST(validtrialIdx);
        instanceinfo(thisInstance).trialInfo.validtrialDigiCode = trialDigiCode(validtrialIdx);
        instanceinfo(thisInstance).trialInfo.validtrialDigiCodeST = trialDigiCodeST(validtrialIdx);
        numValidTrials = sum(validtrialIdx);
        validBeginST = zeros(numValidTrials,1);
        for thisTrial = 1:numValidTrials
            digicode = instanceinfo(thisInstance).trialInfo.validtrialDigiCode{thisTrial};
            digiST = instanceinfo(thisInstance).trialInfo.validtrialDigiCodeST{thisTrial};
            idx = digicode == trialCode;
            validBeginST(thisTrial) = digiST(idx);
        end
        instanceinfo(thisInstance).trialInfo.validBeginST = validBeginST;
        instanceinfo(thisInstance).trialInfo.trialNumberinSerial = trialNumberinSerial;
        instanceinfo(thisInstance).trialInfo.trialAlignST = trialAlign;
        instanceinfo(thisInstance).trialInfo.beginST = TrialTimeStamps';
        instanceinfo(thisInstance).trialInfo.endST = trialNumEncodeST';
        if exist('performance','var')
            instanceinfo(thisInstance).trialInfo.performance = performance';
        end
        if exist('allCurrentLevel','var')
            instanceinfo(thisInstance).trialInfo.allCurrentLevel = allCurrentLevel;
        end
        
        instanceinfo(thisInstance).trialInfo.trialDigiCode = trialDigiCode;
        instanceinfo(thisInstance).trialInfo.trialDigiCodeST = trialDigiCodeST;
        if exist('allFalsePositive','var')
            instanceinfo(thisInstance).trialInfo.AllFalsePositive = AllFalsePositive;
        end
        
%     else
        
%     end
end