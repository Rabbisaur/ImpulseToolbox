function instanceinfo = VnCGetTrialData(instanceinfo,trialParam,reloadflag)

if ispc
    slash = '\';
else
    slash = '/';
end

NumInstances = numel(instanceinfo);

for thisInstance = 1:NumInstances
    disp(['Get trial data working on instance ', num2str(thisInstance), '/', num2str(NumInstances)])
    % build basic trial index
    SampleRate = double(instanceinfo(thisInstance).samplerate);
    TrialTimeIndex = trialParam.startT * SampleRate : trialParam.endT * SampleRate;
    if isfield(instanceinfo(thisInstance).trialInfo,'validCorrectedTrialAlignST')
        validTrialAlignST = instanceinfo(thisInstance).trialInfo.validCorrectedTrialAlignST;
    else
        validTrialAlignST = instanceinfo(thisInstance).trialInfo.validTrialAlignST;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Debugging!
%     validTrialAlignST = instanceinfo(thisInstance).trialInfo.validTrialAlignST;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if isfield(instanceinfo(thisInstance).trialInfo,'validAllCurrentLevel')
%         validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
%         visualidx = validAllCurrentLevel == 0;
%         validTrialAlignST(visualidx) = validTrialAlignST(visualidx) - trialParam.FixT * SampleRate;
%     else
%         % do nothing
%     end
    
    
    numValidTrials = numel(validTrialAlignST);
    numDPinaTrial =numel(TrialTimeIndex);
    
    TrialTimeIndex = repmat(TrialTimeIndex,numValidTrials,1);
    TrialAlignST = repmat(validTrialAlignST,1,numDPinaTrial);
    % build the full trial index
    [~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
    TrialTimeIndex = TrialTimeIndex+TrialAlignST - instanceinfo(thisInstance).trialInfo.Timestamp(trialIdx);
    
    % build trial matrix
    numElec = instanceinfo(thisInstance).numElec;
    
    % load electrode data back and build aligned trial data matrix
    numDP = instanceinfo(thisInstance).trialInfo.NumDP(trialIdx);
    for thisElec = 1:numElec
        idx = strfind(instanceinfo(thisInstance).electrodeCachePath{thisElec},[slash,'tmp',slash]);
        electrodeTrialDataSavePath = [instanceinfo(thisInstance).electrodeCachePath{thisElec}(1:idx),'electrodeTrialDataMatrix.mat'];
        if exist(electrodeTrialDataSavePath,'file') && reloadflag == 0
            continue
        end
        EID = instanceinfo(thisInstance).ElecOrder(thisElec);
        % load data
        elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{thisElec},'r');
        tmpdata = fread(elecfp,numDP,'int16=>double');
        fclose(elecfp);
        trialRawData = tmpdata(TrialTimeIndex);
        % save the raw data to disk
        electrodeTrialData.trialData = trialRawData;
        electrodeTrialData.SampleRate = SampleRate;
        electrodeTrialData.trialParam = trialParam;
        electrodeTrialData.trialInfo = instanceinfo(thisInstance).trialInfo;
        electrodeTrialData.EID = EID;
        electrodeTrialData.ElecIdx = thisElec;
        electrodeTrialData.AlignStamps = validTrialAlignST;
        save(electrodeTrialDataSavePath,'electrodeTrialData')
        instanceinfo(thisInstance).electrodeTrialDataPath{thisElec} = electrodeTrialDataSavePath;
    end
end