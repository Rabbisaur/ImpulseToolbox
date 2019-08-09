function instanceinfo = VnCArtifactRemoval5(instanceinfo,StimulationParam,CurrentLevel,MUAparameters,LFPparameters,trialParam,trialAlignSTfieldname,interpolationLength,DebugFlag)
disp('Removing stimulation artifacts...')

matversion = '-v6';

if ispc
    slash = '\';
else
    slash = '/';
end

% parameters
SupersamplingRatio = 1;
numInstances = numel(instanceinfo);

% find current levels
numTrials = numel(CurrentLevel);
validStimIdx = CurrentLevel > 0;
uniCurrent = unique(CurrentLevel(validStimIdx));
numUniCurrent = numel(uniCurrent);

for thisInstance = 1:numInstances
    fprintf('\n')
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances)])
    numElec = instanceinfo(thisInstance).numElec;
    SampleRate = double(instanceinfo(thisInstance).samplerate);
    
    % remove artifact
    fprintf('\n')
    disp('Removing artifacts')
    StimStartEndTime = instanceinfo(thisInstance).trialInfo.TrialMicroStimStartEndPoint;
    StimStartEndTime = StimStartEndTime - trialParam.startT * instanceinfo(thisInstance).samplerate;
    MStrialIdx = CurrentLevel > 0;
    MSCurrentLevel = CurrentLevel(MStrialIdx);
    % remove artifact electrode by electrode
    for thisElec = 1:numElec
        fprintf('.')
        EID = instanceinfo(thisInstance).ElecOrder(thisElec);
        if EID > 128
            break
        end
        % load in data
        cmd = ['trialAlignST = instanceinfo(thisInstance).trialInfo.',trialAlignSTfieldname,';'];
        eval(cmd);
        electrodeTrialData = VnCGetTrialData2(instanceinfo, thisInstance, EID,trialAlignST,trialParam);
        
        trialData = electrodeTrialData.trialData';
        clearvars electrodeTrialData
        AllArtifactTimepoints = zeros(StimulationParam.numberofStimPulses,numTrials);
        % perform artifact removal for each current level
        for thisCurrentLevel = 1:numUniCurrent
            MSidx = MSCurrentLevel == uniCurrent(thisCurrentLevel);
            CondStimStartEndTime = StimStartEndTime(:,MSidx);
            idx = CurrentLevel == uniCurrent(thisCurrentLevel);
            CondTrialData = trialData(:,idx);
%             DebugFlag = 1;
            [CondTrialData, artifactTimepoints] = VnCRemoveArtifactsSubroutine4(CondTrialData,StimulationParam,CondStimStartEndTime,SampleRate,SupersamplingRatio,DebugFlag);
            AllArtifactTimepoints(:,idx) = artifactTimepoints;
            trialData(:,idx) = CondTrialData;
        end
        % calculate MUA, LFP
        [MUAe, LFP] = GetMUAeLFPwithInterpolationSubroutine3(trialData,SampleRate,MUAparameters,AllArtifactTimepoints,interpolationLength,LFPparameters);
        electrodeArtifactRemovedNeuralData.MUAe = MUAe;
        electrodeArtifactRemovedNeuralData.LFP = LFP;
        electrodeArtifactRemovedNeuralData.samplerate = SampleRate;
        electrodeArtifactRemovedNeuralData.elecID = EID;
        electrodeArtifactRemovedNeuralData.elecIdx = thisElec;
        idx = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,2,'last');
        idx = idx(1);
        savefilename = 'electrodeArtifactRemovedNeuralData.mat';
        savepath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx),savefilename];
        save(savepath,'electrodeArtifactRemovedNeuralData',matversion);
        instanceinfo(thisInstance).electrodeArtifactRemovedNeuralDataPath{thisElec} = savepath;
    end
end
end