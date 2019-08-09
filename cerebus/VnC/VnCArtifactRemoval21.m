function instanceinfo = VnCArtifactRemoval2(instanceinfo,StimulationParam,CurrentLevel,MUAparameters,LFPparameters,trialParam,trialAlignSTfieldname,interpolationLength,DebugFlag)
disp('Removing stimulation artifacts...')

if ispc
    slash = '\';
else
    slash = '/';
end

% parameters
SupersamplingRatio = 16;
numInstances = numel(instanceinfo);

% find current levels
validStimIdx = CurrentLevel > 0;
uniCurrent = unique(CurrentLevel(validStimIdx));
numUniCurrent = numel(uniCurrent);

for thisInstance = 1:numInstances
    fprintf('\n')
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances)])
    numElec = instanceinfo(thisInstance).numElec;
    SampleRate = double(instanceinfo(thisInstance).samplerate);
    
    % detect artifacts
    disp('Calculating average trace')
    for thisElec = 1:numElec
        fprintf('.')
        EID = instanceinfo(thisInstance).ElecOrder(thisElec);
        if EID > 128
            break
        end
        cmd = ['trialAlignST = instanceinfo(thisInstance).trialInfo.',trialAlignSTfieldname,';'];
        eval(cmd);
        electrodeTrialData = VnCGetTrialData2(instanceinfo, thisInstance, EID,trialAlignST,trialParam);
        
        for thisCurrentLevel = 1:numUniCurrent
            idx = CurrentLevel == uniCurrent(thisCurrentLevel);
            if thisElec == 1
                MeanTrace(:,thisCurrentLevel) = squeeze(mean(electrodeTrialData.trialData(idx,:),1))';
            else
                MeanTrace(:,thisCurrentLevel) = MeanTrace(:,thisCurrentLevel) + squeeze(mean(electrodeTrialData.trialData(idx,:),1))';
            end
        end
    end
    % find time using average trace
    MeanTrace = mean(MeanTrace,2)/numElec;
    
    % remove artifact
    fprintf('\n')
    disp('Removing artifacts')
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
        % remove artifact electrode by electrode
        trialData = electrodeTrialData.trialData';
        clearvars electrodeTrialData
        [trialData, artifactTimepoints] = VnCRemoveArtifactsSubroutine(MeanTrace, trialData,StimulationParam,SampleRate,SupersamplingRatio,DebugFlag);
        % calculate MUA, LFP
        [MUAe, LFP] = GetMUAeLFPwithInterpolationSubroutine(trialData,SampleRate,MUAparameters,artifactTimepoints,interpolationLength,LFPparameters);
        electrodeArtifactRemovedNeuralData.MUAe = MUAe;
        electrodeArtifactRemovedNeuralData.LFP = LFP;
        electrodeArtifactRemovedNeuralData.samplerate = SampleRate;
        electrodeArtifactRemovedNeuralData.elecID = EID;
        electrodeArtifactRemovedNeuralData.elecIdx = thisElec;
        idx = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,2,'last');
        idx = idx(1);
        savefilename = 'electrodeArtifactRemovedNeuralData.mat';
        savepath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx),savefilename];
        save(savepath,'electrodeArtifactRemovedNeuralData');
        instanceinfo(thisInstance).electrodeArtifactRemovedNeuralDataPath{thisElec} = savepath;
    end
end
end