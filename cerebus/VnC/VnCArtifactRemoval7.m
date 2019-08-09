function instanceinfo = VnCArtifactRemoval7(instanceinfo,StimulationParam,MUAparameters,trialParam,trialAlignSTfieldname,interpolationLength,DebugFlag)
disp('Removing stimulation artifacts...')

matversion = '-v6';

if ispc
    slash = '\';
else
    slash = '/';
end

% parameters
% SupersamplingRatio = 1;
numInstances = numel(instanceinfo);

% find current levels
% numTrials = numel(CurrentLevel);
% validStimIdx = CurrentLevel > 0;
% uniCurrent = unique(CurrentLevel(validStimIdx));
% numUniCurrent = numel(uniCurrent);


for thisInstance = 1:numInstances
    fprintf('\n')
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances)])
    numElec = instanceinfo(thisInstance).numElec;
    SampleRate = double(instanceinfo(thisInstance).samplerate);

    % remove artifact
    fprintf('\n')
    disp('Removing artifacts')

    StimStartTime = round(-trialParam.startT * SampleRate);
    % remove artifact electrode by electrode
    for thisElec = 1:numElec
        fprintf('.')
        EID = instanceinfo(thisInstance).ElecOrder(thisElec);
        % load in data
        cmd = ['trialAlignST = instanceinfo(thisInstance).trialInfo.',trialAlignSTfieldname,';'];
        eval(cmd);
        electrodeTrialData = VnCGetTrialData2(instanceinfo, thisInstance, EID,trialAlignST,trialParam);
        
        trialData = electrodeTrialData.trialData';
        [trialData, AllArtifactTimepoints] = VnCRemoveArtifactsSubroutine5(trialData,StimulationParam,StimStartTime,SampleRate);
        
        % calculate MUA
        MUAe = GetMUAewithInterpolationSubroutine(trialData,SampleRate,MUAparameters,AllArtifactTimepoints,interpolationLength);
        MUAe.time =  (0:(size(MUAe.data,1)-1))/MUAe.MUAesamplerate + trialParam.startT;
        electrodeArtifactRemovedNeuralData.MUAe = MUAe;
%         electrodeArtifactRemovedNeuralData.LFP = LFP;
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