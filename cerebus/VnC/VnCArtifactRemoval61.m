function instanceinfo = VnCArtifactRemoval6(instanceinfo,StimulationParam,CurrentLevel,MUAparameters,LFPparameters,trialParam,trialAlignSTfieldname,interpolationLength,DebugFlag)
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

% LFP parameters
LFPsamplerate = LFPparameters.LFPsamplingrate;
LFPlowpassFreq = LFPparameters.LFPlowpassFreq;


for thisInstance = 1:numInstances
    fprintf('\n')
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances)])
    numElec = instanceinfo(thisInstance).numElec;
    SampleRate = double(instanceinfo(thisInstance).samplerate);
    % build filters
    Fn = SampleRate/2;
    Fbp=[49,51];
    N  = 2;    % filter order
    [Bbs, Abs] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn],'stop'); % Bandstop
    
    % low pass 150Hz
    Fl = LFPlowpassFreq;
    Fn = SampleRate/2;
    N = 2;
    [BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass
    
    Fl = round(10); % Hz
    Fn = SampleRate*SupersamplingRatio/2;
    N = 2;
    [Bhigh, Ahigh] = butter(N,Fl/Fn,'high'); % Artifact highpass filter
    
    % LFP parameters
    LFPsupersamplingRatio = lcm(SampleRate, LFPsamplerate) / SampleRate;
    LFPdownsamplingRatio = SampleRate * LFPsupersamplingRatio / LFPsamplerate;

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
        % remove 50Hz from data
        % band stop 50Hz
        trialData = filtfilt(Bbs, Abs, trialData);
        
        % calculate LFP
        % LFP filtering
        LFPdata = filtfilt(BLFP, ALFP, trialData);
        % LFP supersampling
        if LFPsupersamplingRatio > 1
            datalength = size(LFPdata,1);
            x = 1:datalength;
            xx = 1:(1/LFPsupersamplingRatio):datalength;
            LFPdata = interp1(x,LFPdata,xx,'linear');
            %     LFPdata = spline(x, LFPdata, xx);
        end
        % LFP downsampling
        LFPdata = downsample(LFPdata,LFPdownsamplingRatio);
        LFP.data = LFPdata;
        LFP.LFPsamplerate = LFPsamplerate;
        LFP.LFPlowpassFreq = LFPlowpassFreq;
        
        trialData = filtfilt(Bhigh, Ahigh,trialData);
        
        
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
        MUAe = GetMUAewithInterpolationSubroutine(trialData,SampleRate,MUAparameters,AllArtifactTimepoints,interpolationLength);
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