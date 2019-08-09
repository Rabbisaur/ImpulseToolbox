function instanceinfo = VnCArtifactRemoval9(instanceinfo,StimulationParam,MUAparameters,LFPparameters,trialParam,trialAlignSTfieldname,interpolationLength,DebugFlag)
disp('Removing stimulation artifacts...')

matversion = '-v6';
interpolationRatio = 1.33;
interpolationLength0 = interpolationLength;

if ispc
    slash = '\';
else
    slash = '/';
end

% parameters
% SupersamplingRatio = 1;
numInstances = numel(instanceinfo);

SampleRate = instanceinfo(1).samplerate;
LFPdownsamplingRatio = SampleRate/LFPparameters.LFPsamplingrate;

% low pass filter 150Hz
Fl=150;
Fn = SampleRate/2;
N = 2;
[BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass

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
        disp(['EID=',num2str(EID)])
        % load in data
        cmd = ['trialAlignST = instanceinfo(thisInstance).trialInfo.',trialAlignSTfieldname,';'];
        eval(cmd);
        electrodeTrialData = VnCGetTrialData2(instanceinfo, thisInstance, EID,trialAlignST,trialParam);
        trialData = electrodeTrialData.trialData';
        [trialData, AllArtifactTimepoints,Interpolateddata,WaveDPbefore,WaveDPafter] = VnCRemoveArtifactsSubroutine7(trialData,StimulationParam,StimStartTime,SampleRate);
        
        % calculate MUA
        interpolationLength = (WaveDPafter + WaveDPbefore) * interpolationRatio / SampleRate;
        if interpolationLength < interpolationLength0 % minimal
            interpolationLength = interpolationLength0;
        end
        MUAe = GetMUAewithInterpolationSubroutine(trialData,SampleRate,MUAparameters,AllArtifactTimepoints,interpolationLength);
        MUAe.time =  (0:(size(MUAe.data,1)-1))/MUAe.MUAesamplerate + trialParam.startT;
        electrodeArtifactRemovedNeuralData.MUAe = MUAe;
        
        % calculate LFP
        LFP.data = filtfilt(BLFP,ALFP,Interpolateddata);
        LFP.data = downsample(LFP.data,LFPdownsamplingRatio);
        LFP.time =(0:(size(LFP.data,1)-1))/LFPparameters.LFPsamplingrate + trialParam.startT;
        electrodeArtifactRemovedNeuralData.LFP = LFP;
        
        if true
            h = figure('position',[0 0 2048 1440]);
            subplot(2,2,1)
            hold on
            plot(mean(trialData,2))
            plot(AllArtifactTimepoints,max(mean(trialData,2))/2*ones(size(AllArtifactTimepoints)),'r.')
            xlim([(abs(trialParam.startT)-0.05)*SampleRate, (abs(trialParam.startT)+0.17)*SampleRate])
            subplot(2,2,2)
            hold on
            plot(mean(Interpolateddata,2))
            plot(AllArtifactTimepoints,max(mean(trialData,2))/2*ones(size(AllArtifactTimepoints)),'r.')
            xlim([(abs(trialParam.startT)-0.05)*SampleRate, (abs(trialParam.startT)+0.17)*SampleRate])
            subplot(2,2,3)
            plot(MUAe.time,mean(MUAe.data,2))
            xlim([-0.05 0.17])
            subplot(2,2,4)
            plot(LFP.time,mean(LFP.data,2))
            xlim([-0.05 0.17])
            
            idx = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,3,'last');
            idx = idx(1);
            dirpath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx),'debugimg'];
            if exist(dirpath,'dir')==7
            else
                mkdir(dirpath)
            end
%             if EID == 49
%                 a = 1
%             end
            savepath = [dirpath,'/',num2str(EID),'.png'];
            saveas(h,savepath)
            savepath = [dirpath,'/',num2str(EID),'.fig'];
            saveas(h,savepath)
            close(h)
        end
        
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