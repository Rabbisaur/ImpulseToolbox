function instanceinfo = VnCArtifactRemovalTTLonly(instanceinfo,StimulationParam,MUAparameters,LFPparameters,trialParam,trialAlignSTfieldname,interpolationLength,ARFineAlignFlag)
disp('Removing stimulation artifacts...')

matversion = '-v6';

interpolationLength0 = interpolationLength;

if ispc
    slash = '\';
else
    slash = '/';
end

% parameters
SampleRate = instanceinfo.samplerate;
interpolationLengthDPupperLimit = round((1/StimulationParam.PulsesFrequency)*3/5 * SampleRate);

LFPdownsamplingRatio = SampleRate/LFPparameters.LFPsamplingrate;
% low pass filter 150Hz
Fl=150;
Fn = SampleRate/2;
N = 2;
[BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass


numElec = instanceinfo.numElec;
SampleRate = double(instanceinfo.samplerate);

% remove artifact
fprintf('\n')
disp('Removing artifacts')

StimStartTime = round(-trialParam.startT * SampleRate);
% remove artifact electrode by electrode
for thisElec = 1:min([numElec,128])
    fprintf('.')
    EID = instanceinfo.ElecOrder(thisElec);
    disp(['EID=',num2str(EID)])
    % load in data
    cmd = ['trialAlignST = instanceinfo.trialInfo.',trialAlignSTfieldname,';'];
    eval(cmd);
    electrodeTrialData = VnCGetTrialData2(instanceinfo, 1, EID,trialAlignST,trialParam);
    trialData = electrodeTrialData.trialData';
    trialData0 = trialData;
    GoodMStrialIdx = true(numel(trialData,2),1);
    GoodTrialIdx = GoodMStrialIdx;
    [trialData, AllArtifactTimepoints,Interpolateddata,WaveDPbefore,WaveDPafter] = VnCRemoveArtifactsSubroutineTTLonly(trialData,StimulationParam,StimStartTime,SampleRate,EID,1,instanceinfo,GoodMStrialIdx,ARFineAlignFlag);
    
    % search for best interpolation length
    baselineMUAe = GetMUAewithInterpolationSubroutine(trialData0(:,GoodTrialIdx),SampleRate,MUAparameters,AllArtifactTimepoints,0);
    if size(baselineMUAe.data,1) == 1 && size(baselineMUAe.data,2)>size(baselineMUAe.data,1)
            baselineMUAe.data = baselineMUAe.data';
    end
    baselineMUAe.time =  (0:(size(baselineMUAe.data,1)-1))/baselineMUAe.MUAesamplerate + trialParam.startT;
    
    counter = 0;
    
    
    numelInterpolationLength = numel(round(WaveDPbefore+WaveDPafter):interpolationLengthDPupperLimit);
    tmpdata = zeros(numel(baselineMUAe.time),numelInterpolationLength);
    for interpolationLengthDP = round(WaveDPbefore+WaveDPafter):interpolationLengthDPupperLimit
        fprintf('.')
        counter = counter + 1;
        interpolationLength = interpolationLengthDP /SampleRate;
        testMUAe = GetMUAewithInterpolationSubroutine(trialData(:,GoodTrialIdx),SampleRate,MUAparameters,AllArtifactTimepoints,interpolationLength);
        if size(testMUAe.data,1) == 1 && size(testMUAe.data,2)>size(testMUAe.data,1)
            testMUAe.data = testMUAe.data';
        end
        testMUAe.time =  (0:(size(testMUAe.data,1)-1))/testMUAe.MUAesamplerate + trialParam.startT;
        
        tmpdata(:,counter) = smooth(mean(testMUAe.data,2),20);
    end
    fprintf('\n')
    
    idx = testMUAe.time < -0.05 & testMUAe.time > testMUAe.time(5);
    meanBase = smooth(mean(baselineMUAe.data,2),20);
    %         testDiff = tmpdata(idx,:)-repmat(meanBase(idx),1,size(tmpdata,2));
    %         testDiff = sum((testDiff.^2),1);
    meanDiff = (mean(tmpdata(idx,:),1)-mean(meanBase(idx))).^2;
    baseLine = meanBase(idx);
    testData = tmpdata(idx,:);
    
    
    metric3 = corr(testData,repmat(baseLine,1,size(testData,2)));
    metric3 = metric3(1:end,1);
    
    meanDiff = (meanDiff-min(meanDiff))/(max(meanDiff)-min(meanDiff));
    idx = find(meanDiff<0.05,1,'first');
    idx = idx(1);
    [~,idx2] = max(metric3);
    idx2 = idx2(1);
    interpolationLengthDP = round(WaveDPbefore+WaveDPafter):interpolationLengthDPupperLimit;
    levelInterpolationLengthDP = interpolationLengthDP(idx);
    interpolationLengthDP = interpolationLengthDP(idx2);
    ARInterpolationLengthDP = round(WaveDPbefore+WaveDPafter);
    
    
    
    h = figure('position',[100 100 900 700],'visible','off');
    subplot(1,2,1)
    hold on
    for thisLength = 1:size(tmpdata,2)
        plot(baselineMUAe.time,tmpdata(:,thisLength),'color',[1,0,1]*thisLength/size(tmpdata,2))
    end
    plot(baselineMUAe.time,meanBase,'black','linewidth',2)
    plot(baselineMUAe.time,tmpdata(:,idx),'red','linewidth',2)
    plot(baselineMUAe.time,tmpdata(:,idx2),'blue','linewidth',2)
    xlabel('Time from stim onset (s)')
    ylabel('MUAe')
    subplot(1,2,2)
    hold on
    xlabel('Interpolation samples')
    %         plot(round(WaveDPbefore+WaveDPafter):interpolationLengthDPupperLimit,(testDiff-min(testDiff))/(max(testDiff)-min(testDiff)),'red')
    h1=plot(round(WaveDPbefore+WaveDPafter):interpolationLengthDPupperLimit,meanDiff,'red');
    h2=plot(round(WaveDPbefore+WaveDPafter):interpolationLengthDPupperLimit,(metric3-min(metric3))/(max(metric3)-min(metric3)),'blue');
    %         legend({'MeanDiff','Correlation','Decision on mean','Decision on correlation'})
    plot([idx idx]+round(WaveDPbefore+WaveDPafter)-1+0.2,[0 1],'color','red','linewidth',2)
    plot([idx2 idx2]+round(WaveDPbefore+WaveDPafter)-1-0.2,[0 1],'color','blue','linewidth',2)
    legend([h1,h2],{'MeanDiff','Correlation'})
    
    idx = find(instanceinfo.electrodeCachePath{EID} == slash,3,'last');
    idx = idx(1);
    dirpath = [instanceinfo.electrodeCachePath{EID}(1:idx),'debugimg'];
    if exist(dirpath,'dir')==7
        
    else
        mkdir(dirpath)
    end
    
    savepath = [dirpath,'/',num2str(EID),'MUAinterpolation.png'];
    saveas(h,savepath)
    %         savepath = [dirpath,'/',num2str(EID),'MUAinterpolation.fig'];
    %         saveas(h,savepath)
    close(h)
    
    disp(['MUA shape preserving interpolation Length = ',num2str(interpolationLengthDP/SampleRate*1000),'ms'])
    
    % calculate MUA
    interpolationLength = interpolationLengthDP /SampleRate;
    if interpolationLength < interpolationLength0 % minimal
        interpolationLength = interpolationLength0;
    end
    
    
    MUAe = GetMUAewithInterpolationSubroutine(trialData,SampleRate,MUAparameters,AllArtifactTimepoints,interpolationLength);
    MUAe.time =  (0:(size(MUAe.data,1)-1))/MUAe.MUAesamplerate + trialParam.startT;
    electrodeArtifactRemovedNeuralData.MUAeMeanLevelInterpolationLength = levelInterpolationLengthDP * 1/SampleRate;
    electrodeArtifactRemovedNeuralData.MUAeCorrelationInterpolationLength = MUAe.interpolationLength;
    electrodeArtifactRemovedNeuralData.MUAeARInterpolationLength = ARInterpolationLengthDP * 1/SampleRate;
    electrodeArtifactRemovedNeuralData.MUAe = MUAe;
    
    % MUA mean level matching interpolation length
    disp(['MUA mean matching interpolation Length = ',num2str(electrodeArtifactRemovedNeuralData.MUAeMeanLevelInterpolationLength*1000),'ms'])
    MUAe = GetMUAewithInterpolationSubroutine(trialData,SampleRate,MUAparameters,AllArtifactTimepoints,electrodeArtifactRemovedNeuralData.MUAeMeanLevelInterpolationLength);
    MUAe.time =  (0:(size(MUAe.data,1)-1))/MUAe.MUAesamplerate + trialParam.startT;
    electrodeArtifactRemovedNeuralData.MUAeMeanLevelInterpolation = MUAe;
    
    % plot shape preserved MUAe
    MUAe = electrodeArtifactRemovedNeuralData.MUAe;
    % calculate LFP
    LFP.data = filtfilt(BLFP,ALFP,Interpolateddata);
    LFP.data = downsample(LFP.data,LFPdownsamplingRatio);
    LFP.lpFilterFreq = Fl;
    LFP.sampleRate = LFPparameters.LFPsamplingrate;
    LFP.time =(0:(size(LFP.data,1)-1))/LFPparameters.LFPsamplingrate + trialParam.startT;
    electrodeArtifactRemovedNeuralData.LFP = LFP;
    
    if true
        h = figure('position',[0 0 2048 1440],'visible','off');
        subplot(2,2,1)
        hold on
        plot(mean(trialData(:,GoodTrialIdx),2))
        plot(mean(trialData(:,GoodMStrialIdx),2))
        plot(mean(trialData(:,(GoodTrialIdx&(~GoodMStrialIdx))),2))
        plot(AllArtifactTimepoints,max(mean(trialData(:,GoodTrialIdx),2))/2*ones(size(AllArtifactTimepoints)),'r.')
        xlim([(abs(trialParam.startT)-0.05)*SampleRate, (abs(trialParam.startT)+0.17)*SampleRate])
        xlabel('Sample points')
        ylabel('Voltage')
        subplot(2,2,2)
        hold on
        plot(mean(Interpolateddata(:,GoodTrialIdx),2))
        plot(mean(Interpolateddata(:,GoodMStrialIdx),2))
        plot(mean(Interpolateddata(:,(GoodTrialIdx&(~GoodMStrialIdx))),2))
        plot(AllArtifactTimepoints,max(mean(trialData,2))/2*ones(size(AllArtifactTimepoints)),'r.')
        xlim([(abs(trialParam.startT)-0.05)*SampleRate, (abs(trialParam.startT)+0.17)*SampleRate])
        xlabel('Sample points')
        ylabel('Voltage')
        subplot(2,2,3)
        hold on
        plot(MUAe.time,mean(MUAe.data(:,GoodTrialIdx),2),'linewidth',1.5)
        plot(MUAe.time,mean(MUAe.data(:,GoodMStrialIdx),2),'linewidth',1.5)
        plot(MUAe.time,mean(MUAe.data(:,(GoodTrialIdx&(~GoodMStrialIdx))),2),'linewidth',1.5)
        legend({'AllTrials','MS trials','Catch trials'})
        xlim([-0.05 0.17])
        xlabel('Time from stim onset (s)')
        ylabel('MUAe')
        subplot(2,2,4)
        hold on
        plot(LFP.time,mean(LFP.data(:,GoodTrialIdx),2),'linewidth',1.5)
        plot(LFP.time,mean(LFP.data(:,GoodMStrialIdx),2),'linewidth',1.5)
        plot(LFP.time,mean(LFP.data(:,(GoodTrialIdx&(~GoodMStrialIdx))),2),'linewidth',1.5)
        xlim([-0.05 0.17])
        legend({'AllTrials','MS trials','Catch trials'})
        xlabel('Time from stim onset (s)')
        ylabel('LFP')
        
        % save figure
        idx = find(instanceinfo.electrodeCachePath{EID} == slash,3,'last');
        idx = idx(1);
        dirpath = [instanceinfo.electrodeCachePath{EID}(1:idx),'debugimg'];
        if exist(dirpath,'dir')==7
        else
            mkdir(dirpath)
        end
        %             if EID == 49
        %                 a = 1
        %             end
        savepath = [dirpath,'/',num2str(EID),'result.png'];
        saveas(h,savepath)
        savepath = [dirpath,'/',num2str(EID),'result.fig'];
        saveas(h,savepath)
        close(h)
    end
    
    electrodeArtifactRemovedNeuralData.samplerate = SampleRate;
    electrodeArtifactRemovedNeuralData.elecID = EID;
    electrodeArtifactRemovedNeuralData.elecIdx = thisElec;
    idx = find(instanceinfo.electrodeCachePath{EID} == slash,2,'last');
    idx = idx(1);
    savefilename = 'electrodeArtifactRemovedNeuralData.mat';
    savepath = [instanceinfo.electrodeCachePath{EID}(1:idx),savefilename];
    save(savepath,'electrodeArtifactRemovedNeuralData',matversion);
    instanceinfo.electrodeArtifactRemovedNeuralDataPath{thisElec} = savepath;
end