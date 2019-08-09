function instanceinfo = VnCMicroStimWaveformFineAlign(instanceinfo,StimulationParam,trialParam,trialAlignSTfieldname,GoodMStrialIdx)
disp('Performing fine aligment...')

if ispc
    slash = '\';
else
    slash = '/';
end

% parameters
SupersamplingRatio = 16;
numInstances = numel(instanceinfo);

for thisInstance = 1:numInstances
    fprintf('\n')
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances)])
    numElec = instanceinfo(thisInstance).numElec;
    SampleRate = double(instanceinfo(thisInstance).samplerate);

    % Fine alignment
    fprintf('\n')
    disp('Performing fine aligment')

    StimStartTime = round(-trialParam.startT * SampleRate);
    % remove artifact electrode by electrode
    for thisElec = 1:min([numElec,128])
        fprintf('.')
        EID = instanceinfo(thisInstance).ElecOrder(thisElec);
        disp(['EID=',num2str(EID)])
        % load in data
        cmd = ['trialAlignST = instanceinfo(thisInstance).trialInfo.',trialAlignSTfieldname,';'];
        eval(cmd);
        electrodeTrialData = VnCGetTrialData2(instanceinfo, thisInstance, EID,trialAlignST,trialParam);
        RawTrace = electrodeTrialData.trialData';
        RawTrace = RawTrace(:,GoodMStrialIdx);
        
        if size(RawTrace,2) > size(RawTrace,1)
            RawTrace = RawTrace';
        end
        % parameters
        StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
        numberofStimPulses = StimulationParam.numberofStimPulses;
        
        Fl = 100;
        Fn = SampleRate/2;
        N = 4;
        [BWave, AWave] = butter(N,Fl/Fn,'high'); % Waveform hi pass
        HipassTrace = filtfilt(BWave, AWave,RawTrace);
        
        
        
        NumTrials = size(RawTrace,2);
%         NumDP = size(RawTrace,1);
        
        StimulationPulseInterval = 1/StimulationPulsesFrequency;
        
        
        
        predictedLoc = (0:(numberofStimPulses-1)) * StimulationPulseInterval * SampleRate;
        predictedLoc = predictedLoc + StimStartTime;
        predictedLoc = predictedLoc + StimulationParam.stimulationwaveformwidth*0.95 * SampleRate;
        % predictedLoc = floor(predictedLoc);
        
        loc0 = predictedLoc;
        
        % Find stimulation exact locations, get waveform
%         numSamplePoints = size(HipassTrace,1);
        
        StimWidth = 1/StimulationPulsesFrequency/2;
        idx = (0:round(StimWidth*SampleRate-1))';
        WaveLength = numel(idx);
        idx = idx - floor(max(idx)/2);
        idx1 = repmat(idx,1,numberofStimPulses);
        idx1 = round(idx1 + repmat(loc0,numel(idx),1));
        
        h = figure('position',[100 100 900 700]);
        subplot(3,1,1)
        hold on
        meanWaveform = zeros(WaveLength,NumTrials);
        for thisTrial = 1:NumTrials
            Waveform = HipassTrace(idx1,thisTrial);
            Waveform = reshape(Waveform,numel(idx),numberofStimPulses);
            meanWaveform(:,thisTrial) = mean(Waveform,2);
        end
        
        % supersampling
        x = 1:WaveLength;
        xx = 1:(1/SupersamplingRatio):WaveLength;
        meanWaveform = spline(x,meanWaveform',xx)';
        
        for thisTrial = 1:NumTrials
            meanWaveform(:,thisTrial) = meanWaveform(:,thisTrial) - mean([meanWaveform(1,thisTrial),meanWaveform(end,thisTrial)]);
        end
        
        plot(meanWaveform)
        
        % alignment
        [~,maxIdx] = max(meanWaveform,[],1);
        [~,minIdx] = min(meanWaveform,[],1);
        zerocrossingpoint = minIdx;
        for thisTrial = 1:NumTrials
            [~,zerocrossingpoint(thisTrial)] = min((meanWaveform(min([maxIdx(thisTrial),minIdx(thisTrial)]):max([maxIdx(thisTrial),minIdx(thisTrial)]),thisTrial)-0).^2);
            zerocrossingpoint(thisTrial) = zerocrossingpoint(thisTrial) + min([maxIdx(thisTrial),minIdx(thisTrial)])-1;
            plot(zerocrossingpoint(thisTrial),meanWaveform(zerocrossingpoint(thisTrial),thisTrial),'rx')
        end
        
        middleIdx = median(zerocrossingpoint);
        plot([middleIdx middleIdx],[min(min(meanWaveform)),max(max(meanWaveform))],'black','linewidth',2)
        
        shiftvalue = zerocrossingpoint - middleIdx;
        
        disp(['Deviation = ', num2str((max(shiftvalue)-min(shiftvalue))/SupersamplingRatio/SampleRate*1000),'ms'])
        
        meanWaveformShifted = meanWaveform;
        for thisTrial = 1:NumTrials
            meanWaveformShifted(:,thisTrial) = circshift(meanWaveform(:,thisTrial),-shiftvalue(thisTrial));
        end
        subplot(3,1,2)
        hold on
        plot(meanWaveformShifted)
        [~,maxIdx] = max(meanWaveformShifted,[],1);
        [~,minIdx] = min(meanWaveformShifted,[],1);
        for thisTrial = 1:NumTrials
            [~,zerocrossingpoint(thisTrial)] = min((meanWaveformShifted(min([maxIdx(thisTrial),minIdx(thisTrial)]):max([maxIdx(thisTrial),minIdx(thisTrial)]),thisTrial)-0).^2);
            zerocrossingpoint(thisTrial) = zerocrossingpoint(thisTrial) + min([maxIdx(thisTrial),minIdx(thisTrial)])-1;
            plot(zerocrossingpoint(thisTrial),meanWaveformShifted(zerocrossingpoint(thisTrial),thisTrial),'go')
        end
        middleIdx2 = median(zerocrossingpoint);
        plot([middleIdx2 middleIdx2],[min(min(meanWaveform)),max(max(meanWaveform))],'black','linewidth',2)
        
        shiftvalue0 = shiftvalue;
        shiftvalue = round(shiftvalue / SupersamplingRatio);
        
        % save result
        correctedtrialAlignST = instanceinfo(thisInstance).trialInfo.correctedtrialAlignST;
        correctedtrialAlignST(GoodMStrialIdx) = correctedtrialAlignST(GoodMStrialIdx) + shiftvalue'; % important! use the TrialAlignPointOffset!
        instanceinfo(thisInstance).trialInfo.correctedtrialAlignST = correctedtrialAlignST;
        instanceinfo(thisInstance).trialInfo.fineAlign.shiftvalue0 = shiftvalue0;
        instanceinfo(thisInstance).trialInfo.fineAlign.SupersamplingRatio = SupersamplingRatio;
        instanceinfo(thisInstance).trialInfo.fineAlign.GoodMSTrialIndex = GoodMStrialIdx;
        
        % validation
        cmd = ['trialAlignST = instanceinfo(thisInstance).trialInfo.',trialAlignSTfieldname,';'];
        eval(cmd);
        electrodeTrialData = VnCGetTrialData2(instanceinfo, thisInstance, EID,trialAlignST,trialParam);
        RawTrace = electrodeTrialData.trialData';
        RawTrace = RawTrace(:,GoodMStrialIdx);
        if size(RawTrace,2) > size(RawTrace,1)
            RawTrace = RawTrace';
        end
        HipassTrace = filtfilt(BWave, AWave,RawTrace);
        meanWaveform = zeros(WaveLength,NumTrials);
        for thisTrial = 1:NumTrials
            Waveform = HipassTrace(idx1,thisTrial);
            Waveform = reshape(Waveform,numel(idx),numberofStimPulses);
            meanWaveform(:,thisTrial) = mean(Waveform,2);
        end
        
        % supersampling
        x = 1:WaveLength;
        xx = 1:(1/SupersamplingRatio):WaveLength;
        meanWaveform = spline(x,meanWaveform',xx)';
        
        for thisTrial = 1:NumTrials
            meanWaveform(:,thisTrial) = meanWaveform(:,thisTrial) - mean([meanWaveform(1,thisTrial),meanWaveform(end,thisTrial)]);
        end
        subplot(3,1,3)
        hold on
        plot(meanWaveform)
        
        
        idx = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,3,'last');
        idx = idx(1);
        dirpath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx),'debugimg'];
        if exist(dirpath,'dir')==7
            
        else
            mkdir(dirpath)
        end
        savepath = [dirpath,'/',num2str(EID),'FineAlignment.png'];
        saveas(h,savepath)
        savepath = [dirpath,'/',num2str(EID),'FineAlignment.fig'];
        saveas(h,savepath)
        close(h)
    end
end
end