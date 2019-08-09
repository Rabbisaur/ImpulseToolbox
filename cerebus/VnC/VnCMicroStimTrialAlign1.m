function instanceinfo = VnCMicroStimTrialAlign(instanceinfo,matdatapath,StimulationParam)

if ispc
    slash = '\';
else
    slash = '/';
end

waveformbeforepoint = 30;
waveformafterpoint = 30;

% stimulation frequency
StimulationFreq = 1/(StimulationParam.stimulationwaveformwidth/1000);
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
FreqRange = 200; %Hz

load(matdatapath);

SampleRate = instanceinfo(1).samplerate;

numberofStimPulses=StimulationParam.numberofStimPulses;

% stimulationwaveformwidthDP=StimulationParam.stimulationwaveformwidthDP; % 0.5ms


%Calculate trial info
% trialstartTime = -initialtrialAlignTime;
% trialendTime = trialstartTime + initialtriallength;
% trialstartIdx = trialstartTime * SampleRate;
% trialendIdx = trialendTime * SampleRate;
% trialidx = trialstartIdx:(trialendIdx-1);
% trialAlignIdx = initialtrialAlignTime * SampleRate;

%calculate stimulation timing profile
% idx = 0:(stimulationwaveformwidthDP-1);
% stimwaveidx = repmat(idx,numberofStimPulses,1);


% build filters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % design the filters
% % try to remove 50Hz from the signal
disp('Building filters...')
d50Hz = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    'DesignMethod','butter','SampleRate',SampleRate);

% low pass filter
lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',StimulationFreq/4,...
    'DesignMethod','butter','SampleRate',SampleRate);

% high pass filter
hiFilt = designfilt('highpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',StimulationFreq/2,...
    'DesignMethod','butter','SampleRate',SampleRate);


% band pass filter
bpFilt2 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',StimulationFreq - FreqRange,'HalfPowerFrequency2',StimulationFreq + FreqRange, ...
    'DesignMethod','butter','SampleRate',SampleRate);

bpFilt3 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',StimulationPulsesFrequency/2,'HalfPowerFrequency2',StimulationPulsesFrequency*2, ...
    'DesignMethod','butter','SampleRate',SampleRate);



% split trials only keep good trials for each electrode
disp('Aligning microstimulation trials and removing stimulation artifacts...')
numInstances = numel(instanceinfo);
skipInstanceFlag = zeros(numInstances,1);
for thisInstance = 1:numInstances
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances),'...'])
    numElec = instanceinfo(thisInstance).numElec;
    numNeuralElec = sum(instanceinfo(thisInstance).ElecOrder<=128);
    % load electrode data back and look for the correct trials and put
    % the trials into a matrix
    numDP = instanceinfo(thisInstance).trialInfo.NumDP;
    trialStartST = instanceinfo(thisInstance).trialInfo.validBeginST;
    trialendST = instanceinfo(thisInstance).trialInfo.validEndST;
    numTrials = numel(trialStartST);
    timeStamp0 = instanceinfo(thisInstance).trialInfo.Timestamp;
    
    SampleRate = double(instanceinfo(thisInstance).samplerate);
    ISIthresholdLow = 1/StimulationPulsesFrequency * SampleRate/2;
    ISIthresholdHi = 1/StimulationPulsesFrequency * SampleRate*3;
    
    % find the best electrode for alignment
    disp('  Searching for the best electrode for artifact alignment...')
    
    validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
    idx = validAllCurrentLevel > 0; % microstimulation trial index
    NumMicroStimTrials = sum(idx);
    
    trialStimEnergy = zeros(NumMicroStimTrials,numElec);
    for thisElec = 1:numElec
        EID = instanceinfo(thisInstance).ElecOrder(thisElec);
        
        % exclude non-neural data
        if EID > 128
            break
        end
        % load data
        elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{thisElec},'r');
        tmpdata = fread(elecfp,numDP,'int16=>double');
        fclose(elecfp);
        
        % create a cell contanning neuronal data in all trials starting
        % from the end of last trial to the end of current trial
        trialrawdata = cell(numTrials,1);
        for thisTrial = 1:numTrials
            if trialStartST(thisTrial) >= trialendST(thisTrial)
                error('Trial start time later than end time or trial has no length')
            end
            idx = (trialStartST(thisTrial):trialendST(thisTrial))-timeStamp0;
            trialrawdata{thisTrial} = tmpdata(idx);
        end
        
        
        % store the raw trial data on disk
        idx = strfind(instanceinfo(thisInstance).electrodeCachePath{thisElec},[slash,'tmp',slash]);
        electroderawTrialDataSavePath = [instanceinfo(thisInstance).electrodeCachePath{thisElec}(1:idx),'electroderawTrialData.mat'];
        electroderawTrialData.samplerate = SampleRate;
        electroderawTrialData.elecID = instanceinfo(thisInstance).ElecOrder(thisElec);
        electroderawTrialData.elecIdx = thisElec;
        electroderawTrialData.trialdata = trialrawdata;
        save(electroderawTrialDataSavePath,'electroderawTrialData');
        
        
        % correct trial alignment for microstimulation trials
        % 1. find the first stimulation artifact
        % 2. align to the center of the first artifact
        
        validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
        idx = validAllCurrentLevel > 0; % microstimulation trial index
        NumMicroStimTrials = sum(idx);
        
        trialrawdata = trialrawdata(idx); % keep only microstimulation trials
        
        for thisTrial = 1:NumMicroStimTrials
            tmpdata = trialrawdata{thisTrial};
            % remove 50 Hz
            tmpdata = filtfilt(d50Hz,tmpdata);
            
            % band pass filter to get the stimulation energy
            bpdata = filtfilt(bpFilt2,tmpdata);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % recalculate the align point based on stimulation artifact
            % 1. Align the trials based on time of the target bit.
            % 2. Remove 50Hz noise from the data
            % 3. Band pass filter 750-5000Hz
            % 4. Find the electrode that has the largest stimulation artifact
            % 5. Find a threshold that can just capture all the 50 artifact pulses within a trial
            % 6. Mark the first transgression event as the point of alignment. This can also be
            %    considered as the begining of the micro-stimulation
            % 7. Re-align all the trials based on the new marker
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find all
            
%             y = fft(bpdata);
%             f = (0:length(bpdata)-1)*SampleRate/length(bpdata);
%             idx = f > StimulationFreq - FreqRange & f < StimulationFreq + FreqRange;
%             y = abs(y);
%             trialStimEnergy(thisTrial,thisElec) = max(y(idx,:),[],1);
            trialStimEnergy(thisTrial,thisElec) = max(bpdata) - min(bpdata);
            %             subplot(1,2,1),plot(f,abs(y))
            %             subplot(1,2,2),plot(tmpdata,'black'),hold on
            %             plot(bpdata,'red')
        end
    end
    %     [~,idxTrialStimEnergy] = sort(trialStimEnergy,2,'descend');
    %     idx = 1:numNeuralElec;
    %     for thisTrial = 1:NumMicroStimTrials
    %         rankStimEnergy(thisTrial,idxTrialStimEnergy(thisTrial,:)) = idx;
    %     end
    %     rankStimEnergy = zeros(NumMicroStimTrials,numNeuralElec);
    trialMeanStimEnergy = min(trialStimEnergy,[],1);
    [~,elecidx] = max(trialMeanStimEnergy);
    
    disp('  Fine tuning alignment...')
    
    % recalculate the alignment point for all microstimulation trials
    % load data
    elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{elecidx},'r');
    tmpdata = fread(elecfp,numDP,'int16=>double');
    fclose(elecfp);
    
    % create a cell containing neuronal data in all trials starting
    % from the end of last trial to the end of current trial
    trialrawdata = cell(numTrials,1);
    for thisTrial = 1:numTrials
        if trialStartST(thisTrial) >= trialendST(thisTrial)
            error('Trial start time later than end time or trial has no length')
        end
        idx = (trialStartST(thisTrial):trialendST(thisTrial))-timeStamp0;
        trialrawdata{thisTrial} = tmpdata(idx);
    end
    
    validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
    
    idx = validAllCurrentLevel > 0; % microstimulation trial index
    DetectableCurrentLevel = validAllCurrentLevel(idx);
    
    NumMicroStimTrials = sum(idx);
    
    trialrawdata = trialrawdata(idx); % keep only microstimulation trials
    TrialAlignOffset = zeros(NumMicroStimTrials,1);
    trialArtifactTimeIdx = zeros(numberofStimPulses,NumMicroStimTrials);
    for thisTrial = 1:NumMicroStimTrials
        tmpdata = trialrawdata{thisTrial};
        
        
        % hi pass to remove low frequency energy
        tmpdata2 = filtfilt(hiFilt,tmpdata);
        
        % remove 50 Hz
        tmpdata2 = filtfilt(d50Hz,tmpdata2);
        
        % band pass to extract the artifact
        tmpdata2 = filtfilt(bpFilt2,tmpdata2);
        % low pass filter
        
        %
        % full wave rectify
%         tmpdata3 = tmpdata.*(tmpdata>0);
        tmpdata2 = abs(tmpdata2);
        tmpdata2 = filtfilt(hiFilt,tmpdata2);
        tmpdata2 = filtfilt(lpFilt,tmpdata2);
        tmpdata2(1:(20*(waveformbeforepoint+waveformafterpoint))) = mean(tmpdata2);
        tmpdata2((end-20*(waveformbeforepoint+waveformafterpoint)):end) = mean(tmpdata2);
        
        threshold = 4* std(tmpdata2);
        
        idx = find(tmpdata2 > threshold,1,'first');
        tmpdata2(1:(idx-numberofStimPulses*waveformbeforepoint)) = tmpdata2(1:(idx-numberofStimPulses*waveformbeforepoint))*0.01;
        idx = find(tmpdata2 > threshold,1,'last');
        tmpdata2((idx+numberofStimPulses*waveformbeforepoint):end) = tmpdata2((idx+numberofStimPulses*waveformbeforepoint):end) * 0.01;
        
        %    tmpdata2 = filtfilt(bpFilt3,tmpdata2);
        tmpdata2 = tmpdata2.^3;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % search for the all stimulation artifact in a trial
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find
%         thresholdhi = max(tmpdata3);
%         thresholdlow = std(tmpdata3);
%         
%         if thresholdhi - thresholdlow < 256
%             tmpdata3 = tmpdata3 * (256 / (thresholdhi - thresholdlow) * 2);
%             thresholdhi = max(tmpdata3);
%             thresholdlow = std(tmpdata3);
%         end
        tmpdata2 = tmpdata2 - mean(tmpdata2);
        tmpdata2 = tmpdata2 / max(tmpdata2) * 256;
        thresholdhi = max(tmpdata2);
        thresholdlow = mean(tmpdata2);
        
        numTransEvent = 0;
        while(numTransEvent < numberofStimPulses && thresholdhi > thresholdlow)
            idx = tmpdata2 >= thresholdhi;
            idx = idx - [idx(2:end);0];
            eventidx = find(idx == -1);
            eventidx = eventidx + 1;
            % remove extra event based on ISI
            numlowISI = 1;
            while(numlowISI > 0)
                ISI = diff(eventidx);
                lowISIidx = find(ISI < ISIthresholdLow,1,'first');
                eventidx(lowISIidx) = []; % remove low ISI spikes;
                numlowISI = numel(lowISIidx);
            end
            
            numHiISI = 1;
            while(numHiISI > 0)
                ISI = diff(eventidx);
                HiISIidx = find(ISI > ISIthresholdHi,1,'first');
                eventidx(HiISIidx) = []; % remove Hi ISI spikes;
                numHiISI = numel(HiISIidx);
            end
            
            numTransEvent = numel(eventidx);
            thresholdhi = thresholdhi-1;
        end
        %         threshold = thresholdhi;
        %         if numTransEvent ~= numberofStimPulses % try search from bottom to top instead
        %             numTransEvent = 0;
        %             while(numTransEvent < numberofStimPulses && thresholdlow < thresholdhi)
        %                 idx = tmpdata2 <= thresholdlow;
        %                 idx = idx - [idx(2:end);0];
        %                 eventidx = find(idx == -1);
        %                 eventidx = eventidx + 1;
        %                 % remove extra event based on ISI
        %                 numlowISI = 1;
        %                 while(numlowISI > 0)
        %                     ISI = diff(eventidx);
        %                     lowISIidx = find(ISI < ISIthreshold,1,'first');
        %                     eventidx(lowISIidx) = []; % remove low ISI spikes;
        %                     numlowISI = numel(lowISIidx);
        %                 end
        %                 numTransEvent = numel(eventidx);
        %                 thresholdlow = thresholdlow+1;
        %             end
        %             threshold = thresholdlow;
        %         end
        
        % if number of TransEvent found still not equal to number of stim pulses delievered
        if numTransEvent ~= numberofStimPulses
            %             %             error('Could not automatically find all stimulation artifacts')
            %             skipInstanceFlag(thisInstance) = 1;
            %             disp(['Could not find correct number of stimulation artifacts on instance ',num2str(thisInstance),'/',num2str(numInstances), '. Skipping...'])
            %             break
            DetectableCurrentLevel(thisTrial) = -1;
            TrialAlignOffset(thisTrial) = 0;
            trialArtifactTimeIdx(:,thisTrial) = zeros(numberofStimPulses,1);
        else
            % find the waveform of the first artifact
            
            for thisEvent = 1:numTransEvent
                wave = tmpdata((eventidx(thisEvent)-waveformbeforepoint) : (eventidx(thisEvent)+waveformafterpoint));
                
                [~,idx1] = max(wave);
                [~,idx2] = min(wave);
                idx1 = idx1(1);
                idx2 = idx2(1);
                idx = floor(mean([idx1,idx2]));
                eventidx(thisEvent) = eventidx(thisEvent) +idx-waveformbeforepoint;
            end
            TrialAlignOffset(thisTrial) = eventidx(1);
            % save in instanceinfo
            trialArtifactTimeIdx(:,thisTrial) = eventidx - eventidx(1);
            %         figure
            %         hold off
            %         plot(tmpdata2)
            %         hold on
            %         plot(eventidx,repmat(threshold,numel(eventidx),1),'ro')
        end
    end
    % remove undetectable trials
    idx = DetectableCurrentLevel == -1;
    TrialAlignOffset(idx) = [];
    trialArtifactTimeIdx(:,idx) = [];
    if sum(DetectableCurrentLevel == -1) == numel(DetectableCurrentLevel)
        skipInstanceFlag(thisInstance) = 1;
        disp(['Could not find correct number of stimulation artifacts on instance ',num2str(thisInstance),'/',num2str(numInstances), '. Skipping...'])
    end
    idx = validAllCurrentLevel > 0;
    DetectableAllCurrentLevel = validAllCurrentLevel;
    DetectableAllCurrentLevel(idx) = DetectableCurrentLevel;
    % update the trial alignment point for micro-stimulation trials.
    idx = DetectableAllCurrentLevel > 0; % microstimulation trial index
    validCorrectedTrialAlignST = instanceinfo(thisInstance).trialInfo.validBeginST;
    validCorrectedTrialAlignST(idx) = validCorrectedTrialAlignST(idx) + TrialAlignOffset;
    instanceinfo(thisInstance).trialInfo.validCorrectedTrialAlignST = validCorrectedTrialAlignST;
    instanceinfo(thisInstance).trialInfo.trialArtifactIdx = trialArtifactTimeIdx;
    instanceinfo(thisInstance).trialInfo.trialAlignOffset = TrialAlignOffset;
    instanceinfo(thisInstance).trialInfo.DetectableAllCurrentLevel = DetectableAllCurrentLevel;
end

% if skip instance = number of instance then we failed
if sum(skipInstanceFlag) == numel(skipInstanceFlag)
    error('We cannot find reliable align of artifact on all instances')
elseif sum(skipInstanceFlag) > 0
    goodIdx = find(skipInstanceFlag == 0, 1,'first');
    badIdx = find(skipInstanceFlag == 1);
    numBadInstance = numel(badIdx);
    % copy answer from the good instance
    for thisBadInstance = 1:numBadInstance
        instancIdx = badIdx(thisBadInstance);
        % copy answer from the good
        TrialAlignOffset = instanceinfo(goodIdx).trialInfo.trialAlignOffset;
        trialArtifactTimeIdx = instanceinfo(goodIdx).trialInfo.trialArtifactIdx;
        %         validAllCurrentLevel = instanceinfo(instancIdx).trialInfo.validAllCurrentLevel;
        %         idx = validAllCurrentLevel > 0; % microstimulation trial index
        DetectableAllCurrentLevel = instanceinfo(goodIdx).trialInfo.DetectableAllCurrentLevel;
        idx = DetectableAllCurrentLevel > 0; % microstimulation trial index
        validCorrectedTrialAlignST = instanceinfo(instancIdx).trialInfo.validBeginST;
        validCorrectedTrialAlignST(idx) = validCorrectedTrialAlignST(idx) + TrialAlignOffset;
        instanceinfo(instancIdx).trialInfo.DetectableAllCurrentLevel = DetectableAllCurrentLevel;
        instanceinfo(instancIdx).trialInfo.validCorrectedTrialAlignST = validCorrectedTrialAlignST;
        instanceinfo(instancIdx).trialInfo.trialArtifactIdx = trialArtifactTimeIdx;
        instanceinfo(instancIdx).trialInfo.trialAlignOffset = TrialAlignOffset;
    end
end