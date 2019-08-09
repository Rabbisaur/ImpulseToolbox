function instanceinfo = VnCMicroStimTrialAlignRef(instanceinfo,matdatapath,StimulationParam,RefCh,UsingTTLflag,UsingMonitorFlag,electrodeFlag,whichmethod)

if ispc
    slash = '\';
else
    slash = '/';
end

debug = 0;

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

% % hipass > 5Hz
% hipass = designfilt('highpassiir','FilterOrder',4, ...
%         'HalfPowerFrequency',5,...
%         'DesignMethod','butter','SampleRate',SampleRate);


% d50Hz = designfilt('bandstopiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
%     'DesignMethod','butter','SampleRate',SampleRate);
%
% d100Hz = designfilt('bandstopiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',98,'HalfPowerFrequency2',102, ...
%     'DesignMethod','butter','SampleRate',SampleRate);

% low pass filter
% lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',StimulationFreq/4,...
%     'DesignMethod','butter','SampleRate',SampleRate);

% high pass filter
hiFilt = designfilt('highpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',StimulationFreq/10,...
    'DesignMethod','butter','SampleRate',SampleRate);

% hiFilt2 = designfilt('highpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',5,...
%     'DesignMethod','butter','SampleRate',SampleRate);


% band pass filter
% bpFilt2 = designfilt('bandpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',StimulationFreq - FreqRange,'HalfPowerFrequency2',StimulationFreq + FreqRange, ...
%     'DesignMethod','butter','SampleRate',SampleRate);
%
% bpFilt3 = designfilt('bandpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',StimulationPulsesFrequency/2,'HalfPowerFrequency2',StimulationPulsesFrequency*2, ...
%     'DesignMethod','butter','SampleRate',SampleRate);



% split trials only keep good trials for each electrode
disp('Aligning microstimulation trials and removing stimulation artifacts...')
numInstances = numel(instanceinfo);
% skipInstanceFlag = zeros(numInstances,1);

% get the microstimulation pulses from the reference channel
thisInstance = RefCh(1,1);
trialStartST = instanceinfo(thisInstance).trialInfo.validBeginST;
trialendST = instanceinfo(thisInstance).trialInfo.validEndST;
[~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
timeStamp0 = instanceinfo(thisInstance).trialInfo.Timestamp(trialIdx);

SampleRate = double(instanceinfo(thisInstance).samplerate);
ISIthresholdLow = 1/StimulationPulsesFrequency * SampleRate/3;
ISIthresholdHi = 1/StimulationPulsesFrequency * SampleRate*3;

numDP = instanceinfo(thisInstance).trialInfo.NumDP;

numTrials = numel(trialStartST);


SampleRate = double(instanceinfo(thisInstance).samplerate);
validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
if iscell(validAllCurrentLevel)
    validAllCurrentLevel = validAllCurrentLevel{1};
end
idx = validAllCurrentLevel > 0; % microstimulation trial index
NumMicroStimTrials = sum(idx);


if UsingTTLflag == 1
    % load data for the trigger channel
    if RefCh(2,1) == 0 % use sum of all channels
        thisInstance = 1;
        trialStartST = instanceinfo(thisInstance).trialInfo.validBeginST;
        trialendST = instanceinfo(thisInstance).trialInfo.validEndST;
        [~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
        timeStamp0 = instanceinfo(thisInstance).trialInfo.Timestamp(trialIdx);
        arrayRefCh = [1;6;7;8;9;10;11;12;13];
        trialrawdata = cell(numTrials,1);
        for thisElec = 1:numel(arrayRefCh)
            EID = instanceinfo(thisInstance).ElecOrder==128+arrayRefCh(thisElec);
            elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{EID},'r');
            tmpdata = fread(elecfp,numDP,'int16=>double');
            fclose(elecfp);
            % create a cell containing neuronal data in all trials starting
            % from the end of last trial to the end of current trial
            
            for thisTrial = 1:numTrials
                if trialStartST(thisTrial) >= trialendST(thisTrial)
                    error('Trial start time later than end time or trial has no length')
                end
                idx = (trialStartST(thisTrial):trialendST(thisTrial))-timeStamp0;
                if thisElec == 1
                    trialrawdata{thisTrial} = zeros(size(idx));
                end
                tmp = tmpdata(idx);
                if size(tmp,1) > size(tmp,2)
                    tmp = tmp';
                end
                trialrawdata{thisTrial} = trialrawdata{thisTrial} + tmp;
                trialrawdata{thisTrial} = trialrawdata{thisTrial} - min(trialrawdata{thisTrial}); % minus baseline
            end
        end
    else
        
        thisInstance = RefCh(2,1);
        trialStartST = instanceinfo(thisInstance).trialInfo.validBeginST;
        trialendST = instanceinfo(thisInstance).trialInfo.validEndST;
        [~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
        timeStamp0 = instanceinfo(thisInstance).trialInfo.Timestamp(trialIdx);
        
        thisElec = instanceinfo(thisInstance).ElecOrder==RefCh(2,2);
        elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{thisElec},'r');
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
    end
    
    
    
    
    
    % find time point of microstimulation pulses marked in the trigger signal for each microstimulation trial
    validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
    if iscell(validAllCurrentLevel)
        validAllCurrentLevel = validAllCurrentLevel{1};
    end
    idx = validAllCurrentLevel > 0; % microstimulation trial index
    
    NumMicroStimTrials = sum(idx);
    
    trialrawdata = trialrawdata(idx); % keep only microstimulation trials
    
    MicroStimStartEndTime = zeros(2,NumMicroStimTrials);
    Threshold = 20000; % threshold for the TTL signal is 2V
    figure
    hold on
    for thisTrial = 1:NumMicroStimTrials
        tmpdata = trialrawdata{thisTrial};
        tmpdata = tmpdata - mean(tmpdata(1:300));
        idx = find(tmpdata > Threshold, 1, 'first');
        MicroStimStartEndTime(1,thisTrial) = idx;
        idx = find(tmpdata > Threshold, 1, 'last');
        MicroStimStartEndTime(2,thisTrial) = idx;
        plot(tmpdata)
    end
end


% Build electrode trial data matrix

if ~UsingMonitorFlag
    figure
    hold on
    % create a cell containing neuronal data in all trials starting
    % from the end of last trial to the end of current trial
    trialrawdata = cell(numTrials,1);
    % load every neuronal channel from every instance
    electrodesCounted = 0;
    for thisInstance = 1:numInstances
        %     thisInstance = RefCh(1,1);
        trialStartST = instanceinfo(thisInstance).trialInfo.validBeginST;
        trialendST = instanceinfo(thisInstance).trialInfo.validEndST;
        [~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
        timeStamp0 = instanceinfo(thisInstance).trialInfo.Timestamp(trialIdx);
        
        
        numElec = sum(instanceinfo(thisInstance).ElecOrder <= 128);
        for thisElec = 1:numElec
            EID = instanceinfo(thisInstance).ElecOrder(thisElec);
            if electrodeFlag(thisInstance,EID)
                %     thisElec = find(instanceinfo(thisInstance).ElecOrder==RefCh(1,2));
                % load data
                elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{EID},'r');
                tmpdata = fread(elecfp,numDP,'int16=>double');
                fclose(elecfp);
                electrodesCounted = electrodesCounted + 1;
                for thisTrial = 1:numTrials
                    if trialStartST(thisTrial) >= trialendST(thisTrial)
                        error('Trial start time later than end time or trial has no length')
                    end
                    
                    if electrodesCounted == 1
                        idx = (trialStartST(thisTrial):trialendST(thisTrial))-timeStamp0;
                        
                        trialrawdata{thisTrial} = tmpdata(idx);
                    else
                        idx = (trialStartST(thisTrial):(trialStartST(thisTrial) + numel(trialrawdata{thisTrial})-1))-timeStamp0;
                        trialrawdata{thisTrial} = trialrawdata{thisTrial} + tmpdata(idx);
                    end
                    plot(trialrawdata{thisTrial})
                end
            end
        end
    end
    
    % get all correct microstimulation trials
%     validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
%     if iscell(validAllCurrentLevel)
%         validAllCurrentLevel = validAllCurrentLevel{1};
%     end
%     idx = validAllCurrentLevel > 0; % microstimulation trial index
%     DetectableCurrentLevel = validAllCurrentLevel(idx);
%     
%     NumMicroStimTrials = sum(idx);
%     
%     trialrawdata = trialrawdata(idx); % keep only microstimulation trials
%     TrialAlignOffset = zeros(NumMicroStimTrials,1);
%     trialArtifactTimeIdx = zeros(numberofStimPulses,NumMicroStimTrials);
else
    thisInstance = RefCh(1,1);
    trialStartST = instanceinfo(thisInstance).trialInfo.validBeginST;
    trialendST = instanceinfo(thisInstance).trialInfo.validEndST;
    [~,trialIdx] = max(instanceinfo(thisInstance).trialInfo.NumDP);
    timeStamp0 = instanceinfo(thisInstance).trialInfo.Timestamp(trialIdx);
    
    thisElec = find(instanceinfo(thisInstance).ElecOrder==RefCh(1,2));
    % load data
    elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{thisElec},'r');
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
    
%     validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
%     
%     idx = validAllCurrentLevel > 0; % microstimulation trial index
%     DetectableCurrentLevel = validAllCurrentLevel(idx);
%     
%     NumMicroStimTrials = sum(idx);
%     
%     trialrawdata = trialrawdata(idx); % keep only microstimulation trials
%     TrialAlignOffset = zeros(NumMicroStimTrials,1);
%     trialArtifactTimeIdx = zeros(numberofStimPulses,NumMicroStimTrials);
end

% get all correct microstimulation trials
    validAllCurrentLevel = instanceinfo(thisInstance).trialInfo.validAllCurrentLevel;
    if iscell(validAllCurrentLevel)
        validAllCurrentLevel = validAllCurrentLevel{1};
    end
    idx = validAllCurrentLevel > 0; % microstimulation trial index
    DetectableCurrentLevel = validAllCurrentLevel(idx);
    
    NumMicroStimTrials = sum(idx);
    
    trialrawdata = trialrawdata(idx); % keep only microstimulation trials
    TrialAlignOffset = zeros(NumMicroStimTrials,1);
    trialArtifactTimeIdx = zeros(numberofStimPulses,NumMicroStimTrials);

% find stimulation artifacts for each trial
% figure
for thisTrial = 1:NumMicroStimTrials
    tmpdata = trialrawdata{thisTrial};
    
    % hi pass to remove low frequency energy
    tmpdata2 = tmpdata;
    tmpdata2 = filtfilt(hiFilt,tmpdata2);

    tmpdata2 = tmpdata2 - mean(tmpdata2);
    if abs(min(tmpdata2))> abs(max(tmpdata2))
        tmpdata2 = -tmpdata2;
    end

    if UsingTTLflag == 1
        tmpdata2(1:((MicroStimStartEndTime(1,thisTrial)-3*waveformbeforepoint))) = mean(tmpdata2);
        tmpdata2((MicroStimStartEndTime(2,thisTrial)+3*waveformbeforepoint):end) = mean(tmpdata2);
    else
        threshold = 5* std(tmpdata2);
        idx1 = find(tmpdata2 > threshold,1,'first');
        tmpdata2(1:(idx1-numberofStimPulses*waveformbeforepoint)) = mean(tmpdata2);
        idx2 = find(tmpdata2 > threshold,1,'last');
        tmpdata2((idx2+numberofStimPulses*waveformbeforepoint):end) = mean(tmpdata2);
    end

    % [~,loc,~,p] = findpeaks(tmpdata2,'WidthReference','halfheight','MinPeakDistance',1/StimulationPulsesFrequency*SampleRate*0.95,'MinPeakHeight',max(tmpdata2)/5);
    
%     [~,loc,~,~] = findpeaks(tmpdata2,'WidthReference','halfheight','MinPeakDistance',1/StimulationPulsesFrequency*SampleRate*0.97,'MinPeakProminence',1.2*threshold);
%     [~,loc,~,~] = findpeaks(tmpdata2,'WidthReference','halfheight','MinPeakDistance',1/StimulationPulsesFrequency*SampleRate*0.95,'SortStr','descend','MinPeakProminence',1* threshold,'Annotate','extents');
    [~,loc,~,~] = findpeaks(tmpdata2,'WidthReference','halfheight','MinPeakDistance',1/StimulationPulsesFrequency*SampleRate*0.95,'SortStr','descend','Annotate','extents');

    loc = loc(1:min(numberofStimPulses,numel(loc)));
    eventidx = sort(loc,'ascend');

    
%     eventidx = loc;
    numTransEvent = numel(eventidx);
    % if number of TransEvent found still not equal to number of stim pulses delievered
    if numTransEvent ~= numberofStimPulses
        disp(['Failed to find all artifact on this trial, current level = ',num2str(DetectableCurrentLevel(thisTrial))])

        DetectableCurrentLevel(thisTrial) = -DetectableCurrentLevel(thisTrial);
        TrialAlignOffset(thisTrial) = 0;
        trialArtifactTimeIdx(:,thisTrial) = zeros(numberofStimPulses,1);
    else
        % find the waveform of the first artifact
        
        for thisEvent = 1:numTransEvent
            wave = tmpdata2((eventidx(thisEvent)-waveformbeforepoint) : (eventidx(thisEvent)+waveformafterpoint));
            x = 1:numel(wave);
            xx = 1:0.1:numel(wave);
            wave2 = spline(x,wave,xx);

            [~,idx1] = max(wave2);
%             wave2(1:idx1) = inf;
            [~,idx2] = min(wave2);
            idx1 = idx1(1);
            idx2 = idx2(1);
            idx = mean([idx1,idx2]);
            idx = round(idx/10);
            eventidx(thisEvent) = eventidx(thisEvent) +idx;
        end
        eventidx = eventidx - waveformbeforepoint;
        eventidx(1) = eventidx(2) - (eventidx(3)-eventidx(2));
        % save in instanceinfo
        if whichmethod == 1 && UsingTTLflag
            TrialAlignOffset(thisTrial) = MicroStimStartEndTime(1,thisTrial);
            trialArtifactTimeIdx(:,thisTrial) = eventidx - MicroStimStartEndTime(1,thisTrial);
        else
            TrialAlignOffset(thisTrial) = eventidx(1);
            trialArtifactTimeIdx(:,thisTrial) = eventidx - eventidx(1);
        end
        
    end
end

% hold off
% plot(tmpdata2)
% hold on
% plot(trialArtifactTimeIdx(:,thisTrial),tmpdata2(trialArtifactTimeIdx(:,thisTrial)),'x')
% pause(0.1)

% remove undetectable trials
idx = DetectableCurrentLevel < 0;
TrialAlignOffset(idx) = [];
trialArtifactTimeIdx(:,idx) = [];

thisInstance = 1;
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


% copy answer from the ref instance
refInstancIdx = 1;
for thisInstance = 1:numInstances
    % copy answer from the good
    TrialAlignOffset = instanceinfo(refInstancIdx).trialInfo.trialAlignOffset;
    trialArtifactTimeIdx = instanceinfo(refInstancIdx).trialInfo.trialArtifactIdx;
    %         validAllCurrentLevel = instanceinfo(instancIdx).trialInfo.validAllCurrentLevel;
    %         idx = validAllCurrentLevel > 0; % microstimulation trial index
    DetectableAllCurrentLevel = instanceinfo(refInstancIdx).trialInfo.DetectableAllCurrentLevel;
    idx = DetectableAllCurrentLevel > 0; % microstimulation trial index
    validCorrectedTrialAlignST = instanceinfo(thisInstance).trialInfo.validBeginST;
    validCorrectedTrialAlignST(idx) = validCorrectedTrialAlignST(idx) + TrialAlignOffset;
    instanceinfo(thisInstance).trialInfo.DetectableAllCurrentLevel = DetectableAllCurrentLevel;
    instanceinfo(thisInstance).trialInfo.validCorrectedTrialAlignST = validCorrectedTrialAlignST;
    instanceinfo(thisInstance).trialInfo.trialArtifactIdx = trialArtifactTimeIdx;
    instanceinfo(thisInstance).trialInfo.trialAlignOffset = TrialAlignOffset;
end

end