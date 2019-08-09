function instanceinfo = VnCMicroStimTrialAlignAndArtifactRemoval(instanceinfo,matdatapath,trialParam,StimulationParam)

if ispc
    slash = '\';
else
    slash = '/';
end

% stimulation frequency
StimulationFreq = 1/(StimulationParam.stimulationwaveformwidth/1000);
FreqRange = 200; %Hz

% internal algorithm parameters
threshold_lowlim = -32767;
threshold_upperlim = 32767;
threshold_stepsize = 1;
threshold = threshold_lowlim:threshold_stepsize:threshold_upperlim;
numberofthresholdsteps = numel(threshold);

load(matdatapath);

SampleRate = instanceinfo(1).samplerate;

initialtriallength = trialParam.initialtriallength; % sec
initialtrialAlignTime = trialParam.initialtrialAlignTime; % sec
endtriallength = trialParam.endtriallength; % sec
endtrialAlignTime = trialParam.endtrialAlignTime; % sec

numberofStimPulses=StimulationParam.numberofStimPulses;
stimulationwaveformwidthDP=StimulationParam.stimulationwaveformwidthDP; % 0.5ms


%Calculate trial info
trialstartTime = -initialtrialAlignTime;
trialendTime = trialstartTime + initialtriallength;
trialstartIdx = trialstartTime * SampleRate;
trialendIdx = trialendTime * SampleRate;
trialidx = trialstartIdx:(trialendIdx-1);
trialAlignIdx = initialtrialAlignTime * SampleRate;

%calculate stimulation timing profile
idx = 0:(stimulationwaveformwidthDP-1);
stimwaveidx = repmat(idx,numberofStimPulses,1);


% build filters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % design the filters
% % try to remove 50Hz from the signal
disp('Building filters...')
d50Hz = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    'DesignMethod','butter','SampleRate',SampleRate);

% band pass filter 300-6000Hz
bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',750,'HalfPowerFrequency2',5000, ...
    'DesignMethod','butter','SampleRate',SampleRate);

% low pass filter < 150 Hz
% lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',150,...
%     'DesignMethod','butter','SampleRate',SampleRate);

% MUAe low pass filter < 500 Hz
% MUAeFilt = designfilt('lowpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',500,...
%     'DesignMethod','butter','SampleRate',SampleRate);


% % High pass filter > 1000 Hz
% hpFilt = designfilt('highpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',1000,...
%     'DesignMethod','butter','SampleRate',SampleRate);


% split trials only keep good trials for each electrode
disp('Aligning microstimulation trials and removing stimulation artifacts...')
numInstances = numel(instanceinfo);
skipinstanceflag = zeros(numInstances,1);
for thisInstance = 1:numInstances
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances),'...'])
    numElec = instanceinfo(thisInstance).numElec;
    numNeuralElec = sum(instanceinfo(thisInstance).ElecOrder<=128);
    numTrials = instanceinfo(thisInstance).numTrials;
    if numTrials == 1
        % load electrode data back and look for the correct trials and put
        % the trials into a matrix
        numDP = instanceinfo(thisInstance).trialInfo(1).NumDP;
        validtrialIdx = instanceinfo(thisInstance).trialInfo.trialAlignST~=-1 & performance' ~= 0;
        numTrials = sum(validtrialIdx);
        trialtimeIdx = repmat(trialidx,numTrials,1);
        trialAlignST = instanceinfo(thisInstance).trialInfo.trialAlignST;
        idx = trialAlignST ~= -1;
        trialAlignST = trialAlignST(idx);
        trialLengthST = numel(trialidx);
        trialAlignST = repmat(trialAlignST,1,trialLengthST);
        trialtimeIdx = trialtimeIdx + trialAlignST;
        trialtimeIdx = trialtimeIdx';
        % precise alignment
        trialtimeIdx = trialtimeIdx - instanceinfo(thisInstance).trialInfo(1).Timestamp;
        clearvars trialAlignST;
        % find all microstimulation trials
        MicroStimTrialsIdx = find(allCurrentLevel(validtrialIdx)>0);
        numberMicrostimTrials = numel(MicroStimTrialsIdx);
        diffvalues = zeros(numberMicrostimTrials,numElec);
        disp('  Performing coarse alignment...')
        numTransEvent = zeros(numberofthresholdsteps,numberMicrostimTrials,numElec);
        trialStimEnergy = zeros(numberMicrostimTrials,numNeuralElec);
        rankStimEnergy = zeros(numberMicrostimTrials,numNeuralElec);
        for thisElec = 1:numElec
            EID = instanceinfo(thisInstance).ElecOrder(thisElec);
%             disp(EID)
            if EID > 128
                break
            end
            % load data
            elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{thisElec},'r');
            tmpdata = fread(elecfp,numDP,'int16=>double');
            fclose(elecfp);
            % do trial alignment
            tmpdata = tmpdata(trialtimeIdx);
            % zero mean
            meandata = mean(tmpdata,1);
            tmpdata = tmpdata - repmat(meandata,trialLengthST,1);
            % remove 50 Hz
            tmpdata = filtfilt(d50Hz,tmpdata);
            % bandpass filter, MUAs, MUAe
            bpdata = filtfilt(bpFilt,tmpdata);
            
            bpdataMicroStim = bpdata(:,MicroStimTrialsIdx);
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
            
            %%%%%%%%%%%%%%%%%%%%%% This algorithms fails when there is very
            %%%%%%%%%%%%%%%%%%%%%% large spikes
            %             maxvalues = max(bpdataMicroStim,[],1);
            %             minvalues = min(bpdataMicroStim,[],1);
            %             diffvalues(:,thisElec) = maxvalues-minvalues;
            
            %%%% we need to perform a search through all trials to find the
            %%%% one electrode that has the most prominent stimulation
            %%%% artifact
            
            % do a threshold search from -65535 to 65535 and calculate
            % the number of threshold crossing events for each
            % threshold and each trial
%             for thisThreshold = 1:numberofthresholdsteps
%                 currentThreshold = threshold(thisThreshold);
%                 if currentThreshold >= 0
%                     idx = bpdataMicroStim >= currentThreshold;
%                 else
%                     idx = bpdataMicroStim <= currentThreshold;
%                 end
%                 idx = idx - [idx(2:end,:);zeros(1,size(idx,2))];
%                 idx = idx == -1;
%                 eventidx = find(idx == -1);
%                 eventidx = eventidx + 1;
%                 numTransEvent(thisThreshold,:,thisElec) = sum(idx,1);
%             end
            y = fft(bpdataMicroStim);     
            f = (0:size(bpdataMicroStim,1)-1)*SampleRate/size(bpdataMicroStim,1);
            idx = f > StimulationFreq - FreqRange & f < StimulationFreq + FreqRange;
            y = abs(y);
            trialStimEnergy(:,thisElec) = max(y(idx,:),[],1);
%             subplot(1,2,1),plot(f,abs(y))
%             subplot(1,2,2),plot(bpdataMicroStim)
            
        end
        [~,idxTrialStimEnergy] = sort(trialStimEnergy,2,'descend');
        idx = 1:numNeuralElec;
        for thisTrial = 1:numberMicrostimTrials
            rankStimEnergy(thisTrial,idxTrialStimEnergy(thisTrial,:)) = idx;
        end
        
        trialMinStimEnergy = min(trialStimEnergy,[],1);
        
        
        disp('  Fine tuning alignment and removing stimulation artifact...')
        % find the electrode with the max minimal difference in all trials
%         elecidx = instanceinfo(thisInstance).ElecOrder >= 1 & instanceinfo(thisInstance).ElecOrder <=128; % only consider electrode on the frontal amplifier
%         [~,idx] = max(min(diffvalues(:,elecidx),[],1));
%         elecID = instanceinfo(thisInstance).ElecOrder(elecidx);
%         elecID = elecID(idx);
%         elecidx = find(instanceinfo(thisInstance).ElecOrder == elecID); % elecidx is the index of the chosen electrode
        %%%%%%%%%%%%%%
        % find the electrode with highest stimulation energy
%         [~,elecidx] = min(sum(rankStimEnergy,1));
        [~,elecidx] = max(trialMinStimEnergy);

        % recalculate the alignment point
        elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{elecidx},'r');
        tmpdata = fread(elecfp,numDP,'int16=>double');
        fclose(elecfp);
        % do trial alignment
        tmpdata = tmpdata(trialtimeIdx);
        % zero mean
        meandata = mean(tmpdata,1);
        tmpdata = tmpdata - repmat(meandata,trialLengthST,1);
        % remove 50 Hz
        tmpdata = filtfilt(d50Hz,tmpdata);
        % bandpass filter, MUAs, MUAe
        bpdata = filtfilt(bpFilt,tmpdata);
        bpdataMicroStim = bpdata(:,MicroStimTrialsIdx);
        
        % set a threshold and search for transgression event from the positive
        % end
        thresholdhi = max(bpdataMicroStim,[],1);
        thresholdlow = min(bpdataMicroStim,[],1);
        %     figure
        
        TrialAlignOffset = zeros(numberMicrostimTrials,1);
        %         figure,hold on
        trialArtifactTimeIdx = zeros(numberofStimPulses,numberMicrostimTrials);
        for thisTrial = 1:numberMicrostimTrials
            
            % find all the points that higher than the threshold
            numTransEvent = 0;
            counter = 0;
            while(numTransEvent < numberofStimPulses && counter < 65536)
                idx = bpdataMicroStim(:,thisTrial) >= thresholdhi(thisTrial);
                idx = idx - [idx(2:end);0];
                eventidx = find(idx == -1);
                eventidx = eventidx + 1;
                numTransEvent = numel(eventidx);
                thresholdhi(thisTrial) = thresholdhi(thisTrial)-threshold_stepsize;
                counter = counter +1;
            end
            threshold = thresholdhi(thisTrial);
            if numTransEvent ~= numberofStimPulses % try search from bottom to top instead
                numTransEvent = 0;
                counter = 0;
                while(numTransEvent < numberofStimPulses && counter < 65536)
                    idx = bpdataMicroStim(:,thisTrial) <= thresholdlow(thisTrial);
                    idx = idx - [idx(2:end);0];
                    eventidx = find(idx == -1);
                    eventidx = eventidx + 1;
                    numTransEvent = numel(eventidx);
                    thresholdlow(thisTrial) = thresholdlow(thisTrial)+threshold_stepsize;
                    counter = counter + 1;
                end
                threshold = thresholdlow(thisTrial);
            end
            % if cannot find correct number of stimulus we can use another
            % instance
            if numel(eventidx)~=numberofStimPulses
                skipinstanceflag(thisInstance) = 1;
                break
            end
            % find the waveform of the first artifact
            waveformbeforepoint = 35;
            waveformafterpoint = 35;
            wave = bpdataMicroStim((eventidx(1)-waveformbeforepoint) : (eventidx(1)+waveformafterpoint),thisTrial);
            
            [~,idx1] = max(wave);
            [~,idx2] = min(wave);
            idx1 = idx1(1);
            idx2 = idx2(1);
            idx = floor(mean([idx1,idx2]));
            %             plot(circshift(wave,35-idx))
            eventidx = eventidx +idx-waveformbeforepoint;
            TrialAlignOffset(thisTrial) = eventidx(1);
            % save in instanceinfo
            trialArtifactTimeIdx(:,thisTrial) = eventidx - eventidx(1);
            
            
            % inspection
%                     hold off
%                     plot(bpdataMicroStim(:,thisTrial))
%                     hold on
%                     plot(eventidx,repmat(threshold,numel(eventidx),1),'ro')
        end
        
        if skipinstanceflag(thisInstance) == 1
            continue;
        end
        
        % plot the data of all trials on that electrode to inspect the result
        %     [~,currentlevelsIdx] = sort(allCurrentLevel(performance~=0),'descend');
        %     offset = sum(allCurrentLevel == 0);
        %     numMicrostimulationTrials = sum(allCurrentLevel ~= 0 & performance ~=0);
        %     figure
        %     subplot(4,1,1), hold on
        %     colors = lines(numMicrostimulationTrials);
        %     for thisTrial = 1:numMicrostimulationTrials
        %         plot(bpdataMicroStim(:,thisTrial),'color',colors(thisTrial,:))
        %     end
        
        % plot the realigned data
        %     subplot(4,1,2), hold on
        %     colors = lines(numMicrostimulationTrials);
        %     for thisTrial = 1:numMicrostimulationTrials
        %         ydata = circshift(bpdataMicroStim(:,thisTrial),2500-TrialAlignOffset(thisTrial));
        % %         ydata = circshift(ydata,-5000);
        %         plot(ydata,'color',colors(thisTrial,:))
        %     end
        %
        %     subplot(4,1,3), hold on
        %     colors = lines(numMicrostimulationTrials);
        %     for thisTrial = 1:numMicrostimulationTrials
        %         plot(bpdataMicroStim(:,thisTrial),'color',colors(thisTrial,:))
        %     end
        %
        %     subplot(4,1,4), hold on
        %     colors = lines(numMicrostimulationTrials);
        %     for thisTrial = 1:numMicrostimulationTrials
        %         ydata = circshift(bpdataMicroStim(:,thisTrial),2500-TrialAlignOffset(thisTrial));
        % %         ydata = circshift(ydata,-5000);
        %         plot(ydata,'color',colors(thisTrial,:))
        %     end
        
        
        % update the trial alignment point for micro-stimulation trials.
        idx = find(validtrialIdx);
        idx = idx(MicroStimTrialsIdx);
        instanceinfo(thisInstance).trialInfo.trialAlignST(idx) = instanceinfo(thisInstance).trialInfo.trialAlignST(idx) - (trialAlignIdx - TrialAlignOffset);
        instanceinfo(thisInstance).trialInfo.trialArtifactIdx = trialArtifactTimeIdx;
        instanceinfo(thisInstance).trialInfo.trialAlignOffset = TrialAlignOffset;
    else
    end
end

if sum(skipinstanceflag) == numInstances
    error('not enough alignment information')
end

% if on some instance(s) we cannot find correct stimulus artifact, we can
% copy answer from a known good instance
if sum(skipinstanceflag) > 0
    badinstanceidx = find(skipinstanceflag);
    numelSkippedInstances = numel(badinstanceidx);
    goodidx = find(skipinstanceflag==0);
    goodidx = goodidx(1);
    for thisInstance = 1:numelSkippedInstances
        instanceIdx = badinstanceidx(thisInstance);
        idx = find(validtrialIdx);
        idx = idx(MicroStimTrialsIdx);
        TrialAlignOffset = instanceinfo(goodidx).trialInfo.trialAlignOffset;
        instanceinfo(instanceIdx).trialInfo.trialAlignST(idx) = instanceinfo(thisInstance).trialInfo.trialAlignST(idx) - (trialAlignIdx - TrialAlignOffset);
        instanceinfo(instanceIdx).trialInfo.trialArtifactIdx = instanceinfo(goodidx).trialInfo.trialArtifactIdx;
        instanceinfo(instanceIdx).trialInfo.trialAlignOffset = trialAlignOffset;
    end
end
        
for thisInstance = 1:numInstances
        % recalculate the trial data based on the new alignment point
        numTrials = instanceinfo(thisInstance).numTrials;
        trialstartTime = -endtrialAlignTime;
        trialendTime = trialstartTime + endtriallength;
        trialstartIdx = trialstartTime * SampleRate;
        trialendIdx = trialendTime * SampleRate;
        trialidx = trialstartIdx:(trialendIdx-1);
        trialAlignIdx = endtrialAlignTime * SampleRate;
        
        trialtimeIdx = repmat(trialidx,numTrials,1);
        trialAlignST = instanceinfo(thisInstance).trialInfo.trialAlignST;
        idx = trialAlignST ~= -1;
        trialAlignST = trialAlignST(idx);
        trialLengthST = numel(trialidx);
        trialAlignST = repmat(trialAlignST,1,trialLengthST);
        trialtimeIdx = trialtimeIdx + trialAlignST;
        trialtimeIdx = trialtimeIdx';
        % precise alignment
        trialtimeIdx = trialtimeIdx - instanceinfo(thisInstance).trialInfo(1).Timestamp;
        clearvars trialAlignST;
        % find all microstimulation trials
        MicroStimTrialsIdx = find(allCurrentLevel(instanceinfo(thisInstance).trialInfo.trialAlignST~=-1)>0);
        numberMicrostimTrials = numel(MicroStimTrialsIdx);
        diffvalues = zeros(numberMicrostimTrials,numElec);
        for thisElec = 1:numElec
            % load data
            elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{thisElec},'r');
            tmpdata = fread(elecfp,numDP,'int16=>double');
            fclose(elecfp);
            % do trial alignment
            tmpdata = tmpdata(trialtimeIdx);
            % zero mean
            meandata = mean(tmpdata,1);
            tmpdata = tmpdata - repmat(meandata,trialLengthST,1);
            
            trialArtifactIdx = instanceinfo(thisInstance).trialInfo.trialArtifactIdx;
            trialArtifactIdx = round(mean(trialArtifactIdx,2));
            trialdata = tmpdata;
            numtrials = size(tmpdata,2);
            for thisTrial = 1:numtrials
%                 trialIdx = MicroStimTrialsIdx(thisTrial);
                tmptrialdata = tmpdata(:,thisTrial);
                triallength = numel(tmptrialdata);
                xindex = 1:triallength;
                % remove stimulus artifact from the xdata
                xindex2 = xindex;
                stimstart = repmat(trialArtifactIdx+trialAlignIdx - floor(stimulationwaveformwidthDP/2)+1,1,stimulationwaveformwidthDP);
                xindex3 = stimwaveidx + stimstart;
                xindex2(xindex3) = [];
                %                 trialdata2 = spline(xindex2,trialdata(xindex2),xindex);
                trialdata(:,thisTrial) = interp1(xindex2,tmptrialdata(xindex2),xindex,'linear');
                %                 figure, hold on
                %                 plot(trialdata,'black')
                %                 plot(trialdata2,'red')
                %                 plot(trialdata3,'blue')
                %                 plot(xindex3,zeros(size(xindex3)),'.')
                %                 plot(trialArtifactIdx(:,thisTrial)+trialAlignIdx,zeros(size(trialArtifactIdx(:,thisTrial))),'rx')
                %                 plot(xindex4,zeros(size(xindex4)),'ro')
            end
            % store the trial data on disk
            electroderawTrialData.samplerate = SampleRate;
            electroderawTrialData.elecID = instanceinfo(thisInstance).ElecOrder(thisElec);
            electroderawTrialData.elecIdx = thisElec;
            electroderawTrialData.trialdata = trialdata;
            electroderawTrialData.trialAlignTime = endtrialAlignTime;
            
            electroderawTrialData.trialIdx = find(instanceinfo(thisInstance).trialInfo.trialAlignST~=-1);
            electroderawTrialData.trialAlignST = instanceinfo(thisInstance).trialInfo.trialAlignST(electroderawTrialData.trialIdx);
            idx = strfind(instanceinfo(thisInstance).electrodeCachePath{thisElec},[slash,'tmp',slash]);
            electroderawTrialDataSaveDir = instanceinfo(thisInstance).electrodeCachePath{thisElec}(1:(idx-1));
            if ~exist(electroderawTrialDataSaveDir,'dir')
                mkdir(electroderawTrialDataSaveDir)
            end
            electroderawTrialDataSavePath = [instanceinfo(thisInstance).electrodeCachePath{thisElec}(1:idx),'electroderawTrialData.mat'];
            instanceinfo(thisInstance).electroderawTrialDataPath{thisElec} = electroderawTrialDataSavePath;
            save(electroderawTrialDataSavePath,'electroderawTrialData')
        end
end