function instanceinfo = VnCArtifactRemoval(instanceinfo,stimParam,MUAparameters,trialParam,whichmethod,loadSpikenMUAsFlag,interpolationLength,DebugFlag)

% DebugFlag = 1;

if DebugFlag
    h=figure('position',[100 100 1024 800]);
end

if ispc
    slash = '\';
else
    slash = '/';
end

% switch whichmethod
%     case 1
%         interpolationLength = 0.0012;% s
%     case 2
%         interpolationLength = 0.0012;% s
%     case 3
%         interpolationLength = 0.0015; % s
% end

spikeTrainEdge = (0:0.001:(trialParam.endT - trialParam.startT)) + trialParam.startT;
psth10Edge = (0:0.01:(trialParam.endT - trialParam.startT)) + trialParam.startT;
waveformLength = 64;
waveformAlignpoint = 0.3;
xEstimation = -3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% super sampling
supersampling = 10; % times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters ends here

MUAeSamplingrate=MUAparameters.MUAeSamplingrate;





numInstances = numel(instanceinfo);
% figure;
disp('Removing artifact...')
for thisInstance = 1:numInstances
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances)])
    numElec = instanceinfo(thisInstance).numElec;
    
    SampleRate = double(instanceinfo(thisInstance).samplerate);
    
    % hipass > 5Hz
    %     hipass = designfilt('highpassiir','FilterOrder',2, ...
    %         'HalfPowerFrequency',5,...
    %         'DesignMethod','butter','SampleRate',SampleRate);
    
    %     d50Hz = designfilt('bandstopiir','FilterOrder',2, ...
    %         'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    %         'DesignMethod','butter','SampleRate',SampleRate*supersampling);
    
    % band pass filter 500-5000Hz
    %     bpFilt = designfilt('bandpassiir','FilterOrder',2, ...
    %         'HalfPowerFrequency1',500,'HalfPowerFrequency2',5000, ...
    %         'DesignMethod','butter','SampleRate',SampleRate);
    
    Fn = SampleRate/2;
    Fbp=[500,5000];
    N  = 2;    % filter order
    [Bbp, Abp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]); % compute filter coefficients
    
    %     % high pass 500Hz
    %     hipass500 = designfilt('highpassiir','FilterOrder',2, ...
    %         'HalfPowerFrequency',700,...
    %         'DesignMethod','butter','SampleRate',SampleRate*supersampling);
    
    % low pass filter < 200 Hz
    %     MUAeFilt = designfilt('lowpassiir','FilterOrder',2, ...
    %         'HalfPowerFrequency',200,...
    %         'DesignMethod','butter','SampleRate',SampleRate*supersampling);
    Fl=200;
    Fn = SampleRate/2;
    N = 2;
    [Blp, Alp] = butter(N,Fl/Fn,'low'); % compute filter coefficients
    
    
    %     AlignIdx = -trialParam.startT * SampleRate;
    ArtifactIdx = instanceinfo(thisInstance).trialInfo.trialArtifactIdx;
    
    switch whichmethod
        case 1
            stimWidthMargin = 32;
            stimWidth = ceil(stimParam.stimulationwaveformwidth/1000 * SampleRate)+stimWidthMargin;
        case 2
            stimWidthMargin = 16;
            stimWidth = ceil(stimParam.stimulationwaveformwidth/1000 * SampleRate)+stimWidthMargin;
        case 3
            stimWidthMargin = 32;
            stimWidth = ceil(stimParam.stimulationwaveformwidth/1000 * SampleRate)+stimWidthMargin;
    end
    stimIdx = 0:stimWidth-1;
    stimIdx = repmat(stimIdx,stimParam.numberofStimPulses,1);
    interpolength = round(interpolationLength*SampleRate);
    stimIdx2 = 0:interpolength-1;
    stimIdx2 = repmat(stimIdx2,stimParam.numberofStimPulses,1);
    
    microstimtrialidx = find(instanceinfo(thisInstance).trialInfo.DetectableAllCurrentLevel > 0);
%     microstimtrialidx = instanceinfo(thisInstance).trialInfo.DetectableAllCurrentLevel > 0;
    numMicroStimTrials = numel(microstimtrialidx);
    MUAeDownSampleRatio = round(SampleRate/MUAeSamplingrate);
    for thisElec = 1:numElec
        %     for thisElec = 1:1
        %         ElecStimulationArtifacts = [];
        %         EID = instanceinfo(thisInstance).ElecOrder(thisElec);
        %         if EID ==128 + 15
        %             a = 1;
        %         end
        %         instanceinfo(thisInstance).electrodeTrialDataPath{thisElec} = strrep(instanceinfo(thisInstance).electrodeTrialDataPath{thisElec},'D:\','..\');
        load(instanceinfo(thisInstance).electrodeTrialDataPath{thisElec});
        trialParam = electrodeTrialData.trialParam;
        AlignIdx = round(-trialParam.startT * SampleRate);
        trialAlignIdx = - electrodeTrialData.trialParam.startT * SampleRate;
        ArtifactIdx2 = ArtifactIdx + trialAlignIdx;
        ArtifactIdx2 = round(mean(ArtifactIdx2,2));
        ArtifactIdx = round(mean(ArtifactIdx,2));
        trialData =  electrodeTrialData.trialData;
        xlimit = [-200,500]/1000;
        if DebugFlag
            subplot(2,3,1) % raw trace
            hold off
            xdata = (1:size(trialData,2)) - trialAlignIdx;
            xdata = xdata/SampleRate;
            plot(xdata,trialData(2:11,:)','blue')
            hold on
            xlim(xlimit)
            %             xlim([trialAlignIdx-300,trialAlignIdx+1000])
            subplot(2,3,2) % after artifact substracting
            hold off
            plot(xdata,trialData(2:11,:)','blue')
            hold on
            xlim(xlimit)
            %             xlim([trialAlignIdx-300,trialAlignIdx+1000])
            
        end
        
        
        if whichmethod == 3 % subtract trial average from the data
            
            StimTrialData = trialData(microstimtrialidx,:);
            % zeromean on each trial
            tmpData = StimTrialData;
            % calculate mean stim trial data for each stimulus conditions
            Condition = instanceinfo(1).trialInfo.Condition;
            UniCond = unique(Condition);
            NumUniCond = numel(UniCond);
            for thisCond = 1:NumUniCond
                trialIdx = Condition(microstimtrialidx) == UniCond(thisCond);
                meanStimTrialData = mean(tmpData(trialIdx,:),1);
                StimTrialData(trialIdx,:) = StimTrialData(trialIdx,:) - repmat(meanStimTrialData,sum(trialIdx),1);
            end
            trialData(microstimtrialidx,:) = StimTrialData;
        else
            
            %load(instanceinfo(thisInstance).electrodeNeuraldataPath{thisElec});
            % load raw trial data
            % perform same operation for all trials?
            % Method 1 is Chris's method
            % find the average shape of the artifact
            % subtract the artifact from the raw trace
            % band filter the trace
            % apply linear interpolation
            % low pass filter to get MUAe
            % down sample
            
            trialArtifactIdx = repmat(ArtifactIdx2,1,stimWidth);
            trialArtifactIdx = stimIdx + trialArtifactIdx - round(stimWidth/2) + 1;
            %             trialArtifactIdx = reshape(trialArtifactIdx',[],1);
            %             StimTrialData = trialData(microstimtrialidx,trialArtifactIdx);
            %             meanArtifacts = mean(StimTrialData,1);
            %             meanArtifacts = reshape(meanArtifacts,stimWidth,stimParam.numberofStimPulses);
            %             for thisStimPulse = 1:size(meanArtifacts,2)
            %                 meanArtifacts(:,thisStimPulse) = meanArtifacts(:,thisStimPulse) - mean([meanArtifacts(1,thisStimPulse),meanArtifacts(end,thisStimPulse)]);
            %             end
            %             meanArtifacts = reshape(meanArtifacts,[],1);
            for thisTrial = 1:numMicroStimTrials
                trialIdx  = microstimtrialidx(thisTrial);
                
                thisTrialData = trialData(trialIdx,:);
                StimulationArtifacts = thisTrialData(trialArtifactIdx);
                numArtifacts = size(StimulationArtifacts,1);
                switch whichmethod
                    case 1
                        %                         calculate the mean artifact shape
                        
                        for thisArtifact = 1:numArtifacts
                            StimulationArtifacts(thisArtifact,:) = StimulationArtifacts(thisArtifact,:) - mean([StimulationArtifacts(thisArtifact,1),StimulationArtifacts(thisArtifact,end)]);
                        end
                        
                        
                        %                         ElecStimulationArtifacts{thisTrial} = StimulationArtifacts;
                        
                        meanStimulationArtifacts = mean(StimulationArtifacts,1);
                        meanStimulationArtifacts = meanStimulationArtifacts - mean([meanStimulationArtifacts(1),meanStimulationArtifacts(end)]);
                        % sustract the mean artifact shape from the trace
                        meanStimulationArtifacts = repmat(meanStimulationArtifacts,stimParam.numberofStimPulses,1);
                        thisTrialData2 = thisTrialData;
                        thisTrialData2(trialArtifactIdx) = thisTrialData2(trialArtifactIdx) - meanStimulationArtifacts;% + baselineLevel;
                        trialData(trialIdx,:) = thisTrialData2;
                        if DebugFlag
                            subplot(2,3,2)
                            plot(((trialArtifactIdx- trialAlignIdx)/SampleRate)',meanStimulationArtifacts','green')
                        end
                        %                           trialData(trialIdx,trialArtifactIdx) = trialData(trialIdx,trialArtifactIdx) - meanArtifacts';
                    case 2
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % find start and stop ponits
                        xx = 1:numel(thisTrialData);
                        x = xx;
                        thisTrialDataTmp = thisTrialData;
                        for thisArtifact = 1:numArtifacts
                            thisArtifactWave = StimulationArtifacts(thisArtifact,:);
                            
                            [~,idx1] = max(thisArtifactWave);
                            [~,idx2] = min(thisArtifactWave);
                            if idx1 < idx2
                                thisArtifactWave = -thisArtifactWave;
                            end
                            
                            [~,idx1] = max(thisArtifactWave);
                            tmp = thisArtifactWave(1:idx1);
                            tmp = abs(tmp-tmp(1));
                            idx1 = find(tmp >= max(tmp)/10,1,'first');
                            if idx1 > 2
                                idx1 = idx1 - 2;
                            end
                            
                            [~,idx2] = min(thisArtifactWave);
                            tmp = thisArtifactWave(idx2:end);
                            
                            tmp = abs(tmp-tmp(end));
                            idx = find(tmp >= max(tmp)/10,1,'last');
                            idx2 = idx + idx2;
                            if length(thisArtifactWave) - idx2 > 2
                                idx2 = idx2 + 2;
                            end
                            
                            if idx1 < 1
                                idx1 = 1;
                            end
                            if idx2 > length(thisArtifactWave)
                                idx2 = length(thisArtifactWave);
                            end
                            
                            idx = trialArtifactIdx(thisArtifact,:);
                            idx = idx(idx1:idx2);
                            x(idx) = nan;
                        end
                        idx = isnan(x);
                        x(idx) = [];
                        thisTrialDataTmp = thisTrialDataTmp(x);
                        thisTrialData = interp1(x,thisTrialDataTmp,xx,'linear');
                        trialData(trialIdx,:) = thisTrialData;
                end
            end
        end
        
        %         trialData3 = trialData;
        trialData = trialData';
        
        if DebugFlag
            subplot(2,3,2) % after artifact substracting
            plot(xdata,trialData(:,2:11),'red')
            subplot(2,3,3) % after high pass filtering
            hold off
            plot(xdata,trialData(:,2:11),'blue')
            hold on
            xlim(xlimit)
            %             xlim([trialAlignIdx-300,trialAlignIdx+1000])
            
        end
        
        % 50Hz removal
        %         x = 1:size(trialData,1);
        %         xx = 1:(1/supersampling):size(trialData,1);
        %         trialData = trialData';
        %         trialData = spline(x,trialData,xx);
        %         trialData = trialData';
        %         trialData = filtfilt(d50Hz,trialData);
        
        
        % band pass 500-5000Hz filter the trace
        %         trialData = filtfilt(bpFilt,trialData);
        trialData = filtfilt(Bbp, Abp, trialData);
        
        %         trialData = downsample(trialData,supersampling);
        
        if DebugFlag
            subplot(2,3,3) % after high pass filtering
            plot(xdata,trialData(:,2:11),'red')
            subplot(2,3,4) % after full wave rectification
            hold off
            plot(xdata,trialData(:,2:11),'blue')
            hold on
            xlim(xlimit)
            %             xlim([trialAlignIdx-300,trialAlignIdx+1000])
        end
        
        
        trialData2 = trialData;
        
        if loadSpikenMUAsFlag
            % Estimate the threshold
            tmpdata = reshape(trialData,[],1);
            elecThreshold = xEstimation * median(abs(tmpdata)/0.6745,1); % 0.6745 is from experience discribed in iterature
            
            % spike detection
            TrialLengthDP = size(trialData,1);
            C_poz = trialData - repmat(elecThreshold,TrialLengthDP,1);
            
            if xEstimation < 0
                C_poz = double(C_poz < 0);
            else
                C_poz = double(C_poz > 0);
            end
            
            C_poz_shifted = [zeros(1,size(C_poz,2)); C_poz(1:end-1,:)];
            C_poz_diff = C_poz - C_poz_shifted;
            clearvars C_poz C_poz_shifted
            
            numtrials = size(trialData,2);
            SpikeSTAligned = cell(numtrials,1);
            SpikeTimeAligned = cell(numtrials,1);
            RawElecSpikes = cell(numtrials,1);
            NumTrialSpikes = zeros(numtrials,1);
            SpikeTrain = zeros(numel(spikeTrainEdge)-1,numtrials);
            Psth10 = zeros(numel(psth10Edge)-1,numtrials);
            for thisTrial = 1: numtrials
                % get spikes
                TrialSpikes = find(C_poz_diff(:,thisTrial) == 1);
                % exclude any spikes within one waveform length from the beginning
                % and one waveform length at the end
                idx = TrialSpikes < 2*waveformLength | size(C_poz_diff,1) - TrialSpikes < 2*waveformLength;
                TrialSpikes(idx) = [];
                RawElecSpikes{thisTrial} = TrialSpikes;
                SpikeSTAligned{thisTrial} = TrialSpikes - AlignIdx;
                SpikeTimeAligned{thisTrial} = SpikeSTAligned{thisTrial} / SampleRate;
                NumTrialSpikes(thisTrial) = numel(TrialSpikes);
                tmp = histc(SpikeTimeAligned{thisTrial},spikeTrainEdge);
                tmp(end) = [];
                SpikeTrain(:,thisTrial) = tmp;
                tmp = histc(SpikeTimeAligned{thisTrial},psth10Edge);
                tmp(end) = [];
                Psth10(:,thisTrial) = tmp;
            end
            
            % get waveform
            elecRawWaveform = zeros(waveformLength,sum(NumTrialSpikes));
            elecRawSpikeStamps = zeros(sum(NumTrialSpikes),1);
            waveidx(1) = 1;
            waveidx(2) = 0;
            for thisTrial = 1: numtrials
                if ~isempty(RawElecSpikes{thisTrial})
                    C_wavestart = RawElecSpikes{thisTrial} - round(waveformLength * waveformAlignpoint); % waveformLength is the length of padded data
                    C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
                    %         align to the minimal or maximal
                    temp = reshape(trialData(C_waveIndex,thisTrial),[],waveformLength);
                    if xEstimation < 0
                        [~,temp] = min(temp,[],2);
                    else
                        [~,temp] = max(temp,[],2);
                    end
                    C_wavestart = C_wavestart + temp - round(waveformLength * waveformAlignpoint);
                    C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
                    ElecTrialWaveform = reshape(trialData(C_waveIndex,thisTrial),[],waveformLength);
                    ElecTrialWaveform = ElecTrialWaveform';
                    waveidx(2) = waveidx(2) + NumTrialSpikes(thisTrial);
                    elecRawWaveform(:,waveidx(1):waveidx(2)) = ElecTrialWaveform;
                    elecRawSpikeStamps(waveidx(1):waveidx(2)) = SpikeSTAligned{thisTrial} + electrodeTrialData.AlignStamps(thisTrial);
                    waveidx(1) = waveidx(2) + 1;
                end
            end
            if size(elecRawWaveform,1) ~= waveformLength
                elecRawWaveform = elecRawWaveform';
            end
        end
        
        % full wave rectification
        trialData2 = abs(trialData2);
        
        if DebugFlag
            subplot(2,3,4) % after full wave rectification
            plot(xdata,trialData2(:,2:11),'red')
            subplot(2,3,5) % after interpolation
            hold off
            plot(xdata,trialData2(:,2:11),'blue')
            hold on
            xlim(xlimit)
            %             xlim([trialAlignIdx-300,trialAlignIdx+1000])
        end
        
        % apply linear interpolation
        
        trialLength = size(trialData2,1);
        trialXXindex = 1:trialLength;
        
        for thisTrial = 1:numMicroStimTrials
            
            trialXindex = trialXXindex;
            trialIdx  = microstimtrialidx(thisTrial);
            % remove samples to be interpolated
            trialArtifactIdx3 = repmat(ArtifactIdx2,1,interpolength);
            trialArtifactIdx3 = stimIdx2 + trialArtifactIdx3 - round(interpolength/2) +1;
            trialXindex(trialArtifactIdx3)=[];
            thisTrialData = trialData2(trialXindex,trialIdx);
            
            
            
            % perform linear interpolation
            trialData2(:,trialIdx) = interp1(trialXindex,thisTrialData,trialXXindex,'linear');
            
        end
        
        if DebugFlag
            subplot(2,3,5) % after interpolation
            hold on
            %             xlim([trialAlignIdx-300,trialAlignIdx+1000])
            plot(xdata,trialData2(:,2:11),'red')
        end
        
        
        % low pass filter to get MUAe
        trialData2 = filtfilt(Blp, Alp, trialData2);
        
        % down sample
        MUAeArtifactRemoved = downsample(trialData2,MUAeDownSampleRatio);
        
        if DebugFlag
            subplot(2,3,6) % MUAe
            hold off
            ydata = trialData2(:,2:11);
            %             ydata = smooth(ydata,200);
            xdata = 1:size(ydata,1);
            xdata = xdata * 1000/SampleRate + trialParam.startT * 1000; % ms
            plot(xdata,ydata)
            hold on
            
            ydata = MUAeArtifactRemoved;
            ydata = mean(ydata,2);
            ydata = smooth(ydata,20);
            
            xdata = 1:numel(ydata);
            xdata = xdata * 1000/MUAeSamplingrate + trialParam.startT * 1000; % ms
            plot(xdata,ydata)
            plot([0 0],[min(ydata),max(ydata)],'black','linewidth',1.5)
            plot([167 167],[min(ydata),max(ydata)],'black','linewidth',1.5)
            xlim([-600 1000])
            %             a = input('Continue?');
            pause(0.001)
        end
        
        
        %save the result
        %         instanceinfo(thisInstance).electrodeNeuraldataPath{thisElec} = strrep(instanceinfo(thisInstance).electrodeNeuraldataPath{thisElec},'D:\','..\');
        if isfield(instanceinfo(thisInstance),'electrodeNeuraldataPath') && thisElec <= numel(instanceinfo(thisInstance).electrodeNeuraldataPath)
            electrodeNeuraldataPath = instanceinfo(thisInstance).electrodeNeuraldataPath{thisElec};
            load(electrodeNeuraldataPath)
        else
            tmppath = instanceinfo(thisInstance).electrodeTrialDataPath{thisElec};
            idx = find(tmppath == slash,1,'last');
            electrodeNeuraldataPath = [tmppath(1:idx),'electrodeNeuraldata.mat'];
            instanceinfo(thisInstance).electrodeNeuraldataPath{thisElec} = electrodeNeuraldataPath;
        end
        %         electrodeTrialdataPath = instanceinfo(thisInstance).electrodeTrialDataPath{thisElec};
        %         load(electrodeTrialdataPath)
        switch whichmethod
            case 1 % Feng
                electrodeNeuraldata.MUAeArtifactRemovedFengMethod = MUAeArtifactRemoved;
                %                 electrodeNeuraldata.RawArtifactsFengMethod = ElecStimulationArtifacts;
                if loadSpikenMUAsFlag
                    electrodeNeuraldata.SpikeTimeArtifactRemovedFengMethod = SpikeTimeAligned;
                    electrodeNeuraldata.PSTH10ArtifactRemovedFengMethod = Psth10;
                    electrodeNeuraldata.elecRawWaveformArtifactRemovedFengMethod = elecRawWaveform;
                end
                
                %                 electrodeTrialData.TrialDataArtifactRemovedChrisMethod = trialData3;
                %                 save(electrodeTrialdataPath,'electrodeTrialData');
            case 2 % Heffer
                %                 electrodeNeuraldata.RawArtifactsHefferMethod = ElecStimulationArtifacts;
                electrodeNeuraldata.MUAeArtifactRemovedHefferMethod = MUAeArtifactRemoved;
                if loadSpikenMUAsFlag
                    electrodeNeuraldata.SpikeTimeArtifactRemovedHefferMethod = SpikeTimeAligned;
                    electrodeNeuraldata.PSTH10ArtifactRemovedHefferMethod = Psth10;
                    electrodeNeuraldata.elecRawWaveformArtifactRemovedHefferMethod = elecRawWaveform;
                end
                %                 electrodeTrialData.TrialDataArtifactRemovedHefferMethod = trialData3;
                %                 save(electrodeTrialdataPath,'electrodeTrialData');
            case 3 % Bram
                %                 electrodeNeuraldata.RawArtifactsBramMethod = ElecStimulationArtifacts;
                electrodeNeuraldata.MUAeArtifactRemovedBramMethod = MUAeArtifactRemoved;
                if loadSpikenMUAsFlag
                    electrodeNeuraldata.SpikeTimeArtifactRemovedBramMethod = SpikeTimeAligned;
                    electrodeNeuraldata.PSTH10ArtifactRemovedBramMethod = Psth10;
                    electrodeNeuraldata.elecRawWaveformArtifactRemovedBramMethod = elecRawWaveform;
                end
                %                 electrodeTrialData.TrialDataArtifactRemovedBramMethod = trialData3;
                %                 save(electrodeTrialdataPath,'electrodeTrialData');
        end
        save(electrodeNeuraldataPath,'electrodeNeuraldata');
    end
end
if DebugFlag
    close(h)
end