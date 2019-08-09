function instanceinfo = VnCgetNeuraldata2(instanceinfo,trialAlignSTfieldname,MUAparameters,LFPparameters, SpikenWaveParameters,trialParam,savefilename,loadSpikenMUAsFlag)

if ispc
    slash = '\';
else
    slash = '/';
end

%Calculate trial info

LFPsamplingrate=LFPparameters.LFPsamplingrate;
MUAeSamplingrate=MUAparameters.MUAeSamplingrate;
MUAshapleySampleRate=MUAparameters.MUAshapleySampleRate;

% spike and waveform parameters
xEstimation=SpikenWaveParameters.xEstimation;
waveformLength=SpikenWaveParameters.waveformLength;
waveformAlignpoint=SpikenWaveParameters.waveformAlignpoint;

spikeTrainEdge = (0:0.001:(trialParam.endT - trialParam.startT)) + trialParam.startT;
psth10Edge = (0:0.01:(trialParam.endT - trialParam.startT)) + trialParam.startT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters ends here

% calculate LFP and MUAs and MUAe



numInstances = numel(instanceinfo);
% figure;
disp('Calculating LFP, MUA, and spikes...')
for thisInstance = 1:numInstances
    fprintf('\n')
    disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances)])
    numElec = instanceinfo(thisInstance).numElec;
    
    SampleRate = double(instanceinfo(thisInstance).samplerate);
    AlignIdx = -trialParam.startT * SampleRate;
    
    % design filters
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % design the filters
    % % try to remove 50Hz from the signal


    Fn = SampleRate/2;
    Fbp=[500,5000];
    N  = 2;    % filter order
    [Bbp, Abp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]); % BandPass
    
    Fn = SampleRate/2;
    Fbs=[49,51];
    N  = 2;    % filter order
    [Bbs, Abs] = butter(N, [min(Fbs)/Fn max(Fbs)/Fn],'stop'); % Bandstop
    

    Fl=200;
    Fn = SampleRate/2;
    N = 2;
    [BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass
    
    Fl=150;
    Fn = SampleRate/2;
    N = 2;
    [BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass
    
    % calculate down sample rate
    LFPdownSampleRatio = round(SampleRate/LFPsamplingrate);
    MUAeDownSampleRatio = round(SampleRate/MUAeSamplingrate);
    
    % get trialAlignST
    cmd = ['trialAlignST = instanceinfo(thisInstance).trialInfo.',trialAlignSTfieldname,';'];
    eval(cmd);
    
    for thisElec = 1:numElec
        fprintf('.')
        EID = instanceinfo(thisInstance).ElecOrder(thisElec);
                    if EID > 128
                        break
                    end
        
        % load electrode trial data from disk cache
        electrodeTrialData = VnCGetTrialData2(instanceinfo, thisInstance, EID,trialAlignST,trialParam);

        tmpdata = electrodeTrialData.trialData';
        AlignStamps = electrodeTrialData.AlignStamps;
        
        clearvars electrodeTrialData
        
        numtrials = size(tmpdata,2);
        %         plot(mean(tmpdata,2),'black')
        % remove 50Hz
        tmpdata = filtfilt(Bbs, Abs,tmpdata);
        %         hold on, plot(mean(tmpdata,2),'red'),hold off
        % remove 100Hz
%         tmpdata = filtfilt(d100Hz,tmpdata);
        %         hold on, plot(mean(tmpdata,2),'blue'),hold off
        % low pass filter and downsample to get LFP
        LFP = filtfilt(BLFP,ALFP,tmpdata);
        LFP = downsample(LFP,LFPdownSampleRatio);
        
        if loadSpikenMUAsFlag
            % MUAshapley
            rawElecMUA = filtfilt(hpFilt,tmpdata);
            %     half wave rectify in the negative direction
            rawElecMUA = rawElecMUA.*(rawElecMUA<0);
            meanElecMUA = mean(rawElecMUA,1);
            stdElecMUA = std(rawElecMUA,0,1);
            % count the number of events more than 3 times std with in 1 ms time
            % bins
            MUAthres = meanElecMUA - 3 * stdElecMUA;
            rawElecMUA = double(rawElecMUA < repmat(MUAthres,size(rawElecMUA,1),1));
            MUAlength = size(rawElecMUA,1);
            NumPointsperBin = round(SampleRate/MUAshapleySampleRate);
            NumBins = floor(MUAlength / NumPointsperBin);
            rawElecMUA = rawElecMUA(1:(NumPointsperBin*NumBins),:);
            MUA = zeros(NumBins,size(tmpdata,2));
            for thispoint = 1:NumPointsperBin
                MUA = MUA+ rawElecMUA(thispoint:NumPointsperBin:(end-(NumPointsperBin-thispoint)),:);
            end
        end
        
        % band pass filter
        bpdata = filtfilt(Bbp,Abp,tmpdata);
        clearvars tmpdata
        
        
        
        % MUAe
        MUAe = abs(bpdata);
        
        if ~loadSpikenMUAsFlag
            clearvars bpdata
        end
        
        MUAe = filtfilt(BMUAe,AMUAe,MUAe);
        MUAe = downsample(MUAe,MUAeDownSampleRatio);
        
        if loadSpikenMUAsFlag
            % Estimate the threshold
            tmpdata = reshape(bpdata,[],1);
            elecThreshold = xEstimation * median(abs(tmpdata)/0.6745,1); % 0.6745 is from experience discribed in iterature
            
            % spike detection
            TrialLengthDP = size(bpdata,1);
            C_poz = bpdata - repmat(elecThreshold,TrialLengthDP,1);
            
            if xEstimation < 0
                C_poz = double(C_poz < 0);
            else
                C_poz = double(C_poz > 0);
            end
            
            C_poz_shifted = [zeros(1,size(C_poz,2)); C_poz(1:end-1,:)];
            C_poz_diff = C_poz - C_poz_shifted;
            clearvars C_poz C_poz_shifted
            
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
                    temp = reshape(bpdata(C_waveIndex,thisTrial),[],waveformLength);
                    if xEstimation < 0
                        [~,temp] = min(temp,[],2);
                    else
                        [~,temp] = max(temp,[],2);
                    end
                    C_wavestart = C_wavestart + temp - round(waveformLength * waveformAlignpoint);
                    C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
                    ElecTrialWaveform = reshape(bpdata(C_waveIndex,thisTrial),[],waveformLength);
                    ElecTrialWaveform = ElecTrialWaveform';
                    waveidx(2) = waveidx(2) + NumTrialSpikes(thisTrial);
                    elecRawWaveform(:,waveidx(1):waveidx(2)) = ElecTrialWaveform;
                    elecRawSpikeStamps(waveidx(1):waveidx(2)) = SpikeSTAligned{thisTrial} + AlignStamps(thisTrial);
                    waveidx(1) = waveidx(2) + 1;
                end
            end
            if size(elecRawWaveform,1) ~= waveformLength
                elecRawWaveform = elecRawWaveform';
            end
        end
        % save data on disk
        electrodeNeuraldata.samplerate = SampleRate;
        electrodeNeuraldata.LFPsamplerate = LFPsamplingrate;
        electrodeNeuraldata.MUAeSamplerate = MUAeSamplingrate;
        electrodeNeuraldata.elecID = EID;
        electrodeNeuraldata.elecIdx = thisElec;
        electrodeNeuraldata.MUAe = MUAe;
        electrodeNeuraldata.LFP = LFP;
        
        if loadSpikenMUAsFlag
            electrodeNeuraldata.MUAshapleySamplerate = MUAshapleySampleRate;
            electrodeNeuraldata.MUAshapley = MUA;
            electrodeNeuraldata.SpikeTime = SpikeTimeAligned;
            electrodeNeuraldata.SpikeST = SpikeSTAligned;
            electrodeNeuraldata.SpikeSTraw = elecRawSpikeStamps;
            electrodeNeuraldata.SpikeWaveform = elecRawWaveform;
            electrodeNeuraldata.SpikeTrain = SpikeTrain;
            electrodeNeuraldata.SpikeTrainEdge = spikeTrainEdge;
            electrodeNeuraldata.Psth10 = Psth10;
            electrodeNeuraldata.Psth10Edge = psth10Edge;
        end
        electrodeNeuraldata.trialParam = trialParam;
%         electrodeNeuraldata.trialInfo = electrodeTrialData.trialInfo;
        idx = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,2,'last');
        idx = idx(1);
        electrodeNeuraldataPath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx),savefilename];
        save(electrodeNeuraldataPath,'electrodeNeuraldata');
        instanceinfo(thisInstance).electrodeNeuraldataPath{thisElec} = electrodeNeuraldataPath;
    end
end
fprintf('\n')