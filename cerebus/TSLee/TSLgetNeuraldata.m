function instanceinfo = TSLgetNeuraldata(instanceinfo,TrialAlignST,MUAparameters,LFPparameters, SpikenWaveParameters,trialParam,savefilename,loadLFPflag, loadMUAflag,loadSpikenMUAsFlag)

if loadLFPflag + loadMUAflag + loadSpikenMUAsFlag ==0
    disp('No load flag is set. Return!')
    return
end

%Calculate trial info
if loadLFPflag
    LFPsamplingrate=LFPparameters.LFPsamplingrate;
end

if loadMUAflag
    MUAeSamplingrate=MUAparameters.MUAeSamplingrate;
end

if loadSpikenMUAsFlag == 2
    MUAshapleySampleRate=MUAparameters.MUAshapleySampleRate;
end

if loadSpikenMUAsFlag
    % spike and waveform parameters
    xEstimation=SpikenWaveParameters.xEstimation;
    waveformLength=SpikenWaveParameters.waveformLength;
    waveformAlignpoint=SpikenWaveParameters.waveformAlignpoint;
    
    spikeTrainEdge = (0:0.001:(trialParam.endT - trialParam.startT)) + trialParam.startT;
    psth10Edge = (0:0.01:(trialParam.endT - trialParam.startT)) + trialParam.startT;
    
    refractionPeriod = 1.5; % ms
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters ends here

% calculate LFP and MUAs and MUAe

numElec = instanceinfo.numElec;

SampleRate = instanceinfo.samplerate;
AlignIdx = -trialParam.startT * SampleRate;

% design filters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % design the filters

% band stop 50Hz
Fn = SampleRate/2;
Fbp=[49,51];
N  = 2;    % filter order
[Bbs, Abs] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn],'stop'); % Bandstop

% band pass 500-5000Hz
Fn = SampleRate/2;
Fbp=[500,5000];
N  = 2;    % filter order
[Bbp, Abp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn]); % BandPass


if loadMUAflag
    % MUA low pass filter 200Hz
    Fl=200;
    Fn = SampleRate/2;
    N = 2;
    [BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass
end

if loadLFPflag
    % LFP low pass filter 150Hz
    Fl=150;
    Fn = SampleRate/2;
    N = 2;
    [BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass
end

% calculate downsample and supersample
if loadLFPflag
    LFPsupersamplingRatio = lcm(SampleRate, LFPsamplingrate) / SampleRate;
    LFPdownSampleRatio = SampleRate * LFPsupersamplingRatio / LFPsamplingrate;
end
if loadMUAflag
    MUAesupersamplingRatio = lcm(SampleRate, MUAeSamplingrate) / SampleRate;
    MUAeDownSampleRatio = SampleRate * MUAesupersamplingRatio / MUAeSamplingrate;
end

for thisElec = 1:numElec
    fprintf('.')
    EID = instanceinfo.ElecOrder(thisElec);
    if EID > 128
        break
    end
    
    electrodeTrialData = TSLGetTrialData(instanceinfo, EID, TrialAlignST,trialParam);
    
    tmpdata = electrodeTrialData.trialData;
    AlignStamps = electrodeTrialData.AlignStamps;
    clearvars electrodeTrialData
    numtrials = size(tmpdata,2);
    
    % remove 50Hz
    tmpdata = filtfilt(Bbs,Abs,tmpdata);
    
    if loadLFPflag
        % low pass filter and downsample to get LFP
        LFP = filtfilt(BLFP,ALFP,tmpdata);
        LFP = downsample(LFP,LFPdownSampleRatio);
    end
    
    if loadSpikenMUAsFlag == 2
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
    MUAe = filtfilt(BMUAe,AMUAe,MUAe);
    
    if MUAesupersamplingRatio > 1
        x = 1:length(MUAe);
        xx = 1:(1/MUAesupersamplingRatio):length(MUAe);
        MUAe = spline(x,MUAe',xx)';
    end
    MUAe = downsample(MUAe,MUAeDownSampleRatio);
    
    if loadSpikenMUAsFlag
        % Estimate the threshold
        tmpdata = reshape(bpdata,[],1);
        elecThreshold = xEstimation * median(abs(tmpdata)/0.6745,1); % 0.6745 is from experience discribed in iterature
        clearvars tmpdata
        
        % spike detection
        TrialLengthDP = size(bpdata,1);
        NumTrials = size(bpdata,2);
        C_poz = bpdata - repmat(elecThreshold,TrialLengthDP,NumTrials);
        
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
            
            % refraction period
            NumSpikes = numel(TrialSpikes);
            if NumSpikes > 0
                
                CurrIdx = 1;
                GoodIdx = false(NumSpikes,1);
                GoodIdx(CurrIdx) = 1;
                for thisSpike = 2:NumSpikes
                    ISI = TrialSpikes(thisSpike) - TrialSpikes(CurrIdx);
                    if ISI <= refractionPeriod * SampleRate / 1000
                        continue
                    else
                        CurrIdx = thisSpike;
                        GoodIdx(CurrIdx) = true;
                    end
                end
                TrialSpikes = TrialSpikes(GoodIdx);
            end
            
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
                %                     temp = reshape(bpdata(C_waveIndex,thisTrial),[],waveformLength);
                %                     if xEstimation < 0
                %                         [~,temp] = min(temp,[],2);
                %                     else
                %                         [~,temp] = max(temp,[],2);
                %                     end
                %                     C_wavestart = C_wavestart + temp - round(waveformLength * waveformAlignpoint);
                %                     C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
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
    %         plot(elecRawWaveform)
    % save data on disk
    electrodeNeuraldata.samplerate = instanceinfo.samplerate;
    
    electrodeNeuraldata.elecID = EID;
    electrodeNeuraldata.elecIdx = thisElec;
    
    if loadMUAflag
        electrodeNeuraldata.MUAeSamplerate = MUAeSamplingrate;
        electrodeNeuraldata.MUAe = MUAe;
        electrodeNeuraldata.MUAeTime = trialParam.startT:(1/MUAeSamplingrate):trialParam.endT;
    end
    
    if loadLFPflag
        electrodeNeuraldata.LFP = LFP;
        electrodeNeuraldata.LFPsamplerate = LFPsamplingrate;
        electrodeNeuraldata.LFPTime = trialParam.startT:(1/LFPsamplingrate):trialParam.endT;
    end
    
    if loadSpikenMUAsFlag == 2
        electrodeNeuraldata.MUAshapleySamplerate = MUAshapleySampleRate;
        electrodeNeuraldata.MUAshapley = MUA;
        electrodeNeuraldata.MUAshapleyTime = trialParam.startT:(1/MUAshapleySampleRate):trialParam.endT;
    end
    
    if loadSpikenMUAsFlag
        electrodeNeuraldata.SpikeTime = SpikeTimeAligned;
        electrodeNeuraldata.SpikeST = SpikeSTAligned;
        electrodeNeuraldata.SpikeSTraw = elecRawSpikeStamps;
        electrodeNeuraldata.SpikeWaveform = elecRawWaveform;
        electrodeNeuraldata.SpikeTrain = SpikeTrain;
        electrodeNeuraldata.SpikeTrainEdge = spikeTrainEdge;
        electrodeNeuraldata.Psth10 = Psth10;
        electrodeNeuraldata.Psth10Time = psth10Edge(1:(end-1));
    end
    electrodeNeuraldata.trialParam = trialParam;

    idx = strfind(instanceinfo.electrodeCachePath{EID},'tmp');
    idx = idx(1);
    electrodeNeuraldataPath = [instanceinfo.electrodeCachePath{EID}(1:(idx-1)),savefilename];
    save(electrodeNeuraldataPath,'electrodeNeuraldata');
    instanceinfo.electrodeNeuraldataPath{EID} = electrodeNeuraldataPath;
end
end