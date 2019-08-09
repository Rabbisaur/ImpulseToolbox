function instanceinfo = VnCgetNeuraldata4(instanceinfo,trialAlignSTfieldname,MUAparameters,LFPparameters, SpikenWaveParameters,trialParam,savefilename,loadSpikenMUAsFlag,linenoiseFreq)
matversion = '-v6';
if ispc
    slash = '\';
else
    slash = '/';
end

%Calculate trial info

LFPsamplingrate=LFPparameters.LFPsamplingrate;
MUAeSamplingrate=MUAparameters.MUAeSamplingrate;

% spike and waveform parameters
xEstimation=SpikenWaveParameters.xEstimation;
waveformLength=SpikenWaveParameters.waveformLength;
waveformAlignpoint=SpikenWaveParameters.waveformAlignpoint;

spikeTrainEdge = (0:0.001:(trialParam.endT - trialParam.startT)) + trialParam.startT;
psth10Edge = (0:0.01:(trialParam.endT - trialParam.startT)) + trialParam.startT;

refractionPeriod = 1.5; % ms

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
    
    % low pass filter 200Hz
    Fl=200;
    Fn = SampleRate/2;
    N = 2;
    [BMUAe, AMUAe] = butter(N,Fl/Fn,'low'); % MUAe low pass
    
    % low pass filter 150Hz
    Fl=150;
    Fn = SampleRate/2;
    N = 2;
    [BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass
    
    % calculate downsample and supersample
    LFPsupersamplingRatio = lcm(SampleRate, LFPsamplingrate) / SampleRate;
    MUAesupersamplingRatio = lcm(SampleRate, MUAeSamplingrate) / SampleRate;
    
    LFPdownSampleRatio = SampleRate * LFPsupersamplingRatio / LFPsamplingrate;
    MUAeDownSampleRatio = SampleRate * MUAesupersamplingRatio / MUAeSamplingrate;
    
    for thisElec = 1:numElec
        fprintf('.')
        EID = instanceinfo(thisInstance).ElecOrder(thisElec);
        if EID > 128
            break
        end
        cmd = ['trialAlignST = instanceinfo(thisInstance).trialInfo.',trialAlignSTfieldname,';'];
        eval(cmd);
        
        electrodeTrialData = VnCGetTrialData2(instanceinfo, thisInstance, EID,trialAlignST,trialParam);

        tmpdata = electrodeTrialData.trialData';
        AlignStamps = electrodeTrialData.AlignStamps;
        clearvars electrodeTrialData
        numtrials = size(tmpdata,2);

        % remove 50Hz
        tmpdata = filtfilt(Bbs,Abs,tmpdata);

        % low pass filter and downsample to get LFP
        LFP = filtfilt(BLFP,ALFP,tmpdata);
        LFP = downsample(LFP,LFPdownSampleRatio);
        
        % band pass filter
        bpdata = filtfilt(Bbp,Abp,tmpdata);
        clearvars tmpdata
        
        % MUAe
        MUAe = abs(bpdata);
        MUAe = filtfilt(BMUAe,AMUAe,MUAe);
        
        if MUAesupersamplingRatio > 1
            x = 1:length(MUAe);
            xx = 1:(1/MUAesupersamplingRatio):size(MUAe,1);
            MUAe = interp1(x,MUAe,xx);
        end
        data = downsample(MUAe,MUAeDownSampleRatio);
        clearvars MUAe;
        MUAe.data = data;
        MUAe.MUAesamplerate = MUAeSamplingrate;
        MUAe.time =  (0:(size(MUAe.data,1)-1))/MUAe.MUAesamplerate + trialParam.startT;
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
        data = LFP;
        clearvars LFP;
        LFP.data = data;
        LFP.samplerate = LFPsamplingrate;
        LFP.time =  (0:(size(LFP.data,1)-1))/LFP.samplerate + trialParam.startT;
        
        % removing 50Hz Xing's method
        muafilt = MUAe.data;
        FsD = MUAe.MUAesamplerate;
        Fn = FsD/2; % Downsampled Nyquist frequency
        for v = [linenoiseFreq linenoiseFreq*2 linenoiseFreq*3]
            Fbp = [v-2,v+2];
            [Blp, Alp] = butter(2, [min(Fbp)/Fn max(Fbp)/Fn],'stop'); % compute filter coefficients
            muafilt = filtfilt(Blp, Alp, muafilt);
        end
        MUAe.data = muafilt;
        
        lfpfilt = LFP.data;
        FsD = LFP.samplerate;
        Fn = FsD/2; % Downsampled Nyquist frequency
        for v = [linenoiseFreq linenoiseFreq*2 linenoiseFreq*3]
            Fbp = [v-2,v+2];
            [Blp, Alp] = butter(2, [min(Fbp)/Fn max(Fbp)/Fn],'stop'); % compute filter coefficients
            lfpfilt = filtfilt(Blp, Alp, lfpfilt);
        end
        LFP.data = lfpfilt;
        
        electrodeNeuraldata.elecID = EID;
        electrodeNeuraldata.elecIdx = thisElec;
        electrodeNeuraldata.MUAe = MUAe;
        electrodeNeuraldata.LFP = LFP;

        if loadSpikenMUAsFlag
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
        idx = find(instanceinfo(thisInstance).electrodeCachePath{EID} == slash,2,'last');
        idx = idx(1);
        electrodeNeuraldataPath = [instanceinfo(thisInstance).electrodeCachePath{EID}(1:idx),savefilename];
        save(electrodeNeuraldataPath,'electrodeNeuraldata',matversion);
        instanceinfo(thisInstance).electrodeNeuraldataPath{thisElec} = electrodeNeuraldataPath;
    end
end
fprintf('\n')