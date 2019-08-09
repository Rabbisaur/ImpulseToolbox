function TSLloadNS6ElecNoSortingParPause2(filepath,xEstimation)
% filepath - full or relative path to the ns5 file
% savewave - 0,1, whether save spike waveform or not, default = 1
% xEstimation - default -5, auto thresholding with x times estimation
%
% LoadNs5 will filter 30KHz Nsx files and find zero phase shifted spike
% time, construct experiment event matrix from nev and save wave form.
%
% wrote and revised by Wang Feng @ BNU, tested and with honour.
%
% latest revision date: Oct. 12, 2016, by Wang Feng @ CMU
%

% verify the validity of input arguments
switch nargin
    case 1
        %         savewave = 1; % spike sorting is very important and mean to be the very next step before any analysis. So if you don't want store spike wave form, state it explicitly.
        xEstimation = -3.5; % auto thresholding with x times estimation
    case 2
        % currently doing nothing
    otherwise
        fprintf('failure: LoadNs5 accepts 3 input arguments.\n'); return;
end

nevmatfilepath = [filepath(1:end-3),'mat'];
if exist(nevmatfilepath,'file') == 2
    delete(nevmatfilepath)
end

if ispc
    slash = '\';
else
    slash= '/';
end
tic;
%% configurations all goes here

datachunksize = 100; % seconds
overlapsize = 1; % seconds

% ref_period = 1; %ms
waveformLength = 64; % how many sample points in a waveform
waveformAlignpoint = 0.33; % align waveform at which point, round(waveformLength * waveformAlignpoint)
matversion = '-v7.3';

discriminationWindow = [0 0.00027];
window1 = [-800 -150];
window2 = [0 300];
LFPSampleRate = 2000; %Hz
MUASampleRate = 1000; %Hz
% TrialTimeBeforeAlign = 0.2; % Sec
LFPTimeBeforeAlign = 0.2; % Sec
LFPTimeAfterAlign = 0.600; % Sec
MUATimeBeforeAlign = 0.2; % Sec
MUATimeAfterAlign = 0.600; % Sec

% NumberPar = 8; % number of channels for parallel processing.
%% construct experiment event matrix from nev file
% check the validity of the nev & mbm files specified
Ns6Path = [filepath(1: end-3) 'ns6'];
idx = find(Ns6Path == slash, 1, 'last');
DataDir = Ns6Path(1:idx-1);
cd(DataDir);
if exist(Ns6Path,'file') == 2
    % create folders for storing mat files
    matdir = filepath(1: end-4);
    if exist(matdir,'dir') == 7
        rmdir(matdir,'s');
    end
    if mkdir(matdir)
    else
        error('matlab:mkdirfailed','failed to make directory for the mat file')
    end
    fileattrib(matdir, '+w');
    nevpath = [filepath(1: end-3) 'nev'];
    if exist(nevpath,'file') == 2
        %         fprintf('nev file found ...\n');
    else
        error('matlab:filenotfound','nev file not found!')
    end
else
    error('matlab:filenotfound','ns6 file not found!');
end
nevpath = [filepath(1:end-4) '.nev'];
expmarkpath = [filepath(1:end-4),slash,'Expmark.mat'];
TSLloadNevPause(nevpath);
load(expmarkpath);

%% load 30KHz continuous signal.
disp('Loading NS6 file')

% process basic header
NsxID = fopen(Ns6Path,'r');
fseek(NsxID, 0, 1);
NsxSize = uint64(ftell(NsxID));
fprintf('NS6 file size is %dMB\n',NsxSize/1024^2)
frewind(NsxID);
NsxBasic.FileType = sprintf('%s', fread(NsxID, 8, '*char'));
% check file type ID
if ~strcmp(NsxBasic.FileType, 'NEURALCD')
    error('NS6 FileType mismatch!')
end
VerNum = fread(NsxID, 2, 'uchar=>uchar');
NsxBasic.FileVersion = ['Spec. ' num2str(VerNum(1)) '.' num2str(VerNum(2))];
NsxBasic.HeaderSize = fread(NsxID, 1, 'uint32=>uint32');
NsxBasic.Label = sprintf('%s', fread(NsxID, 16, 'char=>uint32'));
% NsxBasic.Comment = sprintf('%s', fread(NsxID, 256, '*char'));
fseek(NsxID, 256, 0);
NsxBasic.Period = fread(NsxID, 1, 'uint32=>uint32');
NsxBasic.ClockFs = fread(NsxID, 1, 'uint32=>uint32');
SampleRate = double(NsxBasic.ClockFs);
NsxBasic.Fs = NsxBasic.ClockFs/NsxBasic.Period;
NsxBasic.TimeOrigin = fread(NsxID, [1 8], 'uint16=>uint16');
fprintf('NS6 file created at %d:%d:%d %d-%d-%d\n',NsxBasic.TimeOrigin(5),NsxBasic.TimeOrigin(6),NsxBasic.TimeOrigin(7), ...
    NsxBasic.TimeOrigin(1),NsxBasic.TimeOrigin(2),NsxBasic.TimeOrigin(4));
NsxBasic.NumElec = fread(NsxID, 1, 'uint32=>uint32');
CurrPos = ftell(NsxID);
if CurrPos~=314 % total length of basic header of v2.2
    error('FAILURE: Error in reading basic header');
end
if CurrPos+NsxBasic.NumElec*66~=NsxBasic.HeaderSize
    error('FAILURE: Error in the size of extended headers.');
end
fclose(NsxID);
% process extended headers
NsxID = fopen(Ns6Path);
fseek(NsxID, 314, -1);
NsxExtend = cell(NsxBasic.NumElec, 1);
PacketID = fread(NsxID, [2 NsxBasic.NumElec], '2*char=>char', 66-2);
ElecOrder = zeros(NsxBasic.NumElec, 1);
for thisElec = 1:NsxBasic.NumElec
    CurrStartPos = CurrPos+(thisElec-1)*66+2;  %electrode ID
    fseek(NsxID, CurrStartPos, 'bof');
    NsxExtend{thisElec,1}.PacketID = sprintf('%s', PacketID(:,thisElec));
    switch NsxExtend{thisElec,1}.PacketID
        case 'CC'  %continuous channels
            ElecOrder(thisElec,1) = fread(NsxID, 1, 'uint16');
            NsxExtend{thisElec,1}.ElecID = ElecOrder(thisElec,1);
            NsxExtend{thisElec,1}.ElecLabel = sprintf('%s', fread(NsxID, 16, 'char=>char'));
            NsxExtend{thisElec,1}.Connector = char(64+fread(NsxID, 1, 'char=>char'));
            NsxExtend{thisElec,1}.Pin = fread(NsxID, 1, 'char=>char');
            NsxExtend{thisElec,1}.MinDigitalValue = fread(NsxID, 1, 'int16=>int16');
            NsxExtend{thisElec,1}.MaxDigitalValue = fread(NsxID, 1, 'int16=>int16');
            NsxExtend{thisElec,1}.MinAnalogValue = fread(NsxID, 1, 'int16=>int16');
            NsxExtend{thisElec,1}.MaxAnalogValue = fread(NsxID, 1, 'int16=>int16');
            NsxExtend{thisElec,1}.AnalogUnit = sprintf('%s', fread(NsxID, 16, 'char=>char'));  %mv/uv
            NsxExtend{thisElec,1}.HighCutFreq = fread(NsxID, 1, 'uint32=>uint32');
            NsxExtend{thisElec,1}.HighCutOrder = fread(NsxID, 1, 'uint32=>uint32');  %0 = NONE
            NsxExtend{thisElec,1}.HighCutType = fread(NsxID, 1, 'uint16=>uint16');
            NsxExtend{thisElec,1}.LowCutFreq = fread(NsxID, 1, 'uint32=>uint32');
            NsxExtend{thisElec,1}.LowCutOrder = fread(NsxID, 1, 'uint32=>uint32');  %0 = NONE
            NsxExtend{thisElec,1}.LowCutType = fread(NsxID, 1, 'uint16=>uint16');
        otherwise
            disp('FAILURE: No current PacketID was found.');
            return;
    end
end
if numel(ElecOrder)~=numel(unique(ElecOrder))
    disp('Error in reading extended headers.');
    return;
end
CurrPos = ftell(NsxID);
if CurrPos~=NsxBasic.HeaderSize
    disp('FAILURE: Error in reading extended headers');
    return;
end
clear VerNum PacketID CurrStartPos
fclose(NsxID);

NsxID = fopen(Ns6Path);
fseek(NsxID,NsxBasic.HeaderSize,-1);


NumValidElec = numel(ElecOrder);
fprintf('Got %d valid electrodes.\n',NumValidElec);
disp('Seeking valid trials...')
% Seek for trials
TrialNum = 0;
while ftell(NsxID) < NsxSize
    TrialNum = TrialNum + 1;
    Trial.position(TrialNum) = ftell(NsxID);
    Header = fread(NsxID,1,'int8=>int8');
    if Header ~= 1
        error('Header seek error, abort!');
    end
    Timestamp = fread(NsxID,1,'uint32=>double');
    Trial.Timestamp(TrialNum) = Timestamp;
    NumDP = fread(NsxID,1,'uint32=>double');
    Trial.NumDP(TrialNum) = NumDP;
    fseek(NsxID,Trial.NumDP(TrialNum)*NumValidElec*2,0); % jump to the beginning of the next data package
end

ValidityIdx = RawExpmark(3,:) > 0;

ValidTrial.NumValidTrials = sum(ValidityIdx);
if ValidTrial.NumValidTrials ~= size(ExpmarkTime,2)
    error('Number of valid trials mismatch!');
end

% temporary solotion for Summer's data, don't know why, need further
% investigation.
if numel(Trial.position) ~= numel(ValidityIdx)
    ValidityIdx(1) = [];
end
ValidTrial.position = Trial.position(ValidityIdx);
ValidTrial.Timestamp = Trial.Timestamp(ValidityIdx);
ValidTrial.NumDP = Trial.NumDP(ValidityIdx);

% in order to handle the unfinished experiments, valid trials should be
% truncated to the max number of finished blocks times the number of
% conditions
StimCond = ExpmarkTime(2,:);
UniqueStimCond = unique(StimCond);
NumUniStimCond = numel(UniqueStimCond);
CondNumValidTrials = zeros(NumUniStimCond,1);
for thisCond = 1:NumUniStimCond
    CondNumValidTrials(thisCond) = sum(StimCond == UniqueStimCond(thisCond));
end
MaxNumFinishedBlocks = min(CondNumValidTrials);
NumValidTrials = MaxNumFinishedBlocks * NumUniStimCond;
if ValidTrial.NumValidTrials ~= NumValidTrials
    disp('Experiment unfinished, will only save max number of finished blocks')
    disp(['which are ', num2str(MaxNumFinishedBlocks),' blocks, a total of ', num2str(NumValidTrials), ' trials, instead of '...
        num2str(ValidTrial.NumValidTrials),' trials.'])
    ValidTrial.NumValidTrials = NumValidTrials;
    ValidTrial.position = ValidTrial.position(1:NumValidTrials);
    ValidTrial.Timestamp = ValidTrial.Timestamp(1:NumValidTrials);
    ValidTrial.NumDP = ValidTrial.NumDP(1:NumValidTrials);
    % truncate Expmark as well as nev data
    elecpath = [filepath(1:end-4),slash,'electrodes.mat'];
    load(elecpath)
    spikepath = [filepath(1:end-4),slash,'SpikeTime.mat'];
    load(spikepath)
    spikeSTpath = [filepath(1:end-4),slash,'SpikeTimeStamp.mat'];
    load(spikeSTpath)
    for thisElec = 1:numel(electrodes)
        EID = electrodes(thisElec);
        SpikeTimeAligned{EID} = SpikeTimeAligned{EID}(1:NumValidTrials);
        SpikeTimeStamp{EID} = SpikeTimeStamp{EID}(1:NumValidTrials);
        SpikeTimeStampAligned{EID} = SpikeTimeStampAligned{EID}(1:NumValidTrials);
    end
    save(spikepath,'SpikeTimeAligned');
    save(spikeSTpath,'SpikeTimeStampAligned','SpikeTimeStamp');
    waveformpath = [filepath(1:end-4),slash,'Waveform.mat'];
    load(waveformpath)
    for thisElec = 1:numel(electrodes)
        EID = electrodes(thisElec);
        NumSpikes = 0;
        for thisTrial = 1:NumValidTrials
            NumSpikes = NumSpikes + numel(SpikeTimeAligned{EID}{thisTrial});
        end
        Waveform{EID} = Waveform{EID}(:,1:NumSpikes);
    end
    save(waveformpath,'Waveform')
    clearvars SpikeTimeAligned SpikeTimeStampAligned SpikeTimeStamp Waveform
    ExpmarkST = ExpmarkST(:,1:NumValidTrials);
    ExpmarkTime = ExpmarkTime(:,1:NumValidTrials);
    save(expmarkpath,'ExpmarkST','ExpmarkTime');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% design the filter
% try to remove 60Hz from the signal
d = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',SampleRate);

% band pass filter 300-6000Hz
bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',300,'HalfPowerFrequency2',6000, ...
    'DesignMethod','butter','SampleRate',SampleRate);

% low pass filter < 150 Hz
lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',150,...
    'DesignMethod','butter','SampleRate',SampleRate);

% High pass filter > 1000 Hz
hpFilt = designfilt('highpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',1000,...
    'DesignMethod','butter','SampleRate',SampleRate);

idx = find(filepath == '.', 1, 'last');
savebasedir = filepath(1:(idx-1));
% tmpdir = [savebasedir,slash,'tmp'];
% if exist(tmpdir,'dir') == 7
%     rmdir(tmpdir,'s');
% end
% mkdir(tmpdir);

disp('Making electrode cache files')
elecfp = zeros(NumValidElec,1);
for thisElec = 1:NumValidElec
    EID = ElecOrder(thisElec);
    savedir = [savebasedir,slash,'elec',num2str(EID)];
    if exist(savedir,'dir')==7
        rmdir(savedir,'s');
    end
    mkdir(savedir);
    mkdir([savedir,slash,'tmp']);
    tmpfilepath = [savedir,slash,'tmp',slash,'elecrawdata.tmp'];
    elecfp(thisElec) = fopen(tmpfilepath,'w');
end

% find the shortest valid trial, cut and align trials accordingly
AlignST = ExpmarkST(6,:);
TrialStartST = ValidTrial.Timestamp;
TrialEndST = ExpmarkST(4,:);
TrialLengthbeforAlign = AlignST - TrialStartST;
TrialLengthafterAlign = TrialEndST-AlignST;
minTrialLengthafterAlign = min(TrialLengthafterAlign);
minTrialLengthbeforeAlign = min(TrialLengthbeforAlign);

% split nsx file into electrodes
ValidTrial.tmpposition = zeros(size(ValidTrial.position));
ValidTrial.tmptimestamp = zeros(size(ValidTrial.position));
ValidTrial.tmpnumdp = zeros(size(ValidTrial.position));

for thisTrial = 1:ValidTrial.NumValidTrials
    pos = ValidTrial.position(thisTrial);
    fseek(NsxID,pos,-1);
    
    Header = fread(NsxID,1,'int8=>int8');
    if Header ~= 1
        error('Header seek error, abort!');
    end
    Timestamp = fread(NsxID,1,'uint32=>double');
    
    if ValidTrial.Timestamp(thisTrial) ~= Timestamp
        error('Timestamp mismatch!')
    end
    NumDP = fread(NsxID,1,'uint32=>double');
    if ValidTrial.NumDP(thisTrial) ~= NumDP
        error('NumDP mismatch!')
    end
    
    DataLength = ValidTrial.NumDP(thisTrial);
    
    TrialData = fread(NsxID,[NumValidElec,DataLength],'int16=>double')';
    
    % align and cut data
    TrialStartIdx = AlignST(thisTrial) - minTrialLengthbeforeAlign - Timestamp + 1;
    TrialEndIdx = TrialStartIdx + minTrialLengthafterAlign + minTrialLengthbeforeAlign -1;
    ValidTrial.tmpnumdp(thisTrial) = TrialEndIdx - TrialStartIdx + 1;
    TrialData = TrialData(TrialStartIdx:TrialEndIdx,:);
    ValidTrial.tmptimestamp(thisTrial) = Timestamp + TrialStartIdx - 1;
    
    for thisElec = 1:NumValidElec
        ValidTrial.tmpposition(thisTrial) = ftell(elecfp(thisElec));
        fwrite(elecfp(thisElec),uint8(Header),'uint8'); % write header
        fwrite(elecfp(thisElec),uint32(ValidTrial.tmptimestamp(thisTrial)),'uint32'); % write Timestamp
        fwrite(elecfp(thisElec),uint32(ValidTrial.tmpnumdp(thisTrial)),'uint32'); % write NumDP
        elecData = int16(TrialData(:,thisElec));
        fwrite(elecfp(thisElec),elecData,'int16');
    end
    
end


for thisElec = 1:NumValidElec
    fclose(elecfp(thisElec));
end
fclose(NsxID);

for thisElec = 1:NumValidElec
    EID = ElecOrder(thisElec);
    savedir = [savebasedir,slash,'elec',num2str(EID)];
    tmpfilepath = [savedir,slash,'tmp',slash,'elecrawdata.tmp'];
    elecfp(thisElec) = fopen(tmpfilepath,'r');
end

disp('Reading electrode data back')
for thisElec = 1:NumValidElec
    EID = ElecOrder(thisElec);
    savedir = [savebasedir,slash,'elec',num2str(EID)];
    disp(['Processing electrode #',num2str(EID)])
    % read elec tmp file back
    ElecData = zeros(ValidTrial.tmpnumdp(1),ValidTrial.NumValidTrials);
    for thisTrial = 1:ValidTrial.NumValidTrials
        pos = ValidTrial.tmpposition(thisTrial);
        fseek(elecfp(thisElec),pos,-1);
        
        Header = fread(elecfp(thisElec),1,'int8=>int8');
        if Header ~= 1
            error('Header seek error, abort!');
        end
        Timestamp = fread(elecfp(thisElec),1,'uint32=>double');
        
        if ValidTrial.tmptimestamp(thisTrial) ~= Timestamp
            error('Timestamp mismatch!')
        end
        NumDP = fread(elecfp(thisElec),1,'uint32=>double');
        if ValidTrial.tmpnumdp(thisTrial) ~= NumDP
            error('NumDP mismatch!')
        end
        ElecData(:,thisTrial) = fread(elecfp(thisElec),NumDP,'int16=>double')';
    end
    % remove line noise
    ElecData = filtfilt(d,ElecData);
    % lowpass to get the LFP
    LFP = filtfilt(lpFilt,ElecData);
    % down sampling LFP
    LFP = LFP(1:SampleRate/LFPSampleRate:size(LFP,1),:);
    % align LFP
    LFPStart = round((minTrialLengthbeforeAlign - LFPTimeBeforeAlign * SampleRate) / (SampleRate/LFPSampleRate));
    LFPEnd = round((minTrialLengthbeforeAlign + LFPTimeAfterAlign * SampleRate)/(SampleRate/LFPSampleRate));
    LFP = LFP(LFPStart:LFPEnd,:);
    save([savedir,slash, 'LFP6.mat'],'LFP');
    clearvars LFP;
    
    % hipass to get the MUA
    rawElecMUA = filtfilt(hpFilt,ElecData);
    %     half wave rectify in the negative direction
    rawElecMUA = rawElecMUA.*(rawElecMUA<0);
    meanElecMUA = mean(rawElecMUA,1);
    stdElecMUA = std(rawElecMUA,0,1);
    % count the number of events more than 3 times std with in 1 ms time
    % bins
    MUAthres = meanElecMUA - 3 * stdElecMUA;
    %     edge = 1:(SampleRate/MUASampleRate):(size(rawElecMUA,1)+(SampleRate/MUASampleRate));
    
    rawElecMUA = double(rawElecMUA < repmat(MUAthres,size(rawElecMUA,1),1));
    MUAlength = size(rawElecMUA,1);
    NumPointsperBin = round(SampleRate/MUASampleRate);
    NumBins = floor(MUAlength / NumPointsperBin);
    rawElecMUA = rawElecMUA(1:(NumPointsperBin*NumBins),:);
    MUA = zeros(NumBins,ValidTrial.NumValidTrials);
    for thispoint = 1:NumPointsperBin
        MUA = MUA+ rawElecMUA(thispoint:NumPointsperBin:(end-(NumPointsperBin-thispoint)),:);
    end
    
    %     MUA = zeros(numel(edge),ValidTrial.NumValidTrials);
    %     for thisTrial = 1: ValidTrial.NumValidTrials
    %         temp = find(rawElecMUA(:,thisTrial));
    %         MUA(:,thisTrial) = histc(temp,edge);
    %     end
    
    % align MUA
    MUAStart = round((minTrialLengthbeforeAlign - MUATimeBeforeAlign * SampleRate) / (SampleRate/MUASampleRate));
    MUAEnd = round((minTrialLengthbeforeAlign + MUATimeAfterAlign * SampleRate)/(SampleRate/MUASampleRate));
    MUA = MUA(MUAStart:MUAEnd,:);
    
    clearvars rawElecMUA
    save([savedir,slash, 'MUA6.mat'],'MUA');
    clearvars MUA
    
    % bandpass filter
    ElecData = filtfilt(bpFilt,ElecData);
    
    
    
    
    
    
    
    % Estimate the threshold
    %         disp('Estimating threshold')
    elecThreshold = xEstimation * median(abs(ElecData)/0.6745,1); % 0.6745 is from experience discribed in iterature
    
    % spike detection
    %         disp('Detecting spikes')
    DataChunkLength = size(ElecData,1);
    C_poz = ElecData - repmat(elecThreshold,DataChunkLength,1);
    
    if xEstimation < 0
        C_poz = double(C_poz < 0);
    else
        C_poz = double(C_poz > 0);
    end
    
    C_poz_shifted = [zeros(1,size(C_poz,2)); C_poz(1:end-1,:)];
    C_poz_diff = C_poz - C_poz_shifted;
    clearvars C_poz C_poz_shifted
    
    SpikeSTAligned = cell(ValidTrial.NumValidTrials,1);
    SpikeTimeAligned = cell(ValidTrial.NumValidTrials,1);
    RawElecSpikes = cell(ValidTrial.NumValidTrials,1);
    NumTrialSpikes = zeros(ValidTrial.NumValidTrials,1);
    for thisTrial = 1: ValidTrial.NumValidTrials
        % get spikes
        TrialSpikes = find(C_poz_diff(:,thisTrial) == 1);
        % exclude any spikes within one waveform length from the beginning
        % and one waveform length at the end
        idx = TrialSpikes < 2*waveformLength | size(C_poz_diff,1) - TrialSpikes < 2*waveformLength;
        TrialSpikes(idx) = [];
        RawElecSpikes{thisTrial} = TrialSpikes;
        SpikeSTAligned{thisTrial} = TrialSpikes + ValidTrial.tmptimestamp(thisTrial) - AlignST(thisTrial);
        SpikeTimeAligned{thisTrial} = SpikeSTAligned{thisTrial} / SampleRate;
        NumTrialSpikes(thisTrial) = numel(TrialSpikes);
    end
    % get waveform
    elecRawWaveform = zeros(waveformLength,sum(NumTrialSpikes));
    elecRawSpikeStamps = zeros(sum(NumTrialSpikes),1);
    waveidx(1) = 1;
    waveidx(2) = 0;
    for thisTrial = 1: ValidTrial.NumValidTrials
        if ~isempty(RawElecSpikes{thisTrial})
            C_wavestart = RawElecSpikes{thisTrial} - round(waveformLength * waveformAlignpoint); % waveformLength is the length of padded data
            C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
            %         align to the minimal or maximal
            temp = reshape(ElecData(C_waveIndex,thisTrial),[],waveformLength);
            if xEstimation < 0
                [~,temp] = min(temp,[],2);
            else
                [~,temp] = max(temp,[],2);
            end
            C_wavestart = C_wavestart + temp - round(waveformLength * waveformAlignpoint);
            C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
            ElecTrialWaveform = reshape(ElecData(C_waveIndex,thisTrial),[],waveformLength);
            ElecTrialWaveform = ElecTrialWaveform';
            waveidx(2) = waveidx(2) + NumTrialSpikes(thisTrial);
            elecRawWaveform(:,waveidx(1):waveidx(2)) = ElecTrialWaveform;
            elecRawSpikeStamps(waveidx(1):waveidx(2)) = RawElecSpikes{thisTrial} + ValidTrial.tmptimestamp(thisTrial);
            waveidx(1) = waveidx(2) + 1;
        end
    end
    if size(elecRawWaveform,1) ~= waveformLength;
        elecRawWaveform = elecRawWaveform';
    end
    % save data
    save([savedir,slash, 'spikeST6.mat'],'SpikeSTAligned');
    save([savedir,slash, 'SpikeTime6.mat'],'SpikeTimeAligned');
    save([savedir,slash, 'RawWaveform6.mat'],'elecRawWaveform');
    save([savedir,slash, 'Threshold.mat'],'elecThreshold');
    save([savedir,slash, 'RawSpikeStamps6.mat'],'elecRawSpikeStamps');
    
end
disp('Cleaning up...')
% close files
for thisElec = 1:NumValidElec
    fclose(elecfp(thisElec));
end
% clean up
for thisElec = 1:NumValidElec
    EID = ElecOrder(thisElec);
    savedir = [savebasedir,slash,'elec',num2str(EID)];
    rmdir([savedir,slash,'tmp'],'s');
end
disp('Done!')