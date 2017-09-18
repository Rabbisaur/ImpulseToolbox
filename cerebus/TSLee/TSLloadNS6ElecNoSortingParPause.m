function TSLloadNS6ElecNoSortingParPause(filepath,xEstimation)
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

if ispc
    slash = '\';
else
    slash= '/';
end

nevmatfilepath = [filepath(1:end-3),'mat'];
if exist(nevmatfilepath,'file') == 2
    delete(nevmatfilepath)
end
%% configurations all goes here

datachunksize = 100; % seconds
overlapsize = 1; % seconds

% ref_period = 1; %ms
waveformLength = 64; % how many sample points in a waveform
waveformAlignpoint = 0.25; % align waveform at which point, round(waveformLength * waveformAlignpoint)
matversion = '-v7.3';

discriminationWindow = [0 0.00027];
window1 = [-800 -150];
window2 = [0 300];
LFPSampleRate = 2000; %Hz
MUASampleRate = 1000; %Hz
LFPTimeBeforeAlign = 0.2; % Sec
LFPTimeAfterAlign = 0.600; % Sec
MUATimeBeforeAlign = 0.2; % Sec
MUATimeAfterAlign = 0.600; % Sec

NumberPar = 8; % number of channels for parallel processing.
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
if ValidTrial.NumValidTrials ~= size(ExpmarkTime,2);
    error('Number of valid trials mismatch!');
end
ValidTrial.position = Trial.position(ValidityIdx);
ValidTrial.Timestamp = Trial.Timestamp(ValidityIdx);
ValidTrial.NumDP = Trial.NumDP(ValidityIdx);

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


threshold = zeros(ValidTrial.NumValidTrials,NumValidElec);

idx = find(filepath == '.', 1, 'last');
savebasedir = filepath(1:(idx-1));
tmpdir = [savebasedir,slash,'tmp'];
if exist(tmpdir,'dir') == 7
    rmdir(tmpdir,'s');
end
mkdir(tmpdir);

for thisElec = 1:NumValidElec
    EID = ElecOrder(thisElec);
    savedir = [savebasedir,slash,'elec',num2str(EID)];
    if exist(savedir,'dir')==7
        rmdir(savedir,'s');
    end
    mkdir(savedir);
%     mkdir([savedir,slash,'tmp']);
end

for thisTrial = 1:ValidTrial.NumValidTrials
    pos = ValidTrial.position(thisTrial);
    fseek(NsxID,pos,-1);
    
    Header = fread(NsxID,1,'int8=>int8');
    if Header ~= 1
        error('Header seek error, abort!');
    end
    Timestamp = fread(NsxID,1,'uint32=>double');
    
    if ValidTrial.Timestamp(thisTrial) ~= Timestamp;
        error('Timestamp mismatch!')
    end
    NumDP = fread(NsxID,1,'uint32=>double');
    if ValidTrial.NumDP(thisTrial) ~= NumDP;
        error('NumDP mismatch!')
    end
    
    DataLength = ValidTrial.NumDP(thisTrial);
    
    TrialData = fread(NsxID,[NumValidElec,DataLength],'int16=>double')';
    
    
    %         disp('Filtering the data')
    % remove line noise
    TrialData = filtfilt(d,TrialData);
    % lowpass to get the LFP
    TrialLFP = filtfilt(lpFilt,TrialData);
    % hipass to get the MUA
    rawTrialMUA = filtfilt(hpFilt,TrialData);
    % bandpass filter
    TrialData = filtfilt(bpFilt,TrialData);
    
    % down sampling LFP
    TrialLFP = TrialLFP(1:SampleRate/LFPSampleRate:size(TrialLFP,1),:);
    % half wave rectify in the negative direction
    rawTrialMUA = rawTrialMUA.*(rawTrialMUA<0);
    meanTrialMUA = mean(rawTrialMUA,1);
    stdTrialMUA = std(rawTrialMUA,0,1);
    % count the number of events more than 3 times std with in 1 ms time
    % bins
    MUAthres = meanTrialMUA - 3 * stdTrialMUA;
    edge = 1:(SampleRate/MUASampleRate):(size(rawTrialMUA,1)+(SampleRate/MUASampleRate));
    
    rawTrialMUA = rawTrialMUA < repmat(MUAthres,size(rawTrialMUA,1),1);
    
    TrialMUA = zeros(numel(edge),NumValidElec);
    for thisElec = 1: NumValidElec
        temp = find(rawTrialMUA(:,thisElec));
        TrialMUA(:,thisElec) = histc(temp,edge);
    end
    
    % left/right pad data
    TrialData = [zeros(waveformLength,size(TrialData,2)); TrialData];
    TrialData = [TrialData;zeros(waveformLength,size(TrialData,2))];
    
    
    % Estimate the threshold
    %         disp('Estimating threshold')
    C_Threshold = xEstimation * median(abs(TrialData)/0.6745,1); % 0.6745 is from experience discribed in iterature
    threshold(thisTrial,:) = C_Threshold;
    
    
    % spike detection
    %         disp('Detecting spikes')
    DataChunkLength = size(TrialData,1);
    C_poz = TrialData - repmat(C_Threshold,DataChunkLength,1);
    
    if xEstimation < 0
        C_poz = double(C_poz < 0);
    else
        C_poz = double(C_poz > 0);
    end
    C_poz_shifted = [zeros(1,size(C_poz,2)); C_poz(1:end-1,:)];
    C_poz_diff = C_poz - C_poz_shifted;
    
    % get spikes and waveforms
    %         disp('Preparing waveforms')
    
    TrialSpikeSTAligned = cell(NumValidElec,1);
    TrialSpikeTimeAligned = cell(NumValidElec,1);
    TrialWaveform = cell(NumValidElec,1);
    TrialNumSpikes = zeros(thisElec,1);
    
    TrialAlign = ExpmarkST(6,thisTrial) - ValidTrial.Timestamp(thisTrial);
    for thisElec = 1:NumValidElec
        % Spikestamps
        ElecTrialSpikes = find(C_poz_diff(:,thisElec) == 1);
        % exclude any spikes outside the range of the chunk
        ElecTrialSpikes = ElecTrialSpikes - waveformLength;
        
        ElecDataChunk = TrialData(:,thisElec);
        % get waveforms
        C_wavestart = ElecTrialSpikes - round(waveformLength * waveformAlignpoint) + waveformLength; % waveformLength is the length of padded data
        C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
        %         align to the minimal or maximal
        temp = reshape(ElecDataChunk(C_waveIndex),[],waveformLength);
        if xEstimation < 0
            [~,temp] = min(temp,[],2);
        else
            [~,temp] = max(temp,[],2);
        end
        C_wavestart = C_wavestart + temp - round(waveformLength * waveformAlignpoint);
        C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
        ElecTrialWaveform = reshape(ElecDataChunk(C_waveIndex),[],waveformLength);
        
        
        TrialSpikesAligned = ElecTrialSpikes - TrialAlign;
        TrialSpikeSTAligned{thisElec} = TrialSpikesAligned;
        TrialSpikeTimeAligned{thisElec} = TrialSpikesAligned/SampleRate;

        TrialWaveform{thisElec} = ElecTrialWaveform;
        TrialNumSpikes(thisElec) = numel(ElecTrialSpikes);
    end
    
    % align LFP and MUA
    LFPStart = round((TrialAlign - LFPTimeBeforeAlign * SampleRate) / (SampleRate/LFPSampleRate));
    LFPEnd = round((TrialAlign + LFPTimeAfterAlign * SampleRate)/(SampleRate/LFPSampleRate));
    TrialLFP = TrialLFP(LFPStart:LFPEnd,:);
    TrialLFPlength = size(TrialLFP,1);
    MUAStart = round((TrialAlign - MUATimeBeforeAlign * SampleRate) / (SampleRate/MUASampleRate));
    MUAEnd = round((TrialAlign + MUATimeAfterAlign * SampleRate)/(SampleRate/MUASampleRate));
    TrialMUA = TrialMUA(MUAStart:MUAEnd,:);
    TrialMUAlength = size(TrialMUA,1);
    
    % store spikes, waveform, LFP and MUA
    savedir = [savebasedir,slash,'tmp'];
    save([savedir,slash,'spikenwf',num2str(thisTrial),'.mat'],'TrialWaveform','TrialSpikeSTAligned','TrialSpikeTimeAligned','TrialNumSpikes')
    save([savedir,slash,'LFPMUA',num2str(thisTrial),'.mat'],'TrialLFP','TrialMUA','TrialLFPlength','TrialMUAlength')
end
% very important to close file
fclose(NsxID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split spikes into trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Sorting into trials and saving everything')


% loadData
clearvars -except WaveformCache SpikesCache ValidTrial savebasedir...
    ElecOrder ExpmarkST threshold slash NumValidElec SampleRate...
    LFPTimeBeforeAlign LFPTimeAfterAlign MUATimeBeforeAlign MUATimeAfterAlign...
    LFPSampleRate MUASampleRate TrialLFPlength TrialMUAlength

WaveformCache = cell(NumValidElec,1);
SpikesSTCache = cell(NumValidElec,1);
SpikesTimeCache = cell(NumValidElec,1);
NumSpikes = zeros(ValidTrial.NumValidTrials,NumValidElec);

% spike and wave form
savedir = [savebasedir,slash,'tmp'];
for thisTrial = 1:ValidTrial.NumValidTrials
    matpath = [savedir,slash,'spikenwf',num2str(thisTrial),'.mat'];
    load(matpath);
    for thisElec = 1:NumValidElec
        WaveformCache{thisElec} = [WaveformCache{thisElec};TrialWaveform{thisElec}];
        SpikesSTCache{thisElec} = [SpikesSTCache{thisElec};TrialSpikeSTAligned{thisElec}];
        SpikesTimeCache{thisElec} = [SpikesTimeCache{thisElec};TrialSpikeTimeAligned{thisElec}];
    end
    NumSpikes(thisTrial,:) = TrialNumSpikes;
end



for thisElec = 1:NumValidElec
    EID = ElecOrder(thisElec);
    SpikeSTAligned = mat2cell(SpikesSTCache{thisElec},NumSpikes(:,thisElec),1);
    SpikeTimeAligned = mat2cell(SpikesTimeCache{thisElec},NumSpikes(:,thisElec),1);
    elecRawWaveform = WaveformCache{thisElec};
    elecRawSpikeStamps = SpikesSTCache{thisElec};
    savedir = [savebasedir,slash,'elec',num2str(EID)];
    save([savedir,slash, 'spikeST6.mat'],'SpikeSTAligned');
    save([savedir,slash, 'SpikeTime6.mat'],'SpikeTimeAligned');
    save([savedir,slash, 'RawWaveform6.mat'],'elecRawWaveform');
    elecThreshold = threshold(:,thisElec);
    save([savedir,slash, 'Threshold.mat'],'elecThreshold');
    save([savedir,slash, 'RawSpikeStamps6.mat'],'elecRawSpikeStamps');
end

% % LFP and MUA
clearvars WaveformCache SpikesSTCache SpikesTimeCache
LFPCache = zeros(TrialLFPlength,ValidTrial.NumValidTrials,NumValidElec);
MUACache = zeros(TrialMUAlength,ValidTrial.NumValidTrials,NumValidElec);
savedir = [savebasedir,slash,'tmp'];
for thisTrial = 1:ValidTrial.NumValidTrials
    matpath = [savedir,slash,'LFPMUA',num2str(thisTrial),'.mat'];
    load(matpath);
    LFPCache(:,thisTrial,:) = TrialLFP;
    MUACache(:,thisTrial,:) = TrialMUA;
end

for thisElec = 1:NumValidElec
    EID = ElecOrder(thisElec);
    savedir = [savebasedir,slash,'elec',num2str(EID)];
    LFP = LFPCache(:,:,thisElec);
    MUA = MUACache(:,:,thisElec);
    save([savedir,slash, 'LFP6.mat'],'LFP');
    save([savedir,slash, 'MUA6.mat'],'MUA');
end
tmpdir = [savebasedir,slash,'tmp'];
if exist(tmpdir,'dir') == 7
    rmdir(tmpdir,'s');
end