function TSLloadNS6Elec(filepath,xEstimation)
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
tic;
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
LFPTimeAfterAlign = 1; % Sec
MUATimeBeforeAlign = 0.2; % Sec
MUATimeAfterAlign = 1; % Sec
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
TSLloadNev(nevpath);
load(expmarkpath);

%% load 30KHz continuous signal.
disp('Loading NS6 file')
% NS6 = openNSx(Ns6Path);
%
% SampleRate = double(NS6.MetaTags.SamplingFreq);
% ElecOrder = double(NS6.MetaTags.ChannelID);
% NumValidElec = numel(ElecOrder);
% Data = double(NS6.Data');

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

Header = fread(NsxID,1,'int8=>int8');
if Header ~= 1
    error('Header seek error, abort!');
end
Timestamp = fread(NsxID,1,'uint32=>double');
NumDP = fread(NsxID,1,'uint32=>double');
DataStartFilepos = ftell(NsxID);
% Data = fread(NsxID,[NumValidElec,NumDP],'int16=>double')';

% CurrPos = ftell(NsxID);
% if CurrPos ~= NsxSize
%     disp(CurrPos)
%     disp(NsxSize)
%     error('File not fully read!')
% end


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


DataLength = NumDP;
NumChunks = ceil(DataLength/(datachunksize*SampleRate));
% make it at least two chunks
if NumChunks < 2
    NumChunks = 2;
end
% ChunkSize = floor(DataLength / NumChunks);
ChunkSize = datachunksize * SampleRate;

overlapsize = overlapsize * SampleRate;

% initialize SpikeStamps and waveform
% RawSpikeStamps = cell(128,1);
% RawWaveform = cell(128,1);
% RawUnits = cell(128,1);
% for thisElec = 1:128
%     RawSpikeStamps{thisElec} = [];
%     RawWaveform{thisElec} = [];
%     RawUnits{thisElec} = [];
% end

for thisElec = 1:NsxBasic.NumElec
    fseek(NsxID,DataStartFilepos,-1);
    RawSpikeStamps = [];
    RawWaveform = [];
    % initialize LFP
    rawLFP = [];
    rawMUA = [];
    for thisChunk = 1:NumChunks
        switch thisChunk
            case 1
                startpoz = DataStartFilepos;
                NumDPtoRead = ChunkSize+ overlapsize;
            case NumChunks
                startpoz = DataStartFilepos + ((ChunkSize*(thisChunk-1)) - overlapsize) * NumValidElec * 2;
                NumDPtoRead = DataLength - (NumChunks-1) * ChunkSize + overlapsize;
            otherwise
                startpoz = DataStartFilepos + ((ChunkSize*(thisChunk-1)) - overlapsize) * NumValidElec * 2;
                NumDPtoRead = ChunkSize + 2*overlapsize;
        end
        fseek(NsxID,startpoz,-1);
        DataChunk = fread(NsxID,[NumValidElec,NumDPtoRead],'int16=>double')';
        
        DataChunk = DataChunk(:,thisElec);
        
        disp('Filtering the data')
        % remove line noise
        DataChunk = filtfilt(d,DataChunk);
        % lowpass to get the LFP
        LFPChunk = filtfilt(lpFilt,DataChunk);
        % hipass to get the MUA
        MUAChunk = filtfilt(hpFilt,DataChunk);
        % bandpass filter
        DataChunk = filtfilt(bpFilt,DataChunk);
        
        % cut the LFP data
        switch thisChunk
            case 1
                idx = 1 : ChunkSize;
            case NumChunks
                idx = overlapsize+ 1 : size(DataChunk,1);
            otherwise
                idx = overlapsize+ 1 : overlapsize + ChunkSize;
        end
        LFPChunk = LFPChunk(idx);
        MUAChunk = MUAChunk(idx);
        % down sampling LFP
        rawLFP = [rawLFP;LFPChunk(1:SampleRate/LFPSampleRate:end)];
        % half wave rectify in the negative direction
        MUAChunk = MUAChunk.*(MUAChunk<0);
        meanMUAChunk = mean(MUAChunk);
        stdMUAChunk = std(MUAChunk);
        % count the number of events more than 3 times std with in 1 ms time
        % bins
        MUAthres = meanMUAChunk - 3 * stdMUAChunk;
        edge = 1:(SampleRate/MUASampleRate):(size(MUAChunk,1)+(SampleRate/MUASampleRate));
        if thisChunk ~= NumChunks && floor(size(MUAChunk,1) / (SampleRate/MUASampleRate)) ~= size(MUAChunk,1) / (SampleRate/MUASampleRate)
            error('Please use integar seconds in datachunksize')
        end
        NumMillisec = floor(size(MUAChunk,1) / (SampleRate/MUASampleRate));
        
        MUAChunk = MUAChunk < MUAthres;
        rawMUAElec = find(MUAChunk);
        rawMUAElec = histc(rawMUAElec,edge);
        rawMUAChunk = rawMUAElec(1:NumMillisec);
        
        rawMUA = [rawMUA;rawMUAChunk];
        
        % cut the spike data
        switch thisChunk
            case 1
                idx = 1 : ChunkSize + waveformLength;
            case NumChunks
                idx = overlapsize - waveformLength + 1 : size(DataChunk,1);
            otherwise
                idx = overlapsize - waveformLength + 1 : overlapsize + ChunkSize + waveformLength;
        end
        
        DataChunk = DataChunk(idx,:);
        
        % left pad data for the first chunk and right pad data for the last
        if thisChunk == 1
            DataChunk = [zeros(waveformLength,size(DataChunk,2)); DataChunk];
        end
        if thisChunk == NumChunks
            DataChunk = [DataChunk;zeros(waveformLength,size(DataChunk,2))];
        end
        
        
        % Estimate the threshold
        disp('Estimating threshold')
        C_Threshold = xEstimation * median(abs(DataChunk)/0.6745,1); % 0.6745 is from experience discribed in iterature
        threshold(thisChunk) = C_Threshold;
        
        
        % spike detection
        disp('Detecting spikes')
        DataChunkLength = size(DataChunk,1);
        C_poz = DataChunk - repmat(C_Threshold,DataChunkLength,1);
        
        if xEstimation < 0
            C_poz = double(C_poz < 0);
        else
            C_poz = double(C_poz > 0);
        end
        C_poz_shifted = [zeros(1,size(C_poz,2)); C_poz(1:end-1,:)];
        C_poz_diff = C_poz - C_poz_shifted;
        
        % get spikes and waveforms
        disp('Preparing waveforms')
        
        % Spikestamps
        Spikes = find(C_poz_diff == 1);
        % exclude any spikes outside the range of the chunk
        Spikes(Spikes < waveformLength + 1 | Spikes > ChunkSize) = [];
        Spikes = Spikes - waveformLength;
        
        % get waveforms
        C_wavestart = Spikes - round(waveformLength * waveformAlignpoint) + waveformLength; % waveformLength is the length of padded data
        C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
        %         align to the minimal or maximal
        temp = reshape(DataChunk(C_waveIndex),[],waveformLength);
        if xEstimation < 0
            [~,temp] = min(temp,[],2);
        else
            [~,temp] = max(temp,[],2);
        end
        C_wavestart = C_wavestart + temp - round(waveformLength * waveformAlignpoint);
        C_waveIndex = repmat(C_wavestart,1,waveformLength) + repmat(0:(waveformLength-1),numel(C_wavestart),1);
        tempWF = reshape(DataChunk(C_waveIndex),[],waveformLength);
        
        % store spikes and waveforms
        RawSpikeStamps = [RawSpikeStamps; Spikes + Timestamp + datachunksize* SampleRate * (thisChunk-1)];
        RawWaveform = [RawWaveform;tempWF];
        
    end
    
    CurrPos = ftell(NsxID);
    if CurrPos ~= NsxSize
        disp(CurrPos)
        disp(NsxSize)
        error('File not fully read!')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Split spikes into trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NumTrials = size(ExpmarkST,2);
    % Expmark(1,:) is the start time stamp of a trial
    % Expmark(4,:) is the end time stamp of a trial
    % Expmark(6,:) is where all the spikes aligned to
    disp('Sorting into trials')
    EID = ElecOrder(thisElec);
    SpikeSTAligned = cell(NumTrials,1);
    SpikeTimeAligned = cell(NumTrials,1);
    
    cumidx = zeros(size(RawSpikeStamps));
    NumSpikes = zeros(NumTrials,1);
    for thisTrial = 1:NumTrials
        TrialStart = ExpmarkST(1,thisTrial);
        TrialEnd = ExpmarkST(4,thisTrial);
        TrialAlign = ExpmarkST(6,thisTrial); % aligned to the stimulus onset
        LFPStart = round((TrialAlign - LFPTimeBeforeAlign * SampleRate) / (SampleRate/LFPSampleRate));
        LFPEnd = round((TrialAlign + LFPTimeAfterAlign * SampleRate)/(SampleRate/LFPSampleRate));
        LFP(:,thisTrial) = rawLFP(LFPStart:LFPEnd);
        MUAStart = round((TrialAlign - MUATimeBeforeAlign * SampleRate) / (SampleRate/MUASampleRate));
        MUAEnd = round((TrialAlign + MUATimeAfterAlign * SampleRate)/(SampleRate/MUASampleRate));
        MUA(:,thisTrial) = rawMUA(MUAStart:MUAEnd);
        idx = RawSpikeStamps >=  TrialStart & RawSpikeStamps < TrialEnd;
        TrialSpikes = RawSpikeStamps(idx);
        TrialSpikesAligned = TrialSpikes - TrialAlign;
        if size(TrialSpikesAligned,1) < size(TrialSpikesAligned,2)
            TrialSpikesAligned = TrialSpikesAligned';
        end
        if thisTrial > 108
            a = 1;
        end
        SpikeSTAligned{thisTrial} = TrialSpikesAligned;
        SpikeTimeAligned{thisTrial} = TrialSpikesAligned/SampleRate;
        NumSpikes(thisTrial) = numel(TrialSpikes);
        cumidx = cumidx | idx;
    end
    RawWaveform = RawWaveform(cumidx,:); % exclude spikes that not in a trial
    RawSpikeStamps = RawSpikeStamps(cumidx,:);
    % perform spikesorting
    RawSpikeTime = RawSpikeStamps/ SampleRate;
    %
    % automatical spike sorting with SpikeCluster
    idx = find(filepath == '.', 1, 'last');
    tmpdir = [filepath(1:(idx-1)),slash, 'tmp'];
    if exist(tmpdir,'dir')~=7
        mkdir(tmpdir)
    end
    % save a data file with only one electrode and only one trial in
    % Minggui Chen's format and then call Spike Cluster as a function
    % save spikes waveform ExpMonitor and SpikeBasic as filed by
    % Minggui Chen
    Spike{1} = RawSpikeTime;
    Waveform{1} = RawWaveform;
    ExpMonitor.StartT = Timestamp;
    ExpMonitor.EndT = max(Spike{1})+0.1;
    ExpMonitor.LFPStartT = ExpMonitor.StartT;
    ExpMonitor.LFPEndT = max(Spike{1})+0.1;
    SpikeBasic.WaveformFs = SampleRate;
    
    save([tmpdir,slash,'Elec',num2str(EID),'Spike.mat'],'Spike');
    save([tmpdir,slash,'Elec',num2str(EID),'Waveform.mat'],'Waveform');
    save([tmpdir,slash,'ExpMonitor.mat'],'ExpMonitor');
    save([tmpdir,slash,'SpikeBasic.mat'],'SpikeBasic');
    % direct call SpikeCluster as a function
    disp('Automatically spike sorting')
    FilePath = {
        [tmpdir,slash,'Elec',num2str(EID),'Waveform.mat'];
        };
    ArgIn.FileType = 3;
    ArgIn.FilePath = FilePath;
    SpikeCluster(ArgIn);
    % load unit back
    load([tmpdir,slash,'Elec',num2str(EID),'Unit.mat'])
    RawUnits = double(Unit{1});
    % remove tmp dir
    cd(DataDir);
    if exist(tmpdir,'dir')==7
        rmdir(tmpdir,'s')
    end
    NumSpikesinTrial = zeros(NumTrials,1);
    for thisTrial = 1:NumTrials
        NumSpikesinTrial(thisTrial) = numel(SpikeTimeAligned{thisTrial});
    end
    
    SpikeUnits = mat2cell(Unit{1},NumSpikesinTrial,1);
    
    disp('Saving everything')
    
    idx = find(filepath == '.', 1, 'last');
    savebasedir = filepath(1:(idx-1));
    
    savedir = [savebasedir,slash,'elec',num2str(EID)];
    if exist(savedir,'dir')==7
        rmdir(savedir,'s');
    end
    mkdir(savedir);
    save([savedir,slash, 'LFP6.mat'],'LFP');
    save([savedir,slash, 'MUA6.mat'],'MUA');
    save([savedir,slash, 'SpikeUnits6.mat'],'SpikeUnits');
    % save([filepath(1:(idx-1)),slash, 'waveform6.mat'],'waveform',matversion);
    save([savedir,slash, 'spikeST6.mat'],'SpikeSTAligned');
    save([savedir,slash, 'SpikeTime6.mat'],'SpikeTimeAligned');
    save([savedir,slash, 'RawUnits6.mat'],'RawUnits');
    save([savedir,slash, 'RawWaveform6.mat'],'RawWaveform');
    save([savedir,slash, 'Threshold.mat'],'threshold');
end
% very important to close file
fclose(NsxID);