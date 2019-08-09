function instance = VnCloadcerebusdata(sessionpath)
fclose('all');
NumValidElec = 0;
NsxID = nan;
cacheDir = '';
% switch nargin
%     case 2
%         % do nothing
%     case 1
%         % cacheDir = 'c:\tmpcache';
%         cacheDir = '';
% end
% LFPsamplingrate = 1000; % 1kHz
% MUAeSamplingrate = 2000; % 2kHz


if ~isempty(cacheDir)
    if exist(cacheDir,'dir')==7
        rmdir(cacheDir,'s')
    end
    mkdir(cacheDir)
end

% load the cerebus machine data from one session and save to a matdata
% folder

if ispc
    slash = '\';
else
    slash = '/';
end
NumPhysCores = feature('numcores');

% argin
% sessionpath = 'D:\data\sampledata\barsweeping';

% look for all the data from all instances
filepath = GetFilepath(sessionpath,'*.nev');
NumofNEVfiles = numel(filepath);

% meory cache size = 500M;
memoryCacheSize = 500*10^6;
% read each file in the folder, then combine together

for thisfile = 1:NumofNEVfiles
    try
        % first read all the nev files without reading the wave forms.
        currentNevpath = filepath(thisfile).name;
        % create mat folder for the current file
        idx = find(currentNevpath == '.',1,'last');
        currentMatdir = currentNevpath(1:idx-1);
        if exist(currentMatdir,'dir') == 7
            rmdir(currentMatdir,'s');
        end
        mkdir(currentMatdir);
        nevdata = openNEV(currentNevpath, 'noread', 'report', '8bits', 'nosave',...
            'nomat', 'overwrite');
        Electrodes = double(nevdata.MetaTags.ChannelID);
        SampleRes = double(nevdata.MetaTags.SampleRes);
        Filename = nevdata.MetaTags.Filename;
        
        % store the digi data for trial splitting
        instance(thisfile).nev.digitimestamps = nevdata.Data.SerialDigitalIO.TimeStamp;
        instance(thisfile).nev.digidata = nevdata.Data.SerialDigitalIO.UnparsedData;
        
        % load raw data
        %% load 30KHz continuous signal.
        Ns6Path = [currentNevpath(1:end-3),'ns6'];
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
        instance(thisfile).samplerate = SampleRate;
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
        if numel(ElecOrder) ~= numel(unique(ElecOrder))
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
        % load raw data electrode by electrode and store in cache
        % create cache dir
        
        NsxID = fopen(Ns6Path);
        fseek(NsxID,NsxBasic.HeaderSize,-1);
        
        NumValidElec = numel(ElecOrder);
        fprintf('Got %d valid electrode(s).\n',NumValidElec);
        TrialNum = 0;
        while ftell(NsxID) < NsxSize
            TrialNum = TrialNum + 1;
            Trial.position(TrialNum) = ftell(NsxID);
            Header = fread(NsxID,1,'int8=>int8');
            if Header ~= 1
%                 error('Header seek error, abort!');
                % discard the extradata following the previous trial
                Trial.position(TrialNum) = [];
                TrialNum = TrialNum - 1;
                break
            end
            Timestamp = fread(NsxID,1,'uint32=>double');
            Trial.Timestamp(TrialNum) = Timestamp;
            NumDP = fread(NsxID,1,'uint32=>double');
            Trial.NumDP(TrialNum) = NumDP;
            fseek(NsxID,Trial.NumDP(TrialNum)*NumValidElec*2,0); % jump to the beginning of the next data package
        end
        disp(['Got ',num2str(TrialNum),' trial(s)...'])
        
        % build instance info stucture
        instance(thisfile).numElec = NumValidElec;
        instance(thisfile).ElecOrder = ElecOrder;
        instance(thisfile).numTrials = TrialNum;
        instance(thisfile).trialInfo = Trial;
        
        
        % make electrode cache files on disk
        disp('Making electrode cache files')
        elecfp = zeros(NumValidElec,1);
        for thisElec = 1:NumValidElec
            EID = ElecOrder(thisElec);
            if isempty(cacheDir)
                savedir = [currentMatdir,slash,'elec',num2str(EID)];
            else
                idx = find(currentMatdir == slash,1,'last');
                savedir = [cacheDir,slash,currentMatdir((idx+1):end),slash,'elec',num2str(EID)];
            end
            mkdir([savedir,slash,'tmp']);
            tmpfilepath = [savedir,slash,'tmp',slash,'elecrawdata.tmp'];
            elecfp(thisElec) = fopen(tmpfilepath,'w');
            instance(thisfile).electrodeCachePath{thisElec} = tmpfilepath;
        end
        
        % read all data in the raw form electrodes by datapoints. Single out
        % the data from each electrode and store the data in electrode cache.
        disp(['Working on file ', num2str(thisfile),'/',num2str(NumofNEVfiles)])
        
        % guess the correct trial, the trial with most data

        [~,trialIdx] = max(Trial.NumDP);
        
        for thisTrial = trialIdx
            disp(['  Working on trial ', num2str(thisTrial),'/',num2str(TrialNum)])
            % read a piece of data in the memory cache at a time
            ChunkSize = floor(memoryCacheSize/(NumValidElec*2));
            numChunks = ceil(Trial.NumDP(thisTrial) / ChunkSize);
            ChunkSize = repmat(ChunkSize,numChunks-1,1);
            lastChunkSize = Trial.NumDP(thisTrial)-sum(ChunkSize);
            ChunkSize = [ChunkSize;lastChunkSize];
        end
        % read data from this Chunck
        fseek(NsxID,Trial.position(thisTrial)+9,-1); % jump to the begining of the data package of this trial
        for thisChunk = 1:numChunks
            disp(['    Working on chunk ', num2str(thisChunk), '/', num2str(numChunks)])
            chunck = fread(NsxID,[NumValidElec,ChunkSize(thisChunk)],'int16=>int16')';
            for thisElectrode = 1:NumValidElec
                electrodeData = chunck(:,thisElectrode);
                % write the electrode data to disk cache
                fwrite(elecfp(thisElectrode),electrodeData,'int16');
            end
        end
    catch ME % do the clean ups
        close('all')
        rethrow(ME)
    end
    % close electrode cache files
    for thisElec = 1:NumValidElec % close all disk cache files
        fclose(elecfp(thisElec));
    end
    % close nsx file
    fclose(NsxID);
    clearvars Trial
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % design the filters
% % try to remove 50Hz from the signal
% d = designfilt('bandstopiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
%     'DesignMethod','butter','SampleRate',SampleRate);
%
% % band pass filter 300-6000Hz
% bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',750,'HalfPowerFrequency2',5000, ...
%     'DesignMethod','butter','SampleRate',SampleRate);
%
% % low pass filter < 150 Hz
% lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',150,...
%     'DesignMethod','butter','SampleRate',SampleRate);
%
% % low pass filter < 500 Hz
% MUAeFilt = designfilt('lowpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency',500,...
%     'DesignMethod','butter','SampleRate',SampleRate);
%
% % % High pass filter > 1000 Hz
% % hpFilt = designfilt('highpassiir','FilterOrder',4, ...
% %     'HalfPowerFrequency',1000,...
% %     'DesignMethod','butter','SampleRate',SampleRate);
%
%
% % for each instance read back disk cache and do filtering for MUA, LFP and
% % Spikes (optional right now)
% numInstances = numel(instance);
% NumLogiCores = NumPhysCores * 1;
% disp('Extracting neural activities...')
% for thisInstance = 1:numInstances
%     disp(['  Working on file ', num2str(thisInstance),'/',num2str(numInstances)])
%     numElec = instance(thisInstance).numElec;
%     numTrials = instance(thisInstance).numTrials;
%     NumElecGroups = ceil(numElec / NumLogiCores);
%     remainder = mod(numElec,NumLogiCores);
%     elecIdx = (1:NumElecGroups * NumLogiCores)';
%     devisions = repmat(NumLogiCores,NumElecGroups,1);
%     ElecGroups = mat2cell(elecIdx,devisions,1);
%     if remainder > 0
%         ElecGroups{end} = ElecGroups{end}(1:remainder);
%     end
%     % Open disk cache files for each electrode
%     cacheFP = zeros(numElec,1);
%     for thisElec = 1:numElec
%         cachefilepath = instance(thisInstance).electrodeCachePath{thisElec};
%         cacheFP(thisElec) = fopen(cachefilepath,'r');
%     end
%
%     for thisTrial = 1:numTrials
%         disp(['    Working on trial ', num2str(thisTrial),'/',num2str(numTrials)])
%         numDP = instance(thisInstance).trialInfo.NumDP(thisTrial);
%         % load data in ElecGroups for more efficient filtering
%         for thisElecGroup = 1:NumElecGroups
%             disp(['    Working on electrode group ', num2str(thisElecGroup),'/',num2str(NumElecGroups)])
%             numElecsInGroup = numel(ElecGroups{thisElecGroup});
%             tempData = zeros(numElecsInGroup,numDP);
%             for thisElecInGroup = 1:numElecsInGroup
%                 ElecIdx = ElecGroups{thisElecGroup}(thisElecInGroup);
%                 tempData(thisElecInGroup,:) = fread(cacheFP(ElecIdx),numDP,'int16=>double');
%             end
%             tempData = tempData';
%             % remove 50 Hz
%             tempData = filtfilt(d,tempData);
%             % low pass for LFP
%             LFP = filtfilt(lpFilt,tempData);
%             % down sample LFP
%             LFP = LFP(1:SampleRate/LFPsamplingrate:numDP,:);
%             % band pass for MUA/Spikes
%             BPdata = filtfilt(bpFilt,tempData);
%             % MUAe
%             % starting from BPdata, 1) fullwave rectification, 2) low pass
%             % filter 500Hz
%             MUAe = abs(BPdata);
%             MUAe = filtfilt(MUAeFilt,MUAe);
%             MUAe = MUAe(1:SampleRate/MUAeSamplingrate:numDP,:);
%         end
%     end
% end

