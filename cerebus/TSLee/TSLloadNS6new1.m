function instance = TSLloadNS6new(filepath)
% meory cache size = 4000M;
memoryCacheSize = 4000*10^6;

fclose('all');

% load the cerebus machine data from one session and save to a matdata
% folder

if ispc
    slash = '\';
else
    slash = '/';
end

% first read all the nev files without reading the wave forms.
currentNevpath = [filepath(1:(end-3)),'nev'];
if ~exist(currentNevpath,'file')
    error('Nev file does not exist.')
end
% create mat folder for the current file
idx = find(currentNevpath == '.',1,'last');
currentMatdir = currentNevpath(1:idx-1);
if exist(currentMatdir,'dir') == 7
    % do nothing
else
    mkdir(currentMatdir);
end
% nevdata = openNEV(currentNevpath, 'noread', 'report', '8bits', 'nosave',...
%     'nomat', 'overwrite');
% instance.nev.electrodes = double(nevdata.MetaTags.ChannelID);
% instance.nev.samplerate = double(nevdata.MetaTags.SampleRes);
% instance.nev.filename = nevdata.MetaTags.Filename;
% instance.nev.dateTime = nevdata.MetaTags.DateTime;
% instance.nev.dateTimeRaw = nevdata.MetaTags.DateTimeRaw;
% 
% % store the digi data for trial splitting
% instance.nev.digitimestamps = nevdata.Data.SerialDigitalIO.TimeStamp;
% instance.nev.digidata = nevdata.Data.SerialDigitalIO.UnparsedData;

% load raw data
%% load 30KHz continuous signal.
try
    Ns6Path = [currentNevpath(1:end-3),'ns6'];
    disp('Loading NS6 file')
    
    if ~exist(Ns6Path,'file')
        error('NS6 file does not exist.')
    end
    
    % process basic header
    NsxID = fopen(Ns6Path,'r');
    fseek(NsxID, 0, 1);
    NsxSize = uint64(ftell(NsxID));
    fprintf('NS6 file size is %dMB\n',NsxSize/1024^2)
    
    if NsxSize < memoryCacheSize
        data = openNSx(Ns6Path);

        % write to disk cache
        NumValidElec = numel(data.ElectrodesInfo);
        fprintf('Got %d valid electrode(s).\n',NumValidElec);
        ElecOrder = zeros(NumValidElec,1);
        for thisElec = 1:NumValidElec
            ElecOrder(thisElec) = data.ElectrodesInfo(thisElec).ElectrodeID;
        end
        % make electrode cache files on disk
        disp('Making electrode cache files')
        elecfp = zeros(NumValidElec,1);
        for thisElec = 1:NumValidElec
            EID = ElecOrder(thisElec);
            savedir = [currentMatdir,slash,'elec',num2str(EID)];
            tmpdirpath = [savedir,slash,'tmp'];
            if ~exist(tmpdirpath,'dir')
                mkdir(tmpdirpath);
            end
            tmpfilepath = [savedir,slash,'tmp',slash,'elecrawdata.tmp'];
            elecfp(thisElec) = fopen(tmpfilepath,'w');
            instance.electrodeCachePath{EID} = tmpfilepath;
        end
        
        NumTrials = numel(data.Data);
        
        CurrentTimestamp = 0;
        for thisTrial = 1:NumTrials
            Trial.Timestamp(thisTrial) = data.MetaTags.Timestamp(thisTrial);
            Trial.NumDP(thisTrial) = data.MetaTags.DataPoints(thisTrial);
            disp(['  Working on trial ', num2str(thisTrial),'/',num2str(NumTrials)])
            
            % write in a zero padding to align the data to the NEV file
            NumPadding = data.MetaTags.Timestamp(thisTrial) - CurrentTimestamp;
            
            CurrentTimestamp = data.MetaTags.Timestamp(thisTrial) + data.MetaTags.DataPoints(thisTrial);
            
            if NumPadding > 0
                padding = int16(zeros(NumPadding,1));
                % write padding
                for thisElectrode = 1:NumValidElec
                    % write the electrode data to disk cache
                    fwrite(elecfp(thisElectrode),padding,'int16');
                end
            end
            % read data from this Chunck
            for thisElectrode = 1:NumValidElec
                electrodeData = data.Data{thisTrial}(thisElectrode,:);
                % write the electrode data to disk cache
                fwrite(elecfp(thisElectrode),electrodeData,'int16');
            end
        end
        
        % build instance structure
        instance.samplerate = double(data.MetaTags.SamplingFreq);
        instance.numElec = NumValidElec;
        instance.ElecOrder = ElecOrder;
        instance.numTrials = NumTrials;
        instance.trialInfo = Trial;
    else
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
        instance.samplerate = SampleRate;
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
        
        NumValidElec = numel(ElecOrder);
        fprintf('Got %d valid electrode(s).\n',NumValidElec);
        
        % load raw data electrode by electrode and store in cache
        % create cache dir
        
        NsxID = fopen(Ns6Path);
        fseek(NsxID,NsxBasic.HeaderSize,-1);
        
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
        instance.numElec = NumValidElec;
        instance.ElecOrder = ElecOrder;
        instance.numTrials = TrialNum;
        instance.trialInfo = Trial;
        
        
        % make electrode cache files on disk
        disp('Making electrode cache files')
        elecfp = zeros(NumValidElec,1);
        for thisElec = 1:NumValidElec
            EID = ElecOrder(thisElec);
            savedir = [currentMatdir,slash,'elec',num2str(EID)];
            tmpdirpath = [savedir,slash,'tmp'];
            if ~exist(tmpdirpath,'dir')
                mkdir(tmpdirpath);
            end
            tmpfilepath = [savedir,slash,'tmp',slash,'elecrawdata.tmp'];
            elecfp(thisElec) = fopen(tmpfilepath,'w');
            instance.electrodeCachePath{EID} = tmpfilepath;
        end
        
        % read all data in the raw form electrodes by datapoints. Single out
        % the data from each electrode and store the data in electrode cache.
        
        %         [~,trialIdx] = max(Trial.NumDP);
        CurrentTimestamp = 0;
        for thisTrial = 1:TrialNum
            disp(['  Working on trial ', num2str(thisTrial),'/',num2str(TrialNum)])
            % read a piece of data in the memory cache at a time
            ChunkSize = floor(memoryCacheSize/(NumValidElec*2));
            numChunks = ceil(Trial.NumDP(thisTrial) / ChunkSize);
            
            % write in a zero padding to align the data to the NEV file
            NumPadding = Trial.Timestamp(thisTrial) - CurrentTimestamp;
            
            if numChunks == 0
                % this trial has no data, skip to the next
                disp('This trail has no data, skip to the next.')
                continue
            end
            
            CurrentTimestamp = Trial.Timestamp(thisTrial) + Trial.NumDP(thisTrial);
            
            ChunkSize = repmat(ChunkSize,numChunks-1,1);
            lastChunkSize = Trial.NumDP(thisTrial)-sum(ChunkSize);
            ChunkSize = [ChunkSize;lastChunkSize];
            
            if NumPadding > 0
                padding = int16(zeros(NumPadding,1));
                % write padding
                for thisElectrode = 1:NumValidElec
                    % write the electrode data to disk cache
                    fwrite(elecfp(thisElectrode),padding,'int16');
                end
            end
            % read data from this Chunck
            fseek(NsxID,Trial.position(thisTrial)+9,-1); % jump to the begining of the data package of this trial
            for thisChunk = 1:numChunks
                if numChunks > 1
                    disp(['    Working on chunk ', num2str(thisChunk), '/', num2str(numChunks)])
                end
                chunck = fread(NsxID,[NumValidElec,ChunkSize(thisChunk)],'int16=>int16')';
                for thisElectrode = 1:NumValidElec
                    electrodeData = chunck(:,thisElectrode);
                    % write the electrode data to disk cache
                    fwrite(elecfp(thisElectrode),electrodeData,'int16');
                end
            end
            
        end
        
        % close electrode cache files
        for thisElec = 1:NumValidElec % close all disk cache files
            fclose(elecfp(thisElec));
        end
        % close nsx file
        fclose(NsxID);
        clearvars Trial
    end
catch ME % do the clean ups
    close('all')
    rethrow(ME)
end

end