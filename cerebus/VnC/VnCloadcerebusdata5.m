function instance = VnCloadcerebusdata5(sessionpath,destrootpath,instanceRange,instance)
% if you want to save files in the same folder as the data, give destrootpath
% an empty value

fclose('all');

% load the cerebus machine data from one session and save to a matdata
% folder

if ispc
    slash = '\';
else
    slash = '/';
end

% argin

% look for all the data from all instances
filepath = GetFilepath(sessionpath,'*.nev');
NumofNEVfiles = numel(filepath);

for thisFile = 1:NumofNEVfiles
    % for each nev, we look for instanceN, we only read out instances
    % specified in the input argument
    idx = find(filepath(thisFile).name==slash,1,'last');
    idx2 = find(filepath(thisFile).name=='.',1,'last');
    filename = filepath(thisFile).name((idx+1):(idx2-1));
    idx = strfind(filename,'instance');
    InstanceNumber(thisFile) = str2double(filename((idx+8):end));
end

[~,FileIdx,~] = intersect(InstanceNumber,instanceRange);
NumFiles = numel(FileIdx);


% meory cache size = 2048M;
memoryCacheSize = 2048*10^6;
% read each file in the folder, then combine together

for thisfile = 1:NumFiles
    thisfileIdx = FileIdx(thisfile);
    
%     for thisfile = 6
    try
        % first read all the nev files without reading the wave forms.
        currentNevpath = filepath(thisfileIdx).name;
        % create mat folder for the current file
        idx1 = find(currentNevpath == slash,1,'last');
        idx2 = find(currentNevpath == '.',1,'last');
        currentNevName = currentNevpath((idx1+1):(idx2-1));
        
        if isempty(destrootpath)
            currentMatdir = currentNevpath(1:idx2-1);
        else
            currentMatdir = [destrootpath,slash,currentNevName];
        end
        
        
        if exist(currentMatdir,'dir') == 7
            %rmdir(currentMatdir,'s');
            % do nothing
        else
            mkdir(currentMatdir);
        end
        nevdata = openNEV(currentNevpath, 'noread', 'report', '8bits', 'nosave',...
            'nomat', 'overwrite');
        
        % store the digi data for trial splitting
        instance(thisfile).instanceID = InstanceNumber(thisfileIdx);
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
        NsxID = fopen(Ns6Path,'r');
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
        
        NsxID = fopen(Ns6Path,'r');
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
        instance(thisfile).trialInfo.position = Trial.position;
        instance(thisfile).trialInfo.Timestamp = Trial.Timestamp;
        instance(thisfile).trialInfo.NumDP = Trial.NumDP;
        
        
        % make electrode cache files on disk
        disp('Making electrode cache files')
        elecfp = zeros(NumValidElec,1);
        for thisElec = 1:NumValidElec
            EID = ElecOrder(thisElec);
%             if isempty(cacheDir)
                savedir = [currentMatdir,slash,'elec',num2str(EID)];
%             else
%                 idx = find(currentMatdir == slash,1,'last');
%                 savedir = [cacheDir,slash,currentMatdir((idx+1):end),slash,'elec',num2str(EID)];
%             end
            if ~exist([savedir,slash,'tmp'],'dir')
                mkdir([savedir,slash,'tmp']);
            end
            tmpfilepath = [savedir,slash,'tmp',slash,'elecrawdata.tmp'];
            elecfp(thisElec) = fopen(tmpfilepath,'W');
            relativetmpfilepath = ['../','elec',num2str(EID),'/','tmp'];
            instance(thisfile).relativeElectrodeCachePath{EID} = relativetmpfilepath;
            instance(thisfile).electrodeCachePath{EID} = tmpfilepath;
        end
        
        % read all data in the raw form electrodes by datapoints. Single out
        % the data from each electrode and store the data in electrode cache.
        disp(['Working on file ', num2str(thisfile),'/',num2str(NumofNEVfiles)])
        
        % guess the correct trial, the trial with most data

%         [~,trialIdx] = max(Trial.NumDP);
        CurrentTimestamp = 1;
        for thisTrial = 1:TrialNum
            disp(['  Working on trial ', num2str(thisTrial),'/',num2str(TrialNum)])
            % read a piece of data in the memory cache at a time
            ChunkSize = floor(memoryCacheSize/(NumValidElec*2));
            numChunks = ceil(Trial.NumDP(thisTrial) / ChunkSize);
            
            % check if this trial has data
            if thisTrial < TrialNum
                if Trial.Timestamp(thisTrial) + Trial.NumDP(thisTrial) > Trial.Timestamp(thisTrial + 1)
                    disp('This trial has more data than the time difference between this and next trial, must be wrong!')
                    continue
                end
            end
            

            % write in a zero padding to align the data to the NEV file
            NumPadding = Trial.Timestamp(thisTrial) - CurrentTimestamp;
            if NumPadding < 0
                warning('Num padding less than zero, check data.')
            end

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
            
            % perform a file seek to correct any possible error in previous
            % trials
            for thisElectrode = 1:NumValidElec
                fseek(elecfp(thisElectrode),Trial.Timestamp(thisTrial)*2,'bof'); % int16 takes 2 bits
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