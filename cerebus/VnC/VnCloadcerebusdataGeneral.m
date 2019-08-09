function FileInfo = VnCloadcerebusdataGeneral(filepath)
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

% read each file in the folder, then combine together
% first read all the nev files without reading the wave forms.
% create mat folder for the current file
idx = find(filepath == '.',1,'last');
currentMatdir = filepath(1:idx-1);
mkdir(currentMatdir);
nevdata = openNEV(filepath, 'noread', 'report', '8bits', 'nosave',...
    'nomat', 'overwrite');

% store the digi data for trial splitting
FileInfo.nev.digitimestamps = nevdata.Data.SerialDigitalIO.TimeStamp;
FileInfo.nev.digidata = nevdata.Data.SerialDigitalIO.UnparsedData;

% load raw data
%% load 30KHz continuous signal.

ns6data = openNSx([filepath(1:idx),'ns6'],'read','report','precision','double');

% 
% % build instance info stucture
NumValidElec = double(ns6data.MetaTags.ChannelCount);
ElecOrder = double(ns6data.MetaTags.ChannelID);
FileInfo.numElec = NumValidElec;
FileInfo.ElecOrder = ElecOrder;
% FileInfo(thisfile).numTrials = TrialNum;
% FileInfo(thisfile).trialInfo = Trial;

% % make electrode cache files on disk
% disp('Making electrode cache files')
elecfp = zeros(NumValidElec,1);
for thisElec = 1:NumValidElec
    EID = ElecOrder(thisElec);
    savedir = [currentMatdir,slash,'elec',num2str(EID)];
    mkdir([savedir,slash,'tmp']);
    tmpfilepath = [savedir,slash,'tmp',slash,'elecrawdata.tmp'];
    elecfp(thisElec) = fopen(tmpfilepath,'W');
    FileInfo.electrodeCachePath{EID} = tmpfilepath;
    electrodeData = ns6data.Data(thisElec,:);
    fwrite(elecfp(thisElec),electrodeData,'int16');
    fclose(elecfp(thisElec));
end
