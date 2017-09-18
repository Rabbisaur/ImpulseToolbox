function TSLloadData(machinepath,matpath,reloadflag,spikesortingflag,withpause,threshold)
% This function reads all the raw data recorded with Blackrock Cerebus NSP
% in a given directory and save the output in .MAT file format in a
% destination directory. If the given destination directory is not exist, it will
% automatically create one.
% machinepath: the path to a directory that contains all the raw data that
% you wish to read
% matpath: the path to a directory where you want the processed data to be
% stored
% reloadflag: 0 do not load the session if the result is existed in the
% matpath, 1 discard any result in the matpath and load the data
% spikesortingflag: 0 do not perform automatic offline spike sorting, 1
% perform automatic offline spike sorting using Minggui Chen's method.
% withpause: 0 the data was not recorded with pauses, 1 the data was
% recorded with pauses.

% assuming data recorded with pause and threshold = -3 by default.
switch nargin
    case 4
        withpause = 1;
        threshold = -3;
    case 5
        threhold = -3;
end


if ispc
    slash = '\';
else
    slash = '/';
end

% create the matpath if it is not exist.
if exist(matpath,'dir')~=7
    mkdir(matpath)
end

NS6path = GetFilepath(machinepath,'*.ns6');

NumFiles = numel(NS6path);
for thisFile = 1:NumFiles
    loaddataflag = 0;
    disp(NS6path(thisFile))
    idx = find(NS6path(thisFile).name == '.',1,'last');
    idx2 = find(NS6path(thisFile).name == slash,1,'last');
    sessionname = NS6path(thisFile).name((idx2+1):(idx-1));
    destpath = [matpath,slash,sessionname];
    if reloadflag
        loaddataflag = 1;
    else
        if exist(destpath,'dir')
            disp('Destination exists, skipping.')
        else
            loaddataflag = 1;
        end
    end
    if loaddataflag
        if exist(destpath,'dir')
            disp('Destination exists, removing old files first.')
            rmdir(destpath,'s');
        end
        if withpause
            TSLloadNS6ElecNoSortingParPause2(NS6path(thisFile).name,-3);
        else
            TSLloadNS6ElecNoSortingPar(NS6path(thisFile).name,-3);
        end
        cd(matpath);
        if ispc
            cmd = ['move ' NS6path(thisFile).name(1:idx-1) ' ' matpath];
            dos(cmd);
        else
            cmd = ['mv ' NS6path(thisFile).name(1:idx-1) ' ' matpath];
            unix(cmd);
        end
    end
    if spikesortingflag
        matfilepath = [matpath,NS6path(thisFile).name(idx2:idx-1)];
        electrodespath = [matfilepath,'/electrodes.mat'];
        load(electrodespath)
        electrodes(electrodes==129) = [];
        numElec = numel(electrodes);
        for thisElec = 1:numElec
            EID = electrodes(thisElec);
            elecPath = [matfilepath,'/elec',num2str(EID)];
            elecUnitpath = [elecPath,'/unit.mat'];
            if exist(elecUnitpath,'file')
                continue
            else
                SpikeSortingMGmethod(elecPath);
            end
        end
    end
    
end