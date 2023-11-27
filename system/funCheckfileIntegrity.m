function funCheckfileIntegrity(path1,path2)

% compare the file content in two directories

% path1 = 'Q:\ATLASprobeRF';
% path2 = 'T:\ATLASprobeRF';

% get full list of files recursively
filelist1 = dir(fullfile([path1,'\**\*.*']));
filelist1 = filelist1(~[filelist1.isdir]);
filelist2 = dir(fullfile([path2,'\**\*.*']));
filelist2 = filelist2(~[filelist2.isdir]);


% check file names
numfiles1 = numel(filelist1);
numfiles2 = numel(filelist2);

if numfiles1~=numfiles2
    disp('Warning!!! The number of files in path1 is NOT equal to the number of files in path2.')
end

lengthPath1 = numel(path1);
lengthPath2 = numel(path2);

for thisFile = 1:numfiles1
    filepath1{thisFile} = [filelist1(thisFile).folder((lengthPath1+1):end),'/',filelist1(thisFile).name];
end

for thisFile = 1:numfiles2
    filepath2{thisFile} = [filelist2(thisFile).folder((lengthPath2+1):end),'/',filelist2(thisFile).name];
end

[commonpath,commonpathIdx1, commonpathIdx2] = intersect(filepath1,filepath2);

numCommonFiles = numel(commonpath);

if numCommonFiles ~= numfiles1 || numCommonFiles ~= numfiles2
    disp('Warning!!! File names are differet between path1 and path2.')
    disp('Files are unique on path1:')
    disp(filepath1(~commonpathIdx1))
    disp('Files are unique on path2:')
    disp(filepath2(~commonpathIdx2))
else % file names are equal, calculate and compare checksum
    for thisFile = 1:numfiles1
        filename = [filelist1(thisFile).folder,'/',filelist1(thisFile).name];
        disp(['Checksum working on file: ',filename])
        path1checksums{thisFile} = Simulink.getFileChecksum(filename);
        filename = [filelist2(thisFile).folder,'/',filelist2(thisFile).name];
        path2checksums{thisFile} = Simulink.getFileChecksum(filename);
    end

    for thisFile = 1:numfiles1
        checksum1 = path1checksums{thisFile};
        checksum2 = path2checksums{thisFile};
        if issame(checksum1,checksum2)
            fileflag(thisFile) = 1;
        else
            fileflag(thisFile) = 0;
        end
    end

    if sum(fileflag) == numfiles1
    disp('All files on path1 and path2 are the same.')
    else
%         [commonpath,commonpathIdx1, commonpathIdx2] = intersect(path1checksums,path2checksums);
        idx = fileflag == 0;
        disp('Warning!!! File checksums are differet between path1 and path2.')
        disp('Files are unique on path1:')
        disp(filepath1(idx))
        disp('Files are unique on path2:')
        disp(filepath2(idx))
    end
end

