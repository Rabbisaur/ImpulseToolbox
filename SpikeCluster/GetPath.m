function [FilePath] = GetPath(SearchPath, Pattern)
%Find files named user-specified characters within subfolders
% Usage:
%   [FilePath] = GetPath(SearchPath, Pattern)
% Input:
%   SearchPath     full path of will-be-searching folder
%   Pattern        user-specified characters
% Output:
%   FilePath       paths of the files named user-specified characters
%                       return subfolders when Pattern is empty
% Copyright 2012-2012 Minggui Chen, Beijing Normal University
% Revision: 0.0 Date: 2012/5/15 23:59:59

%% Check SearchPath
SearchPath(SearchPath=='\') = '/';
if SearchPath(end)=='/'
    SearchPath = SearchPath(1:end-1);
end
if exist(SearchPath, 'dir')~=7
    error('FAILURE: SEARCHPATH does not exist.');
end
%% Search for files
FolderPath = genpath(SearchPath);  % Current folder & all subfolders
FolderPath(FolderPath=='\') = '/';
PathSeparator = pathsep;
NumFolder = 0;
FolderList = {};
while 1
    PathSepIdx = strfind(FolderPath, PathSeparator);
    if isempty(PathSepIdx)
        break;
    end
    NumFolder = NumFolder+1;
    FolderList{NumFolder,1} = [FolderPath(1:PathSepIdx(1)-1) '/'];
    FolderPath(1:PathSepIdx(1)) = [];
end
FolderPath = FolderList;
if isempty(Pattern)
    FilePath = FolderPath;
    return
end
NumFile = 0;
for i = 1:numel(FolderPath)	 % Scan folders
    FileList = dir(FolderPath{i});
    for j = 1:numel(FileList)
        if FileList(j).isdir
            continue;
        end
        if ~isempty(Pattern) && ~strcmpi(Pattern, '*')  % All files
            if isempty(strfind(FileList(j).name, Pattern))
                continue;
            end
        end
        NumFile = NumFile+1;
        FilePath{NumFile,1} = [FolderPath{i} FileList(j).name];
        FilePath{NumFile,1}(FilePath{NumFile,1}=='\') = '/';
    end
end