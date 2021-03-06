function filefullpath = GetFilepath(dirpath,searchpattern)
% 
% Use: filefullpath = GetFilepath(dirpath,searchpattern)
% if searchpattern is empty, it will return recursive dirpath under
% search root
%
% GetFilepath searches for filename which contents the search pattern in the
% dirpath and all its subfolders and returns a structure array with 4
% fields: 
%           name - full-path-to-the-file
%           datenum - date of the file in numerical format to help sorting
%           date - date of the file in a human readable format
%           bytes - file size in bytes
%
% Author: Wang Feng (euphoria.wang@gmail.com)
% Date: Nov. 2, 2012
%

% get path delimiter
if isunix
    PathDelimiter = '/';
elseif ispc
    PathDelimiter = '\';
end
% check the avaliability of the dirpath
if (exist(dirpath,'dir')==7)
    Rdirpath = genpath(dirpath); % recursively generate dir paths under root path
    % get individual directory into a cell array
    if isunix
        % path delimiter is : under unix
        cRdirpath = textscan(Rdirpath,'%s','Delimiter',':');
    elseif ispc
        % path delimiter is ; under pc dos
        cRdirpath = textscan(Rdirpath,'%s','Delimiter',';');
    end
    cRdirpath = cRdirpath{1};
    % get wanted file paths under all sub folders
    if exist('searchpattern','var')
        counter = 0; % counter for non-empty folders
        for diridx = 1:numel(cRdirpath)
            tmp_dirpath = cRdirpath{diridx};
            ifilepath = dir([tmp_dirpath, PathDelimiter, searchpattern]);
            if ~isempty(ifilepath)
                counter = counter +1;
                NumFiles = numel(ifilepath);
                for fileidx = 1:NumFiles
                    filefullpath{counter}(fileidx).name = [tmp_dirpath PathDelimiter ifilepath(fileidx).name];
                    filefullpath{counter}(fileidx).datenum = ifilepath(fileidx).datenum;
                    filefullpath{counter}(fileidx).date = ifilepath(fileidx).date;
                    filefullpath{counter}(fileidx).bytes = ifilepath(fileidx).bytes;
                    filefullpath{counter}(fileidx).isdir = ifilepath(fileidx).isdir;
                end
            end
        end
        if counter == 0 % search pattern not found
            filefullpath = [];
        else
            filefullpath = cell2mat(filefullpath);
        end
    else
        for diridx = 1:numel(cRdirpath)-1
            filefullpath(diridx).name = cRdirpath{diridx+1};
        end
    end
else
    error('User input root path not exist!')
end
