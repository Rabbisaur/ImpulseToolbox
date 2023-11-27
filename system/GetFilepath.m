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

if ispc
    PathDelimiter = '\';
    tmpfilepath = 'c:\GestFilepathTmp.txt';
    cmd = ['dir "',dirpath,PathDelimiter,searchpattern,'" /l /s /b /-c > ',tmpfilepath];
%     cmd = ['dir "',dirpath,PathDelimiter,searchpattern,'" /l /s /b /-c'];
    status = dos(cmd);
else
    PathDelimiter = '/';
    disp('Not implemented yet on your OS')
end

% read in the file list
fp = fopen(tmpfilepath,'r');
linenum = 0;
while(~feof(fp))
    linenum = linenum+1;
    thisLine = fgetl(fp);
    filefullpath.path{linenum} = thisLine;
end

fclose(fp);
delete(tmpfilepath)
end