clc, clear,close all
rawmachinedatadir = '/home/wangfeng/projects/SurroundSuppression/machinedata';
ctxrecordeddatadir = '/home/wangfeng/projects/ctxrecordedfiles';
destmachinedatadir = '../leo/surroundsuppression';

% sessionDBpath = '../../SurroundSuppression/db/sessionDB736699.7444.mat';
% load(sessionDBpath);
% get all nev files
nevrawpath = GetFilepath(rawmachinedatadir,'*.nev');

if ispc
    slash = '\';
else
    slash = '/';
end

% for each nev file, find the corresponding data files
% nev
% ns2
% ns6
% ccf
% par
% cnd
% tm
% itm
% set
% lut
% blk
NumNevfiles = numel(nevrawpath);
for thisSession = 1:NumNevfiles
    filename = nevrawpath(thisSession).name;
    
    idx = find(filename == slash,1,'last');
    sessionname = filename(idx+1:end-4);
    idx = find(filename == slash,3,'last');
    expmarkpath = [filename(1:idx(1)),'matdata',slash,sessionname,slash,'Expmark.mat'];
    load(expmarkpath)
    condID = ExpmarkST(2,:);
    unicondID = unique(condID);
    NumCond = numel(unicondID);
    NumRepeats = zeros(NumCond,1);
    for thisCond = 1:NumCond
        NumRepeats(thisCond) = sum(condID == unicondID(thisCond));
    end
    repeats = min(NumRepeats);
    filedir = filename(1:idx(end));
    % two letter monkey name, 8 digit date, last three digits are session
    % number in a day
    sessionnumber = str2double(sessionname(end-2:end));
    monkeyname = sessionname(1:2);
    switch monkeyname
        case 'LE'
            monkeyname = 'leo';
    end
    sessiondate = sessionname(3:3+7);
    if strfind(sessionname,'_x_')
        strtype = 1.1;
        strtypestr = 'RF scan x';
        % find corresponding tm, itm, set, cnd, par file
        tmfilepath = [ctxrecordeddatadir,slash,'GRAT_H.TM'];
        itmfilepath = [ctxrecordeddatadir,slash,'ARR_GRAT.ITM'];
        setfilepath = [ctxrecordeddatadir,slash,'GRAT_H.SET'];
        cndfilepath = [ctxrecordeddatadir,slash,'GRAT_H.CND'];
        parfilepath = [ctxrecordeddatadir,slash,'GRAT_H.PAR'];
    end
    if strfind(sessionname,'_y_')
        strtype = 1.2;
        strtypestr = 'RF scan y';
        % find corresponding tm, itm, set file
        tmfilepath = [ctxrecordeddatadir,slash,'GRAT_V.TM'];
        itmfilepath = [ctxrecordeddatadir,slash,'ARR_GRAT.ITM'];
        setfilepath = [ctxrecordeddatadir,slash,'GRAT_V.SET'];
        cndfilepath = [ctxrecordeddatadir,slash,'GRAT_V.CND'];
        parfilepath = [ctxrecordeddatadir,slash,'GRAT_V.PAR'];
    end
    if strfind(sessionname,'_o_')
        sessiontype = 1.3;
        sessiontypestr = 'RF scan o';
        % find corresponding tm, itm, set file
        tmfilepath = [ctxrecordeddatadir,slash,'GRAT_O.TM'];
        itmfilepath = [ctxrecordeddatadir,slash,'ARR_GRAT.ITM'];
        setfilepath = [ctxrecordeddatadir,slash,'GRAT_O.SET'];
        cndfilepath = [ctxrecordeddatadir,slash,'GRAT_O.CND'];
        parfilepath = [ctxrecordeddatadir,slash,'GRAT_O.PAR'];
    end
    if strfind(sessionname,'_ssexp1_')
        sessiontype = 2;
        sessiontypestr = 'surround suppression experiment 1';
        % find corresponding tm, itm, set file
        tmfilepath = [ctxrecordeddatadir,slash,'SSEXP1.TM'];
        itmfilepath = [ctxrecordeddatadir,slash,'SSEXP1.ITM'];
        setfilepath = [ctxrecordeddatadir,slash,'SSEXP1.SET'];
        cndfilepath = [ctxrecordeddatadir,slash,'SSEXP1.CND'];
        parfilepath = [ctxrecordeddatadir,slash,'SSEXP1.PAR'];
    end
    if strfind(sessionname,'_chess_')
        sessiontype = 3;
        sessiontypestr = 'checker board';
        tmfilepath = [ctxrecordeddatadir,slash,'CHECKER.TM'];
        itmfilepath = [ctxrecordeddatadir,slash,'CHECKER.ITM'];
        setfilepath = [ctxrecordeddatadir,slash,'CHECKER.SET'];
        cndfilepath = [ctxrecordeddatadir,slash,'CHECKER.CND'];
        parfilepath = [ctxrecordeddatadir,slash,'CHECKER.PAR'];
    end
    % lut file
    lutfilepath = [ctxrecordeddatadir,slash,'gamma1.lut'];
    % find blk file
    % blk file has a name of 6 digit date and an extention of the session
    % number
    ctxfilepath = [ctxrecordeddatadir,slash,sessiondate(3:end),'.',num2str(sessionnumber)];
    % cerebus files
    % nev, ns2, ns6, ccf
    nevpath = [filedir,sessionname,'.nev'];
    ns2path = [filedir,sessionname,'.ns2'];
    ns6path = [filedir,sessionname,'.ns6'];
    ccfpath = [filedir,sessionname,'.ccf'];
    % check file availability
    if exist(tmfilepath,'file')~=2
        error([tmfilepath,' does not exist!'])
    end
    if exist(itmfilepath,'file')~=2
        error([itmfilepath,' does not exist!'])
    end
    if exist(setfilepath,'file')~=2
        error([setfilepath,' does not exist!'])
    end
    if exist(cndfilepath,'file')~=2
        error([cndfilepath,' does not exist!'])
    end
    if exist(parfilepath,'file')~=2
        error([parfilepath,' does not exist!'])
    end
    if exist(lutfilepath,'file')~=2
        error([lutfilepath,' does not exist!'])
    end
    if exist(ctxfilepath,'file')~=2
        error([ctxfilepath,' does not exist!'])
    end
    if exist(nevpath,'file')~=2
        error([nevpath,' does not exist!'])
    end
    if exist(ns2path,'file')~=2
        error([ns2path,' does not exist!'])
    end
    if exist(ns6path,'file')~=2
        error([ns6path,' does not exist!'])
    end
    if exist(ccfpath,'file')~=2
        error([ccfpath,' does not exist!'])
    end
    % make new folder structure
    % project name/ date/ 1,2,3 etc.
    if exist(destmachinedatadir,'dir')~=7
        mkdir(destmachinedatadir)
    end
    destdir = [destmachinedatadir,slash,num2str(sessiondate),slash,num2str(sessionnumber+100)];
    if exist(destdir,'dir')==7
        rmdir(destdir,'s')
    end
    mkdir(destdir)
    % make json file
    jsonfilepath = [destdir,slash,'note.json'];
    jsonfp = fopen(jsonfilepath,'w');
    fprintf(jsonfp,'{\n');
    fprintf(jsonfp,'"RF":{"x":"","y":""},\n');
    fprintf(jsonfp,'"blocks":"%d",\n',repeats);
    fprintf(jsonfp,'"notes":"%s",\n',sessiontypestr);
    fprintf(jsonfp,'"sessionname":"%s"\n',sessionname);
    
    fprintf(jsonfp,'}\n');
%     {
%     "RF":{"x":-0.4,"y":-1.4},
%     "blocks":20,
%     "data":"Le_2015_12_17_001.nev",
%     "notes":"100 conditions. F-D2N-D3N-NewNovel"
% }
    % 4 cerebus files, 7 cortex files
    sourcepath{1} = nevpath;
    sourcepath{2} = ns2path;
    sourcepath{3} = ns6path;
    sourcepath{4} = ccfpath;
    sourcepath{5} = tmfilepath;
    sourcepath{6} = itmfilepath;
    sourcepath{7} = setfilepath;
    sourcepath{8} = cndfilepath;
    sourcepath{9} = parfilepath;
    sourcepath{10} = lutfilepath;
    sourcepath{11} = ctxfilepath;
    
    
    destpath{1} = [destdir,slash,'LE_',sessiondate(1:4),'_',sessiondate(5:6),'_',sessiondate(7:8),'_1',sessionname(end-1:end),'.nev'];
    destpath{2} = [destdir,slash,'LE_',sessiondate(1:4),'_',sessiondate(5:6),'_',sessiondate(7:8),'_1',sessionname(end-1:end),'.ns2'];
    destpath{3} = [destdir,slash,'LE_',sessiondate(1:4),'_',sessiondate(5:6),'_',sessiondate(7:8),'_1',sessionname(end-1:end),'.ns6'];
    destpath{4} = [destdir,slash,'LE_',sessiondate(1:4),'_',sessiondate(5:6),'_',sessiondate(7:8),'_1',sessionname(end-1:end),'.ccf'];
    destpath{5} = [destdir,slash];
    destpath{6} = [destdir,slash];
    destpath{7} = [destdir,slash];
    destpath{8} = [destdir,slash];
    destpath{9} = [destdir,slash];
    destpath{10} = [destdir,slash];
    idx1 = find(ctxfilepath==slash,1,'last');
    idx2 = find(ctxfilepath=='.',1,'last');
    destpath{11} = [destdir,slash,ctxfilepath(idx1+1:idx2),num2str(sessionnumber+100)];
    % copy data files to the destination
    Numfiles = numel(sourcepath);
    for thisFile = 1:Numfiles
        if ispc
            cmd = ['copy ', sourcepath{thisFile},' ', destpath{thisFile}];
            dos(cmd);
        else
            cmd = ['cp ', sourcepath{thisFile},' ', destpath{thisFile}];
            unix(cmd);
        end
    end
    destpath{12} = [destdir,slash,'dummy.blk'];
    fp = fopen(destpath{12},'w');
    fclose(fp);
end