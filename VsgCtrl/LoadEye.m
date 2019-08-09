function [EyeData, EyeTrace] = LoadEye(filepath)
EyeHeaderSize = 64;
disp('LoadEye working on file:')
disp(filepath)

% Check File
filepath = [filepath(1:end-4) '.eye'];
% mbmfilepath = [filepath(1:end-4) '.mbm'];
if exist(filepath,'file')
    EyeID = fopen(filepath);
    fseek(EyeID,0,1);
    EyeFileSize = ftell(EyeID);
    fseek(EyeID,0,-1);
else
    error('Eyetrace file not found!')
end
% if exist(mbmfilepath,'file')
%     [MbmNumtrials, MbmInfo, MBMFieldNames, sMbmInfo] = LoadMbm(mbmfilepath);
% else
%     error('Mbm file not found!')
% end
% check header
if fread(EyeID,1,'uint32=>uint32') ~= 1720232635
    error('Eyetrace file found but invalid')
end

header.version = fread(EyeID,1,'single=>double');
% position pointer to eye trace data beginning
fseek(EyeID,EyeHeaderSize,-1);
idx = 0;
while ftell(EyeID) ~= EyeFileSize
    idx = idx + 1;
    EyeTrace(idx) = readsingleeyedata(EyeID);
end
% reshape EyeTrace data structure
EyeData.TrialNum = zeros(idx,1,'double');
EyeData.StimID = zeros(idx,1,'double');
EyeData.interms = zeros(idx,1,'double');
EyeData.xy = cell(idx,1);
EyeData.xydeg = cell(idx,1);
for i = 1:idx
    EyeData.TrialNum(i) = EyeTrace(i).trialID;
    EyeData.StimID(i) = EyeTrace(i).stimID;
    EyeData.xy{i} = EyeTrace(i).xy;
    EyeData.xydeg{i} = EyeTrace(i).xydeg;
    EyeData.interms(i) = 1000 / EyeTrace(i).adfreq;
    EyeData.xypairs(i) = EyeTrace(i).xypairs;
end

% % % save to mat file
% % matversion = '-v7';
% % matpath = [filepath(1:end-4) '/' 'EyeData.mat'];
% % save(matpath,'EyeData',matversion)
% save to file
[savepath, savedir] = GetSavePath(filepath,'EyeData.mat');
if exist(savedir,'dir') == 7
else
    mkdir(savedir)
end
save(savepath,'EyeData','EyeTrace')
% if ispc
%     cmd = ['copy ' filepath ' ' savedir];
%     dos(cmd);
% else
%     cmd = ['cp ' filepath ' ' savedir];
%     unix(cmd);
% end
fclose(EyeID);
end
function data = readsingleeyedata(EyeID)
EyeDataPackteSize = 64;
% avgtimewindow = 100;
% aligntowindow = 1;
if fread(EyeID,1,'uint32=>uint32') ~= 3148515430;
    data = 50;
    return
end
data.trialID = fread(EyeID,1,'uint16=>double');
data.stimID = fread(EyeID,1,'uint16=>double');
data.bitflag = fread(EyeID,1,'uint16=>double');
data.adfreq = fread(EyeID,1,'uint16=>double');
data.gain = fread(EyeID,1,'uint16=>double');
data.winx = fread(EyeID,1,'int16=>double');
data.winy = fread(EyeID,1,'int16=>double');
data.winw = fread(EyeID,1,'int16=>double');
data.winh = fread(EyeID,1,'int16=>double');
data.xypairs = fread(EyeID,1,'int16=>double');
data.fpx = fread(EyeID,1,'single=>double');
data.fpy = fread(EyeID,1,'single=>double');
data.kx = fread(EyeID,1,'single=>double');
data.bx = fread(EyeID,1,'single=>double');
data.ky = fread(EyeID,1,'single=>double');
data.by = fread(EyeID,1,'single=>double');
fseek(EyeID,EyeDataPackteSize - 48,0);
data.xy = fread(EyeID,[2,double(data.xypairs)],'2*int16=>int16');
if ~isempty(data.xy)
    data.xy(1,:) = data.xy(1,:);% - data.winx;%-mean(data.xy(1,1:floor(avgtimewindow/(1000/data.adfreq))),2);
    data.xy(2,:) = data.xy(2,:);% - data.winy;%-mean(data.xy(2,1:floor(avgtimewindow/(1000/data.adfreq))),2);
    xy = double(data.xy);
%     if aligntowindow
%         xy(1,:) = xy(1,:) - data.winx;
%         xy(2,:) = xy(2,:) - data.winy;
%     end
    
    
    xy(1,:) = (xy(1,:) - data.bx)/data.kx; % inverse y = kx + b
    xy(2,:) = (xy(2,:) - data.by)/data.ky;
    
%     if avgtimewindow > 0 && data.xypairs > avgtimewindow
%         xy(1,:) = xy(1,:) - mean(xy(1,1:floor(avgtimewindow/(1000/data.adfreq))));
%         xy(2,:) = xy(2,:) - mean(xy(2,1:floor(avgtimewindow/(1000/data.adfreq))));
%     end
%     if aligntowindow
%         xy(1,:) = xy(1,:) - (single(data.winx)-data.bx)/data.kx;
%         xy(2,:) = xy(2,:) - (single(data.winy)-data.by)/data.ky;
%     end
    data.xydeg = xy;
else
    data.xydeg = -1;
end
end