function [Status] = LoadLFP(LFPPath, ElecID)
%Load NSx and save mat for each electrode in hard disk
% Usage:
%   [Status] = LoadLFP(LFPPath, ElecID, Window)
%   Pad 2^15-1 when no trigger or non-enough recording duration
% Inputs:
%   LFPPath        full path of the NSx/PLX file
%   ElecID         eletrodes to be analyzed
% Outputs:
%   LFPBasic	   basic header information
%   LFPExtend      extended header information
%   Elec*LFP       trial by sample voltage amplitudes, in uv
%   Status         -1/1  failed/succeed
% Copyright 2010-2011 Minggui Chen, Beijing Normal University
% Revision: 9.0 Date: 2011/12/31 23:59:59

%% Check NSx file
Status = -1;
for i = 0:9
    LFPPath = [LFPPath(1:end-4) '.ns' num2str(i)];
    if exist(LFPPath, 'file')==2
        break;
    end
    if i==9
        disp('FAILURE: NSx file does exist.');
        return;
    end
end
disp(['Processing ' LFPPath '...']);
LFPPath(strfind(LFPPath, '\')) = '/';
switch nargin
    case 1
        ElecID = [];
    case 2
    otherwise
        disp('FAILURE: LoadLFP accepts 1 or 2 input arguments.');
        return;
end
%% Create folders for storing mat files
MachineMarker = strfind(lower(LFPPath), '/machinedata/');
if isempty(MachineMarker)
    MatDir = [LFPPath(1:end-4) '/'];
else
    MatDir = [LFPPath(1:MachineMarker) 'MatData/' LFPPath(MachineMarker+13:end-4) '/'];
end
if exist(MatDir, 'dir')~=7
    disp(['FAILURE: ' MatDir ' does NOT exist. Please LoadSpike first.']);
    return;
end
if exist([MatDir 'ExpMonitor.mat'], 'file')~=2
    disp('FAILURE: EXPMONITOR does NOT exist. Please LoadSpike first.');
    return;
end
FileList = dir(MatDir);
for i = 1:numel(FileList)
    if FileList(i).isdir==0 && (strcmpi(FileList(i).name, 'LFPBasic.mat') ||...
            strcmpi(FileList(i).name, 'LFPExtend.mat') ||...
            ~isempty(strfind(FileList(i).name, 'LFP.mat')))
        delete([MatDir FileList(i).name]);
    end
end
%% Read basic header
LFPID = fopen(LFPPath,'r');
fseek(LFPID, 0, 'eof');
LFPSize = ftell(LFPID);
frewind(LFPID);
LFPBasic.FileType = sprintf('%s', fread(LFPID, 8, '*char'));
VerNum = fread(LFPID, 2, 'char');
LFPBasic.FileVersion = ['Spec. ' num2str(VerNum(1)) '.' num2str(VerNum(2))];
LFPBasic.HeaderSize = fread(LFPID, 1, 'uint32');
LFPBasic.Label = sprintf('%s', fread(LFPID, 16, '*char'));
% LFPBasic.Comment = sprintf('%s', fread(LFPID, 256, '*char'));
fseek(LFPID, 256, 'cof');
LFPBasic.Period = fread(LFPID, 1, 'uint32');
LFPBasic.ClockFs = fread(LFPID, 1, 'uint32');
LFPBasic.Fs = LFPBasic.ClockFs/LFPBasic.Period;
LFPBasic.TimeOrigin = fread(LFPID, [1 8], 'uint16');
LFPBasic.NumElec = fread(LFPID, 1, 'uint32');
CurrPos = ftell(LFPID);
if CurrPos~=314
    disp('FAILURE: Error in reading basic header');
    return;
end
if CurrPos+LFPBasic.NumElec*66~=LFPBasic.HeaderSize
    disp('FAILURE: Error in the size of extended headers.');
    return;
end
%% Read extend headers
LFPExtend = cell(LFPBasic.NumElec, 1);
PacketID = fread(LFPID, [2 LFPBasic.NumElec], '2*char=>char', 66-2);
ElecOrder = ones(LFPBasic.NumElec, 1);
for i = 1:LFPBasic.NumElec
    CurrStartPos = CurrPos+(i-1)*66+2;  %electrode ID
    fseek(LFPID, CurrStartPos, 'bof');
    LFPExtend{i,1}.PacketID = sprintf('%s', PacketID(:,i));
    switch LFPExtend{i,1}.PacketID
        case 'CC'  %continuous channels
            ElecOrder(i,1) = fread(LFPID, 1, 'uint16');
            LFPExtend{i,1}.ElecID = ElecOrder(i,1);
            LFPExtend{i,1}.ElecLabel = sprintf('%s', fread(LFPID, 16, '*char'));
            LFPExtend{i,1}.Connector = char(64+fread(LFPID, 1, 'char'));
            LFPExtend{i,1}.Pin = fread(LFPID, 1, 'char');
            LFPExtend{i,1}.MinDigitalValue = fread(LFPID, 1, 'int16');
            LFPExtend{i,1}.MaxDigitalValue = fread(LFPID, 1, 'int16');
            LFPExtend{i,1}.MinAnalogValue = fread(LFPID, 1, 'int16');
            LFPExtend{i,1}.MaxAnalogValue = fread(LFPID, 1, 'int16');
            LFPExtend{i,1}.AnalogUnit = sprintf('%s', fread(LFPID, 16, '*char'));  %mv/uv
            LFPExtend{i,1}.HighCutFreq = fread(LFPID, 1, 'uint32');
            LFPExtend{i,1}.HighCutOrder = fread(LFPID, 1, 'uint32');  %0 = NONE
            LFPExtend{i,1}.HighCutType = fread(LFPID, 1, 'uint16');
            LFPExtend{i,1}.LowCutFreq = fread(LFPID, 1, 'uint32');
            LFPExtend{i,1}.LowCutOrder = fread(LFPID, 1, 'uint32');  %0 = NONE
            LFPExtend{i,1}.LowCutType = fread(LFPID, 1, 'uint16');
        otherwise
            disp('FAILURE: No current PacketID was found.');
            return;
    end
end
if numel(ElecOrder)~=numel(unique(ElecOrder))
    disp('Error in reading extended headers.');
    return;
end
ElecPos = zeros(size(ElecID));
ElecValidity = zeros(size(ElecID));
for i = 1:numel(ElecID)
   if ~isempty(find(ElecOrder==ElecID(i), 1))
       ElecPos(i) = find(ElecOrder==ElecID(i), 1);  %row in RawLFP
       ElecValidity(i) = 1;
   else
       disp(['Elec' num2str(ElecID(i)) ' was not found in extended headers.']);
   end
end
ElecID = ElecID(logical(ElecValidity));
ElecPos = ElecPos(logical(ElecValidity));
CurrPos = ftell(LFPID);
if CurrPos~=LFPBasic.HeaderSize
    disp('FAILURE: Error in reading extended headers');
    return;
end
%% Read timestamp, nunmber of samples and LFP samples
load([MatDir 'ExpMonitor.mat']);
NumTrial = numel(ExpMonitor.StartT);
TimeNumSample = nan(NumTrial, 2, 'double');
NumDP = 0;
NumReadSample = 0;
NumOmitSample = 0;
MissPauseOn = 0;  %set to 1 when missing first PauseOn
while 1
    if ftell(LFPID)==LFPSize  %end of file
        break;
    end
    NumDP = NumDP+1;
    Header = fread(LFPID, 1, '*uint8');
    if Header~=hex2dec('01')
        disp(['FAILURE: The ' num2str(NumDP) 'th header was not set to 0x01.']);
        return;
    end
    TimeNumSample(NumDP,:) = fread(LFPID, 2, '2*uint32');
    TrialLFP = fread(LFPID, [LFPBasic.NumElec TimeNumSample(NumDP,2)], '*int16');
    if size(TrialLFP, 2)~=TimeNumSample(NumDP,2)
        disp('FAILURE: Error in the size of current packet.');
        return;
    end
    if MissPauseOn && NumDP==1
        NumOmitSample = TimeNumSample(1,2);
        CurrStartT = TimeNumSample(1,1)/LFPBasic.ClockFs-ExpMonitor.TTLT(1,1);
        CurrEndT = CurrStartT+(TimeNumSample(1,2)-1)*LFPBasic.Period/LFPBasic.ClockFs;
        MinStartT = ceil((ExpMonitor.StartT(1,1)+ExpMonitor.TTLT(1,1))*LFPBasic.ClockFs/LFPBasic.Period)*LFPBasic.Period/LFPBasic.ClockFs-ExpMonitor.TTLT(1,1);
        MaxEndT = floor((ExpMonitor.EndT(1,1)+ExpMonitor.TTLT(1,1))*LFPBasic.ClockFs/LFPBasic.Period)*LFPBasic.Period/LFPBasic.ClockFs-ExpMonitor.TTLT(1,1);
        StartT = max(CurrStartT, MinStartT);
        EndT = min(CurrEndT, MaxEndT);
        TimeNumSample(1,1) = (StartT+ExpMonitor.TTLT(1,1))*LFPBasic.ClockFs;
        TimeNumSample(1,2) = round((EndT-StartT)*LFPBasic.ClockFs/LFPBasic.Period+1);
        StartID = round((StartT-CurrStartT)*LFPBasic.ClockFs/LFPBasic.Period+1);
        TrialLFP = TrialLFP(:,StartID:StartID+TimeNumSample(1,2)-1);
        NumOmitSample = NumOmitSample-TimeNumSample(1,2);
    end
    if NumDP==1 && exist('RawLFP', 'var')~=1
        if TimeNumSample(NumDP,1)<LFPBasic.Period && ~MissPauseOn
            NumDP = NumDP-1;
            NumSample = ceil((LFPSize-ftell(LFPID)-NumTrial*9)/(LFPBasic.NumElec*2));
            RawLFP = ones(NumSample, numel(ElecID), 'int16');
            continue;
        else
            NumSample = ceil((LFPSize-LFPBasic.HeaderSize-NumTrial*9)/(LFPBasic.NumElec*2))-NumOmitSample;
            RawLFP = ones(NumSample, numel(ElecID), 'int16');
        end
    end
    RawLFP(NumReadSample+(1:size(TrialLFP, 2)),:) = transpose(TrialLFP(ElecPos,:));
    NumReadSample = NumReadSample+size(TrialLFP, 2);
end
if NumDP~=NumTrial
    disp('FAILURE: Error in the number of data packets.');
    return;
end
% TimeNumSample(:,1) = round(TimeNumSample(:,1)/LFPBasic.ClockFs*LFPBasic.Fs)*LFPBasic.ClockFs;  %align to ClockFs
TimeNumSample(:,1) = TimeNumSample(:,1)/LFPBasic.ClockFs-ExpMonitor.TTLT(:,1);
ExpMonitor.LFPStartT = TimeNumSample(:,1);
ExpMonitor.LFPEndT = TimeNumSample(:,1)+(TimeNumSample(:,2)-1)/LFPBasic.Fs;
if sum((ExpMonitor.LFPStartT<ExpMonitor.StartT)+(ExpMonitor.LFPStartT>ExpMonitor.EndT))~=0
    disp('FAILURE: Error in relative timing between LFP and spikes.');
    return;
end
save([MatDir 'LFPBasic.mat'], 'LFPBasic', '-v7');
save([MatDir 'LFPExtend.mat'], 'LFPExtend', '-v7');
save([MatDir 'ExpMonitor.mat'], 'ExpMonitor', '-v7');
%% Orgnize RawLFP by electrode
TrigTrialNum = ~isnan(ExpMonitor.TTLT(:,1));
for i = 1:numel(ElecID)  %col in RawLFP
    LFP = mat2cell(RawLFP(:,i), TimeNumSample(:,2), 1);
    for j = 1:numel(TrigTrialNum)
        if ~TrigTrialNum(j)
            LFP{j,1} = int16.empty([0 1]);
        end
    end
    save([MatDir 'Elec' num2str(ElecID(i)) 'LFP.mat'], 'LFP', '-v7');
end
%% Clean up
fclose(LFPID);
Status = 1;
disp('Done.');