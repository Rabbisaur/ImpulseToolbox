function [Status] = LoadSpike(SpikePath, ElecID, SaveWaveform)
%Load SPIKE and save mat for each electrode in hard disk
% Usage:
%   [Status] = LoadSpike(SpikePath, ElecID, SaveWaveform)
% Inputs:
%   SpikePath      full path of the NEV/PLX file
%   ElecID         eletrodes to be analyzed
%   SaveWaveform   save spike waveform
% Outputs:
%   SpikeBasic	   basic header information
%   SpikeExtend    extended header information
%   Elec*Spike     spike time
%   Elec*Unit      sorted unit of the spike, 0/12345 - unsorted/sorted
%   Elec*Waveform  spike waveform containing 48 samples
%   Status         -1/0/1  failed/user-aborted/succeed
% External Includes:
%   MbmInfo = LoadMbm(MbmPath)
% Copyright 2010-2011 Minggui Chen, Beijing Normal University
% Revision: 7.0 Date: 2011/12/31 23:59:59

%% Check SPIKE and MBM files
Status = -1;
SpikePath = [SpikePath(1:end-3) 'nev'];
disp(['Processing ' SpikePath '...']);
SpikePath(strfind(SpikePath, '\')) = '/';
if exist(SpikePath, 'file')~=2
    disp('FAILURE: SPIKE file does not exist.');
    return;
end
MbmPath = [SpikePath(1:end-3) 'mbm'];
ValidMbm = exist(MbmPath, 'file')==2;
if ~ValidMbm
    MbmReply = input('MBM file does not exist. Press Y to continue? ', 's');
    if ~strcmpi(MbmReply, 'y')
        Status = 0;
        disp('FAILURE: User aborted LoadSpike when no mbm.');
        return;
    end
end
switch nargin
    case 1
        ElecID = [];
        SaveWaveform = 0;
    case 2
        SaveWaveform = 0;
    case 3
    otherwise
        disp('FAILURE: LoadSpike accepts 1 to 3 input arguments.');
        return;
end
%% Read basic header
SpikeID = fopen(SpikePath, 'r');
fseek(SpikeID, 0, 'eof');
SpikeSize = ftell(SpikeID);
frewind(SpikeID);
SpikeBasic.FileType = sprintf('%s', fread(SpikeID, 8, '*char'));
VerNum = fread(SpikeID, 2, 'char');
SpikeBasic.FileVersion = ['Spec. ' num2str(VerNum(1)) '.' num2str(VerNum(2))];
WaveformBit = fread(SpikeID, 16, '*ubit1');
if WaveformBit(1)==1
    SpikeBasic.ExtraFlag = 'all 16-bit wav';
else
    SpikeBasic.ExtraFlag = 'mixed 16-bit wav';
end
SpikeBasic.HeaderSize = fread(SpikeID, 1, 'uint32');
SpikeBasic.DPSize = fread(SpikeID, 1, 'uint32');
SpikeBasic.ClockFs = fread(SpikeID, 1, 'uint32');
SpikeBasic.WaveformFs = fread(SpikeID, 1, 'uint32');
SpikeBasic.TimeOrigin = fread(SpikeID, [1 8], 'uint16');
SpikeBasic.SoftwareVersion = sprintf('%s', fread(SpikeID, 32, '*char'));
% SpikeBasic.Comment = sprintf('%s', fread(SpikeID, 256, '*char'));
fseek(SpikeID, 256, 'cof');
SpikeBasic.NumExtend = fread(SpikeID,1,'uint32');
CurrPos = ftell(SpikeID);
if CurrPos~=336
    disp('FAILURE: Error in reading basic header.');
    return;
end
if SpikeBasic.HeaderSize~=336+SpikeBasic.NumExtend*32
    disp('FAILURE: Error in the size of extended headers.');
    return;
end
NumDP = floor((SpikeSize-SpikeBasic.HeaderSize)/SpikeBasic.DPSize);
if SpikeSize~=SpikeBasic.HeaderSize+NumDP*SpikeBasic.DPSize
    disp('FAILURE: Error in the size of data packets.');
    return;
end
NumWavSamp = (SpikeBasic.DPSize-8)/2;
%% Read extended headers
SpikeExtend = cell(SpikeBasic.NumExtend, 1);
ExtendPacketID = fread(SpikeID, [8 SpikeBasic.NumExtend], '8*char=>char', 32-8);
for i = 1:SpikeBasic.NumExtend
    CurrStartpos = CurrPos+(i-1)*32+8;
    fseek(SpikeID, CurrStartpos, 'bof');
    SpikeExtend{i,1}.PacketID = sprintf('%s',ExtendPacketID(:,i));  
    switch SpikeExtend{i,1}.PacketID
        case 'ARRAYNME'  %array name
            SpikeExtend{i,1}.ArrayName = sprintf('%s', fread(SpikeID, 24, '*char')); 
        case 'ECOMMENT'  %extra comment
            SpikeExtend{i,1}.EComment = sprintf('%s', fread(SpikeID, 24, '*char'));
        case 'CCOMMENT'  %continued comment
            SpikeExtend{i,1}.CComment = sprintf('%s', fread(SpikeID, 24, '*char')); 
        case 'MAPFILE'  %set to MAPFILE+NULL
            SpikeExtend{i,1}.MapFile = sprintf('%s', fread(SpikeID, 24, '*char'));  
        case 'NEUEVWAV'  %neural event waveform   
            SpikeExtend{i,1}.ElecID = fread(SpikeID, 1, 'uint16');
            SpikeExtend{i,1}.Connector = char(64 + fread(SpikeID,1,'char'));
            SpikeExtend{i,1}.Oin = fread(SpikeID, 1, 'char');
            SpikeExtend{i,1}.Scale = fread(SpikeID, 1, 'uint16');  %nV per LSB step 
            SpikeExtend{i,1}.EnergyThre = fread(SpikeID, 1, 'uint16');  %0 = NONE
            SpikeExtend{i,1}.HighThre = fread(SpikeID, 1, 'int16');  %uV
            SpikeExtend{i,1}.LowThre = fread(SpikeID, 1, 'int16');  %uV
            SpikeExtend{i,1}.NumUnit = fread(SpikeID, 1, 'char');  %0 = no sorted unit
            SpikeExtend{i,1}.WaveformSize = max(1, fread(SpikeID, 1, 'char'));  %quatify the length of waveform?
            fseek(SpikeID, 10, 'cof');  %10 bytes reserved
        case 'NEUEVLBL'  %neural event label
            SpikeExtend{i,1}.ElecID = fread(SpikeID, 1, 'uint16');
            SpikeExtend{i,1}.NeuralEventLabel = sprintf('%s', fread(SpikeID, 16, '*char'));
            fseek(SpikeID, 6, 'cof');  %6 bytes reserved
        case 'NEUEVFLT'  %neural event filter
            SpikeExtend{i,1}.ElecID = fread(SpikeID, 1, 'uint16');
            SpikeExtend{i,1}.HighCutFreq = fread(SpikeID, 1, 'uint32');
            SpikeExtend{i,1}.HignCutOrder = fread(SpikeID, 1, 'uint32');  %0 = NONE
            SpikeExtend{i,1}.HighCutType = fread(SpikeID, 1, 'uint16');
            SpikeExtend{i,1}.LowCutFreq = fread(SpikeID, 1, 'uint32');
            SpikeExtend{i,1}.LowCutOrder = fread(SpikeID, 1, 'uint32');  %0 = NONE
            SpikeExtend{i,1}.LowCutType = fread(SpikeID, 1, 'uint16');
            fseek(SpikeID,2,'cof');  %2 bytes reserved
        case 'DIGLABEL'  %digital label
            SpikeExtend{i,1}.DigitalLabel = sprintf('%s', fread(SpikeID, 16, '*char'));
            SpikeExtend{i,1}.DigitalMode = fread(SpikeID, 1, 'uint8');  %o/1-serial/parallel
            fseek(SpikeID, 7, 'cof');  %7 bytes reserved         
        case 'NSASEXEV'  %NSAS experiment information channels
            SpikeExtend{i,1}.PacketFreq = fread(SpikeID, 1, 'uint16');  %0 = NONE
            DigitalInputConfig = fread(SpikeID, 8, '*ubit1');
            ConfigOption = {'Ignored'; 'Informative'};
            SpikeExtend{i,1}.DigitalInputConfig = ConfigOption{DigitalInputConfig(1)+1};
            for k = 1:5
                AnalogConfig = fread(SpikeID, 8, '*ubit1');
                if AnalogConfig(1)==1
                    SpikeExtend{i,1}.(cell2mat(['Analog' num2str(k) 'Config'])) = 'High2Low';
                else
                    SpikeExtend{i,1}.(cell2mat(['Analog' num2str(k) 'Config'])) = 'Low2High';
                end
                SpikeExtend{i,1}.(cell2mat(['Analog' num2str(k) 'EdgeVal'])) = fread(SpikeID, 1, 'int16');  %mV
            end
            fseek(SpikeID, 6, 'cof');  %6 bytes reserved;
        otherwise
            disp('FAILURE: No current PacketID was found.');
            return;
    end
end
CurrPos = ftell(SpikeID);
if CurrPos~=SpikeBasic.HeaderSize
    disp('FAILURE: Error in reading extended headers.');
    return;
end
%check validity of ElecID
for i = 1:numel(SpikeExtend)
    if isfield(SpikeExtend{i,1}, 'ElecID')
        ValidElecID(i) = SpikeExtend{i,1}.ElecID;
    end
end
ValidElecID= unique(ValidElecID);
for i = 1:numel(ElecID)
   if sum(ValidElecID==ElecID(i))==0
       disp(['Elec' num2str(ElecID(i)) ' was not found in extended headers.']);
   end
end
%% Read timestamp & packet ID
NumFreadDP = 10^5;  %number of data packets per fread
RawData = ones([10 NumDP], 'uint8');
for i = 1:ceil(NumDP/NumFreadDP)
    NumReadDP = (i-1)*NumFreadDP;
    CurrNumDP = min(NumFreadDP, NumDP-NumReadDP);
    RawData(:,NumReadDP+(1:CurrNumDP)) = fread(SpikeID, [10 CurrNumDP], '10*uint8=>uint8', SpikeBasic.DPSize-10);
end
RawTime(:,1) = (double(RawData(1,:))*2^0+double(RawData(2,:))*2^8+double(RawData(3,:))*2^16+double(RawData(4,:))*2^24)/SpikeBasic.ClockFs;
RawID(:,1) = uint16(RawData(5,:))*2^0 + uint16(RawData(6,:))*2^8;
RawUnit(:,1) = RawData(7,:);
RawDigInput = RawData(9:10,:);  %9/10 - lower/higher 8 bits
clear RawData ExtendPacketID;
%% Retrieve trial timeline
TTLElec = 129;
InsertDPPos = find(RawID==0);  %position in raw packets
NumStart = 0;
NumPauseOff = 0;
NumStimCnd = 0;
NumVsgTrig = 0;
NumAbort = 0;
NumPauseOn = 0;
NumStop = 0;
for i = 1:numel(InsertDPPos)
    InsertReason = dec2bin(RawUnit(InsertDPPos(i)), 8);
    if str2double(InsertReason(8))==1  %bit 0 set
        switch RawDigInput(1,InsertDPPos(i))  %lower 8 bits
            case 1  %CBEVENT_START
                NumStart = NumStart+1;
            case 2  %CBEVENT_STOP
                NumStop = NumStop+1;
            case 4  %CBEVENT_PAUSEOFF, trial beginning
                NumPauseOff = NumPauseOff+1;
                OffTime(NumPauseOff,1) = RawTime(InsertDPPos(i));
                OffPos(NumPauseOff,1) = InsertDPPos(i);
            case 255  %CBEVENT_STIMCND = 0xff
                NumStimCnd = NumStimCnd+1;
                StimCnd(NumPauseOff, 1) = int16(RawDigInput(2,InsertDPPos(i)));
            case 254  %CBEVENT_STIMCND = 0xfe
                NumStimCnd = NumStimCnd+1;
                StimCnd(NumPauseOff, 1) = int16(256+RawDigInput(2,InsertDPPos(i)));
            case 253  %CBEVENT_STIMCND = 0xfd
                NumStimCnd = NumStimCnd+1;
                StimCnd(NumPauseOff, 1) = int16(512+RawDigInput(2,InsertDPPos(i)));
            case 128  %CBEVENT_VSGTRIG = 0x80 (software TTL, NOT mbm m_fTrigT)
                NumVsgTrig = NumVsgTrig+1;
            case 5  %CBEVENT_ABORT
                NumAbort = NumAbort+1;
                AbortTime(NumPauseOff,1) = RawTime(InsertDPPos(i));
            case 3  %CBEVENT_PAUSEON, trial ending
                NumPauseOn = NumPauseOn+1;
                OnTime(NumPauseOn,1) = RawTime(InsertDPPos(i));
                OnPos(NumPauseOn,1) = InsertDPPos(i);
                if NumPauseOff>=1 && (exist('StimCnd', 'var')~=1 || numel(StimCnd)==NumPauseOff-1)
                    StimCnd(NumPauseOff, 1) = nan;
                end
            otherwise  %CBEVENT_UNDEFINED, should occured?
        end
    end
end
if ~((NumPauseOff==NumPauseOn-1 || NumPauseOff==NumPauseOn) && NumPauseOff>=NumAbort && NumPauseOff>=NumStimCnd)
    disp('FAILURE: Error in the number of PAUSEON/OFF/ABORT.');
%     return;
end
OnPos = OnPos(NumPauseOn-NumPauseOff+1:end );
OnTime = OnTime(NumPauseOn-NumPauseOff+1:end);
if sum(OffTime>OnTime)~=0
    disp('FAILURE: Error in the relative timing between PAUSEOFF and PAUSEON.');
    return;
end
TTLDPPos = find(RawID==TTLElec);
TTLTime = NaN(NumPauseOff, 1000, 'double');
for i = 1:numel(TTLDPPos)
    TrialNum = logical((OffPos<=TTLDPPos(i)).*(OnPos>=TTLDPPos(i)));
    if sum(TrialNum)~=0
        TTLTime(TrialNum,sum(~isnan(TTLTime(TrialNum,:)))+1) = RawTime(TTLDPPos(i));
    end
end
for i = 1:NumPauseOff
    RawTime(OffPos(i):OnPos(i)) = RawTime(OffPos(i):OnPos(i))-TTLTime(i,1);
end
%% Compare MBM and SPIKE
if exist('AbortTime', 'var')
    AbortTime = double(AbortTime);
    AbortTime(AbortTime==0) = NaN;
    AbortTime(end+1:NumPauseOff,1) = NaN;
else
    AbortTime = NaN(NumPauseOff, 1, 'double');
end
if ValidMbm
    MbmInfo = LoadMbm(MbmPath);
    SpikeValidTrial = isnan(AbortTime);
    MbmValidTrial = MbmInfo.RespCode>0;
    if isequal(size(SpikeValidTrial), size(MbmValidTrial)) && ~isequal(SpikeValidTrial, MbmValidTrial)
        AbortTime(SpikeValidTrial~=MbmValidTrial) = 0;  %autocorrect AbortTime  
    end
    SpikeValidTrial = isnan(AbortTime);
    MbmValidTrial = MbmInfo.RespCode>0;
%     MbmValidTrial = logical((MbmInfo.RespCode>0).*(MbmInfo.RespCode~=14));
    if numel(SpikeValidTrial)==numel(MbmValidTrial)  %autocorrection
        AbortTime(SpikeValidTrial~=MbmValidTrial) = 0;
    end
%     AbortTime([260 411]) = 0;
    SpikeValidTrial = isnan(AbortTime);
    if ~isequal(SpikeValidTrial, MbmValidTrial)
        for i = 1:max(numel(SpikeValidTrial), numel(MbmValidTrial))
            if i<=numel(SpikeValidTrial) && i<=numel(MbmValidTrial)
                if SpikeValidTrial(i)==MbmValidTrial(i)
                    continue;
                end
            end
            disp(['FAILURE: The ' num2str(i) 'th RESPCODE in SPIKE and MBM is not matched.']);  %raw trial number
            return;
        end
    end
    StimCnd(isnan(StimCnd)) = MbmInfo.StimID(isnan(StimCnd));
    if numel(StimCnd)==numel(MbmInfo.StimID)
        StimCnd = MbmInfo.StimID;  %autocorrect StimCnd
    end
%     StimCnd([1]) = MbmInfo.StimID([1]);
    if ~isequal(StimCnd, MbmInfo.StimID)
        for i = 1:max(numel(StimCnd), numel(MbmInfo.StimID))
            if i<=numel(StimCnd) && i<=numel(MbmInfo.StimID)
                if SpikeValidTrial(i)==0 || ( SpikeValidTrial(i)==1 && StimCnd(i)==MbmInfo.StimID(i) )
                    continue;
                end
            end
            disp(['FAILURE: The ' num2str(i) 'th STIMCND in SPIKE and MBM is not matched.']);  %raw trial number
            return;
        end
    end
end
%% Generate ExpMoniter
ExpMonitor.StartT = OffTime-TTLTime(:,1);
ExpMonitor.StimCnd = StimCnd;
ExpMonitor.TTLT = TTLTime;
ExpMonitor.AbortT = AbortTime-TTLTime(:,1);
ExpMonitor.AbortT(isnan(TTLTime(:,1))) = -999;
ExpMonitor.EndT = OnTime-TTLTime(:,1);
if ValidMbm
    MbmField = fieldnames(MbmInfo);
    for i = 1:numel(MbmField)
        ExpMonitor.(MbmField{i}) = MbmInfo.(MbmField{i});
    end
end
%% Create folders for storing mat files
MachineMarker = strfind(lower(SpikePath), '/machinedata/');
if isempty(MachineMarker)
    MatDir = [SpikePath(1:end-4) '/'];
else
    MatDir = [SpikePath(1:MachineMarker) 'MatData/' SpikePath(MachineMarker+13:end-4) '/'];
end
if exist(MatDir,'dir')==7
%     rmdir(MatDir,'s');
%     mkdir(MatDir);
else
    mkdir(MatDir);
end
fileattrib(MatDir, '+w');
save([MatDir 'SpikeBasic.mat'], 'SpikeBasic', '-v7');
save([MatDir 'SpikeExtend.mat'], 'SpikeExtend', '-v7');
save([MatDir 'ExpMonitor.mat'], 'ExpMonitor', '-v7');
%% Orgnize spikes by electrode and trial
clearvars -except SpikeBasic NumWavSamp ElecID SaveWaveform NumOff TTLTime OffPos OnPos MatDir SpikeSize SpikeID CurrPos NumDP RawTime RawID RawUnit;
RawTime = single(RawTime);
fseek(SpikeID, CurrPos, 'bof');
if SaveWaveform
    [~, SystemView] = GetMemory;
    MaxArraySize = 0.95*SystemView.PhysicalMemory.Available;
    if SpikeSize<=MaxArraySize
        RawWaveform = fread(SpikeID, [SpikeBasic.DPSize/2 NumDP], '*int16');
    end
end
RawID([1:OffPos(1)-1 OnPos(end)+1:end]) = 0;
TrialMarker = [OffPos; OnPos(end)];
TrigTrialNum = ~isnan(TTLTime(:,1));
for i = 1:numel(ElecID)
    ElecDPPos = find(RawID==ElecID(i));  %position in RawID(TrialDPPos)
    if isempty(ElecDPPos)
        disp(['Elec' num2str(ElecID(i)) ' did NOT show any spikes.']);
        continue;
    end
    NumTrialDP = histc(ElecDPPos, TrialMarker);  %NumTrial+1
    Spike = mat2cell(RawTime(ElecDPPos), NumTrialDP(1:end-1));
    Unit = mat2cell(RawUnit(ElecDPPos), NumTrialDP(1:end-1));
    for j = 1:numel(TrigTrialNum)
        if ~TrigTrialNum(j)
            Spike{j,1} = single.empty([0 1]);
            Unit{j,1} = uint8.empty([0 1]);
        end
    end
    save([MatDir 'Elec' num2str(ElecID(i)) 'Spike.mat'], 'Spike', '-v7');
    save([MatDir 'Elec' num2str(ElecID(i)) 'Unit.mat'], 'Unit', '-v7');
    if SaveWaveform && SpikeSize<=MaxArraySize
        Waveform = mat2cell(transpose(RawWaveform(5:SpikeBasic.DPSize/2, ElecDPPos)), NumTrialDP(1:end-1), NumWavSamp);
        for j = 1:numel(TrigTrialNum)
            if ~TrigTrialNum(j)
                Waveform{j,1} = int16.empty([0 NumWavSamp]);
            end
        end
        save([MatDir 'Elec' num2str(ElecID(i)) 'Waveform.mat'], 'Waveform', '-v7');
    end
end
%% Read waveforms when greater SPIKE
if SaveWaveform && SpikeSize>MaxArraySize
    clear RawTime RawUnit Spike Unit;
    fseek(SpikeID, CurrPos, 'bof');
    [~, SystemView] = GetMemory;
    MaxArraySize = 0.95*SystemView.PhysicalMemory.Available;
    NumFreadDP = floor(MaxArraySize/SpikeBasic.DPSize);
    for i = 1:ceil(NumDP/NumFreadDP)
        VarName = ['Waveform' num2str(i)];
        NumReadDP = (i-1)*NumFreadDP;
        CurrNumDP = min(NumFreadDP, NumDP-NumReadDP);
        RawWaveform = fread(SpikeID, [SpikeBasic.DPSize/2, CurrNumDP], '*int16');
        for j = 1:numel(ElecID)  %orgnize waveforms by electrode
            if exist([MatDir 'Elec' num2str(ElecID(j)) 'Spike.mat'], 'file')==2
                ElecDPPos = find(RawID==ElecID(j))-NumReadDP;  %position in RawWaveform
                CurrStartPos = find(ElecDPPos>=1, 1, 'first');  %position in ElecDPPos
                CurrEndPos = find(ElecDPPos<=CurrNumDP, 1, 'last');  %position in ElecDPPos
                eval(sprintf('%s = transpose(RawWaveform(5:SpikeBasic.DPSize/2, ElecDPPos(CurrStartPos:CurrEndPos)));', VarName));
                MatName = [MatDir 'Elec' num2str(ElecID(j)) 'Waveform.mat'];
                if exist(MatName, 'file')~=2
                    save(MatName, VarName, '-v7');
                else
                    save(MatName, VarName, '-v7', '-append');
                end
            end
        end
        clear RawWaveform;
    end
    for i = 1:numel(ElecID)  %orgnize waveforms by trial        
        MatName = [MatDir 'Elec' num2str(ElecID(i)) 'Waveform.mat'];
        if exist(MatName, 'file')==2
            RawWaveform = load(MatName);
            ElecDPPos = find(RawID==ElecID(i));  %position in RawID
            Waveform = ones([numel(ElecDPPos) NumWavSamp], 'int16');
            WaveformField = fieldnames(RawWaveform);
            NumWave = 0;
            for j = 1:numel(WaveformField)
                CurrNumWave = size(RawWaveform.(WaveformField{j}), 1);
                NumWave = NumWave+CurrNumWave;
                Waveform(NumWave+1+(-CurrNumWave:-1),:) = RawWaveform.(WaveformField{j});
            end
            NumTrialDP = histc(ElecDPPos, TrialMarker);  %NumTrial+1
            Waveform = mat2cell(Waveform, NumTrialDP(1:end-1), NumWavSamp);
%             Waveform(isnan(TTLTime(:,1)),1) = cell(sum(isnan(TTLTime(:,1))), 1);
            for j = 1:numel(TrigTrialNum)
                if ~TrigTrialNum(j)
                    Waveform{j,1} = int16.empty([0 NumWavSamp]);
                end
            end
            save(MatName, 'Waveform', '-v7');
        end
    end
end
%% Clean up
fclose(SpikeID);
Status = 1;
disp('Done.');

function [UserView SystemView] = GetMemory
% Mimic the behavior of function memory on windows system
% Tested on Ubuntu 11.04
% Adapted from Wang Feng@BNU, Jan.4, 2012
if ispc
    [UserView SystemView] = memory;
elseif isunix
    %get MATLAB memory usage values
    [~, PageSize] = unix('getconf PAGE_SIZE');
    Ratio = str2double(PageSize)/1024;
    [~, A] = unix('ps -a|grep MATLAB');
    MatlabPID = textscan(A, '%d', 1);
    Cmd = ['cat /proc/' num2str(MatlabPID{1}) '/statm'];
    [~, A] = unix(Cmd);
    B = textscan(A, '%d');
    MatlabMem = B{1}(2)*Ratio*1024;
    %get vmstat values
    Cmd = 'cat /proc/meminfo';
    [~, B] = unix(Cmd);
    A = textscan(B, '%s');
    TotalMemory = str2double(A{1}(2))*1024;
    FreeMemory = str2double(A{1}(5))*1024;
    TotalSwap = str2double(A{1}(41))*1024;
    FreeSwap = str2double(A{1}(44))*1024;
    %get outputs
    UserView.MaxPossibleArrayBytes = FreeMemory+FreeSwap;
    UserView.MemAvailableAllArrays = FreeMemory+FreeSwap;
    UserView.MemUsedMATLAB = MatlabMem;
    SystemView.VirtualAddressSpace.Available = FreeSwap;
    SystemView.VirtualAddressSpace.Total = TotalSwap;
    SystemView.SystemMemory.Available = FreeMemory;
    SystemView.PhysicalMemory.Available = FreeMemory;
    SystemView.PhysicalMemory.Total = TotalMemory;
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  