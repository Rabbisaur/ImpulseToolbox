function varargout = SpikeCluster(varargin)
% SpikeCluster designed by Minggui Chen @BNU 2/27/2012
% Imitate the GUI of the well-known package WAVE_CLUS
% Without permission, please do NOT distribute this program.

% Coding remained: 1) std estimation in GetThre()
%                  2) tetrode sorting & neuronal distance, programmable GUI
%                  3) overlapping spikes
%                  4) disable Handles.VoltThre
%                  5) hidden Markov model
%                  6) direct file I/O
%                  7) quality metrics, L Ratio etc.

% Initialization, DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpikeCluster_OpeningFcn, ...
                   'gui_OutputFcn',  @SpikeCluster_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% Initialization by input arguments
function SpikeCluster_OpeningFcn(~, ~, Handles, varargin)
% Display initialization info.
% DispInitInfo();
% Start main GUI
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
JavaFrame = get(Handles.SpikeClusterFig,'javaframe');
Handles.Track = randi(10, 1);
% Include files under SpikeCluster
SCPath = mfilename('fullpath');
SCPath(SCPath=='\') = '/';
idx = find(SCPath == '/',1,'last');
SCPath = SCPath(1:(idx-1));
% FolderIdx = strfind(SCPath, '/');
% SCPath = SCPath(1:FolderIdx(end));
% addpath(genpath(SCPath));
% cd(SCPath);
% Load default settings
% FigIcon = javax.swing.ImageIcon([SCPath '/Logo.jpg']);
% JavaFrame.setFigureIcon(FigIcon);
movegui(Handles.SpikeClusterFig, 'center');
ParaPath = [SCPath '/DefaultPara.mat'];
load(ParaPath);
Handles.Para = Para;
Handles.RawPara = Handles.Para;
Handles.Clr = 'kbrgcm';
Handles.ShadingClr = {[.5 .5 .5];[.5 .5 1];[1 .5 .5];[.5 1 .5];[.5 1 1];[1 .5 1];};
Handles.CD = cd;
Handles.IsAutoTemp = 0;
Handles.varargin = varargin;
Handles.JavaFrame = JavaFrame;
guidata(Handles.SpikeClusterFig, Handles);

% Update parameters form Parameters.fig
function ParametersPush_Callback(~, ~, Handles)
Handles.Para = Parameters(Handles.Para);
Handles.RawPara = Handles.Para;
guidata(Handles.SpikeClusterFig, Handles);

% Treat SpikeCluster as a function with input arguments
function SpikeClusterFcn(Handles)
if isempty(Handles.varargin)  % No file selected
    return;
end
Handles.IsAutoTemp = 0;
Handles.Para = Handles.RawPara;
ArgIn = Handles.varargin{1};
Handles.FileType = ArgIn.FileType;
Handles.FilePath = ArgIn.FilePath;
if Handles.FileType>=5
    Handles.ElecID = ArgIn.ElecID;
end
Handles = Port2Core(Handles);
guidata(Handles.SpikeClusterFig, Handles);
close(Handles.SpikeClusterFig);

% Open a file
function OpenPush_Callback(~, ~, Handles)
Handles.IsAutoTemp = 0;
Handles.Para = Handles.RawPara;
% Get file path
[Handles.FileName, Handles.PathName, Handles.FilterType] = uigetfile({...
                                 '*Waveform.mat; *LFP.mat; *WaveformBatch.m; *LFPBatch.m; *SpikeFileBatch.m; *LFPFileBatch.m;','All Supported Files';...
                                 '*Waveform.mat','Waveform File (*Waveform.mat)';...
                                 '*LFP.mat','LFP File (*LFP.mat)';...
                                 '*WaveformBatch.m','Waveform Files (*WaveformBatch.m)';...
                                 '*LFPBatch.m','LFP Files (*LFPBatch.m)';...
                                 '*SpikeFileBatch.m','Spike Files (*SpikeFileBatch.m)';...
                                 '*LFPFileBatch.m','LFP Files (*LFPFileBatch.m)'},...
                                 'Select a file to open');
if isnumeric(Handles.FileName)  % No file selected
    return;
end
Handles.PathName(Handles.PathName=='\') = '/';
% Check file type
CurrFileName = lower(Handles.FileName);
if ~isempty(strfind(CurrFileName, 'waveform.mat'))
    Handles.FileType = 1;
elseif ~isempty(strfind(CurrFileName, 'lfp.mat'))
    Handles.FileType = 2;
elseif ~isempty(strfind(CurrFileName, 'waveformbatch.m'))
    Handles.FileType = 3;
    run([Handles.PathName Handles.FileName]);
    Handles.FilePath = FilePath;
elseif ~isempty(strfind(CurrFileName, 'lfpbatch.m'))
    Handles.FileType = 4;
    run([Handles.PathName Handles.FileName]);
    Handles.FilePath = FilePath;
elseif ~isempty(strfind(CurrFileName, 'spikefilebatch.m'))
    Handles.FileType = 5;
    run([Handles.PathName Handles.FileName]);
    Handles.FilePath = FilePath;
    Handles.ElecID = ElecNum;
elseif ~isempty(strfind(CurrFileName, 'lfpfilebatch.m'))
    Handles.FileType = 6;
    run([Handles.PathName Handles.FileName]);
    Handles.FilePath = FilePath;
    Handles.ElecID = ElecNum;
end
Handles = Port2Core(Handles);
guidata(Handles.SpikeClusterFig, Handles);

% Couple with sorting algorithms
function Handles = Port2Core(Handles)
switch Handles.FileType
    case {1;2;}  % *Waveform.mat; *LFP.mat
        Handles.Progress = '';
        Handles = CentralProcess(Handles);
    case {3;4;}  % *WaveformBatch.m; *LFPBatch.m
        set(Handles.WavPop, 'Value', 2);
        FilePath = Handles.FilePath;
        BatchFileType = Handles.FileType;
        for i = 1:numel(FilePath)
            Handles.Progress = [' {' num2str(i) '/' num2str(numel(FilePath)) ' Files}'];
            % Reset parameters
            Handles.Para = Handles.RawPara;
            Handles.FileType = BatchFileType-2;
            FilePath{i}(FilePath{i}=='\') = '/';
            Mark = strfind(FilePath{i}, '/');
            Handles.PathName = FilePath{i}(1:Mark(end));
            Handles.FileName = FilePath{i}(Mark(end)+1:end);
            if exist([Handles.PathName Handles.FileName], 'file')==2
                Handles = CentralProcess(Handles);
                Handles = ContinuePush_Callback(gcbo, gcbo, Handles);
                SavePush_Callback(gcbo, gcbo, Handles)
            end
        end
        set(Handles.SpikeClusterFig, 'Name', 'Spike Cluster #OK');
        set(Handles.Status, 'String', 'No file in processing');
        drawnow;
    case {5;6;}  % *SpikeFileBatch.m;*LFPFileBatch.m
        set(Handles.WavPop, 'Value', 2);
        FilePath = Handles.FilePath;
        ElecID = Handles.ElecID;
        BatchFileType = Handles.FileType;
        for i = 1:numel(FilePath)
            for j = 1:numel(ElecID{i})
                Handles.Progress = [' {' num2str(i) '/' num2str(numel(FilePath)) ' Files; '...
                                         num2str(j) '/' num2str(numel(ElecID{i})) ' Channels}'];
                % Reset parameters
                Handles.Para = Handles.RawPara;
                % File path
                Mark = strfind(FilePath{i}, '.');
                Handles.PathName = [FilePath{i}(1:Mark(end)-1) '/'];
                Handles.PathName(Handles.PathName=='\') = '/';
                if BatchFileType==5  %Waveform
                    Handles.FileName = ['Elec' num2str(ElecID{i}(j)) 'Waveform.mat'];
                elseif BatchFileType==6  %LFP
                    Handles.FileName = ['Elec' num2str(ElecID{i}(j)) 'LFP.mat'];
                end
                Handles.FileType = BatchFileType-4;
                if exist([Handles.PathName Handles.FileName], 'file')==2
                    Handles = CentralProcess(Handles);
                    Handles = ContinuePush_Callback(gcbo, gcbo, Handles);
                    SavePush_Callback(gcbo, gcbo, Handles)
                end
            end
        end
        set(Handles.SpikeClusterFig, 'Name', 'Spike Cluster #OK');
        set(Handles.Status, 'String', 'No file in processing');
        drawnow;
end
guidata(Handles.SpikeClusterFig, Handles);

function Handles = CentralProcess(Handles)
clc;
cd(Handles.PathName);
set(Handles.Status, 'String', [Handles.PathName Handles.FileName]);
% Reset all axes
FieldName = fieldnames(Handles);
for i = 1:numel(FieldName)
    if ~isempty(strfind(FieldName{i}, 'Axes'))
        cla(Handles.(FieldName{i}), 'reset');
    end
end
for i = 0:3
    set(Handles.(['Unit' num2str(i)]), 'String', '');
    set(Handles.(['Index' num2str(i)]), 'String', '');
end
for i = 1:3
    set(Handles.(['Wav' num2str(i) 'Pop']), 'Value', 1);
end
% Reset all the variables
FieldName = {'NumSpikeTrial';'Wav';'Spike';'Fs';'Coeff';'Clustered';'AssignedRemainder';'AssignedOutlier';...
             'Fate';'Unit';'MinSize';'Model';'Class';'Temp';'ClusterUnit';'Quality';'ExpDuration';...
             'ElecID';'CoeffNum';'VoltThre';'WavNum';'Outlier';'WBSignal';'SpikeSignal';'LFPSignal';...
             'TSpikeSignal';'TrialID';'RecordedWav';'RecordedSpike';'Rethreshold';'RecordedNumSpikeTrial';};
for i = 1:numel(FieldName)
    if isfield(Handles, FieldName{i})
        Handles = rmfield(Handles, FieldName{i});
    end
end
% Reset parameters
Handles.Para = Handles.RawPara;
% Load
set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Loading' Handles.Progress]);
drawnow;
Handles = LoadData(Handles);
% Filter
if Handles.Para.Filter
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Filtering'  Handles.Progress]);
    drawnow;
    for i = 1:numel(Handles.WBSignal)
        SpikePara = Handles.Para;
        SpikePara.HipassFreq = SpikePara.HipassFreq(1);
        SpikePara.LopassFreq = SpikePara.LopassFreq(1);
        SpikePara.HipassOrder = SpikePara.HipassOrder(1);
        SpikePara.LopassOrder = SpikePara.LopassOrder(1);
        Handles.SpikeSignal{i,1} = Filtering(Handles.WBSignal{i}, SpikePara);
        if numel(Handles.Para.HipassFreq)==2 && numel(Handles.Para.HipassFreq)==2 && ...
                numel(Handles.Para.HipassOrder)==2 && numel(Handles.Para.HipassOrder)==2
            LFPPara = Handles.Para;
            LFPPara.HipassFreq = LFPPara.HipassFreq(2);
            LFPPara.LopassFreq = LFPPara.LopassFreq(2);
            LFPPara.HipassOrder = LFPPara.HipassOrder(2);
            LFPPara.LopassOrder = LFPPara.LopassOrder(2);
            Handles.LFPSignal{i,1} = Filtering(Handles.WBSignal{i}, LFPPara);
            %remove line noise
%             LN.tapers = [3 5]; LN.Fs = Handles.Fs; LN.fpass = [0 500];
%             Handles.LFPSignal{i,1} = rmlinesc(Handles.LFPSignal{i,1}, LN, [], 'n', 50);
        end
    end
end
% Transform spike signals
if Handles.Para.Det && Handles.FileType~=1  %not waveform
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Transforming' Handles.Progress]);
    drawnow;
    switch Handles.Para.DetMethod
        case 1  %Raw Voltage
            Handles.TSpikeSignal = Handles.SpikeSignal;
        case 2  %Nonlinear Energy Operator; Mukhopadhyay IEE98
            for i = 1:numel(Handles.SpikeSignal)
                CurrSig = Handles.SpikeSignal{i};
                if isempty(CurrSig)
                    Handles.TSpikeSignal{i,1} = [];
                else
                    NLE = CurrSig.^2-CurrSig([1 1:end-1]).*CurrSig([2:end end]);
                    %NLE = NLE+CurrSig.^2-CurrSig([1 1 1:end-2]).*CurrSig([3:end end end]);
                    %NLE = NLE+CurrSig.^2-CurrSig([1 1 1 1:end-3]).*CurrSig([4:end end end end]);
                    %NLE = NLE+CurrSig.^2-CurrSig([1 1 1 1 1:end-4]).*CurrSig([5:end end end end end]);
                    Handles.TSpikeSignal{i,1} = NLE;%conv(NLE, bartlett(round(6/30000*Handles.Fs)), 'same');
                end
            end
        case 3  %Mathematical Morphonogy; Geng et al., Neurocomputing 2010
            SpikeSignal = cell2mat(Handles.SpikeSignal);
            SpikeSignal = SpikeSignal(1:min(10000, numel(SpikeSignal)));
            Gauss = gmdistribution.fit(SpikeSignal', 3);
            Amp = max(Gauss.mu)-min(Gauss.mu);
            N = round(0.875*Handles.Fs/1000);
            Element = Amp*(0.54-0.46*cos(2*pi*(0:N)/N));
            for i = 1:numel(Handles.SpikeSignal)
                CurrSig = Handles.SpikeSignal{i};
                OpenClose = MathMorph(MathMorph(CurrSig, Element, 'opening'), Element, 'closing');
                CloseOpen = MathMorph(MathMorph(CurrSig, Element, 'closing'), Element, 'opening');
                PVE = CurrSig-0.5*(OpenClose+CloseOpen)';  % Peak-to-valley energy
                Handles.TSpikeSignal{i,1} = PVE;
            end
    end
    Handles.TrialID = 1;
    NumTrial = numel(Handles.WBSignal);
    warning('off');
    set(Handles.TrialIDSlider, 'min', 1, 'max', NumTrial, 'sliderstep', [1 1]/NumTrial, 'Value', Handles.TrialID);
    set(Handles.TrialIDEdit, 'String', Handles.TrialID);
    Handles.Thre = GetThre(Handles);
    PlotContSignal(Handles);
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Detecting' Handles.Progress]);
    drawnow;
    for i = 1:numel(Handles.SpikeSignal)
        Para = Handles.Para;
        Para.Thre = Handles.Thre(i,:);
        [Spike{i,1}, Waveform{i,1}] = Detecting(Handles.TSpikeSignal{i}, Handles.SpikeSignal{i}, Para);
    end
    Handles.Wav = single(cell2mat(Waveform));
    Handles.Spike = cell2mat(Spike);
    for i = numel(Waveform):-1:1
        Handles.NumSpikeTrial(i,1) = size(Waveform{i}, 1);
    end
    NumSpike = 0;
    for i = 1:size(Handles.StartT, 1)
%         if ~isnan(Handles.TTLT(i))
%             Handles.Spike(NumSpike+1:NumSpike+Handles.NumSpikeTrial(i)) = ...
%                 Handles.Spike(NumSpike+1:NumSpike+Handles.NumSpikeTrial(i))+Handles.TTLT(i)+Handles.StartT(i);
%         end
        NumSpike = NumSpike+Handles.NumSpikeTrial(i);
    end
    Handles.RawSpike = Handles.Spike;  %for outlier
    AlignID = round(Handles.Para.AlignTime*Handles.Fs/1000)+1;
    Handles.VoltThre = max(Handles.Wav(:,AlignID));
end
Handles.RecordedWav = Handles.Wav;
Handles.RecordedSpike = Handles.RawSpike;
Handles.RecordedNumSpikeTrial = Handles.NumSpikeTrial;
if Handles.Para.Det && Handles.FileType==1
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Rethresholding' Handles.Progress]);
    drawnow;
    Handles.Rethreshold = Rethresholding(Handles);
else
    Handles.Rethreshold = true([size(Handles.Wav, 1) 1]);
end
Handles.Wav = Handles.Wav(Handles.Rethreshold,:);
Handles.Spike = Handles.Spike(Handles.Rethreshold,:);
IsSpike = mat2cell(Handles.Rethreshold, Handles.NumSpikeTrial);
for i = 1:numel(IsSpike)
    Handles.NumSpikeTrial(i,1) = sum(IsSpike{i});
end
NumRow = size(Handles.Wav, 1);
Handles.Unit = zeros([NumRow 1], 'int32')-2;
Handles.Outlier = false([NumRow 1]);
if NumRow==0
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #No Spikes'  Handles.Progress]);
    return;
end
% Reverse filter to restore the real waveforms
if Handles.Para.ReverseFilterWav==2 && Handles.FileType==1  %recorded waveform
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Reverse Filtering' Handles.Progress]);
    drawnow;
    ReversePara = Handles.Para;
    ReversePara.HipassFreq = ReversePara.HipassFreq(1);
    ReversePara.LopassFreq = ReversePara.LopassFreq(1);
    ReversePara.HipassOrder = ReversePara.HipassOrder(1);
    ReversePara.LopassOrder = ReversePara.LopassOrder(1);
    ReversePara.FilterDirection = 3;  %Reverse Only
    Handles.Wav = Filtering(Handles.Wav', ReversePara)';
end
% Align
if Handles.Para.Align
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Aligning' Handles.Progress]);
    drawnow;
    InterpMethod = Handles.Para.InterpMethod;
    switch InterpMethod
        case 1  %Raw Waveform
        case {2;3;4;5;}
            Str = {'cubic';'v5cubic';'linear';'nearest'};
            Handles.Wav = Resampling(Handles.Wav, Handles.Fs, Handles.Para.InterpReso, Str{InterpMethod});
    end
    [Handles.Wav, Handles.AlignedID] = Aligning(Handles.Wav, Handles.Para.AlignMethod);
else
    Handles.AlignedID = false(size(Handles.Unit));
end
% Extract features
set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Extracting' Handles.Progress]);
drawnow;
if NumRow>10^4  
    set(Handles.WavPop, 'Value', min(2, get(Handles.WavPop, 'value')));
end
NotAlignID = 1:NumRow;
NotAlignID = NotAlignID(~Handles.AlignedID);
MaxNum = min(numel(NotAlignID), Handles.Para.NumClusteredSpike);
if MaxNum<numel(NotAlignID)
    Handles.Para.AssignRemainder = min(5, Handles.Para.AssignRemainder);
end
Handles.WavNum = NotAlignID(randperm(numel(NotAlignID), MaxNum));  % Waveforms for extracting features
if Handles.Para.Extract
    Method = Handles.Para.ExtractMethod;
    if Method==6
        WaveletType = Handles.Para.MultiwaveletTypeOpt{Handles.Para.WaveletType};
    else
        WaveletType = Handles.Para.WaveletTypeOpt{Handles.Para.WaveletType};
    end
    WaveletScale = Handles.Para.WaveletScale;
    Handles.Unit(Handles.WavNum) = 1;
    Coeff = ExtractFeature(Handles.Wav(Handles.WavNum,:), Method, WaveletType, WaveletScale);
else
    Coeff = Handles.Wav(Handles.WavNum,:);
    Handles.Unit(Handles.WavNum) = 1;
end
% Select coefficients
if Handles.Para.Select
    Method = Handles.Para.SelectMethod;
    VarRatio = Handles.Para.VarRatio;
    NumCoeff = Handles.Para.NumCoeff;
    try
        Handles.CoeffNum = SortFeature(Coeff, Method, VarRatio, NumCoeff);
    catch
        Handles.CoeffNum = SortFeature(Coeff, 3, VarRatio, NumCoeff);  %KS test
    end
    Handles.Coeff = Coeff(:,Handles.CoeffNum);
else
    Handles.Coeff = Coeff;
end
clear Coeff;
% Set processing state, for Continue pushbotton
Handles.Clustered = 0;
Handles.AssignedRemainder = 0;
Handles.AssignedOutlier = 0;
Handles.Fate = [1 1 1];  % Unit1/2/3Pop
% Minimum number of waveforms to do clustering
Handles.MinNumWav = 30;
if size(Handles.Wav, 1)<Handles.MinNumWav
    Handles.Unit(Handles.Unit~=1) = 1;
    PlotAllWav(Handles);
    PlotUnitWav(Handles);
    PlotIndex(Handles);
    guidata(Handles.SpikeClusterFig, Handles);
    return;
end
% Cluster
if Handles.Para.Cluster
    OutlierRemoval = Handles.Para.OutlierRemoval;
    ClusterMethod = Handles.Para.ClusterMethod;
    NumCluster = Handles.Para.NumCluster;
    Handles.Unit(Handles.WavNum) = 1;
    if (OutlierRemoval==3 || OutlierRemoval==4)  % Do a pre-cluster, then remove outliers in each cluster
        set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Pre-clustering'  Handles.Progress]);
        drawnow;
        if OutlierRemoval==4  % Remove global outliers
            Handles.Unit(Handles.WavNum) = RemoveOutlier(Handles.Coeff, Handles.Unit(Handles.WavNum), 1e-4);
            WavIdx = Handles.Unit(Handles.WavNum)>=0;
            Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
            Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
        end
        WavIdx = Handles.Unit(Handles.WavNum)>=0;
        [Model Class] = Clustering(Handles.Coeff(WavIdx,:), ClusterMethod, NumCluster);
        Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
        if ClusterMethod==1  % Single-linkage, do a mandatory removal of outliers in dendrogram
            [Model Class] = Linkage2SPC(Model, Handles.MinSize);
            Handles.Temp = SelectTemp(Model, Handles.MinSize);
            Class = Temp2Class(Class, Handles.Temp, Handles.MinSize);
            Handles.Unit(Handles.WavNum(WavIdx)) = ChkClusterSize(Class, Handles.MinSize);
            Handles.Unit(Handles.Unit==0) = -1;
            Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 0.01);
%             Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
            Handles.Unit(Handles.WavNum) = SortCluster(Handles.Wav(Handles.WavNum,:), Handles.Coeff, Handles.Unit(Handles.WavNum));
            WavIdx = Handles.Unit(Handles.WavNum)>=0;
            [Model, ~] = Clustering(Handles.Coeff(WavIdx,:), ClusterMethod, NumCluster);
            Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
            [Model Class] = Linkage2SPC(Model, Handles.MinSize);
        end
        Handles.Model = Model;
        Handles.Class = Class;
        if ClusterMethod==1 || ClusterMethod==5  % Single-linkage/SPC
            Handles.Temp = SelectTemp(Model, Handles.MinSize);
            Class = Temp2Class(Class, Handles.Temp, Handles.MinSize);
        end
        Handles.Unit(Handles.WavNum(WavIdx)) = ChkClusterSize(Class, Handles.MinSize);
        Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 0.01, 0);
        Handles.Unit = MergeCluster(Handles.Wav, Handles.Unit);
%         Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
        Handles.Unit(Handles.WavNum) = SortCluster(Handles.Wav(Handles.WavNum,:), Handles.Coeff, Handles.Unit(Handles.WavNum));
        Handles.Unit = SortSNR(Handles.Wav, Handles.Unit, Handles.Outlier);
        Handles.ClusterUnit = Handles.Unit;
        % Display results
        PlotModel(Handles);
        PlotAllWav(Handles);
        PlotUnitWav(Handles);
        PlotIndex(Handles);
    else  % Do final clustering
        Handles = ContinuePush_Callback(gcbo, gcbo, Handles);
        Handles.Clustered = 1;
    end
    Handles.Track = TrackCD(Handles.Track);
end
% Hidden Markov Model
if Handles.Para.HMM
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #HMM' Handles.Progress]);
    drawnow;
end
% Save
guidata(Handles.SpikeClusterFig, Handles);
set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster' Handles.Progress]);
drawnow;

function Handles = ContinuePush_Callback(~, ~, Handles)
if size(Handles.Wav, 1)<Handles.MinNumWav
    return;
end
if ~Handles.Clustered  % Do final clustering
    OutlierRemoval = Handles.Para.OutlierRemoval;
    ClusterMethod = Handles.Para.ClusterMethod;
    NumCluster = Handles.Para.NumCluster;
    if Handles.Para.Cluster
        set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Clustering' Handles.Progress]);
        drawnow;
        Unit = Handles.Unit(Handles.WavNum);
        Handles.Unit = zeros([size(Handles.Unit, 1) 1], 'int32')-2;
        Handles.Unit(Handles.WavNum) = Unit;
        if OutlierRemoval==2 || OutlierRemoval==4  % Remove global outliers
            Handles.Unit(Handles.WavNum) = RemoveOutlier(Handles.Coeff, Handles.Unit(Handles.WavNum), 1e-4);
            WavIdx = Handles.Unit(Handles.WavNum)>=0;
            Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
            Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
        end
        WavIdx = Handles.Unit(Handles.WavNum)>=0;
        Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
        [Model Class] = Clustering(Handles.Coeff(WavIdx,:), ClusterMethod, NumCluster);
        if ClusterMethod==1  % Single-linkage, do a mandatory removal of outliers in dendrogram
            [Model Class] = Linkage2SPC(Model, Handles.MinSize);
            Temp = SelectTemp(Model, Handles.MinSize);
            Class = Temp2Class(Class, Temp, Handles.MinSize);
            Handles.Unit(Handles.WavNum(WavIdx)) = ChkClusterSize(Class, Handles.MinSize);
            Handles.Unit(Handles.Unit==0) = -1;
            Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 0.01);
            WavIdx = Handles.Unit(Handles.WavNum)>=0;
            Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
%             Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
            WavIdx = Handles.Unit(Handles.WavNum)>=0;
            Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
            [Model, ~] = Clustering(Handles.Coeff(WavIdx,:), ClusterMethod, NumCluster);
            [Model Class] = Linkage2SPC(Model, Handles.MinSize);
        end
        Handles.Model = Model;
        Handles.Class = Class;
        if ClusterMethod==1 || ClusterMethod==5  % Single-linkage/SPC
            Handles.Temp = SelectTemp(Model, Handles.MinSize);
            Class = Temp2Class(Class, Handles.Temp, Handles.MinSize);
        end
        Handles.Unit(Handles.WavNum(WavIdx)) = ChkClusterSize(Class, Handles.MinSize);
        Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 0.05, 0);
        Handles.Unit = MergeCluster(Handles.Wav, Handles.Unit);
%         Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
        Handles.Unit(Handles.WavNum) = SortCluster(Handles.Wav(Handles.WavNum,:), Handles.Coeff, Handles.Unit(Handles.WavNum));
        Handles.Unit = SortSNR(Handles.Wav, Handles.Unit, Handles.Outlier);
        Handles.ClusterUnit = Handles.Unit;
    end
    PlotModel(Handles);
    PlotAllWav(Handles);
    PlotUnitWav(Handles);
    PlotIndex(Handles);
    Handles.Clustered = 1;
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
else
    % Clustering procedure only identifies one cluster!!!
    if Handles.Para.AutoTemplate~=1 && sum(unique(Handles.Unit)>0)==1
        Handles = AutoTemplatePush_Callback(gcbo, gcbo, Handles);
    end
    % Assign other uncalssified spikes
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Assigning Remainders'  Handles.Progress]);
    drawnow;
%     if ~Handles.AssignedRemainder
        Handles = AssignRemainderPush_Callback(gcbo, gcbo, Handles);
        for i = 1:5
            if sum(Handles.Unit<=0)/numel(Handles.Unit)>=0.01 || sum(Handles.Unit<=0)>1000
                Handles = AssignRemainderPush_Callback(gcbo, gcbo, Handles);
            end
        end
%     end
%     if ~Handles.AssignedOutlier  % Very time-consuming
        for i = 1:3
            if sum(Handles.Unit<=0)/numel(Handles.Unit)>=10^-3 && sum(Handles.Unit<=0)<=1000
                Handles = AssignOutlierPush_Callback(gcbo, gcbo, Handles);
            end
        end
%     end
    for i = 1:2
        Handles = ExcludeLFPPush_Callback(gcbo, gcbo, Handles);
    end
    Handles = MergePush_Callback(gcbo, gcbo, Handles);
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
    set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster'  Handles.Progress]);
    drawnow;
end
if isfield(Handles, 'Quality')
    Handles = rmfield(Handles, 'Quality');
end
guidata(Handles.SpikeClusterFig, Handles);

% Autocreate templates for assigning
function Handles = AutoTemplatePush_Callback(~, ~, Handles)
% Autocreate templates only when outputed a single cluster with very low SNR
if size(Handles.Wav, 1)<Handles.MinNumWav
    return;
end
ClassSNR = GetSNR(Handles.Wav, Handles.Unit);
if ~isempty(ClassSNR) && ClassSNR(1)>=Handles.Para.AutotempSNRThre;
    return;
end
if Handles.Para.AssignRemainder>0
    Handles.IsAutoTemp = 1;
end
Handles.Para.AssignRemainder = 3;
NumRow = size(Handles.Wav, 1);
Handles.Unit = zeros([NumRow 1], 'int32');
OutlierRemoval = Handles.Para.OutlierRemoval;
if OutlierRemoval==2 || OutlierRemoval==4  % Remove global outliers
    Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit+1, 1e-4);
    Handles.Unit(Handles.Unit==1) = 0;
end
NotAlignID = 1:NumRow;
NotAlignID = NotAlignID(~logical(Handles.AlignedID+(Handles.Unit<0)));
[NumWav, NumPoint] = size(Handles.Wav(NotAlignID,:));
NumKeyWav = ceil(NumWav*0.005);
AutoTemplate = max(2, Handles.Para.AutoTemplate);
NumRound = round(0.15*Handles.Fs/1000);
MedianWav = median(Handles.Wav(NotAlignID,:), 1);
AlignTime = Handles.Para.AlignTime;
AlignToID = round(AlignTime*Handles.Fs/1000+1);
PreID = max(1, AlignToID-round(0.3*Handles.Fs/1000+1));
PostID = min(NumPoint, AlignToID+round(0.3*Handles.Fs/1000+1));
[~, MaxIdx] = max(MedianWav(PreID:PostID));
[~, MinIdx] = min(MedianWav(PreID:PostID));
MaxIdx = MaxIdx+(-NumRound:NumRound);
if MaxIdx(end)>NumPoint
    MaxIdx = MinIdx+(1:2*NumRound);
end
MaxIdx = max(1, MaxIdx(1)):min(NumPoint, MaxIdx(end));
if MaxIdx(1)>MaxIdx(end)
    MaxIdx = 1:NumPoint;
end
MinIdx = MinIdx+(-NumRound:NumRound);
if MinIdx(1)<1
    MinIdx = round(0.3*Handles.Fs/1000)+(-NumRound:NumRound);
end
MinIdx = max(1, MinIdx(1)):min(NumPoint, MinIdx(end));
if MinIdx(1)>MinIdx(end)
    MinIdx = 1:NumPoint;
end
if AutoTemplate==2  % Peak-Valley Maximum
    PVDiff = max(Handles.Wav(NotAlignID,MaxIdx), [], 2)-min(Handles.Wav(NotAlignID,MinIdx), [], 2);
    [~, Idx] = sort(PVDiff, 'descend');
    Handles.Unit(NotAlignID(Idx(1:NumKeyWav))) = 1;
    Handles.Unit(NotAlignID(Idx(end-NumKeyWav:end))) = 2;
end
if AutoTemplate==3  % Mathematical Morphonogy
    Peak = sort(max(abs(Handles.Wav(NotAlignID,:)), [], 2), 'descend');
    Element = mean(Peak(round(numel(Peak)*0.05)))*(0.54-0.46*cos(2*pi*(1:25)/26));
    PaddedWav = [repmat(Handles.Wav(NotAlignID,1), [1 NumPoint]) Handles.Wav(NotAlignID,:) repmat(Handles.Wav(NotAlignID,end), [1 NumPoint])];
    OpenClose = MathMorph(MathMorph(PaddedWav, Element, 'opening'), Element, 'closing');
    CloseOpen = MathMorph(MathMorph(PaddedWav, Element, 'closing'), Element, 'opening');
    clear PaddedWav;
    % Peak-to-valley energy
    PVE = 0.5*(OpenClose+CloseOpen);
    PVE = Handles.Wav(NotAlignID,:)-PVE(:,NumPoint+1:2*NumPoint);
    [~, Idx] = sort(max(PVE(:,MaxIdx), [], 2), 'descend');
    Handles.Unit(NotAlignID(Idx(1:NumKeyWav))) = 1;
    Handles.Unit(NotAlignID(Idx(end-NumKeyWav:end))) = 2;
end
if AutoTemplate==4  % Peak-Valley Maximum
    NumKeyWav = ceil(numel(Handles.WavNum)*0.01);
    PVDiff = max(Handles.Wav(Handles.WavNum,MaxIdx), [], 2)-min(Handles.Wav(Handles.WavNum,MinIdx), [], 2);
    [~, Idx] = sort(PVDiff, 'descend');
    Handles.Unit(Handles.WavNum(Idx(1:NumKeyWav))) = 1;
    Handles.Unit(Handles.WavNum(Idx(end-NumKeyWav:end))) = 2;
end
if AutoTemplate==5  % Mathematical Morphonogy
    NumKeyWav = ceil(numel(Handles.WavNum)*0.01);
    Peak = sort(max(abs(Handles.Wav(Handles.WavNum,:)), [], 2), 'descend');
    Element = mean(Peak(round(numel(Peak)*0.05)))*(0.54-0.46*cos(2*pi*(1:25)/26));
    PaddedWav = [repmat(Handles.Wav(Handles.WavNum,1), [1 NumPoint]) ...
                 Handles.Wav(Handles.WavNum,:)...
                 repmat(Handles.Wav(Handles.WavNum,end), [1 NumPoint])];
    OpenClose = MathMorph(MathMorph(PaddedWav, Element, 'opening'), Element, 'closing');
    CloseOpen = MathMorph(MathMorph(PaddedWav, Element, 'closing'), Element, 'opening');
    clear PaddedWav;
    % Peak-to-valley energy
    PVE = 0.5*(OpenClose+CloseOpen);
    PVE = Handles.Wav(Handles.WavNum,:)-PVE(:,NumPoint+1:2*NumPoint);
    [~, Idx] = sort(max(PVE(:,MaxIdx), [], 2), 'descend');
    Handles.Unit(Handles.WavNum(Idx(1:NumKeyWav))) = 1;
    Handles.Unit(Handles.WavNum(Idx(end-NumKeyWav:end))) = 2;
end
if AutoTemplate==4 || AutoTemplate==5  % Do clustering
    ClusterMethod = Handles.Para.ClusterMethod;
    NumCluster = Handles.Para.NumCluster;
    WavIdx = Handles.Unit(Handles.WavNum)>=0;
    [Model Class] = Clustering(Handles.Coeff(WavIdx,:), ClusterMethod, NumCluster);
    Handles.MinSize = max(2, NumKeyWav*Handles.Para.MinClusterRatio*0.01);
    if ClusterMethod==1  % Single-linkage, do a mandatory removal of outliers in dendrogram
        [Model Class] = Linkage2SPC(Model, Handles.MinSize);
        Temp = SelectTemp(Model, Handles.MinSize);
        Class = Temp2Class(Class, Temp, Handles.MinSize);
        Handles.Unit(Handles.WavNum(WavIdx)) = ChkClusterSize(Class, Handles.MinSize);
        Handles.Unit(Handles.Unit==0) = -1;
        Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 0.01, -1, 0);
%         Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
        WavIdx = Handles.Unit(Handles.WavNum)>=0;
        Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
        [Model, ~] = Clustering(Handles.Coeff(WavIdx,:), ClusterMethod, NumCluster);
        [Model Class] = Linkage2SPC(Model, Handles.MinSize);
    end
    Handles.Model = Model;
    Handles.Class = Class;
    if ClusterMethod==1 || ClusterMethod==5  % Single-linkage/SPC
        Handles.Temp = SelectTemp(Model, Handles.MinSize);
        Class = Temp2Class(Class, Handles.Temp, Handles.MinSize);
    end
    Handles.Unit(Handles.WavNum(WavIdx)) = ChkClusterSize(Class, Handles.MinSize);
    Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 0.05, 0, 0);
    Handles.Unit = MergeCluster(Handles.Wav, Handles.Unit);
%     Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
    Handles.Unit(Handles.WavNum) = SortCluster(Handles.Wav(Handles.WavNum,:), Handles.Coeff, Handles.Unit(Handles.WavNum));
    Handles.ClusterUnit = Handles.Unit;
end
Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 1e-2);
WavIdx = Handles.Unit(Handles.WavNum)>0;
Handles.MinSize = max(2, sum(WavIdx)*Handles.Para.MinClusterRatio*0.01);
Class = ChkClusterSize(Handles.Unit, Handles.MinSize);
if sum(unique(Class)>0)>=2
    Handles.Unit = Class;
end
Handles.Unit = SortSNR(Handles.Wav, Handles.Unit, Handles.Outlier);
PlotModel(Handles);
PlotAllWav(Handles);
PlotUnitWav(Handles);
PlotIndex(Handles);
Handles.Fate(3) = 3;  % Leave only two clusters
Handles.AssignedRemainder = 0;
Handles.AssignedOutlier = 0;
if isfield(Handles, 'Quality')
    Handles = rmfield(Handles, 'Quality');
end
guidata(Handles.SpikeClusterFig, Handles);

% Plot all the axes
function ContSignalPop_Callback(~, ~, Handles)
if Handles.FileType~=2  %LFP
    return;
end
PlotContSignal(Handles);
function TrialIDEdit_Callback(~, ~, Handles)
if Handles.FileType~=2  %LFP
    return;
end
NumTrial = numel(Handles.WBSignal);
Handles.TrialID = min(NumTrial, round(eval(get(gcbo, 'String'))));
Handles.TrialID = max(1, Handles.TrialID);
set(Handles.TrialIDEdit, 'String', num2str(Handles.TrialID));
SliderMax = get(Handles.TrialIDSlider, 'max');
set(Handles.TrialIDSlider, 'Value', min(SliderMax, Handles.TrialID(1)));
PlotContSignal(Handles);
guidata(Handles.SpikeClusterFig, Handles);
function TrialIDSlider_Callback(~, ~, Handles)
if Handles.FileType~=2  %LFP
    return;
end
Handles.TrialID = round(get(gcbo, 'Value'));
set(Handles.TrialIDEdit, 'String', Handles.TrialID(1));
PlotContSignal(Handles);
guidata(Handles.SpikeClusterFig, Handles);

function WavPop_Callback(~, ~, Handles)
PlotAllWav(Handles);
PlotUnitWav(Handles);
function IndexPop_Callback(~, ~, Handles)
PlotIndex(Handles);

% Unit fate
function Wav1Pop_Callback(~, ~, Handles)
if sum(unique(Handles.Unit)==1)>0
    Handles.Fate(1) = get(Handles.Wav1Pop, 'Value');
    if Handles.Fate(1)==4  % Rejected & Refresh
        Handles.Unit(Handles.Unit==1) = 0;
        Handles.Unit(Handles.Unit>1) = Handles.Unit(Handles.Unit>1)-1;
        Handles.Fate = [Handles.Fate(2:end) 1];
        PlotAllWav(Handles);
        PlotUnitWav(Handles);
        PlotIndex(Handles);
    end
    for i = 1:numel(Handles.Fate)
        set(Handles.(['Wav' num2str(i) 'Pop']), 'Value', Handles.Fate(i));
    end
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
    guidata(Handles.SpikeClusterFig, Handles);
else
    set(Handles.Wav1Pop, 'Value', 1);
end

function Wav2Pop_Callback(~, ~, Handles)
if sum(unique(Handles.Unit)==2)>0
    Handles.Fate(2) = get(Handles.Wav2Pop, 'Value');
    if Handles.Fate(2)==4  % Rejected & Refresh
        Handles.Unit(Handles.Unit==2) = 0;
        Handles.Unit(Handles.Unit>2) = Handles.Unit(Handles.Unit>2)-1;
        Handles.Fate = [Handles.Fate(1) Handles.Fate(3:end) 1];
        PlotAllWav(Handles);
        PlotUnitWav(Handles);
        PlotIndex(Handles);
    end
    for i = 1:numel(Handles.Fate)
        set(Handles.(['Wav' num2str(i) 'Pop']), 'Value', Handles.Fate(i));
    end
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
    guidata(Handles.SpikeClusterFig, Handles);
else
    set(Handles.Wav2Pop, 'Value', 1);
end

function Wav3Pop_Callback(~, ~, Handles)
if sum(unique(Handles.Unit)==3)>0
    Handles.Fate(3) = get(Handles.Wav3Pop, 'Value');
    if Handles.Fate(3)==4  % Rejected & Refresh
        Handles.Unit(Handles.Unit==3) = 0;
        Handles.Unit(Handles.Unit>3) = Handles.Unit(Handles.Unit>3)-1;
        Handles.Fate = [Handles.Fate(1:2) Handles.Fate(4:end) 1];
        PlotAllWav(Handles);
        PlotUnitWav(Handles);
        PlotIndex(Handles);
    end
    for i = 1:numel(Handles.Fate)
        set(Handles.(['Wav' num2str(i) 'Pop']), 'Value', Handles.Fate(i));
    end
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
    guidata(Handles.SpikeClusterFig, Handles);
else
    set(Handles.Wav3Pop, 'Value', 1);
end

% Change unit number
function Unit1NumPop_Callback(~, ~, Handles)
UnitNum = get(Handles.Unit1NumPop, 'Value')-1;
if UnitNum>0 && UnitNum~=1 && sum(Handles.Unit==UnitNum)>0
    Unit = Handles.Unit;
    Handles.Unit(Unit==1) = UnitNum;
    Handles.Unit(Unit==UnitNum) = 1;
    PlotAllWav(Handles);
    PlotUnitWav(Handles);
    PlotIndex(Handles);
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
end
set(Handles.Unit1NumPop, 'Value', 1);
guidata(Handles.SpikeClusterFig, Handles);
function Unit2NumPop_Callback(~, ~, Handles)
UnitNum = get(Handles.Unit2NumPop, 'Value')-1;
if UnitNum>0 && UnitNum~=2 && sum(Handles.Unit==UnitNum)>0
    Unit = Handles.Unit;
    Handles.Unit(Unit==2) = UnitNum;
    Handles.Unit(Unit==UnitNum) = 2;
    PlotAllWav(Handles);
    PlotUnitWav(Handles);
    PlotIndex(Handles);
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
end
set(Handles.Unit2NumPop, 'Value', 1);
guidata(Handles.SpikeClusterFig, Handles);
function Unit3NumPop_Callback(~, ~, Handles)
UnitNum = get(Handles.Unit3NumPop, 'Value')-1;
if UnitNum>0 && UnitNum~=3 && sum(Handles.Unit==UnitNum)>0
    Unit = Handles.Unit;
    Handles.Unit(Unit==3) = UnitNum;
    Handles.Unit(Unit==UnitNum) = 3;
    PlotAllWav(Handles);
    PlotUnitWav(Handles);
    PlotIndex(Handles);
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
end
set(Handles.Unit3NumPop, 'Value', 1);
guidata(Handles.SpikeClusterFig, Handles);

function AllPojectionsPush_Callback(~, ~, Handles)
PlotProjections(Handles);

function AdjClusterPush_Callback(~, ~, Handles)
if size(Handles.Wav, 1)<Handles.MinNumWav
    return;
end
switch Handles.Para.ClusterMethod
    case {1;5;}  % Single-linkage/SPC
        % Index for template waveforms
        TempUnit = find(Handles.Fate==2);  % Template
        TempIdx = cell([0 0]);
        for i = 1:numel(TempUnit)
            TempIdx{i} = Handles.Unit==TempUnit(i);
        end
        % Retrieve new class
        Handles.Unit = Handles.ClusterUnit;
        axes(Handles.IndexAxes);
        hold off;
        [XVal YVal]= ginput(1);
        TempDiff = abs(Handles.Model(:,2)-XVal);
        [~, Handles.Temp] = min(TempDiff);
        YLim = get(gca, 'YLim');
        Handles.MinSize = max(2, min(max(YLim), max(min(YLim), YVal)));
        Class = Temp2Class(Handles.Class, Handles.Temp, Handles.MinSize);
        WavIdx = Handles.Unit(Handles.WavNum)>=0;
        Handles.Unit(Handles.WavNum(WavIdx)) = ChkClusterSize(Class, Handles.MinSize);
        Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 0.05, 0);
%         Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
        Handles.Unit(Handles.WavNum) = SortCluster(Handles.Wav(Handles.WavNum,:), Handles.Coeff, Handles.Unit(Handles.WavNum));
        % Assgin unit from the resultant new class and fixed template
        Idx = false(size(Handles.Unit));
        for i = 1:numel(TempIdx)
            Handles.Unit(TempIdx{i}) = i;
            Idx = Idx+TempIdx{i};
        end
        NonTempIdx = ~logical(Idx);  % Non-template waveforms
        UnitPool = sort(unique(Handles.Unit(NonTempIdx)), 'ascend');
        Unit = numel(TempUnit);
        for i = 1:numel(UnitPool)
            if UnitPool(i)>0
                Unit = Unit+1;
                Idx = logical(NonTempIdx.*(Handles.Unit==UnitPool(i)));
                Handles.Unit(Idx) = Unit;
            end
        end
        Handles.Unit = SortSNR(Handles.Wav, Handles.Unit, Handles.Outlier);
        % Display results
        PlotModel(Handles);
        PlotAllWav(Handles);
        PlotUnitWav(Handles);
        PlotIndex(Handles);
        Handles.Fate = [Handles.Fate(TempUnit) ones([1 max(0, numel(Handles.Fate)-numel(TempUnit))])];
        for i = 1:numel(Handles.Fate)
            set(Handles.(['Wav' num2str(i) 'Pop']), 'Value', Handles.Fate(i));
        end
        Handles.AssignedRemainder = 0;
        Handles.AssignedOutlier = 0;
        if isfield(Handles, 'Quality')
            Handles = rmfield(Handles, 'Quality');
        end
        guidata(Handles.SpikeClusterFig, Handles);
end

function Handles = AssignRemainderPush_Callback(~, ~, Handles)
TempUnit = find(Handles.Fate<=2);  % Auto Fate/Template
Idx = zeros(size(Handles.Unit), 'uint8');
for i = 1:numel(TempUnit)
    Idx = Idx+uint8(Handles.Unit==TempUnit(i));
    Handles.Unit(Handles.Unit==TempUnit(i)) = i;
end
Handles.Unit(~logical(Idx)) = 0;
Handles.Fate = [Handles.Fate(TempUnit) ones([1 max(0, numel(Handles.Fate)-numel(TempUnit))])];
AssignMethod = {'nn';'center';'ml';'mahal'};
NumRow = size(Handles.Wav, 1);
MaxNum = min(NumRow, Handles.Para.NumClusteredSpike);
if MaxNum<NumRow
    Handles.Para.AssignRemainder = min(5, Handles.Para.AssignRemainder);
end
if sum([6:9]==Handles.Para.AssignRemainder)>0
    Handles.Unit(Handles.Unit<=0) = Matching(Handles.Coeff(Handles.Unit(Handles.WavNum)<=0,:), ...
                                             Handles.Coeff(Handles.Unit(Handles.WavNum)>0,:),...
                                             Handles.Unit(Handles.Unit>0),...
                                             AssignMethod{Handles.Para.AssignRemainder-5}, 0.01, 10, 10);
elseif sum([2:5]==Handles.Para.AssignRemainder)>0
    Handles.Unit(Handles.Unit<=0) = Matching(Handles.Wav(Handles.Unit<=0,:), Handles.Wav(Handles.Unit>0,:),...
                                             Handles.Unit(Handles.Unit>0),...
                                             AssignMethod{Handles.Para.AssignRemainder-1}, 0.01, 10, 10);
end
Handles.AssignedRemainder = 1;
Handles.Unit = RemoveOutlier(Handles.Wav, Handles.Unit, 1e-3, -1, 1);
PlotAllWav(Handles);
PlotUnitWav(Handles);
PlotIndex(Handles);
for i = 1:numel(Handles.Fate)
    set(Handles.(['Wav' num2str(i) 'Pop']), 'Value', Handles.Fate(i));
end
Handles.AssignedRemainder = 1;
Handles.Outlier = Handles.Unit<=0;
if isfield(Handles, 'Quality')
    Handles = rmfield(Handles, 'Quality');
end
guidata(Handles.SpikeClusterFig, Handles);

function Handles = AssignOutlierPush_Callback(~, ~, Handles)
if Handles.IsAutoTemp
    return;
end
set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster #Assigning Outliers' Handles.Progress]);
drawnow;
switch Handles.Para.AssignOutlier
    case 1  % Never
    case 2  % Segment matching, use spike waveform only
        if ~Handles.AssignedRemainder
            Handles = AssignRemainderPush_Callback(gcbo, gcbo, Handles);
        end
        if sum(Handles.Unit>0)==0
            return;
        end
        NumCol = size(Handles.Wav, 2);
        WavOut = Handles.Wav(Handles.Unit<=0,:);
        ClassOut = zeros([size(WavOut, 1) 1]);
        Spike = Handles.RawSpike(Handles.Unit<=0);
        [~, MinOutIdx] = min(WavOut, [], 2);
        [~, MaxOutIdx] = max(WavOut, [], 2);
%         MuWav = mean(Handles.Wav(Handles.Unit>0,:), 1);
%         [~, MinInIdx] = min(MuWav);
%         [~, MaxInIdx] = max(MuWav);
        [~, MinInIdx] = min(Handles.Wav(Handles.Unit>0,:), [], 2);
        [~, MaxInIdx] = max(Handles.Wav(Handles.Unit>0,:), [], 2);
        MinInIdx = min(NumCol, max(1, ceil(median(MinInIdx))));
        MaxInIdx = min(NumCol, max(1, ceil(median(MaxInIdx))));
        if MinInIdx<=MaxInIdx  % Extracellular
            PreNum = min(MinOutIdx-1, MinInIdx-1);
            PostNum = min(NumCol-MinOutIdx, NumCol-MinInIdx);
            OutLoIdx = MinOutIdx-PreNum;
            OutHiIdx = MinOutIdx+PostNum;
            InLoIdx = MinInIdx-PreNum;
            InHiIdx = MinInIdx+PostNum;
        else  % Intracellular
            PreNum = min(MaxOutIdx-1, MaxInIdx-1);
            PostNum = min(NumCol-MaxOutIdx, NumCol-MaxInIdx);
            OutLoIdx = MaxOutIdx-PreNum;
            OutHiIdx = MaxOutIdx+PostNum;
            InLoIdx = MaxInIdx-PreNum;
            InHiIdx = MaxInIdx+PostNum;
        end
        for i = 1:numel(ClassOut)
            AssignMethod = {'nn';'center';'ml';'mahal'};
            ClassOut(i) = Matching(WavOut(i,OutLoIdx(i):OutHiIdx(i)), ...
                                   Handles.Wav(Handles.Unit>0,InLoIdx(i):InHiIdx(i)), ...
                                   Handles.Unit(Handles.Unit>0),...
                                   AssignMethod{Handles.Para.AssignRemainder-1}, 0.01, 10, 10);
            if ClassOut(i)>0  %adjust spike time
                if MinInIdx<=MaxInIdx  % Extracellular
                    Spike(i) = Spike(i)+(MinOutIdx(i)-10)/Handles.Fs;
                else  % Intracellular
                    Spike(i) = Spike(i)+(MaxOutIdx(i)-10)/Handles.Fs;
                end
            end
            Handles.Unit = ChkClusterSize(Handles.Unit, Handles.MinSize);
        end
        Handles.Spike(Handles.Unit<=0) = Spike;
        Handles.Unit(Handles.Unit<=0) = ClassOut;
        PlotAllWav(Handles);
        PlotUnitWav(Handles);
        PlotIndex(Handles);
        Handles.AssignedOutlier = 1;
        guidata(Handles.SpikeClusterFig, Handles);
end
set(Handles.SpikeClusterFig, 'Name', ['Spike Cluster' Handles.Progress]);
drawnow;

function Handles = ExcludeLFPPush_Callback(~, ~, Handles)
if sum(unique(Handles.Unit)>0)>1
    Handles.Unit = LFPCluster(Handles.Wav, Handles.Wav, Handles.Unit);
    PlotAllWav(Handles);
    PlotUnitWav(Handles);
    PlotIndex(Handles);
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
    guidata(Handles.SpikeClusterFig, Handles);
end

function Handles = MergePush_Callback(~, ~, Handles)
if sum(unique(Handles.Unit)>0)>1
    for i = 1:2
        Handles.Unit = MergeCluster(Handles.Wav, Handles.Unit);
    end
    Handles.Unit = SortSNR(Handles.Wav, Handles.Unit, Handles.Outlier);
    PlotAllWav(Handles);
    PlotUnitWav(Handles);
    PlotIndex(Handles);
    if isfield(Handles, 'Quality')
        Handles = rmfield(Handles, 'Quality');
    end
    guidata(Handles.SpikeClusterFig, Handles);
end

function Handles = QualityPush_Callback(~, ~, Handles)
% Not outlier & transformed spikes
QualityWav = false(size(Handles.Outlier));
QualityWav(Handles.WavNum) = 1;
QualityWav(Handles.Outlier) = 0;
QualityWav = QualityWav(Handles.WavNum);  %sequence in transormed spikes
% Sort by SNR
if sum(unique(Handles.Unit)>0)>1
    Handles.Unit = SortSNR(Handles.Wav, Handles.Unit, Handles.Outlier);
    PlotAllWav(Handles);
    PlotUnitWav(Handles);
    PlotIndex(Handles);
end
% SNR
if isfield(Handles, 'Quality')
    Handles = rmfield(Handles, 'Quality');
end
UnitPool = unique(Handles.Unit);
Idx = find(UnitPool>0);
if numel(Idx)>0 && Handles.Para.SNR
    for i = numel(Idx):-1:1
        ID = (Handles.Unit==UnitPool(Idx(i)))&(~Handles.Outlier);
        Quality.SNR(i) = SNR(Handles.Wav(ID,:), 1);
    end
end
% Joshua Quality, 2007
if numel(Idx)>0 && (Handles.Para.IsolationScore || Handles.Para.FalseNeg || Handles.Para.FalsePos)
    % Estimate the std of noise
    Wav = Handles.Wav(Handles.WavNum(QualityWav),:);
    Unit = Handles.Unit(Handles.WavNum(QualityWav),:);
    RefWav = Wav(Unit==1,:);
    MuWav = mean(RefWav, 1);
    Noise = RefWav-ones([size(RefWav, 1) 1])*MuWav;  % Detrended noise
    SD = std(Noise(:));
    % Simulate noise waveforms & transform them into coefficients
    NoiseID = ~Handles.AlignedID & Handles.Unit==0;
    NoiseWav = Handles.Coeff((QualityWav(NoiseID(Handles.WavNum))),:);
    if isempty(NoiseWav)
        NoiseWav = normrnd(0, SD, [250 size(RefWav, 2)]);
        if Handles.Para.Extract
            Method = Handles.Para.ExtractMethod;
            if Method==6
                WaveletType = Handles.Para.MultiwaveletTypeOpt{Handles.Para.WaveletType};
            else
                WaveletType = Handles.Para.WaveletTypeOpt{Handles.Para.WaveletType};
            end
            WaveletScale = Handles.Para.WaveletScale;
            NoiseWav = ExtractFeature(NoiseWav, Method, WaveletType, WaveletScale);
        end
        if Handles.Para.Select
            NoiseWav = NoiseWav(:,Handles.CoeffNum);
        end
    end
    Wav = Handles.Coeff(QualityWav,:);
    for i = numel(Idx):-1:1
        [Quality.FalseNegNoise(i), Quality.FalsePosNoise(i), Quality.IsolationScore(i)] = ...
                 JoshuaQuality(Wav(Unit==UnitPool(Idx(i)),:), NoiseWav);
    end
end
% Single vs. MultiUnit, 2009
if numel(Idx)>0 && (Handles.Para.SingleMultiUnit)
    Par.w_pre = round(Handles.Fs*Handles.Para.AlignTime/1000)+1;
    for i = 1:max(Handles.Unit)
        ID = (Handles.Unit==i)&(~Handles.Outlier);
        CurrUnit = Handles.Unit(ID);
        CurrSpike = Handles.Spike(ID)*1000;
        CurrWav = Handles.Wav(ID,:);
        CurrWav = Aligning(CurrWav, 2);  %aligning waveforms to global minimum, assuming extracellular recording
        [~, Par.w_pre] = min(mean(CurrWav));
        Quality.SingleMultiUnit(i) = eval_class_quality([CurrUnit CurrSpike], Par, CurrWav, 0);
    end
end
% David Kleinfeld 2011
if numel(Idx)>0 && Handles.Para.FalseNegDet
    for i = numel(Idx):-1:1
        VoltThre = Handles.VoltThre(Handles.VoltThre~=0);
        ID = (Handles.Unit==UnitPool(Idx(i)))&(~Handles.Outlier);
        Quality.FalseNegDet(i) = Undetected(Handles.Wav(ID,:), VoltThre, 'auto');
    end
end
if numel(Idx)>0 && Handles.Para.FalsePosRP
    for i = numel(Idx):-1:1
        UnitIdx = (Handles.Unit==UnitPool(Idx(i)))&(~Handles.Outlier);
        NumWav = sum(UnitIdx);
        ISI = zeros(size(Handles.Spike));
        ISI(UnitIdx) = [diff(Handles.Spike(UnitIdx)); 0]*1000;
        InvalidIdx = cumsum(Handles.NumSpikeTrial(Handles.NumSpikeTrial~=0));
        UnitIdx(InvalidIdx) = 0;
        CurrISI = ISI(UnitIdx);
        NumInRP = sum(CurrISI<3);
        RPCensor = 0.003-size(Handles.Wav, 2)/Handles.Fs;
        Quality.FalsePosRP(i) = RPVContamination(NumWav, Handles.ExpDuration, RPCensor, NumInRP);  %two independent neurons
        Quality.FalsePosRP(i) = NumInRP/NumWav;
    end
end
if numel(Idx)>0 && Handles.Para.FalseNegCensor
    for i = numel(Idx):-1:1
        NumWav = sum(Handles.Unit~=UnitPool(Idx(i)));
        CensorT = size(Handles.Wav, 2)/Handles.Fs;
        Quality.FalseNegCensor(i) = Censored(CensorT, NumWav, Handles.ExpDuration);
    end
end
if numel(Idx)>0 && (Handles.Para.FalseNegCluster || Handles.Para.FalsePosCluster)
    if numel(Idx)>=1
        FalseNeg = zeros([numel(Idx) numel(Idx)]);  % Row -- Unit
        FalsePos = zeros([numel(Idx) numel(Idx)]);
        Comb = nchoosek(0:numel(Idx), 2);
        %Wav = zscore(Handles.Wav(~Handles.Outlier,:));
        Wav = zscore(Handles.Coeff(QualityWav,:));
        Unit = Handles.Unit(Handles.WavNum(QualityWav));
        for i = 1:size(Comb, 1)
            if sum(Unit==Comb(i,1))==0 || sum(Unit==Comb(i,2))==0
                Confusion = zeros(2);
            else
                Confusion = GaussianOverlap(Wav(Unit==Comb(i,1),:), Wav(Unit==Comb(i,2),:));
            end
            FalsePos(Comb(i,1)+1,Comb(i,2)+1) = Confusion(1,1);
            FalsePos(Comb(i,2)+1,Comb(i,1)+1) = Confusion(2,2);
            FalseNeg(Comb(i,1)+1,Comb(i,2)+1) = Confusion(1,2);
            FalseNeg(Comb(i,2)+1,Comb(i,1)+1) = Confusion(2,1);
        end
        Quality.FalseNegNoiseCluster = FalseNeg(2:end,1)';
        Quality.FalseNegCluster = transpose(sum(FalseNeg, 2));
        Quality.FalsePosNoiseCluster = FalsePos(2:end,1)';
        Quality.FalsePosCluster = transpose(sum(FalsePos, 2));
        Quality.FalseNegCluster = Quality.FalseNegCluster(2:end);
        Quality.FalsePosCluster = Quality.FalsePosCluster(2:end);
    else
        Quality.FalseNegNoiseCluster = [];
        Quality.FalseNegCluster = [];
        Quality.FalsePosNOiseCluster = [];
        Quality.FalsePosCluster = [];
    end
end
if numel(Idx)==0
    Quality.SNR = [];
    Quality.FalsePosNoise = [];
    Quality.IsolationScore = [];
    Quality.FalseNegDet = [];
    Quality.FalsePosRP = [];
    Quality.FalseNegCensor = [];
    Quality.FalseNegNoiseCluster = [];
    Quality.FalseNegCluster = [];
    Quality.FalsePosNOiseCluster = [];
    Quality.FalsePosCluster = [];
    Quality.FalsePosComposite = [];
    Quality.FalseNegComposite = [];
end
if isfield(Quality, 'FalsePosRP') && isfield(Quality, 'FalsePosCluster') && ...
        isfield(Quality, 'FalseNegDet') && isfield(Quality, 'FalseNegCensor') && ...
        isfield(Quality, 'FalseNegCluster')
    for i = numel(Idx):-1:1
        Quality.FalsePosComposite(i) = max(Quality.FalsePosRP(i), Quality.FalsePosCluster(i));
        Quality.FalseNegComposite(i) = 1-(1-Quality.FalseNegDet(i))*(1-Quality.FalseNegCensor(i))+Quality.FalseNegCluster(i);
    end
end
Handles.Quality = Quality;
disp(Quality);
guidata(Handles.SpikeClusterFig, Handles);

function SummaryPush_Callback(~, ~, Handles)
if ~isfield(Handles, 'Quality')
    Handles = QualityPush_Callback(gcbo, gcbo, Handles);
end
% Clustering, not outlier & transformed spikes
QualityWav = false(size(Handles.Outlier));
QualityWav(Handles.WavNum) = 1;
QualityWav(Handles.Outlier) = 0;
QualityWav = QualityWav(Handles.WavNum);  %sequence in transormed spikes
Wav = Handles.Coeff(QualityWav,:);
Unit = Handles.Unit(Handles.WavNum(QualityWav));
UnitPool = unique(Unit);
UnitPool = UnitPool(UnitPool>0);
figure;
subplot(221);
for i = 1:numel(UnitPool)
    CurrCoeff = Wav(Unit==UnitPool(i),:);
    ClrIdx = mod(UnitPool(i)-1, 5)+2;
    if UnitPool(i)==0
        ClrIdx = 1;
    end
    hold on;
    scatter(CurrCoeff(:,1), CurrCoeff(:,2), 10, Handles.Clr(ClrIdx), 'fill', 'o');
    %APDensity(CurrCoeff(:,1:2), Handles.Clr(ClrIdx));
end
axis square; box off; grid off;
%view([0 0]);
set(gca, 'LineWidth', .5, 'TickLength', [.01 .025]);
xlabel('Coefficient 1'); ylabel('Coefficient 2'); zlabel('Coefficient 3');
title(['Elec' num2str(Handles.ElecNum) ' Clustering'], 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'FontName', 'Arial', 'FontSize', 10);
% Waveform distributions
Wav = Handles.Wav(~logical(Handles.Outlier+Handles.AlignedID),:);
Unit = Handles.Unit(~logical(Handles.Outlier+Handles.AlignedID));
T = (0:size(Wav, 2)-1)/Handles.Fs*1000;
XLim = [min(T)-T(2) max(T)+T(2)];
MaxWav = max(Wav, [], 2);
MinWav = min(Wav, [], 2);
YLim = [mean(MinWav)-7*std(single(MinWav))  mean(MaxWav)+7*std(single(MaxWav))];
for i = 1:max(Unit)  % 0/1/2/3/4 -- Unknown/1/2/3...
    if sum(Unit==i)==0
        continue;
    end
    subplot(2, 2, double(i+1));
%     clf;
    WavDensity(Wav(Unit==i,:), 1/200, Handles.Fs, YLim, Handles.Clr(i+1));
    axis([XLim YLim]);
    axis square; box off; grid off;
    title(['Unit ' num2str(i)], 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'LineWidth', .5, 'TickLength', [.01 .025]);
    xlabel('Time (ms)'); ylabel('Voltage (uv)');
    set(gca, 'FontName', 'Arial', 'FontSize', 10);
%     if i>0
%         XLim = get(gca, 'xlim');
%         YLim = get(gca, 'ylim');
%         text(XLim(2), YLim(1),...
%                {['SNR: ' num2str(Handles.Quality.SNR(i), '%.2f')];...
%                ['Joshua: ' num2str(Handles.Quality.IsolationScore(i), '%.2f') '   '...
%                'SUMU: ' num2str(Handles.Quality.SingleMultiUnit(i), '%.2f')];...
%                ['FPos: ' num2str(Handles.Quality.FalsePosComposite(i), '%.2f') '       '...
%                'FNeg: ' num2str(Handles.Quality.FalseNegComposite(i), '%.2f')];},...
%                'FontName', 'Arial', 'FontSize', 8, ...
%                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
%     end
end

function SavePush_Callback(~, ~, Handles)
if size(Handles.Wav, 1)<Handles.MinNumWav
    return;
end
set(Handles.IndexPop, 'Value', 2);
PlotIndex(Handles);
FileType = {'Waveform';'LFP'};
switch Handles.FileType
    case {1;2;}  % Waveform.mat/LFP.mat
        Handles = QualityPush_Callback(gcbo, gcbo, Handles);
        RawWavPop = get(Handles.WavPop, 'Value');
        set(Handles.WavPop, 'Value', 3);
        PlotAllWav(Handles);
        Quality = Handles.Quality;
        WavMark = strfind(Handles.FileName, FileType{Handles.FileType});
        Spike = Handles.RecordedSpike;
        %Spike(Handles.Rethreshold) = Handles.Spike;
        NumSpikeTrial = Handles.RecordedNumSpikeTrial;
        NumSpike = 0;
        for i = 1:size(Handles.StartT, 1)
%             if ~isnan(Handles.TTLT(i))
%                 Spike(NumSpike+1:NumSpike+NumSpikeTrial(i)) = ...
%                     Spike(NumSpike+1:NumSpike+NumSpikeTrial(i))-Handles.TTLT(i);
%             end
            NumSpike = NumSpike+NumSpikeTrial(i);
        end
        if size(Spike,1)<size(Spike,2)
            Spike = Spike';
        end
        if size(Spike,1) ~= sum(NumSpikeTrial)
            a = 1;
        end
        Spike = mat2cell(Spike, NumSpikeTrial);
        if Handles.FileType==2
            save([Handles.PathName Handles.FileName(1:WavMark(end)-1) 'Spike.mat'], 'Spike', '-v7');
            Waveform = mat2cell(Handles.RecordedWav, NumSpikeTrial);
            save([Handles.PathName Handles.FileName(1:WavMark(end)-1) 'Waveform.mat'], 'Waveform', '-v7');
        end
        Unit = zeros([NumSpike 1])+255;
        Unit(Handles.Rethreshold) = Handles.Unit;
        Unit = mat2cell(uint8(Unit), NumSpikeTrial);
        save([Handles.PathName Handles.FileName(1:WavMark(end)-1) 'Unit.mat'], 'Unit', '-v7');
        print(Handles.SpikeClusterFig, [Handles.PathName Handles.FileName(1:WavMark(end)-1) 'Cluster.png'], '-dpng', '-cmyk');
        save([Handles.PathName Handles.FileName(1:WavMark(end)-1) 'Quality.mat'], 'Quality', '-v7');
        set(Handles.WavPop, 'Value', RawWavPop);
end

function SpikeCluster_GenHistory()
try
    ComputerName = char(java.net.InetAddress.getLocalHost);
catch
    [Status, ComputerName] = system('hostname');
    if Status~=0
        if ispc
            ComputerName = getenv('COMPUTERNAME');
        else
            ComputerName = getenv('HOSTNAME');
        end
    end
end
% ComputerName = strtrim(ComputerName);
% DateTime = clock;
% DateTime = [date ' ' num2str(DateTime(4)) ':' num2str(DateTime(5), '%2u')];
% SCPath = what('Spike Cluster');
% CopyrightPath = [SCPath(1).path '\Copyright.txt'];
% CopyrightPath(CopyrightPath=='\') = '/';
% if ~exist(CopyrightPath, 'file')
%     CopyrightID = fopen(CopyrightPath, 'a');
%     fprintf(CopyrightID, ['--------Statements for SpikeCluster----------\n'...
%         '\nSpikeCluster is a software for extracting and classifying action potentials '...
%         '\nfrom different neurons. This utility is developed by Minggui Chen at Beijing '...
%         '\nNormal University. The first version was released on Mar 20th, 2012. It is '...
%         '\nfor academic use only. Anyone who want to get SpikeCluster should contact '...
%         '\nthe author directly (minggui.chen@gmail.com). Without permission, please do '...
%         '\nNOT distribute the program.\n\nMinggui Chen\n\n\n\nUser history:\n']);
% else
%     CopyrightID = fopen(CopyrightPath, 'a');
% end
% fprintf(CopyrightID, ['\nVer 1.2, initiated by ' ComputerName ' at ' DateTime]);
% fclose(CopyrightID);

function SpikeCluster_OutputFcn(~, ~, Handles)
% cd(Handles.CD);
SpikeClusterFcn(Handles);
