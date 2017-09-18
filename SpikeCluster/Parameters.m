function varargout = Parameters(varargin)
% Initialization, DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @Parameters_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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
function Parameters_OpeningFcn(~, ~, Handles, varargin)
Pos = get(gcf, 'OuterPosition');
ScrSize = get(0, 'screensize');
Pos(1) = ScrSize(3)/2-Pos(3)/2;
Pos(2) = ScrSize(4)/2-Pos(4)/2;
set(gcf, 'OuterPosition', Pos);
% Set raw parameters
if nargin==4
    Handles.RawPara = varargin{1};
else
    ParaPath = what('SpikeCluster');
    ParaPath(ParaPath=='\') = '/';
    ParaPath = [ParaPath.path '/DefaultPara.mat'];
    load(ParaPath);
    Handles.RawPara = Para;
end
Handles.RawPara.Load = 1;
DisplayPara(Handles, Handles.RawPara);
Handles.Para = Handles.RawPara;
switch get(Handles.ExtractMethodPop, 'Value')
    case {4;5;}  %wavelet/wavelet packets
        set(Handles.WaveletTypePop, 'String', Handles.Para.WaveletTypeOpt);
    case {1;2;3;6;}  %multiwavelet
        set(Handles.WaveletTypePop, 'String', Handles.Para.MultiwaveletTypeOpt);
end
switch get(Handles.ClusterMethodPop, 'Value')
    case {1;4;5;}  %single linkage/t-dist EM/SPC
        set(Handles.NumClusterPop, 'String', Handles.Para.NumClusterOpt(1));
        Handles.Para.NumCluster = 1;
    case {2;3;}  %k-mean/EM
        set(Handles.NumClusterPop, 'String', Handles.Para.NumClusterOpt(2:end));
        Handles.Para.NumCluster = max(2, Handles.Para.NumCluster);
end
DisplayPara(Handles, Handles.Para);
guidata(gcf, Handles);
uiwait(Handles.ParametersFig);

% Load parameters
function LoadToggle_Callback(~, ~, Handles)
set(gcbo, 'Value', 1);
Handles.Para.Load = get(gcbo, 'Value');
guidata(gcbo, Handles);
function ElecNumEdit_Callback(~, ~, Handles)
Handles.Para.ElecNum = eval(get(gcbo, 'string'));
guidata(gcbo, Handles);
function FsEdit_Callback(~, ~, Handles)
Handles.Para.Fs = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);

% Filter parameters
function FilterToggle_Callback(~, ~, Handles)
Handles.Para.Filter = get(gcbo, 'Value');
guidata(gcbo, Handles);
function HipassFreqEdit_Callback(~, ~, Handles)
Handles.Para.HipassFreq = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function LopassFreqEdit_Callback(~, ~, Handles)
Handles.Para.LopassFreq = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function FilterDirectionPop_Callback(~, ~, Handles)
Handles.Para.FilterDirection = get(gcbo, 'Value');
guidata(gcbo, Handles);
function HipassTypePop_Callback(~, ~, Handles)
Handles.Para.HipassType = get(gcbo, 'Value');
guidata(gcbo, Handles);
function HipassOrderEdit_Callback(~, ~, Handles)
Handles.Para.HipassOrder = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function LopassTypePop_Callback(~, ~, Handles)
Handles.Para.LopassType = get(gcbo, 'Value');
guidata(gcbo, Handles);
function LopassOrderEdit_Callback(~, ~, Handles)
Handles.Para.LopassOrder = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);

% Detection parameters
function DetToggle_Callback(~, ~, Handles)
Handles.Para.Det = get(gcbo, 'Value');
guidata(gcbo, Handles);
function DetMethodPop_Callback(~, ~, Handles)
Handles.Para.DetMethod = get(gcbo, 'Value');
guidata(gcbo, Handles);
function SDMethodPop_Callback(~, ~, Handles)
Handles.Para.SDMethod = get(gcbo, 'Value');
guidata(gcbo, Handles);
function SDDurationEdit_Callback(~, ~, Handles)
Handles.Para.SDDuration = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function DetThreEdit_Callback(~, ~, Handles)
Handles.Para.DetThre = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function DetPolarPop_Callback(~, ~, Handles)
Handles.Para.DetPolar = get(gcbo, 'Value');
guidata(gcbo, Handles);
function SpikeDurationEdit_Callback(~, ~, Handles)
Handles.Para.SpikeDuration = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function RefractoryPeriodEdit_Callback(~, ~, Handles)
Handles.Para.RefractoryPeriod = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);

% Alignment parameters
function AlignToggle_Callback(~, ~, Handles)
Handles.Para.Align = get(gcbo, 'Value');
guidata(gcbo, Handles);
function ReverseFilterWavPop_Callback(~, ~, Handles)
Handles.Para.ReverseFilterWav = get(gcbo, 'Value');
guidata(gcbo, Handles);
function InterpMethodPop_Callback(~, ~, Handles)
Handles.Para.InterpMethod = get(gcbo, 'Value');
guidata(gcbo, Handles);
function InterpResoEdit_Callback(~, ~, Handles)
Handles.Para.InterpReso = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function AlignMethodPop_Callback(~, ~, Handles)
Handles.Para.AlignMethod = get(gcbo, 'Value');
guidata(gcbo, Handles);
function AlignTimeEdit_Callback(~, ~, Handles)
Handles.Para.AlignTime = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);

% Extraction parameters
function ExtractToggle_Callback(~, ~, Handles)
Handles.Para.Extract = get(gcbo, 'Value');
guidata(gcbo, Handles);
function NumClusteredSpikeEdit_Callback(~, ~, Handles)
Handles.Para.NumClusteredSpike = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function ExtractMethodPop_Callback(~, ~, Handles)
Handles.Para.ExtractMethod = get(gcbo, 'Value');
switch get(Handles.ExtractMethodPop, 'Value')
    case {4;5;}  %wavelet/wavelet packets
        set(Handles.WaveletTypePop, 'String', Handles.Para.WaveletTypeOpt);
    case {1;2;3;6;}  %multiwavelet
        set(Handles.WaveletTypePop, 'String', Handles.Para.MultiwaveletTypeOpt);
end
guidata(gcbo, Handles);
function WaveletTypePop_Callback(~, ~, Handles)
Handles.Para.WaveletType = get(gcbo, 'Value');
guidata(gcbo, Handles);
function WaveletScaleEdit_Callback(~, ~, Handles)
Handles.Para.WaveletScale = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);

%Selection parameters
function SelectToggle_Callback(~, ~, Handles)
Handles.Para.Select = get(gcbo, 'Value');
guidata(gcbo, Handles);
function SelectMethodPop_Callback(~, ~, Handles)
Handles.Para.SelectMethod = get(gcbo, 'Value');
guidata(gcbo, Handles);
function VarRatioEdit_Callback(~, ~, Handles)
Handles.Para.VarRatio = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function NumCoeffEdit_Callback(~, ~, Handles)
Handles.Para.NumCoeff = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);

% Cluster parameters
function ClusterToggle_Callback(~, ~, Handles)
Handles.Para.Cluster = get(gcbo, 'Value');
guidata(gcbo, Handles);
function OutlierRemovalPop_Callback(~, ~, Handles)
Handles.Para.OutlierRemoval = get(gcbo, 'Value');
guidata(gcbo, Handles);
function ClusterMethodPop_Callback(~, ~, Handles)
Handles.Para.ClusterMethod = get(gcbo, 'Value');
switch get(Handles.ClusterMethodPop, 'Value')
    case {1;4;5;}  %single linkage/t-dist EM/SPC
        set(Handles.NumClusterPop, 'String', Handles.Para.NumClusterOpt(1));
        Handles.Para.NumCluster = 1;
    case {2;3;}  %k-mean/EM
        set(Handles.NumClusterPop, 'String', Handles.Para.NumClusterOpt(2:end));
        Handles.Para.NumCluster = max(2, Handles.Para.NumCluster);
end
DisplayPara(Handles, Handles.Para);
guidata(gcbo, Handles);
function NumClusterPop_Callback(~, ~, Handles)
Handles.Para.NumCluster = get(gcbo, 'Value');
guidata(gcbo, Handles);
function MinClusterRatioEdit_Callback(~, ~, Handles)
Handles.Para.MinClusterRatio = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function AutoTemplatePop_Callback(~, ~, Handles)
Handles.Para.AutoTemplate = get(gcbo, 'Value');
guidata(gcbo, Handles);
function AutotempSNRThreEdit_Callback(~, ~, Handles)
Handles.Para.AutotempSNRThre = eval(get(gcbo, 'String'));
guidata(gcbo, Handles);
function AssignRemainderPop_Callback(~, ~, Handles)
Handles.Para.AssignRemainder = get(gcbo, 'Value');
guidata(gcbo, Handles);
function AssignOutlierPop_Callback(~, ~, Handles)
Handles.Para.AssignOutlier = get(gcbo, 'Value');
guidata(gcbo, Handles);

% Hidden Markov Model parameters
function HMMToggle_Callback(~, ~, Handles)
Handles.Para.HMM = get(gcbo, 'Value');
guidata(gcbo, Handles);

% Quality Metrics parameters
function SNRChk_Callback(~, ~, Handles)
Handles.Para.SNR = get(gcbo, 'Value');
guidata(gcbo, Handles);
function IsolationScoreChk_Callback(~, ~, Handles)
Handles.Para.IsolationScore = get(gcbo, 'Value');
guidata(gcbo, Handles);
function FalseNegChk_Callback(~, ~, Handles)
Handles.Para.FalseNeg = get(gcbo, 'Value');
guidata(gcbo, Handles);
function FalsePosChk_Callback(~, ~, Handles)
Handles.Para.FalsePos = get(gcbo, 'Value');
guidata(gcbo, Handles);
function SingleMultiUnitChk_Callback(~, ~, Handles)
Handles.Para.SingleMultiUnit = get(gcbo, 'Value');
guidata(gcbo, Handles);
function FalseNegDetChk_Callback(~, ~, Handles)
Handles.Para.FalseNegDet = get(gcbo, 'Value');
guidata(gcbo, Handles);
function FalsePosRPChk_Callback(~, ~, Handles)
Handles.Para.FalsePosRP = get(gcbo, 'Value');
guidata(gcbo, Handles);
function FalseNegClusterChk_Callback(~, ~, Handles)
Handles.Para.FalseNegCluster = get(gcbo, 'Value');
guidata(gcbo, Handles);
function FalsePosClusterChk_Callback(~, ~, Handles)
Handles.Para.FalsePosCluster = get(gcbo, 'Value');
guidata(gcbo, Handles);
function FalseNegCensorChk_Callback(~, ~, Handles)
Handles.Para.FalseNegCensor = get(gcbo, 'Value');
guidata(gcbo, Handles);
function LRatioChk_Callback(~, ~, Handles)
Handles.Para.LRatio = get(gcbo, 'Value');
guidata(gcbo, Handles);
function IsolationDistChk_Callback(~, ~, Handles)
Handles.Para.IsolationDist = get(gcbo, 'Value');
guidata(gcbo, Handles);
function IsolationInfoChk_Callback(~, ~, Handles)
Handles.Para.IsolationInfo = get(gcbo, 'Value');
guidata(gcbo, Handles);

% Ouput parameters & close Paramters.fig
function OKPush_Callback(~, ~, Handles)
uiresume(Handles.ParametersFig);

% Use FVTool to check filter settings
function VisFilterPush_Callback(~, ~, Handles)
[~, HipassFilter, LopassFilter] = Filtering([], Handles.Para);
if ~isempty(HipassFilter) && ~isempty(LopassFilter)
    HVisFilter = fvtool(HipassFilter.B, HipassFilter.A, LopassFilter.B, LopassFilter.A);
    legend(HVisFilter, 'High-pass', 'Low-pass');
elseif ~isempty(HipassFilter)
    HVisFilter = fvtool(HipassFilter.B, HipassFilter.A);
    legend(HVisFilter, 'High-pass');
elseif ~isempty(LopassFilter)
    HVisFilter = fvtool(LopassFilter.B, LopassFilter.A);
    legend(HVisFilter, 'Low-pass');
end
if exist('HVisFilter', 'var')
    set(HVisFilter, 'analysis', 'freq', 'Fs', Handles.Para.Fs, 'PhaseUnits', 'degrees', 'Legend', 'on');
    set(HVisFilter, 'FrequencyScale', 'log');
end

% Restore parameters to default
function RestoreDefaultPush_Callback(~, ~, Handles)
ParaPath = what('Spike Cluster');
ParaPath = [ParaPath.path '/DefaultPara.mat'];
load(ParaPath);
Handles.Para = Para;
DisplayPara(Handles, Handles.Para);
switch get(Handles.ExtractMethodPop, 'Value')
    case {4;5;}  %wavelet/wavelet packets
        set(Handles.WaveletTypePop, 'String', Handles.Para.WaveletTypeOpt);
    case {1;2;3;6;}  %multiwavelet
        set(Handles.WaveletTypePop, 'String', Handles.Para.MultiwaveletTypeOpt);
end
switch get(Handles.ClusterMethodPop, 'Value')
    case {1;4;5;}  %single linkage/t-dist EM/SPC
        set(Handles.NumClusterPop, 'String', Handles.Para.NumClusterOpt(1));
        Handles.Para.NumCluster = 1;
    case {2;3;}  %k-mean/EM
        set(Handles.NumClusterPop, 'String', Handles.Para.NumClusterOpt(2:end));
        Handles.Para.NumCluster = max(2, Handles.Para.NumCluster);
end
DisplayPara(Handles, Handles.Para);
guidata(gcbo, Handles);

% Save current parameters as default
function SaveDefaultPush_Callback(~, ~, Handles)
ParaPath = what('Spike Cluster');
ParaPath = [ParaPath.path '/DefaultPara.mat'];
Para = Handles.Para;
save(ParaPath, 'Para', '-v7');

% Undo all the changes & Close Parameters.fig
function CancelPush_Callback(~, ~, Handles)
DisplayPara(Handles, Handles.RawPara);
Handles.Para = Handles.RawPara;
switch get(Handles.ExtractMethodPop, 'Value')
    case {4;5;}  %wavelet/wavelet packets
        set(Handles.WaveletTypePop, 'String', Handles.Para.WaveletTypeOpt);
    case {1;2;3;6;}  %multiwavelet
        set(Handles.WaveletTypePop, 'String', Handles.Para.MultiwaveletTypeOpt);
end
switch get(Handles.ClusterMethodPop, 'Value')
    case {1;4;5;}  %single linkage/t-dist EM/SPC
        set(Handles.NumClusterPop, 'String', Handles.Para.NumClusterOpt(1));
        Handles.Para.NumCluster = 1;
    case {2;3;}  %k-mean/EM
        set(Handles.NumClusterPop, 'String', Handles.Para.NumClusterOpt(2:end));
        Handles.Para.NumCluster = max(2, Handles.Para.NumCluster);
end
DisplayPara(Handles, Handles.Para);
guidata(gcbo, Handles);
uiresume(Handles.ParametersFig);

% Set all parameters according to input
function DisplayPara(Handles, Para)
set(Handles.LoadToggle, 'Value', Para.Load);
set(Handles.ElecNumEdit, 'String', num2str(Para.ElecNum));
set(Handles.FsEdit, 'String', Para.Fs);
set(Handles.FilterToggle, 'Value', Para.Filter);
set(Handles.HipassFreqEdit, 'String', num2str(Para.HipassFreq));
set(Handles.LopassFreqEdit, 'String', num2str(Para.LopassFreq));
set(Handles.FilterDirectionPop, 'Value', Para.FilterDirection);
set(Handles.HipassTypePop, 'Value', Para.HipassType);
set(Handles.HipassOrderEdit, 'String', num2str(Para.HipassOrder));
set(Handles.LopassTypePop, 'Value', Para.LopassType);
set(Handles.LopassOrderEdit, 'String', num2str(Para.LopassOrder));
set(Handles.DetToggle, 'Value', Para.Det);
set(Handles.DetMethodPop, 'Value', Para.DetMethod);
set(Handles.SDMethodPop, 'Value', Para.SDMethod);
set(Handles.SDDurationEdit, 'String', num2str(Para.SDDuration));
set(Handles.DetThreEdit, 'String', num2str(Para.DetThre));
set(Handles.DetPolarPop, 'Value', Para.DetPolar);
set(Handles.SpikeDurationEdit, 'String', Para.SpikeDuration);
set(Handles.RefractoryPeriodEdit, 'String', Para.RefractoryPeriod);
set(Handles.AlignToggle, 'Value', Para.Align);
set(Handles.ReverseFilterWavPop, 'Value', Para.ReverseFilterWav);
set(Handles.InterpMethodPop, 'Value', Para.InterpMethod);
set(Handles.InterpResoEdit, 'String', Para.InterpReso);
set(Handles.AlignMethodPop, 'Value', Para.AlignMethod);
set(Handles.AlignTimeEdit, 'String', num2str(Para.AlignTime));
set(Handles.ExtractToggle, 'Value', Para.Extract);
set(Handles.NumClusteredSpikeEdit, 'String', Para.NumClusteredSpike);
set(Handles.ExtractMethodPop, 'Value', Para.ExtractMethod);
set(Handles.WaveletTypePop, 'Value', Para.WaveletType);
set(Handles.WaveletScaleEdit, 'String', Para.WaveletScale);
set(Handles.SelectToggle, 'Value', Para.Select);
set(Handles.SelectMethodPop, 'Value', Para.SelectMethod);
set(Handles.VarRatioEdit, 'String', Para.VarRatio);
set(Handles.NumCoeffEdit, 'String', Para.NumCoeff);
set(Handles.ClusterToggle, 'Value', Para.Cluster);
set(Handles.OutlierRemovalPop, 'Value', Para.OutlierRemoval);
set(Handles.ClusterMethodPop, 'Value', Para.ClusterMethod);
set(Handles.NumClusterPop, 'Value', Para.NumCluster);
set(Handles.MinClusterRatioEdit, 'String', Para.MinClusterRatio);
set(Handles.AutoTemplatePop, 'Value', Para.AutoTemplate);
set(Handles.AutotempSNRThreEdit, 'String', Para.AutotempSNRThre);
set(Handles.AssignRemainderPop, 'Value', Para.AssignRemainder);
set(Handles.AssignOutlierPop, 'Value', Para.AssignOutlier);
set(Handles.HMMToggle, 'Value', Para.HMM);
set(Handles.SNRChk, 'Value', Para.SNR);
set(Handles.IsolationScoreChk, 'Value', Para.IsolationScore);
set(Handles.FalseNegChk, 'Value', Para.FalseNeg);
set(Handles.FalsePosChk, 'Value', Para.FalsePos);
set(Handles.SingleMultiUnitChk, 'Value', Para.SingleMultiUnit);
set(Handles.FalseNegDetChk, 'Value', Para.FalseNegDet);
set(Handles.FalsePosRPChk, 'Value', Para.FalsePosRP);
set(Handles.FalseNegClusterChk, 'Value', Para.FalseNegCluster);
set(Handles.FalsePosClusterChk, 'Value', Para.FalsePosCluster);
set(Handles.FalseNegCensorChk, 'Value', Para.FalseNegCensor);
set(Handles.HipassTypePop, 'String', Para.HipassOpt);
set(Handles.LopassTypePop, 'String', Para.LopassOpt);
set(Handles.DetMethodPop, 'String', Para.DetMethodOpt);
set(Handles.SDMethodPop, 'String', Para.SDMethodOpt);
set(Handles.DetPolarPop, 'String', Para.DetPolarOpt);
set(Handles.InterpMethodPop, 'String', Para.InterpMethodOpt);
set(Handles.AlignMethodPop, 'String', Para.AlignMethodOpt);
set(Handles.ExtractMethodPop, 'String', Para.ExtractMethodOpt);
set(Handles.LRatioChk, 'Value', Para.LRatio);
set(Handles.IsolationDistChk, 'Value', Para.IsolationDist);
set(Handles.IsolationInfoChk, 'Value', Para.IsolationInfo);
switch get(Handles.ExtractMethodPop, 'Value')
    case {4;5;}  %wavelet/wavelet packets
        set(Handles.WaveletTypePop, 'String', Para.WaveletTypeOpt);
    case {1;2;3;6;}  %multiwavelet
        set(Handles.WaveletTypePop, 'String', Para.MultiwaveletTypeOpt);
end
set(Handles.SelectMethodPop, 'String', Para.SelectMethodOpt);
set(Handles.OutlierRemovalPop, 'String', Para.OutlierRemovalOpt);
set(Handles.ClusterMethodPop, 'String', Para.ClusterMethodOpt);
switch get(Handles.ClusterMethodPop, 'Value')
    case {1;4;5;6;}  %single linkage/t-dist EM/SPC/Auto-EM
        set(Handles.NumClusterPop, 'String', Para.NumClusterOpt(1));
        Handles.Para.NumCluster = 1;
    case {2;3;}  %k-means/EM
        set(Handles.NumClusterPop, 'String', Para.NumClusterOpt(2:end));
        Handles.Para.NumCluster = max(2, Para.NumCluster);
end
set(Handles.AutoTemplatePop, 'String', Para.AutoTemplateOpt);
set(Handles.AssignRemainderPop, 'String', Para.AssignRemainderOpt);
set(Handles.AssignOutlierPop, 'String', Para.AssignOutlierOpt);

% Output
function varargout = Parameters_OutputFcn(~, ~, Handles)
varargout{1} = Handles.Para;
close(Handles.ParametersFig);
