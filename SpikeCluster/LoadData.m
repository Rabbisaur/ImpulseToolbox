function Handles = LoadData(Handles)
% Load data and save them
%   ContSignal: ContSignal.mat as UserData in Handle.ContSignalAxes
%   Waveforms: Wav.mat as UserData in Handle.WavAxes
%   Spike Time: Spike.mat as UserData in Handle.Wav1Axes
switch Handles.FileType
    case 1  % *Waveform.mat
        ElecMark = strfind(Handles.FileName, 'Elec');
        WavMark = strfind(Handles.FileName, 'Waveform');
        Handles.ElecNum = str2double(Handles.FileName(ElecMark(end)+4:WavMark(end)-1));
        Handles.Para.ElecNum = Handles.ElecNum;
        load([Handles.PathName Handles.FileName]);
        for i = numel(Waveform):-1:1
            Handles.NumSpikeTrial(i,1) = size(Waveform{i}, 1);
        end
        Handles.Wav = single(cell2mat(Waveform));
        load([Handles.PathName Handles.FileName(1:WavMark(end)-1) 'Spike.mat']);
        Handles.Spike = cell2mat(Spike);
        load([Handles.PathName 'ExpMonitor.mat']);
%         Handles.TTLT = ExpMonitor.TTLT(:,1);
        Handles.StartT = ExpMonitor.StartT;
        NumSpike = 0;
        for i = 1:size(Handles.StartT, 1)
%             if ~isnan(Handles.TTLT(i))
%                 Handles.Spike(NumSpike+1:NumSpike+Handles.NumSpikeTrial(i)) = ...
%                     Handles.Spike(NumSpike+1:NumSpike+Handles.NumSpikeTrial(i))+Handles.TTLT(i);
%             end
            NumSpike = NumSpike+Handles.NumSpikeTrial(i);
        end
        Handles.RawSpike = Handles.Spike;  %for outlier
        ExpDuration = ExpMonitor.EndT-ExpMonitor.StartT;
        Handles.ExpDuration = sum(ExpDuration(~isnan(ExpDuration)));
        load([Handles.PathName 'SpikeBasic.mat']);
        Handles.Fs = SpikeBasic.WaveformFs;
        MuWav = mean(Handles.Wav, 1);
%         load([Handles.PathName 'SpikeExtend.mat']);
%         for i = 1:numel(SpikeExtend)
%             if strcmpi(SpikeExtend{i,1}.PacketID, 'NEUEVWAV') && SpikeExtend{i,1}.ElecID==Handles.ElecNum
%                 Handles.VoltThre = [SpikeExtend{i,1}.EnergyThre SpikeExtend{i,1}.LowThre SpikeExtend{i,1}.HighThre];
%             end
%         end
%         AlignID = find( abs(diff(MuWav<=Handles.VoltThre(Handles.VoltThre~=0))), 1, 'first')+1;
        if ~exist('AlignID', 'var') || isempty(AlignID)
            [~, AlignID] = min(MuWav);
        end
        Handles.VoltThre = max(Handles.Wav(:,AlignID));
        Handles.Para.AlignTime = (AlignID-1)/Handles.Fs*1000;
        Handles.Para.Filter = 0;
    case 2  % *LFP.mat; could be better if sampling rate is larger than 20K, as many spikes jump to peak in .1-.2 ms
        ElecMark = strfind(Handles.FileName, 'Elec');
        WavMark = strfind(Handles.FileName, 'LFP');
        Handles.ElecNum = str2double(Handles.FileName(ElecMark(end)+4:WavMark(end)-1));
        Handles.Para.ElecNum = Handles.ElecNum;
        load([Handles.PathName Handles.FileName]);
        Handles.WBSignal = LFP;
        load([Handles.PathName 'LFPBasic.mat']);
        Handles.Fs = LFPBasic.Fs;
        load([Handles.PathName 'ExpMonitor.mat']);
        Handles.StartT = ExpMonitor.LFPStartT;
        Handles.EndT = ExpMonitor.LFPEndT;
%         Handles.TTLT = ExpMonitor.TTLT(:,1);
        ExpDuration = ExpMonitor.EndT-ExpMonitor.StartT;
        Handles.ExpDuration = sum(ExpDuration(~isnan(ExpDuration)));
        Handles.Para.Filter = 1;
        Handles.Para.Det = 1;
    otherwise
        errordlg([FilePath ' is not supported.']);
end
Handles.Para.Fs = Handles.Fs;