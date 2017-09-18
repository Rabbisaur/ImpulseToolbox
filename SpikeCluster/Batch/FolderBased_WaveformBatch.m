% Batch file for SpikeCluster;
% This file must be named in *WaveformBatch.m;
% FilePath -- string cell indicating the full path of waveform.mat;

% Edit file paths here ***************************************
SCPath = what('SpikeCluster');
SCPath = [SCPath.path '\SampleData'];
FilePath = GetPath(SCPath, 'Waveform.mat');
SortedID = [];
for i = 1:numel(FilePath)
    Mark = strfind(FilePath{i}, 'Waveform.mat');
    if exist([FilePath{i}(1:Mark(end)-1) 'Cluster.png'], 'file')==2
        SortedID = [SortedID i];
    end
end
% FilePath(SortedID) = [];