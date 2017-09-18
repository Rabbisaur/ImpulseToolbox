% Batch file for SpikeCluster;
% This file must be named in *LFPFileBatch.m;
% FilePath -- string cell indicating the full path;
% ElecNum -- numeric cell indicating the electrode

% Edit file paths & electrodes here ***************************************
FilePath = {
    'D:\Research\Program\Included\NSD\SpikeCluster\SampleData.nev';
    'D:\Research\Program\Included\NSD\SpikeCluster\SampleData.nev';
};
ElecNum = {
    [12 12];
    [12 12];
};
% *************************************************************************

% Correct the full path of waveforms, Do not edit below!!!
for i = 1:numel(FilePath)
    Mark = strfind(lower(FilePath{i}), 'machine');
    if ~isempty(Mark)
        FilePath{i} = [FilePath{i}(1:Mark-1) 'Mat' FilePath{i}(Mark+7:end)];
    end
end
assert(numel(FilePath)==numel(ElecNum), 'FilePath and ElecNum must be in same size.');