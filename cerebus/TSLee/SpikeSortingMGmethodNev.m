function SpikeSortingMGmethodNev(filepath)
% This function perform the automatic spike sorting using Minggui Chen's
% method.
% filepath is the path of the processed raw data in Feng Wang's format
% It currently will only handle data recorded with Blackrock Cerebus NSP
% and processed by TSLloadNS6ElecNoSorting*

SampleRate = 30000;
if ispc
    slash = '\';
else
    slash = '/';
end

tmppath = [filepath,slash,'tmp'];
if exist(tmppath,'dir') == 7
    rmdir(tmppath,'s');
end
mkdir(tmppath);

RawSpikeStampspath = [filepath,slash,'rawSpikeST.mat'];
RawWaveformpath = [filepath,slash,'rawWaveform.mat'];
electrodePath = [filepath,slash,'electrodes.mat'];
load(RawSpikeStampspath)
load(RawWaveformpath)
load(electrodePath)

NumElec = numel(electrodes);

for thisElec = 1:NumElec
    EID = electrodes(thisElec);
    elecRawSpikeStamps = rawSpikeTimeStamp{EID};
    elecRawWaveform = rawWaveform{EID};
    Spike{1} = elecRawSpikeStamps / SampleRate;
    if size(Spike{1},1)<size(Spike{1},2)
        Spike{1} = Spike{1}';
    end
    if size(elecRawWaveform,1) == numel(Spike{1})
        Waveform{1} = elecRawWaveform;
    elseif size(elecRawWaveform,2) == numel(Spike{1})
        Waveform{1} = elecRawWaveform';
    else
        error('Number of waveforms does not match number of spikes')
    end
    
    
    ExpMonitor.StartT = min(Spike{1})-0.1;
    ExpMonitor.EndT = max(Spike{1})+0.1;
    ExpMonitor.LFPStartT = ExpMonitor.StartT;
    ExpMonitor.LFPEndT = ExpMonitor.EndT;
    SpikeBasic.WaveformFs = SampleRate;
    
    save([tmppath,slash,'Elec1','Spike.mat'],'Spike');
    save([tmppath,slash,'Elec1','Waveform.mat'],'Waveform');
    save([tmppath,slash,'ExpMonitor.mat'],'ExpMonitor');
    save([tmppath,slash,'SpikeBasic.mat'],'SpikeBasic');
    
    disp('Automatically spike sorting')
    FilePath = {
        [tmppath,slash,'Elec1','Waveform.mat'];
        };
    ArgIn.FileType = 3;
    ArgIn.FilePath = FilePath;
    SpikeCluster(ArgIn);
    
    UnitPath = [tmppath,slash,'Elec1Unit.mat'];
    load(UnitPath)
    nevunit{EID} = Unit{1};
end
save([filepath,slash,'nevunit.mat'],'nevunit');
cd ..
rmdir(tmppath,'s')