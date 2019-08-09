function nevdata = TSLSpikeSortingMGmethodNEV(nevdata,filepath)
% This function perform the automatic spike sorting using Minggui Chen's
% method.
% filepath is the path of the processed raw data in Feng Wang's format
% It currently will only handle data recorded with Blackrock Cerebus NSP
% and processed by TSLloadNS6ElecNoSorting*

samplerate = nevdata.samplerate;
electrodes = nevdata.electrodes;
numElec = numel(electrodes);

for thisElec = 1:numElec
    disp(['Spike sorting electrode ', num2str(thisElec),'/',num2str(numElec)])
    EID = electrodes(thisElec);
    tmppath = [filepath,'/tmp'];
    if exist(tmppath,'dir')
        % do nothing
    else
        mkdir(tmppath);
    end
    
    Spike{1} = nevdata.rawspikeTimeStamp{EID}/samplerate;
    if size(Spike{1},1)<size(Spike{1},2)
        Spike{1} = Spike{1}';
    end

    if size(nevdata.rawWaveform{EID},1) == numel(Spike{1})
        Waveform{1} = nevdata.rawWaveform{EID};
    elseif size(nevdata.rawWaveform{EID},2) == numel(Spike{1})
        Waveform{1} = nevdata.rawWaveform{EID}';
    else
        error('Number of waveforms does not match number of spikes')
    end
    
    
    ExpMonitor.StartT = min(Spike{1})-0.1;
    ExpMonitor.EndT = max(Spike{1})+0.1;
    ExpMonitor.LFPStartT = ExpMonitor.StartT;
    ExpMonitor.LFPEndT = ExpMonitor.EndT;
    SpikeBasic.WaveformFs = samplerate;
    
    save([tmppath,'/','Elec1','Spike.mat'],'Spike');
    save([tmppath,'/','Elec1','Waveform.mat'],'Waveform');
    save([tmppath,'/','ExpMonitor.mat'],'ExpMonitor');
    save([tmppath,'/','SpikeBasic.mat'],'SpikeBasic');
    
    disp('Automatically spike sorting')
    FilePath = {
        [tmppath,'/','Elec1','Waveform.mat'];
        };
    ArgIn.FileType = 3;
    ArgIn.FilePath = FilePath;
    SpikeCluster(ArgIn);
    
    UnitPath = [tmppath,'/','Elec1Unit.mat'];
    load(UnitPath)
    nevdata.SpikeUnit{EID} = Unit{1};

    QualityPath = [tmppath,'/','Elec1Quality.mat'];
    load(QualityPath)
    nevdata.SpikeSoringQuality{EID} = Quality;
    
    cd ..
    rmdir(tmppath,'s')
end