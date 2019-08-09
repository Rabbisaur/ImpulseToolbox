function TSLSpikeSortingMGmethodNS6(instanceinfo)
% This function perform the automatic spike sorting using Minggui Chen's
% method.
% filepath is the path of the processed raw data in Feng Wang's format
% It currently will only handle data recorded with Blackrock Cerebus NSP
% and processed by TSLloadNS6ElecNoSorting*

samplerate = instanceinfo.samplerate;
numElec = instanceinfo.numElec;

for thisElec = 1:numElec
    EID = instanceinfo.ElecOrder(thisElec);
    tmppath = instanceinfo.electrodeCachePath{EID};
    idx = strfind(tmppath,'tmp');
    idx = idx(1);
    tmppath = tmppath(1:(idx+2));
    if exist(tmppath,'dir')
        % do nothing
    else
        mkdir(tmppath);
    end
    
    load(instanceinfo.electrodeNeuraldataPath{EID})
    
    Spike{1} = electrodeNeuraldata.SpikeSTraw/samplerate;
    if size(Spike{1},1)<size(Spike{1},2)
        Spike{1} = Spike{1}';
    end

    if size(electrodeNeuraldata.SpikeWaveform,1) == numel(Spike{1})
        Waveform{1} = electrodeNeuraldata.SpikeWaveform;
    elseif size(electrodeNeuraldata.SpikeWaveform,2) == numel(Spike{1})
        Waveform{1} = electrodeNeuraldata.SpikeWaveform';
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
    electrodeNeuraldata.SpikeUnit = Unit{1};

    QualityPath = [tmppath,'/','Elec1Quality.mat'];
    load(QualityPath)
    electrodeNeuraldata.SpikeSoringQuality = Quality;
    
    save(instanceinfo.electrodeNeuraldataPath{EID},'electrodeNeuraldata')
    
    cd ..
    rmdir(tmppath,'s')
end