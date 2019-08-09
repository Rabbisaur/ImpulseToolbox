function ConvertToCDTtable(matfilepath,targetpath)
% load('/home/wangfeng/Downloads/121715_1.mat')

% matfilepath = '/home/wangfeng/storage/projects/SurroundSuppression/ssexp2/matdata/LE20170426_x_-3-3_001';
% targetpath = '/home/wangfeng/Downloads/LE20170426_x_-3-3_001.mat';

% load Expmark
ExpmarkPath = [matfilepath,'/','Expmark.mat'];
load(ExpmarkPath)
NumTrials = size(ExpmarkTime,2);

ElectrodesPath = [matfilepath,'/','electrodes.mat'];
load(ElectrodesPath)
NumElec = numel(electrodes);

condition = ExpmarkTime(2,:)';
starttime = ExpmarkTime(6,:) - ExpmarkTime(6,:);
stoptime = ExpmarkTime(11,:) - ExpmarkTime(6,:);

spikeTimes = cell(NumTrials,1);
spikeUnits = cell(NumTrials,1);
spikeElectrode = cell(NumTrials,1);

% create cache for spike times for each trial
SpikeCache = cell(NumElec,1);
SpikeUnitCache = cell(NumElec,1);
for thisElec = 1:NumElec
    EID = electrodes(thisElec);
    % load elec data
    elecspikepath = [matfilepath,'/','elec',num2str(EID),'/','SpikeTime6.mat'];
    load(elecspikepath);
    SpikeCache{thisElec} = SpikeTimeAligned;
    % load elecunit
    elecunitpath = [matfilepath,'/','elec',num2str(EID),'/','unit.mat'];
    load(elecunitpath)
    spikecount = zeros(NumTrials,1);
    for thisTrial = 1:NumTrials
        spikecount(thisTrial) = numel(SpikeTimeAligned{thisTrial});
    end
    elecunit = mat2cell(elecunit,spikecount,1);
    SpikeUnitCache{thisElec} = elecunit;
end



% populate spikeTimes for each trial
for thisTrial = 1:NumTrials
    % get valid electrodes
    ValidElecIndex = true(NumElec,1);
    for thisElec = 1:NumElec
        if isempty(SpikeCache{thisElec}{thisTrial})
            ValidElecIndex(thisElec) = false;
        end
    end
    EID = electrodes(ValidElecIndex)';
    spikeElectrode{thisTrial} = EID;
    spikeTimes{thisTrial} = cell(sum(ValidElecIndex),1);
    spikeUnits{thisTrial} = cell(sum(ValidElecIndex),1);
    ValidElecIndex = find(ValidElecIndex);
    for thisElec = 1:numel(ValidElecIndex)
        spikeTimes{thisTrial}{thisElec} = SpikeCache{ValidElecIndex(thisElec)}{thisTrial};
        spikeUnits{thisTrial}{thisElec} = SpikeUnitCache{ValidElecIndex(thisElec)}{thisTrial};
    end
end

CDTTables{1}.condition = condition;
CDTTables{1}.starttime = starttime';
CDTTables{1}.stoptime = stoptime';
CDTTables{1}.spikeElectrode = spikeElectrode;
CDTTables{1}.spikeTimes = spikeTimes;
CDTTables{1}.spikeUnits = spikeUnits;
CDTTables{1}.enentCodes = [];
CDTTables{1}.eventtimes = [];

save(targetpath,'CDTTables');