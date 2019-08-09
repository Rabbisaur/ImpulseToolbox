function PSTH = TSLgetPSTHelec(sessionpath,edges,UnitSelect,unitfilename)
% load data
if ispc
    slash = '\';
else
    slash = '/';
end

% % load spikes
% SpikePath = [sessionpath,slash,'SpikeTime6.mat'];
% load(SpikePath);
% % load unit
% UnitPath = [sessionpath,slash,'SpikeUnits6.mat'];
% load(UnitPath);

% load Expmark
ExpmarkPath = [sessionpath,slash,'Expmark.mat'];
load(ExpmarkPath);

% load elec
elecpath = [sessionpath, slash, 'electrodes.mat'];
load(elecpath);

electrodes(electrodes>128) = [];

NumElec = numel(electrodes);

CondID = ExpmarkTime(2,:);
UniCondID = unique(CondID);
NumConditions = numel(UniCondID);
NumTrials = size(ExpmarkTime,2);
NumRepeats = sum(CondID==UniCondID(1));

PSTH = zeros(numel(edges)-1,NumRepeats,NumConditions,NumElec);
edges = repmat(edges,NumTrials,1);
Cedges = mat2cell(edges,ones(size(edges,1),1),size(edges,2));
SpikeCount = zeros(NumRepeats,NumConditions,NumElec);
for thisElec = 1:NumElec
    EID = electrodes(thisElec);
    elecpath = [sessionpath,slash,'elec',num2str(EID)];
    elecspikepath = [elecpath,slash,'SpikeTime6.mat'];
    if UnitSelect ~= 'n'
        unitspikepath = [elecpath,slash,unitfilename];
        load(unitspikepath)
        switch unitfilename
            case 'userunit.mat'
                elecUnit = userelecunit;
            case 'unit.mat'
                elecUnit = elecunit;
            otherwise
                error('Unrecognized spike unit file name.')
        end
        if size(elecUnit,1) < size(elecUnit,2)
            elecUnit = elecUnit';
        end
        
    end
    
    load(elecspikepath)
    
    
    elecSpike = SpikeTimeAligned;
    for thisTrial = 1:NumTrials
        NumSpikesinTrial(thisTrial) = numel(elecSpike{thisTrial});
    end
    
    if UnitSelect ~= 'n'
        elecUnit = mat2cell(elecUnit,NumSpikesinTrial,1);
    end
    
    
    
    
    for thisTrial = 1:NumTrials
        % Unit select
        switch UnitSelect
            case '0'
                elecSpike{thisTrial} = elecSpike{thisTrial}(elecUnit{thisTrial} == 0);
            case '1'
                elecSpike{thisTrial} = elecSpike{thisTrial}(elecUnit{thisTrial} == 1);
            case ['a','A']
                % do nothing
        end
    end
    % calculate PSTH
    elecPSTH = cellfun(@histc,elecSpike,Cedges,'UniformOutput',false);
    % correct vector direction
    for thisTrial = 1:NumTrials
        if size(elecPSTH{thisTrial},1) < size(elecPSTH{thisTrial},2)
            elecPSTH{thisTrial} = elecPSTH{thisTrial}';
        end
    end
    elecPSTH = cell2mat(elecPSTH');
    %     thisElec
    %     size(elecPSTH)
    % some trials may not have spikes
    if size(elecPSTH,2) ~=NumTrials
        elecPSTH = [elecPSTH,zeros(size(elecPSTH,1),NumTrials-size(elecPSTH,2))];
    end
    elecPSTH(end,:) = [];
    for thisCond = 1:NumConditions
        idx = CondID == UniCondID(thisCond);
        PSTH(:,:,thisCond,thisElec) = elecPSTH(:,idx);
        trials = elecSpike(idx);
        for thisRepeat = 1:numel(trials)
            SpikeCount(thisRepeat,thisCond,thisElec) = numel(trials{thisRepeat});
        end
    end
end

meanPSTH = squeeze(mean(PSTH,2));
NormalizedPSTH = meanPSTH;
% Temporary solution, need to be revisited!
if NumElec > 1 && NumConditions > 1
    NormalizationFactor = max(max(meanPSTH,[],2),[],1); % dim 2 is conditions dim 1 is time bins
    for thisElec = 1:NumElec
        NormalizedPSTH(:,:,thisElec) = meanPSTH(:,:,thisElec) ./ NormalizationFactor(thisElec);
    end
elseif NumElec == 1 && NumConditions > 1
    NormalizationFactor = max(meanPSTH,[],1);
    NormalizedPSTH(:,thisElec) = meanPSTH(:,thisElec) ./ NormalizationFactor;
end
temp.PSTH = PSTH;
temp.NormalizedPSTH = NormalizedPSTH;
temp.meanPSTH = meanPSTH;
temp.SpikeCount = squeeze(sum(PSTH,1));
temp.meanSpikeCount = squeeze(mean(temp.SpikeCount,1));
temp.semSpikeCount = squeeze(std(temp.SpikeCount,0,1)/sqrt(size(temp.SpikeCount,1)));
temp.edges = edges;
temp.NumRepeats = NumRepeats;
temp.NumConditions = NumConditions;
temp.NumElec = NumElec;
temp.electrodes = electrodes;
PSTH = temp;