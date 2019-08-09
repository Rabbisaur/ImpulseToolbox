function TSLloadNevNew(filepath)
% clc, clear, close all
switch nargin
    case 1
        % do nothing
    otherwise
        error(' ')
end
if exist(filepath,'file') ~= 2
    error('NEV file not exist, please check')
else
    nevdata = openNEV(filepath);
end

% get digital port input
digicode = nevdata.Data.SerialDigitalIO.UnparsedData;
digihi = floor(double(digicode) / bin2dec('11111111'));
digilow = mod(double(digicode),bin2dec('1111111100000000'));

% split data into trials
% a trial should start with event 15 (START_PRE_TRIAL)
% and end with event 18 (END_POST_TRIAL)
start_pre_trial = 15;
end_post_trial = 18;
idx_start_pre_trial = find(digilow == start_pre_trial);
idx_end_post_trial = find(digilow == end_post_trial);

% if numel(idx_start_pre_trial) ~= numel(idx_end_post_trial)
%     error('Number of start_pre_trial events does not equal to number of end_post_trial events')
% end



% construct trial event list
numberEventsinEachTrial = [idx_start_pre_trial(1)-1;...
    diff(idx_start_pre_trial);numel(digilow)-idx_start_pre_trial(end)+1];
trialeventslow = mat2cell(digilow, numberEventsinEachTrial,1);
trialeventshigh = mat2cell(digihi, numberEventsinEachTrial,1);
trialeventstime = mat2cell(double(nevdata.Data.SerialDigitalIO.TimeStamp)',numberEventsinEachTrial,1);

NumberTotalTrials = numel(trialeventslow);

% exclude extra data in every trial beyond the first end_post_trial
for thisTrial = 1:NumberTotalTrials
    idx = find(trialeventslow{thisTrial} == end_post_trial,1,'first');
    trialeventslow{thisTrial} = trialeventslow{thisTrial}(1:idx);
    trialeventshigh{thisTrial} = trialeventshigh{thisTrial}(1:idx);
    trialeventstime{thisTrial} = trialeventstime{thisTrial}(1:idx);
end


% look into each trial to see whether it is valid
% a valid trial must have at least a number > 140 as the condition ID
% and a 23 which is TURN_TEST0_ON
% and a 24 which is TURN_TEST0_OFF
% and a 17 which is START_POST_TRIAL
% and a 96 which is REWARD_GIVEN
turn_test0_on = 23;
turn_test0_off = 24;
start_post_trial = 17;
reward_given = 96;

criteria = [turn_test0_on,turn_test0_off,start_post_trial,reward_given];

validity = zeros(NumberTotalTrials,1);
for thisTrial = 1:NumberTotalTrials
    validity(thisTrial) = sum(ismember(criteria,trialeventslow{thisTrial}))==numel(criteria);
end
validity = logical(validity);
% exclude invalid trials
trialeventslow = trialeventslow(validity);
trialeventshigh = trialeventshigh(validity);
trialeventstime = trialeventstime(validity);
numbervalidtrials = numel(trialeventslow);

turn_testnp1_on = [23,25,27,29,31];
turn_testnp1_off = [24,26,28,30,32];

% exclude any turnon/ turn off events before first turn_test0_on, caused by
% bugs in timing file
% for thisTrial = 1:numbervalidtrials
%     temp = trialeventslow{thisTrial};
%     idx = find(temp == turn_test0_on,1,'first');
%     [~,idx2,~] = intersect(temp,[turn_testnp1_on,turn_testnp1_off],'stable');
%     idx2 = idx2(idx2<idx);
%     trialeventslow{thisTrial}(idx2) = [];
%     trialeventshigh{thisTrial}(idx2) = [];
%     trialeventstime{thisTrial}(idx2) = [];
% end


% construct exp events matrix
Expmark = zeros(15,numbervalidtrials);
% Expmark (1,:) is the time of START_PRE_TRIAL
% Expmark (2,:) is the condition ID
% Expmark (3,:) is the reserved
% Expmark (4,:) is the time of END_POST_TRIAL
% Expmark (5,:) is the time of TURN_TEST0_ON
% Expmark (6,:) is the time of TURN_TEST1_ON
% Expmark (7,:) is the time of TURN_TEST2_ON
% Expmark (8,:) is the time of TURN_TEST3_ON
% Expmark (9,:) is the time of TURN_TEST4_ON
% Expmark (10,:) is the time of TURN_TEST0_OFF
% Expmark (11,:) is the time of TURN_TEST1_OFF
% Expmark (12,:) is the time of TURN_TEST2_OFF
% Expmark (13,:) is the time of TURN_TEST3_OFF
% Expmark (14,:) is the time of TURN_TEST4_OFF

% turn_testnp1_on = [23,25,27,29,31];
% turn_testnp1_off = [24,26,28,30,32];

for thisTrial = 1:numbervalidtrials
    % Expmark (1,:) is the time of START_PRE_TRIAL
    eventsinthisTrial = trialeventslow{thisTrial};
    timeinthisTrial = trialeventstime{thisTrial};
    idx = eventsinthisTrial == start_pre_trial;
    Expmark(1,thisTrial) = timeinthisTrial(idx);
    % Expmark (2,:) is the condition ID
    idx = find(eventsinthisTrial >= 192);
    if numel(idx) <2
        error(['wrong number of condition ID in trial #', num2str(thisTrial)])
    end
    CondID = (double(eventsinthisTrial(idx(1)))-192) * 64 + (double(eventsinthisTrial(idx(2))) - 192);
%     idx = find(eventsinthisTrial >=140 & eventsinthisTrial <192);
%     CondID2 = double(eventsinthisTrial(idx(1))) - 140;
    Expmark(2,thisTrial) = CondID;
%     Expmark(3,thisTrial) = CondID2;
    % Expmark (4,:) is the time of END_POST_TRIAL
    idx = eventsinthisTrial == end_post_trial;
    Expmark(4,thisTrial) = timeinthisTrial(idx);
    % Expmark (5-9,:) is the time of TURN_TESTN_ON
    for testnp1 = 1:numel(turn_testnp1_on)
        if ismember(turn_testnp1_on(testnp1),eventsinthisTrial)
            idx = eventsinthisTrial == turn_testnp1_on(testnp1);
            if sum(idx) > 1
                idx = find(idx);
                idx = idx(1);
            end
            Expmark(testnp1+4,thisTrial) = timeinthisTrial(idx);
        end
    end
    % Expmark (10-14,:) is the time of TURN_TESTN_OFF
    for testnp1 = 1:numel(turn_testnp1_off)
        if ismember(turn_testnp1_off(testnp1),eventsinthisTrial)
            idx = eventsinthisTrial == turn_testnp1_off(testnp1);
            idx = find(idx);
            idx = idx(1);
            Expmark(testnp1+9,thisTrial) = timeinthisTrial(idx);
        end
    end
end
% check conditions
conditionsID = unique(Expmark(2,:));
numtrialsCond = zeros(numel(conditionsID),1);
for thisCond = 1:numel(conditionsID)
    numtrialsCond(thisCond) = sum(Expmark(2,:) == conditionsID(thisCond));
end

% make sure each condition ID has equal number of trials
if numel(unique(numtrialsCond)) ~= 1
    disp('!!WARNING!!! not all conditions have equal number of trials')
    disp(numtrialsCond)
end

% split spike into electrode and trials
% exclude spikes that before Expmark(1,1) or after Expmark(4,end)
SpikeData = nevdata.Data.Spikes;
idx = SpikeData.TimeStamp >= Expmark(1,1) & SpikeData.TimeStamp <= Expmark(4,end);
SpikeData.Electrode = double(SpikeData.Electrode(idx));
SpikeData.TimeStamp = double(SpikeData.TimeStamp(idx));
SpikeData.Unit = double(SpikeData.Unit(idx));
SpikeData.Waveform = double(SpikeData.Waveform(:,idx));

% determine electrodes in the data
electrodes = unique(SpikeData.Electrode);
numberElectrodes = numel(electrodes);
SpikeTimeStamp = cell(128,1);
for thisEID = 1: numberElectrodes
    idx = SpikeData.Electrode == electrodes(thisEID);
    SpikeTimeStamp{electrodes(thisEID)} = SpikeData.TimeStamp(idx);
    rawSpikeTimeStamp{electrodes(thisEID)} = SpikeData.TimeStamp(idx);
    rawUnit{electrodes(thisEID)} = SpikeData.Unit(idx);
    rawWaveform{electrodes(thisEID)} = SpikeData.Waveform(:,idx);
end

% split data into trials
edgelower = [Expmark(1,:),Expmark(4,end)];
edgehigher = [Expmark(1,1),Expmark(4,:)];
interalign = @(x,y,z) intersect(x,y) - z;
f = @(x) x/double(nevdata.MetaTags.TimeRes);
SpikeTimeStampAligned = cell(128,1);
SpikeTimeAligned = cell(128,1);
for thisEID = 1:numberElectrodes
    spts = SpikeTimeStamp{electrodes(thisEID)};
    sptslower = spts';
    sptshigher = spts';
    spcontlower = histc(sptslower,edgelower);
    if spcontlower(end)
        sptslower(end) = [];
    end
    spcontlower(end) = [];
    
    spconthigher = histc(sptshigher,edgehigher);
    if spconthigher(end)
       sptshigher(end) = [];
    end
    spconthigher(end) = [];
    
    sptslower = mat2cell(sptslower,spcontlower,1);
    sptshigher = mat2cell(sptshigher,spconthigher,1);
    SpikeTimeStamp{electrodes(thisEID)} = cellfun(@intersect,sptslower,sptshigher,'UniformOutput',false);
    temp = cellfun(interalign,sptslower,sptshigher,num2cell(Expmark(6,:))','UniformOutput',false); % must be aligned to the stimulus onset
    SpikeTimeStampAligned{electrodes(thisEID)} = temp;
    SpikeTimeAligned{electrodes(thisEID)} = cellfun(f,temp,'UniformOutput',false);
end

Waveform = cell(128,1);
for thisEID = 1:numberElectrodes
    idx = cellfun(@isempty,SpikeTimeStamp{electrodes(thisEID)});
    ElecSPST = cell2mat(SpikeTimeStamp{electrodes(thisEID)}(~idx));
    eidx = SpikeData.Electrode == electrodes(thisEID);
    [~,idx,~] = intersect(SpikeData.TimeStamp(eidx),ElecSPST,'stable');
    if numel(idx) ~= numel(ElecSPST)
        error('lost spike during waveform assignment!')
    end
    temp = SpikeData.Waveform(:,eidx);
    Waveform{electrodes(thisEID)} = temp(:,idx);
end

% convert ExpmarkST to ExpmarkTime
ExpmarkST = double(Expmark);
idx = [1,4:size(Expmark,1)];
ExpmarkTime = ExpmarkST;
ExpmarkTime(idx,:) = double(ExpmarkST(idx,:)) / double(nevdata.MetaTags.TimeRes);


% convert Expmark to a human readable format
% Expmark (1,:) is the time of START_PRE_TRIAL
% Expmark (2,:) is the condition ID
% Expmark (3,:) is the reserved
% Expmark (4,:) is the time of END_POST_TRIAL
% Expmark (5,:) is the time of TURN_TEST0_ON
% Expmark (6,:) is the time of TURN_TEST1_ON
% Expmark (7,:) is the time of TURN_TEST2_ON
% Expmark (8,:) is the time of TURN_TEST3_ON
% Expmark (9,:) is the time of TURN_TEST4_ON
% Expmark (10,:) is the time of TURN_TEST0_OFF
% Expmark (11,:) is the time of TURN_TEST1_OFF
% Expmark (12,:) is the time of TURN_TEST2_OFF
% Expmark (13,:) is the time of TURN_TEST3_OFF
% Expmark (14,:) is the time of TURN_TEST4_OFF
ExpmarkStr.machineTime.START_PRE_TRIAL = ExpmarkST(1,:);
ExpmarkStr.machineTime.conditionID = ExpmarkST(2,:);
ExpmarkStr.machineTime.END_POST_TRIAL = ExpmarkST(4,:);
ExpmarkStr.machineTime.TURN_TEST0_ON = ExpmarkST(5,:);
ExpmarkStr.machineTime.TURN_TEST1_ON = ExpmarkST(6,:);
ExpmarkStr.machineTime.TURN_TEST2_ON = ExpmarkST(7,:);
ExpmarkStr.machineTime.TURN_TEST3_ON = ExpmarkST(8,:);
ExpmarkStr.machineTime.TURN_TEST4_ON = ExpmarkST(9,:);
ExpmarkStr.machineTime.TURN_TEST0_OFF = ExpmarkST(10,:);
ExpmarkStr.machineTime.TURN_TEST1_OFF = ExpmarkST(11,:);
ExpmarkStr.machineTime.TURN_TEST2_OFF = ExpmarkST(12,:);
ExpmarkStr.machineTime.TURN_TEST3_OFF = ExpmarkST(13,:);
ExpmarkStr.machineTime.TURN_TEST4_OFF = ExpmarkST(14,:);

ExpmarkStr.wallTime.START_PRE_TRIAL = ExpmarkTime(1,:);
ExpmarkStr.wallTime.conditionID = ExpmarkTime(2,:);
ExpmarkStr.wallTime.END_POST_TRIAL = ExpmarkTime(4,:);
ExpmarkStr.wallTime.TURN_TEST0_ON = ExpmarkTime(5,:);
ExpmarkStr.wallTime.TURN_TEST1_ON = ExpmarkTime(6,:);
ExpmarkStr.wallTime.TURN_TEST2_ON = ExpmarkTime(7,:);
ExpmarkStr.wallTime.TURN_TEST3_ON = ExpmarkTime(8,:);
ExpmarkStr.wallTime.TURN_TEST4_ON = ExpmarkTime(9,:);
ExpmarkStr.wallTime.TURN_TEST0_OFF = ExpmarkTime(10,:);
ExpmarkStr.wallTime.TURN_TEST1_OFF = ExpmarkTime(11,:);
ExpmarkStr.wallTime.TURN_TEST2_OFF = ExpmarkTime(12,:);
ExpmarkStr.wallTime.TURN_TEST3_OFF = ExpmarkTime(13,:);
ExpmarkStr.wallTime.TURN_TEST4_OFF = ExpmarkTime(14,:);


% save to disk
idx = find(filepath == '.', 1,'last');
dirpath = filepath(1:idx-1);
if ~exist(dirpath,'dir')
    mkdir(dirpath)
end
if ispc
    slash = '\';
else
    slash = '/';
end


nev.spikeTimeStamp = SpikeTimeStampAligned;
nev.expmark = ExpmarkStr;
nev.expmarkST = ExpmarkST;
nev.expmarkTime = ExpmarkTime;
nev.spikeTime = SpikeTimeAligned;
nev.electrodes = electrodes;
nev.rawspikeTimeStamp = rawSpikeTimeStamp;
nev.rawWaveform = rawWaveform;
nev.rawUnit = rawUnit;
nev.samplerate = double(nevdata.MetaTags.SampleRes);
nev.filename = nevdata.MetaTags.Filename;
nev.dateTime = nevdata.MetaTags.DateTime;
nev.dateTimeRaw = double(nevdata.MetaTags.DateTimeRaw);
nev.digitimestamps = double(nevdata.Data.SerialDigitalIO.TimeStamp);
if size(nev.digitimestamps,2) > size(nev.digitimestamps,1)
    nev.digitimestamps = nev.digitimestamps';
end
nev.digidata = double(nevdata.Data.SerialDigitalIO.UnparsedData);
if size(nev.digidata,2) > size(nev.digidata,1)
    nev.digidata = nev.digidata';
end

nevdata = nev;

matpath = [dirpath, slash, 'nevdata.mat'];
save(matpath,'nevdata')

% matpath = [dirpath, slash, 'SpikeTimeStamp.mat'];
% save(matpath,'SpikeTimeStamp','SpikeTimeStampAligned')
% matpath = [dirpath, slash, 'SpikeTime.mat'];
% save(matpath,'SpikeTimeAligned')
% matpath = [dirpath,slash,'Expmark.mat'];
% save(matpath,'ExpmarkST','ExpmarkTime')
% matpath = [dirpath,slash,'Waveform.mat'];
% save(matpath,'Waveform')
% matpath = [dirpath,slash,'electrodes.mat'];
% save(matpath,'electrodes')
% 
% matpath = [dirpath,slash,'rawspikeST.mat'];
% save(matpath,'rawSpikeTimeStamp')
% matpath = [dirpath,slash,'rawWaveform.mat'];
% save(matpath,'rawWaveform')
% matpath = [dirpath,slash,'rawUnit.mat'];
% save(matpath,'rawUnit')


end