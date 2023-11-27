function alignCereStimTTL(basepath,blockNumber,TTLchannel,instanceRange,TrialAlignTimeStamp,trialParam)
debug = 0;
% Align microstimulation trials using the TTL signal from CereStim

% there are 7-8 Cerestim 96 connected to analog channel 7-13(14)

Threshold = 2e4;
NumTTLchannel = numel(TTLchannel);
sumdata = [];

if debug
h = figure;
end

trialAlignSTfieldname = 'trialAlignST';
instanceInfoPath = [basepath,'/','instance',num2str(1),'_B',num2str(blockNumber),'_instanceinfo.mat'];
load(instanceInfoPath)
for thisTTLchannel = 1:NumTTLchannel
    fprintf('.')
    EID = TTLchannel(thisTTLchannel);
    % load in data
    cmd = ['trialAlignST = instanceinfo.trialInfo.',trialAlignSTfieldname,';'];
    eval(cmd);
    electrodeRawTrialData = VnCGetTrialData3(instanceinfo, EID,trialAlignST,trialParam,basepath);
    %             electrodeRawTrialData.trialData = electrodeRawTrialData.trialData';
    savedir = [basepath,'/','instance',num2str(1),'_B',num2str(blockNumber),'/','elec',num2str(EID)];
    savepath = [savedir,'/','electrodeRawTrialData.bmat'];
%     comment = instance.filename;
    comment = '';
    SaveMatFast(savepath,int16(electrodeRawTrialData.trialData),comment)
    instanceinfo.electrodeRawTrialDataPath{thisTTLchannel} = savepath;
end
save(instanceInfoPath,'instanceinfo')


for thisTTLchannel = 1:NumTTLchannel


tmp = LoadMatFast([basepath,'\','instance1_B',num2str(blockNumber),'\','elec',num2str(TTLchannel(thisTTLchannel)),'\','electrodeRawTrialData.bmat']);
if thisTTLchannel == 1
    sumdata = tmp;
else
    sumdata = sumdata + tmp;
    
end
end

idx = sumdata>Threshold;
if debug
    subplot(2,2,1)
    hold off, plot([0 max(size(tmp))],[Threshold Threshold]),hold on, 
    plot(sumdata(1:10,:)'),xlim([5000 10000])
end


idx = diff(idx,1,2);


NumTrials = size(sumdata,1);
StimOnsetIdx = zeros(NumTrials,1);
StimOffsetIdx = StimOnsetIdx;
TrialShiftAmount = zeros(NumTrials,1);
for thisTrial = 1:NumTrials
    if sum(idx(thisTrial,(TrialAlignTimeStamp-0.1*30000):(TrialAlignTimeStamp+0.1*30000))>0) > 0 % if there is TTL in this trial
        StimOnsetIdx(thisTrial) = find(idx(thisTrial,(TrialAlignTimeStamp-0.1*30000):(TrialAlignTimeStamp+0.1*30000)) == 1,1,'first');
        StimOffsetIdx(thisTrial) = find(idx(thisTrial,(TrialAlignTimeStamp-0.1*30000):end) == -1,1,'first');
        TrialShiftAmount(thisTrial) = TrialAlignTimeStamp - (StimOnsetIdx(thisTrial)+(TrialAlignTimeStamp-0.1*30000)-1);
    else
        TrialShiftAmount(thisTrial) = 0;
    end
end
RealStimParam.StimOnsetIdx = StimOnsetIdx;
RealStimParam.StimOffsetIdx = StimOffsetIdx;
RealStimParam.StimDuration = StimOffsetIdx - StimOnsetIdx;
RealStimParam.StimOnsetIdxCorrected = StimOnsetIdx + TrialShiftAmount;
RealStimParam.StimOffsetIdxCorrected = StimOffsetIdx + TrialShiftAmount;
RealStimParam.Threshold = Threshold;
path = [basepath,'\','instance1_B',num2str(blockNumber),'\','RealStimParam.mat'];
save(path,'RealStimParam')

if debug
for thisTrial = 1:NumTrials
            if TrialShiftAmount(thisTrial) ~= 0
                sumdata(thisTrial,:) = circshift(sumdata(thisTrial,:),TrialShiftAmount(thisTrial));
            end
end


    subplot(2,2,2)
    hold off, plot([0 max(size(tmp))],[Threshold Threshold]),hold on, 
    plot(sumdata(1:10,:)'),xlim([5000 10000])
end

for thisInstance = instanceRange
    instanceInfoPath = [basepath,'\','instance',num2str(thisInstance),'_B',num2str(blockNumber),'_instanceinfo.mat'];
    load(instanceInfoPath)
    instanceinfo.trialInfo.StimTTLTrialShiftAmount = TrialShiftAmount;
    save(instanceInfoPath,'instanceinfo')
end
end