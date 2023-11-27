function alignPhotoDiode(basepath,blockNumber,PDchannel,instanceRange)
% align visual trials based on photo diode

% photodiode was recorded on instance 1. On instance 1 for every trial in the cut data,
% go to the channel where photodiode was recorded, calculate the time shift
% between the photodiode and the StimB and apply the same amount anti-shift
% to all data recorded.

% for alignment only, do not change
TrialAlignTimeStamp = 0.2*30000;
trialParam.startT = -0.2; %s
trialParam.endT = 0.2; %s
ExpectedMisalignmentRange = 0.15*30000;

trialAlignSTfieldname = 'trialAlignST';
instanceInfoPath = [basepath,'/','instance',num2str(1),'_B',num2str(blockNumber),'_instanceinfo.mat'];
load(instanceInfoPath);
for thisElec = 1:numel(PDchannel)
    fprintf('.')
    EID = PDchannel(thisElec);
    % load in data
    cmd = ['trialAlignST = instanceinfo.trialInfo.',trialAlignSTfieldname,';'];
    eval(cmd);
    electrodeRawTrialData = VnCGetTrialData3(instanceinfo, EID,trialAlignST,trialParam,basepath);
    savedir = [basepath,'/','instance',num2str(1),'_B',num2str(blockNumber),'/','elec',num2str(EID)];
    savepath = [savedir,'/','electrodeRawTrialData.bmat'];
    comment = '';
    SaveMatFast(savepath,int16(electrodeRawTrialData.trialData),comment)
    instanceinfo.electrodeRawTrialDataPath{thisElec} = savepath;
end
save(instanceInfoPath,'instanceinfo')

PDdata = LoadMatFast([basepath,'\','instance1_B',num2str(blockNumber),'\','elec',num2str(PDchannel),'\','electrodeRawTrialData.bmat']);

threshold = -3.1e4;
idx = PDdata<threshold;

debug = 0;

idx = diff(idx,1,2);

NumTrials = size(PDdata,1);
StimOnsetIdx = zeros(NumTrials,1);
TrialShiftAmount = zeros(NumTrials,1);

if debug
    figure
end

for thisTrial = 1:NumTrials
    idxThisTrial = find(idx(thisTrial,(TrialAlignTimeStamp-ExpectedMisalignmentRange) : (TrialAlignTimeStamp+ExpectedMisalignmentRange)) == 1,1,'first');
    if isempty(idxThisTrial)
        % do nothing
        disp(['Warning: no diode signal found in trial #',num2str(thisTrial)])
    else
        StimOnsetIdx(thisTrial) = idxThisTrial;
        TrialShiftAmount(thisTrial) = TrialAlignTimeStamp - (StimOnsetIdx(thisTrial)+TrialAlignTimeStamp-ExpectedMisalignmentRange-1);
    end
    if debug
        subplot(1,2,1)
        hold off
        plot(PDdata(thisTrial,:))
        hold on
        plot([1,size(PDdata,2)],[threshold threshold])
        plot([TrialAlignTimeStamp-TrialShiftAmount(thisTrial),TrialAlignTimeStamp-TrialShiftAmount(thisTrial)],threshold,'o')
        subplot(1,2,2)
        hold off
        ydata = PDdata(thisTrial,:);
        ydata = circshift(ydata,TrialShiftAmount(thisTrial));
        plot(ydata)
        hold on
        plot([1,size(PDdata,2)],[threshold threshold])
        plot([TrialAlignTimeStamp,TrialAlignTimeStamp],threshold,'o')
    end
end

for thisInstance = instanceRange
    instanceInfoPath = [basepath,'\','instance',num2str(thisInstance),'_B',num2str(blockNumber),'_instanceinfo.mat'];
    %
    load(instanceInfoPath);
    %
    instanceinfo.trialInfo.DiodeTrialShiftAmount = TrialShiftAmount;
    save(instanceInfoPath,'instanceinfo')
end
end