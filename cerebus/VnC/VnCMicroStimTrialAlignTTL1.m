function instanceinfo = VnCMicroStimTrialAlignTTL(instanceinfo,CurrentLevel,arrayRefCh,StimulationParam)

% parameter
supersamplingRatio = 16;

disp('Aligning microstimulation trials...')
% load(matdatapath);

% validIdx = performance ~= 0; % not fixation broken;
% microstimIdx = allCurrentLevel > 0;

% Fl = 150; % Hz
% Fn = instanceinfo(1).samplerate/2;
% N = 2;
% [Bhp, Ahp] = butter(N,Fl/Fn,'high'); % pulse low pass filter



% load data for the trigger channel
NumRefCh = size(arrayRefCh,1);
MSTrialIdx = find(CurrentLevel>0);
NumMicrostimTrials = numel(MSTrialIdx);
for thisElec = 1:NumRefCh
    thisInstance = arrayRefCh(thisElec,1);
    EID = arrayRefCh(thisElec,2);

    electrodeTrialData = VnCGetTrialDataCell(instanceinfo, thisInstance, EID);
    if thisElec == 1
        AllElecData = electrodeTrialData.trialData(MSTrialIdx);
    else
        for thisTrial = 1:NumMicrostimTrials
            AllElecData{thisTrial} = AllElecData{thisTrial} + electrodeTrialData.trialData{MSTrialIdx(thisTrial)};
        end
    end
end

% only look at valid microstimulation trials

% AllElecData = filtfilt(Bhp, Ahp,AllElecData')';

% find start and end time for each trial
% thisInstance = arrayRefCh(1,1);
SampleRate = instanceinfo(thisInstance).samplerate;
% DPbeforeAlignPoint = (0-trialParam.startT) * SampleRate;
beginST = instanceinfo(thisInstance).trialInfo.beginST;
beginST = beginST(MSTrialIdx);
% TrialMaxvalue = cellfun(@max,AllElecData);
% TrialThresholdValue = TrialMaxvalue *2/3;
TrialAlignPointOffset = zeros(NumMicrostimTrials,1);
TrialMicroStimStartEndPoint = zeros(2,NumMicrostimTrials);
trialAlignST = instanceinfo(thisInstance).trialInfo.trialAlignST(MSTrialIdx);

for thisTrial = 1:NumMicrostimTrials
    data = interp(AllElecData{thisTrial},supersamplingRatio);
%     data = filtfilt(Bhp, Ahp,data);
    data = data - mean(data(1:(90*supersamplingRatio)));
    TrialMaxvalue = max(data);
    TrialThresholdValue = TrialMaxvalue *2/3;
    idx = data>TrialThresholdValue;
    sumIdx = cumsum(idx);
    [~,maxIdx] = max(sumIdx);
    idx2 = find(diff(idx) == 1);
    idx3 = idx2 - maxIdx;
    idx3(idx3 > 0) = [];
    idx2 = idx2(numel(idx3));
    idx = [idx2,maxIdx];
    
%     while 1
%         idx2 = find(diff(idx) > 1);
%         if isempty(idx2)
%             break
%         end
%         idx(idx2) = [];
%     end
    
    MStrialsAlignST = beginST(thisTrial) + round(idx(1)/supersamplingRatio)-1;
    TrialAlignPointOffset(thisTrial) = MStrialsAlignST - trialAlignST(thisTrial);
    TrialMicroStimStartEndPoint(1,thisTrial) = 0;

    TrialMicroStimStartEndPoint(2,thisTrial) = round(idx(2)/supersamplingRatio) - round(idx(1)/supersamplingRatio)+1;
    
%     pID = mod(thisTrial,9);
%     if pID == 1
%         figure
%     end
%     if pID == 0
%         pID = 9;
%     end
%     subplot(3,3,pID),hold on
%     plot(data)
%     plot([1,numel(data)],[TrialThresholdValue TrialThresholdValue])
%     plot(idx,[TrialThresholdValue TrialThresholdValue],'ro')
end

% check data
TrialMicroStimStartEndPointDiff = TrialMicroStimStartEndPoint(2,:) - TrialMicroStimStartEndPoint(1,:);
TrialMicroStimStartEndPointDiff = TrialMicroStimStartEndPointDiff/SampleRate;
StimLengh = StimulationParam.numberofStimPulses * 1/StimulationParam.PulsesFrequency;

idx = TrialMicroStimStartEndPointDiff > StimLengh * 1.05 | TrialMicroStimStartEndPointDiff < StimLengh * 0.95;
if sum(idx) > 0
    disp('Warning! Detected stimulation period is significantly different from the predicted stimulation length')
end


% save result to instance info
numInstances = numel(instanceinfo);
trialidx = CurrentLevel > 0;
for thisInstance = 1:numInstances
    correctedtrialAlignST = instanceinfo(thisInstance).trialInfo.trialAlignST;
    correctedtrialAlignST(trialidx) = correctedtrialAlignST(trialidx) + TrialAlignPointOffset; % important! use the TrialAlignPointOffset!
    instanceinfo(thisInstance).trialInfo.correctedtrialAlignST = correctedtrialAlignST;
    instanceinfo(thisInstance).trialInfo.TrialMicroStimStartEndPoint = TrialMicroStimStartEndPoint;
end

end