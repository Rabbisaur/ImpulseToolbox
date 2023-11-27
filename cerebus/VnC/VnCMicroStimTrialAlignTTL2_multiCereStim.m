function [instanceinfo,data] = VnCMicroStimTrialAlignTTL2_multiCereStim(instanceinfo,CurrentLevel,arrayRefCh,StimulationParam,debugflag)

% parameter
supersamplingRatio = 16;
if ispc
    slash = '\';
else
    slash = '/';
end

disp('Aligning microstimulation trials...')
% load(matdatapath);


% load data for the trigger channel
NumRefCh = size(arrayRefCh,1);
MSTrialIdx = find(CurrentLevel>0);
NumMicrostimTrials = numel(MSTrialIdx);
StimLengh = 1;
StimLength = StimulationParam.numberofStimPulses * 1/StimulationParam.PulsesFrequency;

Checkdata = VnCGetTrialDataCell(instanceinfo, 1, 33);
Checkdata = Checkdata.trialData(MSTrialIdx);

if NumMicrostimTrials == 0
    % do nothing
    TrialAlignPointOffset = 0;
    TrialMicroStimStartEndPoint = 0;
else
    AllStimStartST = cell(NumMicrostimTrials,1);
    AllStimStopST = cell(NumMicrostimTrials,1);
    for thisTrial = 1:NumMicrostimTrials
        AllStimStartST{thisTrial} = [];
        AllStimStopST{thisTrial} = [];
    end
    counter = 0;
    ValidRefChflag = true(NumRefCh,1);
    ValidRefTrial = true(NumMicrostimTrials,NumRefCh);
    StimStartST = nan(NumMicrostimTrials,NumRefCh);
    StimStopST = nan(NumMicrostimTrials,NumRefCh);
    for thisElec = 1:NumRefCh
        thisInstance = arrayRefCh(thisElec,1);
        EID = arrayRefCh(thisElec,2);
        if isempty(instanceinfo(thisInstance).electrodeCachePath{EID})
            % skip
            ValidRefChflag(thisElec) = false;
            continue
        else
            counter = counter + 1;
            electrodeTrialData = VnCGetTrialDataCell(instanceinfo, thisInstance, EID);
            elecData = electrodeTrialData.trialData(MSTrialIdx);
            % check if there is TTL signal on this channel
            
        end
        
        % find TTL on this channel
        % find start and end time for each trial
%         SampleRate = instanceinfo(thisInstance).samplerate;
        beginST = instanceinfo(thisInstance).trialInfo.beginST;
        beginST = beginST(MSTrialIdx);
        TrialAlignPointOffset = zeros(NumMicrostimTrials,1);
        TrialMicroStimStartEndPoint = zeros(2,NumMicrostimTrials);
%         trialAlignST = instanceinfo(thisInstance).trialInfo.trialAlignST(MSTrialIdx);

        for thisTrial = 1:NumMicrostimTrials
            data = interp(elecData{thisTrial},supersamplingRatio);
            %     data = filtfilt(Bhp, Ahp,data);
            data = data - mean(data(1:(90*supersamplingRatio)));
            TrialMaxvalue = max(data);
            TrialThresholdValue = max([1.5e4,TrialMaxvalue *2/3]);
            idx = data>TrialThresholdValue;
            idx1 = find(diff(idx) == 1); % all rising phase
            idx2 = find(diff(idx) == -1); % all falling phase
            if isempty(idx1) || isempty(idx2)
                ValidRefTrial(thisTrial,thisElec) = false;
                continue
            end
            tmp = find(idx1 < idx2(end), 1 ,'last');
            if isempty(tmp)
                ValidRefTrial(thisTrial,thisElec) = false;
                continue
            end
            idx1 = idx1(1:tmp);
            tmp = find(idx2 > idx1(1), 1 ,'first');
            if isempty(tmp)
                ValidRefTrial(thisTrial,thisElec) = false;
                continue
            end
            idx2 = idx2(tmp:end);
            
            TTLlength = idx2-idx1;
%             [tmp1,tmp2] = min(abs(TTLlength/supersamplingRatio - StimLength*instanceinfo(1).samplerate));
            tmp1 = abs(TTLlength/supersamplingRatio - StimLength*instanceinfo(1).samplerate) < StimLength * instanceinfo(1).samplerate * 0.05;
          
            if isempty(tmp1)
                ValidRefTrial(thisTrial,thisElec) = false;
%                 ValidRefChflag(thisElec) = false;
%                 warning('Error! Detected stimulation period is significantly different from the predicted stimulation length')
%                 break
            else
%                 idx = [idx1(tmp1), idx2(tmp1)];
                
                StimStartST = beginST(thisTrial) + round(idx1(tmp1)/supersamplingRatio)-1;
                StimStopST = beginST(thisTrial) + round(idx2(tmp1)/supersamplingRatio)-1;
                AllStimStartST{thisTrial} = union(AllStimStartST{thisTrial},StimStartST);
                AllStimStopST{thisTrial} = union(AllStimStopST{thisTrial},StimStopST);
            end
            
            
%             if debugflag
%                 %h = figure('visible','off');hold on
%                 h = figure;hold on
%                 plot(data)
%                 plot([1,numel(data)],[TrialThresholdValue TrialThresholdValue])
%                 plot(idx,[TrialThresholdValue TrialThresholdValue],'ro')
%                 tmpY = interp(Checkdata{thisTrial},supersamplingRatio);
%                 plot(tmpY)
%                 %             xlim([idx(1)-100 idx(2)+100])
%                 %             ylim([min(tmpY) max(tmpY)])
%                 idx = find(instanceinfo(1).electrodeCachePath{1}==slash,1,'last');
%                 savepath = instanceinfo(1).electrodeCachePath{1}(1:idx);
%                 %             saveas(h,[savepath,num2str(thisTrial),'.fig'])
% %                 saveas(h,[savepath,num2str(thisTrial),'.png'])
%                 close(h)
%             end
%             if tmp1 > StimLength * instanceinfo(1).samplerate * 0.05
%                 ValidRefTrial(thisTrial,thisElec) = false;
% %                 ValidRefChflag(thisElec) = false;
% %                 warning('Error! Detected stimulation period is significantly different from the predicted stimulation length')
% %                 break
%             end
            
            
            
        end
        if sum(ValidRefTrial(:,thisElec)) == 0
            ValidRefChflag(thisElec) = false;
%             warning('Error! Detected stimulation period is significantly different from the predicted stimulation length')
        end
    end
    
    StimStartST(~ValidRefTrial) = nan;
    StimStopST(~ValidRefTrial) = nan;
    
    MStrialsAlignST = nan(NumMicrostimTrials,1);
    
    for thisTrial = 1:NumMicrostimTrials
        AllStimStartST{thisTrial} = unique(AllStimStartST{thisTrial});
        AllStimStopST{thisTrial} = unique(AllStimStopST{thisTrial});
        minST = min(AllStimStartST{thisTrial});
        MStrialsAlignST(thisTrial) = minST(1);
    end
    
    TrialAlignPointOffset = MStrialsAlignST - beginST;
    
    if counter == 0 || sum(ValidRefChflag) == 0
        error('No reference channel recorded in this file!')
    end
end

% save result to instance info
numInstances = numel(instanceinfo);
trialidx = CurrentLevel > 0;
for thisInstance = 1:numInstances
    
    correctedtrialAlignST = instanceinfo(thisInstance).trialInfo.trialAlignST;
    if NumMicrostimTrials == 0
        % do nothing
    else
        correctedtrialAlignST(trialidx) = instanceinfo(thisInstance).trialInfo.beginST(trialidx) + TrialAlignPointOffset; % important! use the TrialAlignPointOffset!
        instanceinfo(thisInstance).trialInfo.TrialMicroStimStartPoint = AllStimStartST;
        instanceinfo(thisInstance).trialInfo.TrialMicroStimStopPoint = AllStimStopST;
    end
    instanceinfo(thisInstance).trialInfo.correctedtrialAlignST = correctedtrialAlignST;
end


end