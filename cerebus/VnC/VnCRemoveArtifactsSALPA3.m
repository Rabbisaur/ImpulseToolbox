function [ARdata] = VnCRemoveArtifactsSALPA3(RawTrace,StimulationParam,StimStartTime,SampleRate,AllStimulationPulsesFrequency)
% SALPA on full trial, dynamic interpolation length factor

% parameters
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
% stimulationwaveformwidth = StimulationParam.stimulationwaveformwidth*SampleRate;
StimulationPulseInterval = 1/StimulationPulsesFrequency;
tau = 0.001; % sec
% tau = 1/StimulationPulsesFrequency; % sec

% inputs
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end


NumTrials = size(RawTrace,2);
ARdata = RawTrace;
% for thisTrial = 1:NumTrials
%     y = RawTrace(:,thisTrial);
%     % perform Salpa here
%     y = salpa(y,'tau',round(tau*SampleRate));
%     y(isnan(y)) = 0;
%     ARdata(:,thisTrial) = y;
% end

% perform averaging
% find all segments after a stimulation pulse
for thisTrial = 1:NumTrials
    NumStartPoint = numel(StimStartTime{thisTrial});
    trialStartTime = StimStartTime{thisTrial};
    for thisStim = 1:NumStartPoint
        % get predicted stimulation locations
        numberofStimPulses = StimulationParam.numberofStimPulses;
        
        FirstPulseTimepoint = trialStartTime(thisStim);
        
        stimFreq = AllStimulationPulsesFrequency(thisTrial,thisStim);
        StimulationPulseInterval = 1/stimFreq;
        
        predictedLoc = (0:(numberofStimPulses-1)) * StimulationPulseInterval * SampleRate;
        predictedLoc = predictedLoc + FirstPulseTimepoint;
        loc0 = round(predictedLoc);
        
        
        % determine AR interpolation length as 1.5 times pulse width

        % calculate artefect locations
%         ArtefactLocation = zeros(numberofStimPulses,NumTrials);
        idx1 = 0 : mean(diff(loc0));
        
        % numinterpDP = numel(idx1);
        if size(idx1,2)>size(idx1,1) % column
            idx1 = idx1';
        end
        if size(loc0,1)>size(loc0,2) % row
            loc0 = loc0';
        end
        idx1 = idx1 + loc0;
        
        y = ARdata(:,thisTrial);

        segments = y(idx1);
        avgSegments = mean(segments,2);
        avgSegments = avgSegments - mean(avgSegments);
        y(idx1) = y(idx1) - avgSegments;
        
        ARdata(:,thisTrial)=y;
    end
end

% salpa second pass
for thisTrial = 1:NumTrials
    y = RawTrace(:,thisTrial);
    % perform Salpa here
    y = salpa(y,'tau',round(tau*SampleRate));
    y(isnan(y)) = 0;
    ARdata(:,thisTrial) = y;
end

end