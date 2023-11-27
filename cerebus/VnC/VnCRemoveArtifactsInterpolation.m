function [Interpolateddata,ArtefactLocation] = VnCRemoveArtifactsInterpolation(RawTrace,StimulationParam,StimStartTime,SampleRate)
% SALPA on full trial, dynamic interpolation length factor

% parameters
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
stimulationwaveformwidth = StimulationParam.stimulationwaveformwidth*SampleRate;
StimulationPulseInterval = 1/StimulationPulsesFrequency;
% tau = 0.002; % sec
tau = 1/StimulationPulsesFrequency; % sec

% inputs
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
NumTrials = size(RawTrace,2);

% get predicted stimulation locations 1
% numberofStimPulses = StimulationParam.numberofStimPulses;
% FirstPulseTimepoint = StimStartTime;
% predictedLoc = (0:(numberofStimPulses-1)) * StimulationPulseInterval * SampleRate;
% predictedLoc = round(predictedLoc + FirstPulseTimepoint);
% WaveDPbefore = 0;
% WaveDPafter1 = round(StimulationPulseInterval * SampleRate);
% idx = WaveDPbefore : WaveDPafter1;
% idx1 = bsxfun(@plus,idx,predictedLoc');

% dynamically determine interpolation factors
interpolationLengthFactor = 2.2;




% padding
NumDPpadding = 0;
if NumDPpadding > 0
    RawTrace = [zeros(NumDPpadding,NumTrials);RawTrace;zeros(NumDPpadding,NumTrials)];
end
NumDP = size(RawTrace,1);



Interpolateddata = RawTrace;

for thisTrial = 1:NumTrials
    NumStartPoint = numel(StimStartTime{thisTrial});
    trialStartTime = StimStartTime{thisTrial};
    for thisStim = 1:NumStartPoint
        % get predicted stimulation locations
        trialStartTime(thisStim) = trialStartTime(thisStim) + NumDPpadding;
        NumPulsesBefore = floor(trialStartTime(thisStim)/SampleRate / StimulationPulseInterval);
        NumPulsesAfter = floor((NumDP - trialStartTime(thisStim))/SampleRate / StimulationPulseInterval);
        numberofStimPulses = NumPulsesBefore + NumPulsesAfter + 1;
        
        FirstPulseTimepoint = floor(trialStartTime(thisStim) - NumPulsesBefore * StimulationPulseInterval* SampleRate);
        
        predictedLoc = (0:(numberofStimPulses-1)) * StimulationPulseInterval * SampleRate;
        predictedLoc = predictedLoc + FirstPulseTimepoint;
        loc0 = round(predictedLoc);
        
        WaveDPbefore = 0;
        WaveDPafter1 = stimulationwaveformwidth*interpolationLengthFactor(1);
        
        % determine AR interpolation length as 1.5 times pulse width

        % calculate artefect locations
        ArtefactLocation = zeros(numberofStimPulses,NumTrials);
        idx1 = round(-WaveDPbefore) : round(WaveDPafter1);
        
        % numinterpDP = numel(idx1);
        if size(idx1,2)>size(idx1,1) % column
            idx1 = idx1';
        end
        if size(loc0,1)>size(loc0,2) % row
            loc0 = loc0';
        end
        idx1 = idx1 + loc0;
        
        
        
        ArtefactLocation(:,thisTrial) = loc0;
        y = Interpolateddata(:,thisTrial);
        
        idx1(idx1<1) = 1;
        idx1(idx1>numel(y)) =numel(y);
        
        % perform linear interpolation 1
        x = (1:numel(y))';
        xx = x;
        x(idx1) = [];
        y(idx1) = [];
        y = interp1(x,y,xx,'linear');
        
        Interpolateddata(:,thisTrial)=y;
    end
end
end