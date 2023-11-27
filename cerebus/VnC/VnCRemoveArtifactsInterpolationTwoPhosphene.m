function [Interpolateddata,AllStimulationPulsesFrequency,AllInterpolationPoints] = VnCRemoveArtifactsInterpolationTwoPhosphene(RawTrace,StimulationParam,StimStartTime,SampleRate)
% SALPA on full trial, dynamic interpolation length factor

% parameters
% StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
stimulationwaveformwidth = StimulationParam.stimulationwaveformwidth*SampleRate;
% StimulationPulseInterval = 1/StimulationPulsesFrequency;
% tau = 0.002; % sec
% tau = 1/StimulationPulsesFrequency; % sec

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
interpolationLengthFactor = 2.5;




% padding
NumDPpadding = 0;
if NumDPpadding > 0
    RawTrace = [zeros(NumDPpadding,NumTrials);RawTrace;zeros(NumDPpadding,NumTrials)];
end
NumDP = size(RawTrace,1);



Interpolateddata = RawTrace;
AllInterpolationPoints = cell(NumTrials,1);
for thisTrial = 1:NumTrials
    AllInterpolationPoints{thisTrial} = false(size(Interpolateddata,1),1);
    NumStartPoint = numel(StimStartTime{thisTrial});
    trialStartTime = StimStartTime{thisTrial};
    for thisStim = 1:NumStartPoint
        StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
        StimulationPulseInterval = 1/StimulationPulsesFrequency;
        y = Interpolateddata(:,thisTrial);
        HW = 0;
        stimLength = round(StimulationParam.numberofStimPulses * StimulationPulseInterval * SampleRate);
        data = y(trialStartTime(thisStim):(trialStartTime(thisStim)+stimLength));
        x = 1:numel(data);
        xx = 1:0.1:numel(data);
        data = spline(x,data,xx);
        [f,fftY] = ez_powermeasure(data',SampleRate,HW);
        idx = f > StimulationPulsesFrequency - 10 & f < StimulationPulsesFrequency + 10;
        f = f(idx);
        fftY = fftY(idx);
        x = 1:0.01:numel(f);
        f = spline(1:numel(f),f,x);
        fftY = spline(1:numel(fftY),fftY,x);
        idx = f > StimulationPulsesFrequency - 5 & f < StimulationPulsesFrequency + 5;
        [~,idx2] = max(fftY(idx));
        idx = find(idx);
        StimulationPulsesFrequency = f(idx(idx2));
        if abs(StimulationPulsesFrequency - StimulationParam.PulsesFrequency) > 2
            disp(['StimulationPulsesFrequency = ', num2str(StimulationPulsesFrequency)])
            warning('Detected stimulation pulse frequency significantly different from commanded frequency')
            disp('Using commanded stimulation pulse frequency')
            StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
        end
%         StimulationPulsesFrequency= StimulationPulsesFrequency + 0.72;
        StimulationPulseInterval = 1/StimulationPulsesFrequency;
        AllStimulationPulsesFrequency(thisTrial,thisStim) = StimulationPulsesFrequency;
        
        % substract average
        % get predicted stimulation locations
        numberofStimPulses = StimulationParam.numberofStimPulses;
        
        FirstPulseTimepoint = trialStartTime(thisStim);
        
        predictedLoc = (0:(numberofStimPulses-1)) * StimulationPulseInterval * SampleRate;
        predictedLoc = predictedLoc + FirstPulseTimepoint;
        loc0 = round(predictedLoc);
        
        
        % determine AR interpolation length as x times pulse width

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
        
        % make the trace level1
%         figure
%         plot(y),hold on
        for thisArtifact = 1:size(idx1,2)
            d = y(idx1(1,thisArtifact) + ceil(stimulationwaveformwidth*interpolationLengthFactor(1))) - y(idx1(1,thisArtifact));
            y(idx1(1,thisArtifact)+ ceil(stimulationwaveformwidth*interpolationLengthFactor(1)):end) = y(idx1(1,thisArtifact)+ ceil(stimulationwaveformwidth*interpolationLengthFactor(1)):end) - d;
%             k = (y(idx1(end,thisArtifact)) - y(idx1(1,thisArtifact)+ ceil(stimulationwaveformwidth*interpolationLengthFactor(1)))) / (size(idx1,1)- ceil(stimulationwaveformwidth*interpolationLengthFactor(1)));
%             x = 0:((size(idx1,1)- ceil(stimulationwaveformwidth*interpolationLengthFactor(1)))-1);
%             yy = x*k;
%             y(idx1(1,thisArtifact)+ ceil(stimulationwaveformwidth*interpolationLengthFactor(1)):idx1(end,thisArtifact)) = y(idx1(1,thisArtifact)+ ceil(stimulationwaveformwidth*interpolationLengthFactor(1)):idx1(end,thisArtifact)) - yy';
        end
        d = y(idx1(end,thisArtifact)+1)-y(idx1(end,thisArtifact));
        y(idx1(end,thisArtifact)+1:end) = y(idx1(end,thisArtifact)+1:end) - d;
        
        % substract average before interpolation
        segments = y(idx1);
        avgSegments = mean(segments,2);
        avgSegments = avgSegments - mean(avgSegments);
        y(idx1) = y(idx1) - avgSegments;
        
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
        
        % determine AR interpolation length as x times pulse width

        % calculate artefect locations
%         ArtefactLocation = zeros(numberofStimPulses,NumTrials);
        idx1 = round(-WaveDPbefore) : round(WaveDPafter1);
        
        % numinterpDP = numel(idx1);
        if size(idx1,2)>size(idx1,1) % column
            idx1 = idx1';
        end
        if size(loc0,1)>size(loc0,2) % row
            loc0 = loc0';
        end
        idx1 = idx1 + loc0;
        
        
        
%         ArtefactLocation(:,thisTrial) = loc0;
        
        
        
        idx1(idx1<1) = 1;
        idx1(idx1>numel(y)) =numel(y);
        AllInterpolationPoints{thisTrial}(sort(reshape(idx1,[],1))) = true;
        
        % make the trace level
        for thisArtifact = 1:size(idx1,2)
            if thisArtifact == 1
                if idx1(1,thisArtifact)-1 > 0
                    
                    d = y(idx1(end,thisArtifact)+1) - y(idx1(1,thisArtifact)-1);% + (y(idx1(end,thisArtifact)) - y(idx1(1,thisArtifact+1)))/(idx1(1,thisArtifact+1) - idx1(end,thisArtifact))*size(idx1,1);
                    y(idx1(end,thisArtifact):end) = y(idx1(end,thisArtifact):end) - d;
                end
            else
                if idx1(end,thisArtifact)+1 <= numel(y)
                    d = y(idx1(end,thisArtifact)+1) - y(idx1(1,thisArtifact)-1);% + (y(idx1(end,thisArtifact-1)) - y(idx1(1,thisArtifact)))/(idx1(1,thisArtifact) - idx1(end,thisArtifact-1))*size(idx1,1);
                    y(idx1(end,thisArtifact):end) = y(idx1(end,thisArtifact):end) - d;
                end
            end
        end
        if idx1(end,thisArtifact)+1 <= numel(y)
            d = y(idx1(end,thisArtifact)+1)-y(idx1(end,thisArtifact));
            y(idx1(end,thisArtifact)+1:end) = y(idx1(end,thisArtifact)+1:end) - d;
        end
        
        % perform linear interpolation 1
        x = (1:numel(y))';
        xx = x;
        x(idx1) = [];
        y(idx1) = [];
        
        % SALPA
        tau = 0.001; % sec
        y = salpa(y,'tau',round(tau*SampleRate));
        
        y = interp1(x,y,xx,'linear');
        
        Interpolateddata(:,thisTrial)=y;
    end
end
end