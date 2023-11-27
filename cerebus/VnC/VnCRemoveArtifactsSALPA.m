function [ARdata] = VnCRemoveArtifactsSALPA(RawTrace,StimulationParam,SampleRate)
% SALPA on full trial, dynamic interpolation length factor

% parameters
StimulationPulsesFrequency = StimulationParam.PulsesFrequency;
tau = 0.002; % sec
% tau = 1/StimulationPulsesFrequency; % sec

% inputs
if size(RawTrace,2) > size(RawTrace,1)
    RawTrace = RawTrace';
end
NumTrials = size(RawTrace,2);
ARdata = RawTrace;
for thisTrial = 1:NumTrials
    y = RawTrace(:,thisTrial);
    % perform Salpa here
    y = salpa(y,'tau',round(tau*SampleRate));
    y(isnan(y)) = 0;
    ARdata(:,thisTrial) = y;
end

end