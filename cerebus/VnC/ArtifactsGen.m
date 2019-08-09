function trialData = ArtifactsGen(StimulationParam,trialParam,SampleRate)

% SampleRate = 30000;
% StimulationParam.numberofStimPulses = 50;
% StimulationParam.PulsesFrequency = 300.6942; %Hz
% StimulationParam.stimulationwaveformwidth = 0.0004; % sec
% StimulationParam.Strength = 5000;
% trialParam.startT = -400/1000; % a trial starts at -200 ms relative to the trigger
% trialParam.endT = 500/1000;


Length = trialParam.endT - trialParam.startT;% sec
supersamplingRatio = 16;
StimulationParam.Var = 0.05; % 5 percent
tau = 0.002; % sec

spreadRatio = 1;
TrialStartTime = -round(trialParam.startT * SampleRate * supersamplingRatio);

% calulate artifact time point
predictedLoc = (0:(StimulationParam.numberofStimPulses-1)) * 1/StimulationParam.PulsesFrequency * SampleRate * supersamplingRatio;
predictedLoc = predictedLoc + TrialStartTime;
predictedLoc = floor(predictedLoc);

% total length
totalLength = Length*SampleRate*supersamplingRatio;
trialData = zeros(Length*SampleRate*supersamplingRatio,1);

for thisArtifact = 1:StimulationParam.numberofStimPulses

artifact = zeros(Length*SampleRate*supersamplingRatio,1);

% sine part
x = (-StimulationParam.stimulationwaveformwidth*SampleRate*supersamplingRatio*spreadRatio):(StimulationParam.stimulationwaveformwidth*SampleRate*supersamplingRatio*spreadRatio);
sineLength = numel(x);
a = StimulationParam.Strength + StimulationParam.Strength * rand(1) * StimulationParam.Var;
b = 2*pi/sineLength;
y1 = @(t) a * sin(b*t);
artifact(1:numel(x)) = y1(x);

% exponential decay part
r  = 0.88 + 0.05 * rand(1);
idx = round(sineLength * r);
N0 = 0.3*artifact(idx);
lamda = 1/((tau + tau * rand(1) * 0.05) *SampleRate* supersamplingRatio);
y2 = @(t) N0*exp(-lamda*t);
x = 0:(totalLength-idx);
tmp = y2(x)';
idx2 = (artifact(idx:end) - tmp) < 0;
tmp2 = artifact(idx:end);
tmp2(idx2) = tmp(idx2);
artifact(idx:end) = tmp2;
% artifact(idx:end) = y2(x)';

artifactLength = totalLength - predictedLoc(thisArtifact) + 1;
trialData(predictedLoc(thisArtifact):end) = trialData(predictedLoc(thisArtifact):end) + artifact(1:artifactLength);
end

trialData = circshift(trialData,round(rand(1)*supersamplingRatio-supersamplingRatio/2));
trialData = downsample(trialData,supersamplingRatio);

plot(trialData)