function [ARdata, MUAe, AllArtifactTimepoints] = VnCArtifactRemovalIndie(RawData,StimulationParam,MUAparameters,trialParam,SampleRate,interpolationLength)

StimStartTime = round(-trialParam.startT * SampleRate);
[ARdata, AllArtifactTimepoints] = VnCRemoveArtifactsSubroutine5(RawData,StimulationParam,StimStartTime,SampleRate);
% calculate MUA
MUAe = GetMUAewithInterpolationSubroutine(ARdata,SampleRate,MUAparameters,AllArtifactTimepoints,interpolationLength);
MUAe.time =  (0:(size(MUAe.data,1)-1))/MUAe.MUAesamplerate + trialParam.startT;
end