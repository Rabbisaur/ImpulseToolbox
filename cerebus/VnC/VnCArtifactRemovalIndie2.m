function [ARdata, AllArtifactTimepoints] = VnCArtifactRemovalIndie2(RawData,StimulationParam,trialParam,SampleRate)
StimStartTime = round(-trialParam.startT * SampleRate);
[ARdata, AllArtifactTimepoints] = VnCRemoveArtifactsSubroutine14(RawData,StimulationParam,StimStartTime,SampleRate);
end