function [Output, InterpT] = Resampling(Input, Fs, Resolution, Method)
%Resample the waveform by cubic spline
% Usage:
%   [Output] = Resampling(Input, Resolution)
% Inputs:
%   Input         a waveform in each row
%   Fs            sampling frequency
%   Resolution    time resolution for interpolation, in ms
%   Method        interpolation method
% Outputs:
%   Output	   interpolated waveforms
% Copyright 2012-2012 Minggui Chen, Beijing Normal University
% Revision: 1.0 Date: 2012/2/9 23:59:59

NumSample = size(Input, 2);
T = (0:NumSample-1)/Fs*1000;
if Resolution==0
    Output = Input;
    return;
end
InterpT = T(1):Resolution:T(end);
if numel(InterpT)~=NumSample
    Output = interp1(T, Input', InterpT', Method)';
else
    Output = Input;
end