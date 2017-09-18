function [Output, HipassFilter, LopassFilter] = Filtering(Input, Para)
%Detect spikes in filtered signals
% Usage:
%   Output = Filtering(Input, Fs, Para)
% Inputs:
%   Input      filtered signals with evenly spaced samples
%   Fs         sampling frequency
%   Para       parameters for detection
%                Para.Fs          sampling frequency
%                Para.FilterDirection  Forward & Reverse/Forward Only/Reverse Only
%                Para.HipassType  None/Butter/Cheby1/Cheby2/Ellip
%                Para.HipassOrder 3
%                Para.HipassFreq  250 Hz
%                Para.LopassType  None/Butter/Cheby1/Cheby2/Ellip
%                Para.LopassOrder 4
%                Para.LopassFreq  7500 Hz
% Outputs:
%   Output      filtered high-frequency component
% Copyright 2013-2013 Minggui Chen, Beijing Normal University
% Revision: 1.0 Date: 2013/10/23 23:59:59

Input = double(Input);
if isrow(Input)
    Input = Input';
end
Fs = Para.Fs;
FilterDir = Para.FilterDirection;
HipassType = Para.HipassType;
HipassOrder = Para.HipassOrder;
HipassFreq = Para.HipassFreq;
LopassType = Para.LopassType;
LopassOrder = Para.LopassOrder;
LopassFreq = Para.LopassFreq;

% High-pass the input
switch HipassType
    case 2  %Butterworth
        [HipassB, HipassA] = butter(HipassOrder, HipassFreq*2/Fs, 'high');
    case 3  %Cheby1
        [HipassB, HipassA] = cheby1(HipassOrder, .5, HipassFreq*2/Fs, 'high');
    case 4  %Cheby2
        [HipassB, HipassA] = cheby2(HipassOrder, .5, HipassFreq*2/Fs, 'high');
    case 5  %Ellip
        [HipassB, HipassA] = ellip(HipassOrder, .5, 20, HipassFreq*2/Fs, 'high');
end
HipassFilter = [];
if HipassType~=1
    HipassFilter.B = HipassB; HipassFilter.A = HipassA;
    Input = DirectedFilter(FilterDir, HipassB, HipassA, Input);
end

% Low-pass the input
LopassFreq = min(LopassFreq, 0.9999*Fs/2);
switch LopassType
    case 2  %Butterworth
        [LopassB, LopassA] = butter(LopassOrder, LopassFreq*2/Fs, 'low');
    case 3  %Cheby1
        [LopassB, LopassA] = cheby1(LopassOrder, .5, LopassFreq*2/Fs, 'low');
    case 4  %Cheby2
        [LopassB, LopassA] = cheby2(LopassOrder, .5, LopassFreq*2/Fs, 'low');
    case 5  %Ellip
        [LopassB, LopassA] = ellip(LopassOrder, .5, 20, LopassFreq*2/Fs, 'low');
end
LopassFilter = [];
if LopassType~=1
    LopassFilter.B = LopassB; LopassFilter.A = LopassA;
    Input = DirectedFilter(FilterDir, LopassB, LopassA, Input);  
end
Output = Input;



%%%%%%%%%%%%%%%%%%%Called functions
function Output = DirectedFilter(FilterDir, B, A, Input)
% Filter signal in terms of direction settings
switch FilterDir
    case 1  %Forward & Reverse
        Output = filtfilt(B, A, Input);
    case 2  %Forward Only
        Output = filter(B, A, Input);
    case 3  %Reverse Only
        Input = flipud(Input);
        Output = filter(B, A, Input);
        Output = flipud(Output);
end