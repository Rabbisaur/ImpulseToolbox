function [Spike, Wav] = Detecting(Input, Signal, Para)
%Detect spikes in filtered signals
% Usage:
%   [Spike, Wav] = Detecting(Input, Fs, Para)
% Inputs:
%   Input      transformed signals with evenly spaced samples
%   Signal     filtered signals
%   Para       parameters for detection
%                Para.Thre
%                Para.Fs         sampling frequency
%                Para.DetPolar   Negative/Positive/Either/Both
%                Para.SpikeDuration
%                Para.AlignTime  Set pre-threshold duration
%                Para.RefractoryPeriod   in ms
% Outputs:
%   Spike      spike times
%   Wav        detected waveforms
% Copyright 2013-2013 Minggui Chen, Beijing Normal University
% Revision: 1.0 Date: 2013/10/23 23:59:59

if ~isnumeric(Input) || sum(size(Input)>1)>1
    error('INPUT must be a numeric vector.');
end
NumPoint = numel(Input);
NumSample = round(Para.SpikeDuration*Para.Fs/1000);
RPSample = round(Para.RefractoryPeriod*Para.Fs/1000);
switch Para.DetPolar
    case 1  %Negative
        IsSpike = Input<=Para.Thre(1);
    case 2  %Positive
        IsSpike = Input>=Para.Thre(1);
    case 3  %Either
        IsSpike = logical((Input<=Para.Thre(1))+(Input>=Para.Thre(2)));
    case 4  %Both
        NegSpike = Input<=Para.Thre(1);
        PosSpike = Input>=Para.Thre(2);
        IsSpike = false(size(NegSpike));
        for i = 1:NumSample-1
            %Neg is earlier than Pos
            IsSpike(1:end-i) = IsSpike(1:end-i) +NegSpike(1:end-i).*PosSpike(i+1:end);
            %Neg is later than Pos
            IsSpike(i+1:end) = IsSpike(i+1:end) +NegSpike(i+1:end).*PosSpike(1:end-i);
        end
        IsSpike = logical(IsSpike);
end
if Para.DetMethod==2 %NLE
    IsSpike = Input>=Para.Thre(1);
    switch Para.DetPolar
        case 1  %Negative
            IsSpike = logical(IsSpike.*(Signal<0));
        case 2  %Positive
            IsSpike = logical(IsSpike.*(Signal>0));
    end
end
%eliminate pseudospikes
i = 1;
while 1
    if i>NumPoint
        break;
    end
    if IsSpike(i)
        IsSpike( i+1:min(i+RPSample, NumPoint) ) = 0;
        i = i+RPSample;
    end
    i = i+1;
end
Spike = 1:NumPoint; Spike = Spike(IsSpike)';
AlignID = round(Para.AlignTime*Para.Fs/1000);
SpikeID = repmat((1:NumSample)-AlignID, [numel(Spike) 1]);
if numel(Spike)>0
    SpikeID = SpikeID+repmat(Spike, [1 NumSample]);
end
RowID = sum((SpikeID<=0)+(SpikeID>NumPoint), 2)==0;
SpikeID = SpikeID(RowID,:);
Spike = Spike(RowID);
Wav = reshape( Signal(SpikeID(:)), size(SpikeID) );
Spike = Spike/Para.Fs;