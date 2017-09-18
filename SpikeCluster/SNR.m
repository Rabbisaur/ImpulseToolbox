function [SNR] = SNR(Wav, Method)
% Estimate the signal-to-noise ratio, for one unit only
% Method   1/2 -- detrended noise/sampling points after crossing threshold

if nargin==1
    Method = 1;
end
MuWav = mean(Wav, 1);
switch Method
    case 1  % Detrended noise
        Noise = single(Wav)-repmat(MuWav, [size(Wav, 1) 1]);
        SD = sqrt(25)*std(Noise(:));
%         Noise = abs(Noise)/0.6745;
%         SD = sqrt(25)*median(Noise(:));
    case 2  % Sampling points before crossing threshold
        Diff = diff(MuWav);
        DiffSD = std(Diff(1:5));
        Idx = find(Diff>Diff(1)+3*DiffSD, 1, 'first');
        Noise = single(Wav(:,1:Idx))-ones([size(Wav, 1) 1])*MuWav(:,1:Idx);
        SD = sqrt(25)*std(Noise(:));
%         Noise = abs(Noise)/0.6745;
%         SD = sqrt(25)*median(Noise(:));
end
[MaxWav, MaxIdx] = max(MuWav);
[MinWav, MinIdx] = min(MuWav);
SNR = (MaxWav-MinWav)/SD;