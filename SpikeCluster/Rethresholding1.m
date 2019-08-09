function IsSpike = Rethresholding(Handles)
% Rethreshold the waveforms based on original waveforms

% Transforming
Wav = Handles.Wav;
DetMethod = Handles.Para.DetMethod;
switch DetMethod
    case 1  %Raw Voltage
    case 2  %Nonlinear Energy Operator; Mukhopadhyay IEE98
        NLE = Wav.^2-Wav(:,[1 1:end-1]).*Wav(:,[2:end end]);
        NLE = NLE+Wav.^2-Wav(:,[1 1 1:end-2]).*Wav(:,[3:end end end]);
        NLE = NLE+Wav.^2-Wav(:,[1 1 1 1:end-3]).*Wav(:,[4:end end end end]);
        NLE = NLE+Wav.^2-Wav(:,[1 1 1 1 1:end-4]).*Wav(:,[5:end end end end end]);
        Wav = mean(NLE, 2);%conv(NLE, bartlett(round(6/30000*Handles.Fs)), 'same');
    case 3  %Mathematical Morphonogy; Geng et al., Neurocomputing 2010
        NPerm = min(10000, numel(Wav));
        PermID = randperm(numel(Wav), NPerm);
        PermWav = Wav(PermID);
        Gauss = gmdistribution.fit(PermWav', 3);
        Amp = max(Gauss.mu)-min(Gauss.mu);
        N = round(0.875/2*Handles.Fs/1000);
        Element = Amp*(0.54-0.46*cos(2*pi*(0:N)/N));
        OpenClose = MathMorph(MathMorph(Wav, Element, 'opening'), Element, 'closing');
        CloseOpen = MathMorph(MathMorph(Wav, Element, 'closing'), Element, 'opening');
        PVE = Wav-0.5*(OpenClose+CloseOpen);  % Peak-to-valley energy
        Wav = PVE;
end

% Re-thresholding
Spike = Handles.Spike;
[NumSpike, NumSample] = size(Wav);
SDDuration = Handles.Para.SDDuration;
SDMethod = Handles.Para.SDMethod;
DetThre = Handles.Para.DetThre;
DetPolar = Handles.Para.DetPolar;
for i = 1:NumSpike
    EndID = find( (Spike-Spike(i))>SDDuration, 1)-1;
    EndID = max(i, min([NumSpike EndID]));
    StartID = find( (Spike-Spike(EndID))<-SDDuration, 1, 'last')+1;
    StartID = min(i, max([1 StartID]));
    Input = Wav(StartID:EndID,:);
    NSpikes = size(Input, 1);
    Mu = repmat(median(Input), [NSpikes 1]);
    CoarseSD = repmat(std(Input), [NSpikes 1]);
    Nonoutlier = Input((Input>Mu-2*CoarseSD)&(Input<Mu+2*CoarseSD));
    %SD(i,1) = median(abs(Nonoutlier-mean(Nonoutlier)))/0.6745;
    SD(i,1) = min(std(Input));
    if i==1 && StartID==1 && EndID==NumSpike
        SD = repmat(SD, [NumSpike 1]);
        break;
    end
end
switch DetPolar
    case 1  %Negative
        IsSpike = sum(Wav<=repmat(-DetThre*SD, [1 NumSample]), 2)>0;
    case 2  %Positive
        IsSpike = sum(Wav>=repmat(DetThre*SD, [1 NumSample]), 2)>0;
    case 3  %Either
        IsSpike1 = sum(Wav<=repmat(-DetThre(1)*SD, [1 NumSample]), 2)>0;
        IsSpike2 = sum(Wav>=repmat(DetThre(2)*SD, [1 NumSample]), 2)>0;
        IsSpike = logical(IsSpike1+IsSpike2);
    case 4  %Both
        IsSpike1 = sum(Wav<=repmat(-DetThre(1)*SD, [1 NumSample]), 2)>0;
        IsSpike2 = sum(Wav>=repmat(DetThre(2)*SD, [1 NumSample]), 2)>0;
        IsSpike = logical(IsSpike1.*IsSpike2);
end
if DetMethod>=2  %NLE
    IsSpike = sum(Wav>=repmat(DetThre*SD, [1 NumSample]), 2)>0;
end