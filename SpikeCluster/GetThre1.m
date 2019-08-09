function Thre = GetThre(Handles)

Signal = Handles.TSpikeSignal;
NumTrial = numel(Signal);
StartT = Handles.StartT;
% if isfield(Handles, 'TTLT')
%     TTLT = Handles.TTLT;
%     StartT = StartT+TTLT;
% end
SDDuration = Handles.Para.SDDuration;
SDMethod = Handles.Para.SDMethod;
DetThre = Handles.Para.DetThre;
DetPolar = Handles.Para.DetPolar;
for i = 1:NumTrial
    EndID = find( (StartT-StartT(i))>SDDuration, 1)-1;
    EndID = max(i, min([NumTrial EndID]));
    StartID = find( (StartT-StartT(EndID))<-SDDuration, 1, 'last')+1;
    StartID = min(i, max([1 StartID]));
    Input = cell2mat(Signal(StartID:EndID));
    switch SDMethod  %Chan et al., J Neurosci Methods 2008
        case 1  %Root mean square
            Mu = median(Input);
            CoarseSD = std(Input);
            Nonoutlier = Input((Input>Mu-2*CoarseSD)&(Input<Mu+2*CoarseSD));
            SD(i,1) = std(Nonoutlier);
        case 2  %Median absolute deviation
            Mu = median(Input);
            CoarseSD = std(Input);
            Nonoutlier = Input((Input>Mu-2*CoarseSD)&(Input<Mu+2*CoarseSD));
            SD(i,1) = median(abs(Nonoutlier-mean(Nonoutlier)))/0.6745;
        case 3  %Cap fitting
        case 4  %Duty-cycle keeping
        case 5  %Max-min spread
    end
    if i==1 && StartID==1 && EndID==NumTrial
        SD = repmat(SD, [NumTrial 1]);
        break;
    end
end
switch DetPolar
    case 1  %Negative
        Thre(:,1) = -DetThre*SD;
    case 2  %Positive
        Thre(:,1) = DetThre*SD;
    case {3;4}  %Either/Both
        Thre(:,1) = -DetThre(1)*SD;
        Thre(:,2) = DetThre(2)*SD;
end
if Handles.Para.DetMethod==2  %NLE
    Thre = abs(Thre);
end