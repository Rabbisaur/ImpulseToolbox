function [Class] = LFPCluster(Wav, Coeff, Class)
% Decide whether a cluster is LFP activity
% Adapted from SAC_ver_2_15 tdistEM  (Shy Shoham, JNeurosciMehtods2003)

Coeff = single(Coeff);
ClassPool = unique(Class);
Idx = find(ClassPool>0);
NumCluster = numel(Idx);
for i = NumCluster:-1:1
    Center(i,:) = mean(Coeff(Class==ClassPool(Idx(i)),:));
    Covar{i} = cov(Coeff(Class==ClassPool(Idx(i)),:));
end
NumCol = size(Coeff, 2);
Loading = eye(NumCol);
NewClass=[];
if NumCluster>1  % Decide which clusters contain local field potentials or overlaps
   for i = 1:NumCluster
      sig1 = Loading(:,1:2)'*(Covar{i}*Loading(:,1:2));
      f(i) = prod(diag(sig1));
   end
%    Garbage = find(f>10*median(f));  % Large covariance indicates a garbage collector
%    NonGarbage = setdiff(1:NumCluster, Garbage);
   NonGarbage = 1:NumCluster;
   [~, NewOrder] = sort(sum(detrend(Center(NonGarbage,:)').^2));  %sort by ascending energy
   NewOrder = [NewOrder(NewOrder~=1) NewOrder(NewOrder==1)];
   for i = 1:numel(NewOrder)
       IsLFP(i) = 0;
       CurrWav = single(Wav(Class==NewOrder(i),:));
       LFPWav = mean(CurrWav, 1);
       LFPWavSD = std(CurrWav, 1);
       Bound = [LFPWav+1*LFPWavSD; LFPWav-1*LFPWavSD];
       Bound = Bound(:,1:5);
       WavDiff = diff(LFPWav);
       UniqueWav = unique(LFPWav);
       MinIdx = find(LFPWav==UniqueWav(1), 1, 'first');
       MaxIdx = find(LFPWav==UniqueWav(end), 1, 'last');
       IsLowLFP(i) = ((UniqueWav(1)<min(min(Bound)))+(UniqueWav(end)>max(max(Bound))))<2;
       UniqueWav = unique(LFPWav(MinIdx:end));
       LastMaxIdx = find(LFPWav(MinIdx:end)==UniqueWav(end), 1, 'last');
       LastMaxIdx = LastMaxIdx+MinIdx-1;
       ExtremRange = 5:numel(LFPWav)-5;
       if MinIdx==MaxIdx || sum(ExtremRange==MinIdx)==0 || sum(ExtremRange==MaxIdx)==0
           IsLFP(i) = 1;
       elseif MinIdx<MaxIdx  % Extracellular
           if sum(WavDiff(MinIdx:MaxIdx-1)<0)>1% || sum(WavDiff(MaxIdx:end)>0)>1
               IsLFP(i) = 1;
           end
       elseif MinIdx>MaxIdx  % Intracellular
          if sum(WavDiff(MaxIdx:MinIdx-1)>0)>1% || sum(WavDiff(MinIdx:LastMaxIdx-1)<0)>1
               IsLFP(i) = 1;
           end
           % Extracellular with local first peaked
%            if LastMaxIdx<(numel(LFPWav)-5) && sum(WavDiff(LastMaxIdx:end)>0)>1
%                IsLFP(i) = 1;
%            end
       end
   end
   if numel(NewOrder)==2 && sum(IsLFP)>=1  %keep the largest SNR
       NewOrder = [2 1];
       IsLFP = [1 0];
       IsLowLFP = [1 0];
   end
   NewOrder = [NewOrder(logical(IsLFP)) NewOrder(~logical(IsLFP))];
   IsLFP = [IsLFP(logical(IsLFP)) IsLFP(~logical(IsLFP))];
   IsLowLFP = [IsLowLFP(logical(IsLFP)) IsLowLFP(~logical(IsLFP))];
   sig1 = Loading(:,1:2)'*(Covar{NonGarbage(NewOrder(1))}*Loading(:,1:2));
   mu1 = Center(NonGarbage(NewOrder(1)),:)*Loading(:,1:2);
   if sum(mu1.^2)<4*sum(diag(sig1)) && (IsLFP(1) || IsLowLFP(1))  % LFP decision
        NumUnit = length(NewOrder)-1;
        NewClass(NonGarbage(NewOrder)) = 0:NumUnit;
        Unit = Class;
        for i = 1:numel(NewClass)
            Class(Unit==i) = NewClass(i);
        end
        Class = SortCluster(Wav, Coeff, Class);
   else
      NumUnit = length(NewOrder);
      NewClass(NonGarbage) = 1:NumUnit;
   end
%    NewClass(Garbage) = 255;
end