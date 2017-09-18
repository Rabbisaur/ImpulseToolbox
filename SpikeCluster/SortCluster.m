function Class = SortCluster(Wav, Coeff, Class)
% Sort well-identified classes maximally separated

ClassPool = unique(Class);
for thisClass = 1:numel(ClassPool)
    idx{thisClass} = Class == ClassPool(thisClass);
end
for thisClass = 1:numel(ClassPool)
    Class(idx{thisClass}) = thisClass-1;
end
ClassPool = unique(Class);
ValidIdx = find(ClassPool>0);
% Check whether the calss is 1:n
if max(ClassPool(ValidIdx))~=numel(ValidIdx)
    error('Class in SortCluster must be 1:n.');
end
if numel(ValidIdx)>1
    NewClass = zeros([numel(ValidIdx) 1]);
    Temp = zeros([numel(ValidIdx) size(Coeff, 2)]);
    for i = 1:numel(ValidIdx)
        Temp(i,:) = mean(Coeff(Class==ClassPool(ValidIdx(i)),:), 1);
        TempWav(i,:) = mean(Wav(Class==ClassPool(ValidIdx(i)),:), 1);
    end
    Dist = pdist(Temp, 'euclidean');
    TempIdx = fliplr(nchoosek(1:numel(ValidIdx), 2));
    MaxMinDiff = max(TempWav, [], 2)-min(TempWav, [], 2);
    % The two most separated classes
    [~, MaxDistIdx] = sort(Dist, 'descend');
    MaxDistDiff = MaxMinDiff(TempIdx(MaxDistIdx(1),:));
    [~, Idx] = sort(MaxDistDiff, 'descend');  % Larger max-min difference in waveform
    NewClass(1:2) = TempIdx(MaxDistIdx(1),Idx);
    % The third class with maximal joint distance relative to the above two
    if numel(ValidIdx)>2
        Idx = 1:numel(ValidIdx);
        Idx(Idx==NewClass(1)) = [];
        Idx(Idx==NewClass(2)) = [];
        Idx1 = [NewClass(1) Idx];
        Dist1 = pdist(Temp(Idx1), 'euclidean');
        Dist1 = Dist(1:numel(ValidIdx)-2);
        Idx2 = [NewClass(2) Idx];
        Dist2 = pdist(Temp(Idx2), 'euclidean');
        Dist2 = Dist(1:numel(ValidIdx)-2);
        [~, MaxDistIdx] = sort(sqrt(Dist1)+sqrt(Dist2), 'descend');
        NewClass(3) = Idx(MaxDistIdx(1));
    end
    % Sort the calsses remained by number of spikes, that is, raw orders
    if numel(ValidIdx)>3
        Idx(Idx==NewClass(3)) = [];
        NewClass(4:end) = Idx;
    end
    % Change the value of Class
    Clus = Class;
    for i = 1:numel(ValidIdx)
        Class(Clus==NewClass(i)) = i;
    end
end