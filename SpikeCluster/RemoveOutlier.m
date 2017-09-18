function [Class] = RemoveOutlier(Input, Class, Percentile, NewClass, ChkPoint)
% Remove outliers in each cluster & assign them to calss '-1'
% Usage: [Class] = RemoveOutlier(Input, Class, Percentile)

if nargin==3
    NewClass = -1;
    ChkPoint = 1;
end
if nargin==4
    ChkPoint = 1;
end
ClassPool = unique(Class);
Thre = norminv(1-Percentile/2, 0, 1);
for i = 1:numel(ClassPool)
    if ClassPool(i)<=0  %outliers/unclassified
        continue;
    end
    Index = find(Class==ClassPool(i));
    Data = single(Input(Index,:));
    [NumRow NumCol] = size(Data);
    % Distances between every sampling points
    Mu = median(Data, 1);
    SD = std(Data, 0, 1);
    Dist = Data-repmat(Mu, [NumRow 1]);
    if ChkPoint
        if numel(ClassPool)>1
            Bound = median(Dist, 1)-Thre*SD;
            Class(Index(sum(Dist<repmat(Bound, [NumRow 1]), 2)>0)) = NewClass;
            Bound = median(Dist, 1)+Thre*SD;
            Class(Index(sum(Dist>repmat(Bound, [NumRow 1]), 2)>0)) = NewClass;
        end
    end
    % Multidimensional distance
    Dist = sqrt(sum(Dist.^2, 2));
    Class(Index(Dist<median(Dist)-Thre*std(Dist))) = NewClass;
    Class(Index(Dist>median(Dist)+Thre*std(Dist))) = NewClass;
end
% Change the value of Class
ClassPool = unique(Class);
ValidIdx = find(ClassPool>0);
Clus = Class;
for i = 1:numel(ValidIdx)
    Class(Clus==ClassPool(ValidIdx(i))) = i;
end
% Check whether the calss is 1:n
ClassPool = unique(Class);
ValidIdx = find(ClassPool>0);
if max(ClassPool(ValidIdx))~=numel(ValidIdx)
    for i = 1:numel(ValidIdx)
        CurrClass = ClassPool(ValidIdx(i));
        Class(Class==CurrClass) = i;
    end
end