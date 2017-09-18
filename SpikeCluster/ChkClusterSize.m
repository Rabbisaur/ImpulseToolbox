function [Class] = ChkClusterSize(Class, MinSize)
ClassPool = unique(Class);
for i = 1:numel(ClassPool)
    if ClassPool(i)>0  %outliers/unclassified
        if sum(Class==ClassPool(i))<MinSize
            Class(Class==ClassPool(i)) = 0;
        end
    end
end
% Change the value of Class
ClassPool = unique(Class);
ValidIdx = find(ClassPool>0);
Clus = Class;
for i = 1:numel(ValidIdx)
    Class(Clus==ClassPool(ValidIdx(i))) = i;
end