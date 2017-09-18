function [Class] = MergeCluster(Input, Class)
% Merge clusters that every point is statistically unseparable

Input = single(Input);
ClassPool = unique(Class);
ClassPool = ClassPool(ClassPool>0);
ReservedClass = ClassPool;
NumClass = numel(ClassPool>0);
for i = 1:NumClass
    if sum(ReservedClass==ClassPool(i))==0
        continue;
    end
    for j = i+1:NumClass
        if sum(ReservedClass==ClassPool(j))==0
            continue;
        end
        %RejNull = ttest2(Input(Class==ClassPool(i),:), Input(Class==ClassPool(j),:), 1e-5);
        for k = 1:size(Input, 2)
            [~, RejNull(k)] = ranksum(Input(Class==ClassPool(i),k), Input(Class==ClassPool(j),k), 'alpha', 1e-5);
        end
        if sum(RejNull)==0
            ReservedClass(ReservedClass==ClassPool(j)) = [];
            Class(Class==ClassPool(j)) = ClassPool(i);
        end
    end
end