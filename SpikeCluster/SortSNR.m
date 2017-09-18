function [Class, ClassSNR] = SortSNR(Wav, Class, Outlier)
% Sort classes by SNR
ClassSNR = [];
ClassPool = unique(Class(~Outlier));
Idx = find(ClassPool>0);
if numel(Idx)>=1
    for i = numel(Idx):-1:1
        CurrWav = Wav((Class==ClassPool(Idx(i)))&(~Outlier),:);
        ClassSNR(i) = SNR(CurrWav, 1);
        MedianWav = median(CurrWav, 1);
        WeightedSNR(i) = ClassSNR(i)*(max(MedianWav)-min(MedianWav));
    end
    [~, NewClass] = sort(WeightedSNR, 'descend');
    ClassSNR = ClassSNR(NewClass);
    Clus = Class;
    for i = 1:numel(NewClass)
        Class(Clus==ClassPool(Idx(NewClass(i)))) = i;
    end
end