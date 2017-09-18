function ClassSNR = GetSNR(Wav, Class)
% Get SNRs
ClassSNR = [];
ClassPool = unique(Class);
Idx = find(ClassPool>0);
if numel(Idx)>=1
    for i = numel(Idx):-1:1
        CurrWav = Wav(Class==ClassPool(Idx(i)),:);
        ClassSNR(i) = SNR(CurrWav, 1);
    end
end