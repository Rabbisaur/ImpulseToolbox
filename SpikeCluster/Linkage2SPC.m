function [Model Class] = Linkage2SPC(LinkageModel, MinSize)
% Usage: [Model Class] = Linkage2SPC(LinkageModel) for select and get class
% Horizontally cut dendrogram generated by Linkage, then return the
% statistic information for each cutoff
NumRow = size(LinkageModel, 1);
Iter = 0;
HiBoundOK = 0;
NumCluster = 2*NumRow;
StepVal = ceil(NumRow*0.01);
while 1
    NumCluster = min(NumRow, max(1, NumCluster-StepVal));
    Clus = cluster(LinkageModel, 'maxclust', NumCluster);
    Table = tabulate(Clus);
    Iter = Iter+1;
    if HiBoundOK==0 && max(Table(:,2))>=MinSize  %locate the upper bound
        HiBoundOK = 1;
        Iter = 0;
        StepVal = -ceil(StepVal*0.1);
        if NumCluster==NumRow
            break;
        end
    end
    if HiBoundOK==1  %locate the lower bound
        if NumCluster==NumRow || max(Table(:,2))<MinSize*0.2
            break;
        end
    end
end
% 2 SPC
% NumCluster = 1:ceil((NumCluster-1)/50):NumCluster;
NumCluster = ceil(logspace(log10(1), log10(NumCluster), 100));
NumCluster = unique(NumCluster);
Clus = cluster(LinkageModel, 'maxclust', NumCluster);
[NumRow NumCol] = size(Clus);
NewClus = zeros([NumRow NumCol]);
% Model
Model = zeros([NumCol NumRow+4]);
Model(:,1) = 0:NumCol-1;  %order
Model(:,2) = NumCluster;  %number of cluster
Model(:,3) = zeros([NumCol 1]);  %0
Model(:,4) = zeros([NumCol 1]);  %0
for i = 1:NumCol
   Table = tabulate(Clus(:,i));
   [~, Index] = sort(Table(:,2), 'descend');
   Table = Table(Index,:);
   [Row Col] = size(Table);
   Model(i,5:5+Row-1) = Table(:,2);
   for j = 1:Row
       Index = Clus(:,i)==Table(j,1);
       NewClus(Index,i) = j-1;
   end
end
% Class
Class = ones([NumCol NumRow+2]);
Class(:,1) = 0:NumCol-1;  %order
Class(:,2) = NumCluster;  %number of cluster
Class(:,3:end) = NewClus';