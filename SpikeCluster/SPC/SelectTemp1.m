function [Temprature] = SelectTemp(Model, MinSize)

% Looks for the most stable clusters as many as possible
for i = 1:size(Model, 1)
    if sum(Model(i,5:end)>MinSize)>0
        NCluster(i) = sum(Model(i,5:end)>MinSize);
        ID = find(Model(i,5:end)>MinSize);
        MinCluster(i) = min(Model(i,4+ID));
    end
end
MaxNCID = NCluster==max(NCluster);
for i = 1:numel(MaxNCID)
    ConvID = conv(single(MaxNCID), ones([1 i]));
    ConvID = ConvID(i:end-i+1);
    if sum(ConvID==i)>0
        StableID(i) = find(ConvID==i, 1);
    else
        StableID(i) = nan;
    end
end
StableNTemp = find(~isnan(StableID), 1, 'last');
StableID = StableID(StableNTemp);
MinCluster(1:StableID-1) = 0;
MinCluster(StableID+StableNTemp:end) = 0;
[~, Temprature] = max(MinCluster);

% If the second cluster is too small
if (Temprature==1 && Model(Temprature,6)<MinSize)
    Temprature = 2;
end