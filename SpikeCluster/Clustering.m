function [Model Class] = Clustering(Input, Method, varargin)
% Usage: [Model Class] = Clustering(Input, Method, varargin)
Input = single(Input);
NumRow = size(Input, 1);
Class = zeros([NumRow 1]);
NumCluster = varargin{1};
switch Method
    case 1  %Single Linkage
        Model = linkage(Input, 'single', 'euclidean');
        Class = Model;
    case 2  %K-Means
        [Class Model] = kmeans(Input, NumCluster);
    case 3  %Expectation-maximization
        [Class Model] = EMGM(Input, NumCluster);
    case 4  %t-Distribution Expectation-maximization
        Input = double(Input);
        try
            [Class Model] = tdistEM(Input, 1);
        catch
            Class = ones([size(Input, 1) 1]);
            Model.Center = median(Input, 1);
            Model.Sigma{1} = cov(Input);
        end
    case 5  %Superparamagnetic
        Temprature = 0.01:0.01:0.26;
        SWCycles = 100;
        KNearNeighb = 11;
        [Class Model] = SPC(Input, Temprature, SWCycles, KNearNeighb);
    case 6  %Automatic Expectation-maximization
        for i = 1:10
            warning('off');
            for j = 1:100
                try
                    Gauss{i} = gmdistribution.fit(Input, i); clc;
                    break;
                end
                if j==100
                    errordlg('Current data can NOT be clustered using EM.');
                end
            end
            AIC(i) = Gauss{i}.AIC;
            BIC(i) = Gauss{i}.BIC;
        end
        [~, MinID] = min(BIC);
        if MinID==numel(BIC)
            ID = (BIC>=(min(BIC)+0.01*(BIC(1)-min(BIC))));
            NCluster = sum(ID(1:MinID));
        else
            NCluster = MinID;
        end
        Model = Gauss{NCluster};
        Class = cluster(Model, Input);
    otherwise
        errordlg('No current clustering algorithm was found.');
end