function [Index] = SortFeature(Input, Method, VarRatio, NumCoeff)
[~, NumCol] = size(Input);
switch Method
    case 1  %no further selection
        Index = 1:NumCol;
        return;
    case 2  %maximum variance
        Input = single(Input);
        VarCol = var(Input, 1);
        [Score, Index] = sort(VarCol, 'descend');
    case 3  %KS test for normality
        KSStat = zeros([1 NumCol]);
        for i = 1:NumCol
            CurrVar = single(Input(:,i));
            Mu = mean(CurrVar);
            SD = std(CurrVar);
            LoBound = Mu-SD*3;
            HiBound = Mu+SD*3;
            Nonoutlier = Input(logical((CurrVar>LoBound).*(CurrVar<HiBound)),i);
            if numel(Nonoutlier)>10
                Nonoutlier = single(Nonoutlier);
                [~, ~, KSStat(i)] = kstest(zscore(Nonoutlier));
                %KSStat = KSStatistics(zscore(double(Nonoutlier)));
            end
        end
        [Score, Index] = sort(KSStat, 'descend');
    case 4  %squre sum of the difference between emprical and normal CDFs
        SSCDF = zeros([1 NumCol]);
        for i = 1:NumCol
            CurrVar = single(Input(:,i));
            Mu = mean(CurrVar);
            SD = std(CurrVar);
            LoBound = Mu-SD*3;
            HiBound = Mu+SD*3;
            Nonoutlier = Input(logical((CurrVar>LoBound).*(CurrVar<HiBound)),i);
            if numel(Nonoutlier)>10
                [CurrProb, Volt] = ecdf(CurrVar);
                TheoProb = cdf('norm', Volt, Mu, SD);
                SSCDF(i) = sum(abs(CurrProb-TheoProb).^2)/sum(CurrProb.^2);
                %SSCDF(i) = sum(abs(CurrProb-TheoProb));
                %SSCDF(i) = KLDiv(Volt,CurrProb/sum(CurrProb),TheoProb/sum(TheoProb), 'js');
            end
        end
        [Score, Index] = sort(SSCDF, 'descend');
    case 5  %goodness of fit for 1 vs. multiple Gaussians
        ImprovedGOF = zeros([1 NumCol]);
        for i = 1:NumCol
            CurrVar = single(Input(:,i));
            Mu = mean(CurrVar);
            SD = std(CurrVar);
            LoBound = Mu-SD*3;
            HiBound = Mu+SD*3;
            Nonoutlier = Input(logical((CurrVar>LoBound).*(CurrVar<HiBound)),i);
            if numel(Nonoutlier)>10
                [CurrProb, Volt] = ksdensity(CurrVar);
                for j = 1:3
                    warning('off');
                    try
                        Gauss = gmdistribution.fit(CurrVar, j); clc;
                    catch
                        Gauss = gmdistribution.fit(CurrVar, j); clc;
                    end
                    Gauss = pdf(Gauss, reshape(Volt, [numel(Volt) 1]))';
                    Gauss = Gauss/sum(Gauss); CurrProb = CurrProb/sum(CurrProb);
                    UnexpVar(j) = sum(abs(Gauss(:)-CurrProb(:)).^2)/sum(CurrProb.^2);
                    %UnexpVar(j) = sum(abs(Gauss(:)-CurrProb(:)));
                    %UnexpVar(j) = KLDiv(Volt,Gauss,CurrProb, 'js');
                end
                ImprovedGOF(1,i) = UnexpVar(1)-min(UnexpVar);
            end
        end
        [Score, Index] = sort(ImprovedGOF(1,:), 'descend');
    otherwise
        error('No current Method was found.');
end
ScoreRatio = cumsum(Score)/sum(Score);
NumCol = find(ScoreRatio>=VarRatio*0.01, 1, 'first');
NumCoeff = min(size(Input, 2), min([NumCoeff NumCol]));
Index = Index(1:NumCoeff);
if isempty(Index)
    Index = 1:min(2, size(Input,2));
end