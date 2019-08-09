function [ClassOut] = Matching(CoeffOut, CoeffIn, ClassIn, Method, varargin)
% Adapted from wave_clus by Quian Quiroga, Neural Computation2004
% Usage: [Class] = Matching(Input, Class, Method, varargin)
% Given classified points, try to classify new points via template matching
% Method  - nn, center, ml, mahal
% Percentile - max radius of cluster, in percentile.
% k     - # of nearest neighbors
% k_min - min # of nn for vote

CoeffOut = single(CoeffOut);
CoeffIn = single(CoeffIn);
sdnum = norminv(1-varargin{1}/2, 0, 1);
k = varargin{2};
k_min = varargin{3};
nspk = size(CoeffOut,1);
ClassOut = zeros(1,size(CoeffOut,1));
switch  lower(Method)
    case 'nn'
        sd    = sqrt(sum(var(CoeffIn,1)))*ones(1,size(CoeffIn,1));
        for i=1:nspk,
            nn = nearest_neighbor(CoeffOut(i,:),CoeffIn,sdnum*sd,Inf*ones(size(CoeffIn)),Inf,k);
            if( nn )
                Count = tabulate(ClassIn(nn));
                [~, MaxID] = max(Count(:,2));
                ClassOut(i) = Count(MaxID,1);
            else
                ClassOut(i) = 0;
            end
        end
    case 'center'
        [centers, sd, pd] = build_templates(ClassIn,CoeffIn); % we are going to ignore pd
        for i=1:nspk,
            ClassOut(i) = nearest_neighbor(CoeffOut(i,:),centers,sdnum*sd);        
        end
    case 'ml'
        [mu sigma] = fit_gaussian(CoeffIn,ClassIn);
        for i=1:nspk,
            ClassOut(i) = ML_gaussian(CoeffOut(i,:),mu,sigma);
        end
    case 'mahal'
        ClassPool = unique(ClassIn);
        for i = 1:numel(ClassPool)
            Dist(:,i) = mahaldist(CoeffOut, CoeffIn(ClassIn==ClassPool(i),:));            
        end
        [~, ClassOut] = min(Dist, [], 2);
        ClassOut = ClassPool(ClassOut);
%         [mu sigma] = fit_gaussian(CoeffIn,ClassIn);
%         for i=1:nspk,
%             ClassOut(i) = nearest_mahal(CoeffOut(i,:),mu,sigma);
%         end
    otherwise
        errordlg('No current Template Matching method was found.');
end

function index = nearest_neighbor(x,vectors,maxdist,varargin)
% function index = nearest_neigbor(x,vectors,maxdist,pointdist*,pointlimit*,k*)
% x is a row vector
% pointdist (optional) - vector of standard deviations
% pointlimit (optional) - upper bound on number of points outside pointdist
% k (optional) - number of points used for nearest neighbor
% Find the distance to all neighbors. Consider only those neighbors where
% the point falls in the radius of possibility for that point. Find the
% nearest possible neighbor.
% Return 0 if there is no possible nearest neighbor.
distances = sqrt(sum((ones(size(vectors,1),1)*x - vectors).^2,2)');
conforming = find(distances < maxdist);
if( length(varargin) > 0 )
    pointdist = varargin{1};
    if( length(varargin) > 1 )
        pointlimit = varargin{2};
    else
        pointlimit = Inf;
    end
    pointwise_conforming = [];
    for i=1:size(vectors,1),
        if( sum( abs(x-vectors(i,:)) > pointdist(i,:) ) < pointlimit )  % number of deviations from pointdist allowed.
            pointwise_conforming = [pointwise_conforming i];
        end
    end
    conforming = intersect(conforming, pointwise_conforming);
end
if( length( conforming ) == 0 )
    index = 0;
else
    if( length(varargin) > 2 )
        k = varargin{3};
        [y i] = sort(distances(conforming)); % k-nearest neighbors
        i = i(1:min(length(i),k));
    else
        [y i] = min(distances(conforming));   
    end
    index = conforming(i);
end

function [templates, maxdist, pointdist] = build_templates(classes,features)
max_class = max(classes);
feature_dim = size(features,2);
templates = zeros(max_class, feature_dim);
maxdist   = zeros(1,max_class);
pointdist = zeros(max_class,feature_dim);
for i=1:max_class,
    fi = features(find(classes==i),:);
    templates(i,:) = mean(fi);
    maxdist(i)     = sqrt(sum(var(fi,1)));   % the 1 means that we want sum(x-m)^2/N, not N-1
                                             % maxdist is the std dev of
                                             % the euclidean distance from
                                             % mean.
    pointdist(i,:)   = sqrt(var(fi,1));      % the std dev of the variation along each dimension.
    
end

function [mu, sigma] = fit_gaussian(x,class)
N = max(class);
mu = zeros(N,size(x,2));
sigma = zeros(size(x,2),size(x,2),N);

for i=1:N,
    mu(i,:) = mean(x(class==i,:));
    sigma(:,:,i) = cov(x(class==i,:));
end

function index = ML_gaussian(x,mu,sigma)
% function index = ML_gaussian(x,mu,sigma)
% x is a vector drawn from some multivariate gaussian
% mu(i,:) is the mean of the ith Gaussian
% sigma(:,:,i) is the covariance of the ith Gaussian
% 
% Returns the index of the Gaussian with the highest value of p(x).
N = size(mu,1);  % number of Gaussians

if( N == 0 )
    index = 0;
else
    for i=1:N,
        % leave out factor of 1/(2*pi)^(N/2) since it doesn't affect argmax
        p(i) = 1/sqrt(det(sigma(:,:,i)))*exp(-0.5*(x-mu(i,:))*inv(sigma(:,:,i))*(x-mu(i,:))');
    end
    [m index] = max(p);
end

function index = nearest_mahal(x,mu,sigma)
% function index = nearest_mahal(x,mu,sigma)
% x is a vector
% mu(i,:) is the mean of the ith Gaussian
% sigma(:,:,i) is the covariance of the ith Gaussian
% 
% Returns the index of the Gaussian closest (by the Mahalanobis distance)
% to x.
N = size(mu,1);  % number of Gaussians
d = [];
if( N == 0 )
    index = 0;
else
    for i=1:N,
        d(i) = (x-mu(i,:))*inv(sigma(:,:,i))*(x-mu(i,:))';
    end
    [m index] = min(d);
end

function d = mahaldist(Y,X)
%MAHAL Mahalanobis distance.
[rx,cx] = size(X);
[ry,cy] = size(Y);
m = mean(X,1);
M = m(ones(ry,1),:);
C = X - m(ones(rx,1),:);
[Q,R] = qr(C,0);
warning('off');
ri = R'\(Y-M)';
d = sum(ri.*ri,1)'*(rx-1);