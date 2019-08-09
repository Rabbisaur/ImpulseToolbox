function pattern = cartesiangratings2(sizeX,sizeY,meanLum,contrast,radius,orientations,cyclesperimage,phases,supersampling,bgLum,maskflag,maskRadius,maskstd,maskonedgeflag,edgemaskRadius,edgemaskstd,edgemasklum)

x = linspace(0,1,max(sizeX,sizeY)*supersampling);
x = x-mean(x);
y = linspace(0,1,max(sizeX,sizeY)*supersampling);
y = y-mean(y);

phase = phases(1);
% orientation = orientations(1) - pi/2;
orientation = -orientations(1);
L1 = @(x,y) cos(pi.*cyclesperimage.*(x.*cos(orientation)-y.*sin(orientation))+phase);

xx = linspace(-1,1,max(sizeX,sizeY)*supersampling);
xx = repmat(xx,max(sizeX,sizeY)*supersampling,1);
yy = xx';
rr = (xx.^2 + yy.^2).^0.5;

pattern = rr;
radius = radius/max(sizeX,sizeY)*supersampling;
idx = rr <= radius;
pattern(idx) = L1(xx(idx),yy(idx));

if numel(orientations) == 2 && numel(phases) == 2
    phase = phases(2);
    orientation = -orientations(2);
    L2 = @(x,y) cos(pi.*cyclesperimage.*(x.*cos(orientation)-y.*sin(orientation))+phase);
    idx = rr > radius;
    pattern(idx) = L2(xx(idx),yy(idx));
end


% orientation = orientations(1) - pi/2;
% phase = phases(1);
% L1 = @(x,y) cos(2 *pi*cyclesperimage*(x*cos(orientation)-y*sin(orientation))+phase);
% orientation = orientations(2) - pi/2;
% phase = phases(2);
% L2 = @(x,y) cos(2 *pi*cyclesperimage*(x*cos(orientation)-y*sin(orientation))+phase);
%
% pattern = zeros(max(sizeX,sizeY)*supersampling,max(sizeX,sizeY)*supersampling);
%
% for xi = 1:max(sizeX,sizeY)*supersampling
%     for yi = 1:max(sizeX,sizeY)*supersampling
%         if ((xi-sizeX/2*supersampling)^2+(yi-sizeY/2*supersampling)^2)^0.5 <= radius * supersampling
%             pattern(xi,yi) = L1(x(xi),y(yi));
%         else
%             pattern(xi,yi) = L2(x(xi),y(yi));
%         end
%     end
% end

% normalization
pattern = pattern-min(min(pattern));
pattern = pattern./max(max(pattern));
maxvalue = (contrast+1) * meanLum;
minvalue = (1-contrast) * meanLum;
pattern = pattern*(maxvalue-minvalue) + minvalue;

for xi = 1:max(sizeX,sizeY)*supersampling
    for yi = 1:max(sizeX,sizeY)*supersampling
        if abs(x(xi)) > sizeX/max(sizeX,sizeY)/2 || abs(y(yi)) > sizeY/max(sizeX,sizeY)/2
            pattern(xi,yi) = bgLum;
        end
    end
end

% outter mask
if maskflag
    maskstd = maskstd/max(sizeX,sizeY);
    maskRadius = maskRadius / max(sizeX,sizeY);
    g = @(x)exp(-(x-maskRadius).^2./(2.*maskstd.^2));
%     filter1dprofile = ones(max(sizeX,sizeY)*supersampling,1);
    xx = linspace(0,1,max(sizeX,sizeY)*supersampling) - 0.5;
    maskpattern = zeros(size(pattern,1),size(pattern,2));
    for i = 1:size(maskpattern,1)
        for j = 1:size(maskpattern,2)
            maskpattern(i,j) = (xx(i) ^ 2 + xx(j) ^2) ^ 0.5;
        end
    end
    idx = maskpattern < maskRadius;
    maskpattern = g(maskpattern);
    maskpattern(idx) = 1;
    idx = maskpattern < 0.01;
    maskpattern(idx) = 0;
    background = ones(size(pattern,1),size(pattern,2));
    background = background * bgLum;
    pattern = pattern .* maskpattern + background .* (1-maskpattern);
end

% mask on edge
if maskonedgeflag == 1
    edgemaskstd = 1/3/ max(sizeX,sizeY)*edgemaskstd;
    maskRadius = edgemaskRadius(1) / max(sizeX,sizeY);
    g = @(x)exp(-(x-maskRadius)^2/(2*edgemaskstd^2));
    
    filter1dprofile1 = ones(max(sizeX,sizeY)*supersampling,1);
    xx = linspace(0,1,max(sizeX,sizeY)*supersampling);
    for xi = 1:numel(filter1dprofile1)
        if xx(xi) > maskRadius
            filter1dprofile1(xi) = g(xx(xi));
        end
    end
    maskRadius = edgemaskRadius(2) / max(sizeX,sizeY);
    g = @(x)exp(-(x-maskRadius)^2/(2*edgemaskstd^2));
    
    filter1dprofile2 = ones(max(sizeX,sizeY)*supersampling,1);
    xx = linspace(0,1,max(sizeX,sizeY)*supersampling);
    for xi = 1:numel(filter1dprofile2)
        if xx(xi) < maskRadius
            filter1dprofile2(xi) = g(xx(xi));
        end
    end
    
    filter1dprofile = max(filter1dprofile1,filter1dprofile2);
    
    maskpattern = zeros(size(pattern,1),size(pattern,2));
    for xi = 1:max(sizeX,sizeY)*supersampling
        for yi = 1:max(sizeX,sizeY)*supersampling
            [~,idx] = min(abs(xx - sqrt(x(xi)^2+y(yi)^2)));
            maskpattern(xi,yi) = filter1dprofile(idx);
        end
    end
    background = ones(size(pattern,1),size(pattern,2));
    background = background * edgemasklum;
    pattern = pattern .* maskpattern + background .* (1-maskpattern);
elseif maskonedgeflag == 2
    edgemaskRadius = (edgemaskRadius+1) / max(sizeX,sizeY)*supersampling;
    xx = linspace(-1,1,max(sizeX,sizeY)*supersampling);
    xx = repmat(xx,max(sizeX,sizeY)*supersampling,1);
    %     yy = linspace(-1,1,max(sizeX,sizeY)*supersampling);
    %     yy = repmat(yy',1,max(sizeX,sizeY)*supersampling);
    yy = xx';
    rr = (xx.^2 + yy.^2).^0.5;
    rr = rr > edgemaskRadius(1) & rr < edgemaskRadius(2);
    pattern(rr) = edgemasklum;
elseif maskonedgeflag == 3 % center only
    edgemaskRadius = (edgemaskRadius+1) / max(sizeX,sizeY)*supersampling;
    xx = linspace(-1,1,max(sizeX,sizeY)*supersampling);
    xx = repmat(xx,max(sizeX,sizeY)*supersampling,1);
    %     yy = linspace(-1,1,max(sizeX,sizeY)*supersampling);
    %     yy = repmat(yy',1,max(sizeX,sizeY)*supersampling);
    yy = xx';
    rr = (xx.^2 + yy.^2).^0.5;
    rr = rr > edgemaskRadius(1);
    pattern(rr) = edgemasklum;
elseif maskonedgeflag == 4 % surround only
    %     edgemaskRadius = edgemaskRadius / max(sizeX,sizeY)*supersampling;
    edgemaskRadius = (edgemaskRadius+1) / max(sizeX,sizeY)*supersampling;
    xx = linspace(-1,1,max(sizeX,sizeY)*supersampling);
    xx = repmat(xx,max(sizeX,sizeY)*supersampling,1);
    %     yy = linspace(-1,1,max(sizeX,sizeY)*supersampling);
    %     yy = repmat(yy',1,max(sizeX,sizeY)*supersampling);
    yy = xx';
    rr = (xx.^2 + yy.^2).^0.5;
    rr = rr <= edgemaskRadius(1);
    pattern(rr) = edgemasklum;
end

pattern = imresize(pattern,1/supersampling,'bilinear');

pattern = rot90(pattern);
pattern = fliplr(pattern);