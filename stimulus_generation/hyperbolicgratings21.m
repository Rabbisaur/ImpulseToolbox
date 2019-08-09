function pattern = hyperbolicgratings2(sizeX,sizeY,meanLum,contrast,radius,frequency,orientation,phase,supersampling,bgLum,maskflag,maskRadius,maskstd,maskonedgeflag,edgemaskRadius,edgemaskstd,edgemasklum)

x = linspace(0,1,max(sizeX,sizeY)*supersampling);
x = x-mean(x);
y = linspace(0,1,max(sizeX,sizeY)*supersampling);
y = y-mean(y);


L1 = @(x,y)cos(2.*pi.*frequency(1).*abs(((x.*cos(orientation(1))-y.*sin(orientation(1))).*(x.*sin(orientation(1))+y.*cos(orientation(1)))).^0.5)+phase(1));

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

xx = linspace(-1,1,max(sizeX,sizeY)*supersampling);
xx = repmat(xx,max(sizeX,sizeY)*supersampling,1);
yy = xx';
rr = (xx.^2 + yy.^2).^0.5;

pattern = rr;
radius = radius/max(sizeX,sizeY)*supersampling;
idx = rr <= radius;
pattern(idx) = L1(xx(idx),yy(idx));

if numel(phase) == 2
    L2 = @(x,y)cos(2.*pi.*frequency(2).*abs(((x.*cos(orientation(2))-y.*sin(orientation(2))).*(x.*sin(orientation(2))+y.*cos(orientation(2)))).^0.5)+phase(2));
    idx = rr > radius;
    pattern(idx) = L2(xx(idx),yy(idx));
end


pattern = pattern-min(min(pattern));
pattern = pattern./max(max(pattern));
maxvalue = (contrast+1) * meanLum;
minvalue = (1-contrast) * meanLum;
pattern = pattern*(maxvalue-minvalue) + minvalue;

x = linspace(0,1,max(sizeX,sizeY)*supersampling);
x = x-mean(x);
y = linspace(0,1,max(sizeX,sizeY)*supersampling);
y = y-mean(y);

for xi = 1:max(sizeX,sizeY)*supersampling
    for yi = 1:max(sizeX,sizeY)*supersampling
        if abs(x(xi)) > sizeX/max(sizeX,sizeY)/2 || abs(y(yi)) > sizeY/max(sizeX,sizeY)/2
            pattern(xi,yi) = bgLum;
        end
    end
end

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
if maskonedgeflag
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
end

pattern = imresize(pattern,1/supersampling,'bilinear');
pattern = rot90(pattern);
