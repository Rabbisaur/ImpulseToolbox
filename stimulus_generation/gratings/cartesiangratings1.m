function pattern = cartesiangratings(sizeX,sizeY,meanLum,contrast,orientation,cyclesperimage,phase,supersampling,bgLum,maskflag,maskRadius,maskstd)

x = linspace(0,1,max(sizeX,sizeY)*supersampling);
x = x-mean(x);
y = linspace(0,1,max(sizeX,sizeY)*supersampling);
y = y-mean(y);


orientation = orientation - pi/2;

L = @(x,y) cos(2 *pi*cyclesperimage*(x*cos(orientation)-y*sin(orientation))+phase);

pattern = zeros(max(sizeX,sizeY)*supersampling,max(sizeX,sizeY)*supersampling);

for xi = 1:max(sizeX,sizeY)*supersampling
    for yi = 1:max(sizeX,sizeY)*supersampling
        pattern(xi,yi) = L(x(xi),y(yi));
    end
end

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

if maskflag
    maskstd = 1/3/ max(sizeX,sizeY)*maskstd;
    maskRadius = maskRadius / max(sizeX,sizeY);
    g = @(x)exp(-(x-maskRadius)^2/(2*maskstd^2));
    filter1dprofile = ones(max(sizeX,sizeY)*supersampling,1);
    xx = linspace(0,1,max(sizeX,sizeY)*supersampling);
    for xi = 1:numel(filter1dprofile)
        if xx(xi) > maskRadius
            filter1dprofile(xi) = g(xx(xi));
        end
    end
    maskpattern = zeros(size(pattern,1),size(pattern,2));
    for xi = 1:max(sizeX,sizeY)*supersampling
        for yi = 1:max(sizeX,sizeY)*supersampling
            [~,idx] = min(abs(xx - sqrt(x(xi)^2+y(yi)^2)));
            maskpattern(xi,yi) = filter1dprofile(idx);
        end
    end
    background = ones(size(pattern,1),size(pattern,2));
    background = background * bgLum;
    pattern = pattern .* maskpattern + background .* (1-maskpattern);
end

pattern = imresize(pattern,1/supersampling,'bilinear');

pattern = rot90(pattern);
pattern = fliplr(pattern);