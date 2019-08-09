% common parameters

patternsize = 6;
drawsize = 8;

% sizeX = 8;
% sizeY = 8;
maskstd = 0.3;
pixelperdegree = 21;

sizeX = drawsize * pixelperdegree;
sizeY = drawsize * pixelperdegree;
maskRadius = patternsize * pixelperdegree/2;
maskstd = maskstd * pixelperdegree;

meanLum = 0.5;
supersampling = 16;
bgLum = 0.5;
maskflag = 1;
notes = '';
dmns=[8, sizeX, sizeY, 1];
contrastpct = 30;
contrast = contrastpct /100;
%% Cartesian gratings
spatialfreq = 2;
cyclesperimage = spatialfreq * drawsize;
phase = 0;
NumOri = 8;
for thisOri = 1:NumOri
    orientation = (thisOri-1)*pi/NumOri;
    pattern = cartesiangratings(sizeX,sizeY,meanLum,contrast,orientation,cyclesperimage,phase,supersampling,bgLum,maskflag,maskRadius,maskstd);
    imwrite(pattern,['./image/cart',num2str(thisOri),'_',num2str(contrastpct),'.bmp'],'bmp')
    filename = ['./ctx/cart',num2str(thisOri),'_',num2str(contrastpct),'.ctx'];
    pattern = round(pattern*127)+128;
    pattern = uint8(pattern);
    savecx(filename, notes, dmns, pattern)
end

%% polar

% radial
concentricspatialfreq = 0;
spatialfreq = 2;
radialspatialfreq = round(spatialfreq*pi*patternsize/2);
phase = 0;
pattern = polargratings(sizeX,sizeY,meanLum,contrast,concentricspatialfreq,radialspatialfreq,phase,supersampling,bgLum,maskflag,maskRadius,maskstd);
imwrite(pattern,['./image/radial0','_',num2str(contrastpct),'.bmp'],'bmp')
filename = ['./ctx/radial0','_',num2str(contrastpct),'.ctx'];
pattern = round(pattern*127)+128;
pattern = uint8(pattern);
savecx(filename, notes, dmns, pattern)

% target
radialspatialfreq = 0;
spatialfreq = 2;
concentricspatialfreq = round(spatialfreq*drawsize);
phase = 0;
pattern = polargratings(sizeX,sizeY,meanLum,contrast,concentricspatialfreq,radialspatialfreq,phase,supersampling,bgLum,maskflag,maskRadius,maskstd);
imwrite(pattern,['./image/target0','_',num2str(contrastpct),'.bmp'],'bmp')
filename = ['./ctx/target0','_',num2str(contrastpct),'.ctx'];
pattern = round(pattern*127)+128;
pattern = uint8(pattern);
savecx(filename, notes, dmns, pattern)


% spiral
spatialfreq = 2;
radialspatialfreq = round(spatialfreq*pi*patternsize/2);
concentricspatialfreq = round(spatialfreq*drawsize);
phase = 0;
pattern = polargratings(sizeX,sizeY,meanLum,contrast,concentricspatialfreq,radialspatialfreq,phase,supersampling,bgLum,maskflag,maskRadius,maskstd);
imwrite(pattern,['./image/spiral0','_',num2str(contrastpct),'.bmp'],'bmp')
filename = ['./ctx/spiral0','_',num2str(contrastpct),'.ctx'];
pattern = round(pattern*127)+128;
pattern = uint8(pattern);
savecx(filename, notes, dmns, pattern)
% mirror spiral
concentricspatialfreq = -concentricspatialfreq;
pattern = polargratings(sizeX,sizeY,meanLum,contrast,concentricspatialfreq,radialspatialfreq,phase,supersampling,bgLum,maskflag,maskRadius,maskstd);
imwrite(pattern,['./image/spiral1','_',num2str(contrastpct),'.bmp'],'bmp')
filename = ['./ctx/spiral1','_',num2str(contrastpct),'.ctx'];
pattern = round(pattern*127)+128;
pattern = uint8(pattern);
savecx(filename, notes, dmns, pattern)

%% hyperbolic
spatialfreq = 2;
frequency = spatialfreq * drawsize+6;
orientation = pi/6;
phase = pi;
NumOri = 4;
for thisOri = 1:NumOri
    orientation = (thisOri-1) * pi / 2 / NumOri;
    pattern = hyperbolicgratings(sizeX,sizeY,meanLum,contrast,frequency,orientation,phase,supersampling,bgLum,maskflag,maskRadius,maskstd);
    imwrite(pattern,['./image/hyper',num2str(thisOri),'_',num2str(contrastpct),'.bmp'],'bmp')
    filename = ['./ctx/hyper',num2str(thisOri),'_',num2str(contrastpct),'.ctx'];
    pattern = round(pattern*127)+128;
    pattern = uint8(pattern);
    savecx(filename, notes, dmns, pattern)
end