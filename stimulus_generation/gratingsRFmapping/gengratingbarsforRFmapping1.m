clear,clc, close all
sizeH = 8; % degrees
sizeV = 8; % degrees
meanLum = 0.5;
contrast = 0.99;
orientation = pi/4;
spatialfreqency = 2;
phase = 0;
supersampling = 4;
bgLum = 0.5;
maskflag=0;
maskRad=0;
maskstd=0;
barwidth = 0.3;
barlength = 8;

pixelsperdegree = 21;
barwidthP = round(barwidth*pixelsperdegree);
barlengthP = round(barlength*pixelsperdegree);
sizeX = sizeH * pixelsperdegree;
sizeY = sizeV * pixelsperdegree;
cyclesperimage = spatialfreqency * sizeH;
notes = '';

pattern = cartesiangratings(sizeX,sizeY,meanLum,contrast,orientation,cyclesperimage,phase,supersampling,bgLum,maskflag,maskRad,maskstd);

vertbar1 = pattern(1:barlengthP,1:barwidthP);
vertbar2 = fliplr(vertbar1);
horzbar2 = rot90(vertbar1);
horzbar1 = flipud(horzbar2);

imwrite(vertbar1,'./image/vertbar1.bmp','bmp');
vertbar1 = round(vertbar1*127)+128;
vertbar1 = uint8(vertbar1);
dmns=[8, barwidthP, barlengthP, 1];
savecx('./ctx/vertbar1.ctx', notes, dmns, vertbar1)

imwrite(vertbar2,'./image/vertbar2.bmp','bmp');
vertbar2 = round(vertbar2*127)+128;
vertbar2 = uint8(vertbar2);
dmns=[8, barwidthP, barlengthP, 1];
savecx('./ctx/vertbar2.ctx', notes, dmns, vertbar2)

imwrite(horzbar1,'./image/horzbar1.bmp','bmp');
horzbar1 = round(horzbar1*127)+128;
horzbar1 = uint8(horzbar1);
dmns=[8, barlengthP, barwidthP, 1];
savecx('./ctx/horzbar1.ctx', notes, dmns, horzbar1)

imwrite(horzbar2,'./image/horzbar2.bmp','bmp');
horzbar2 = round(horzbar2*127)+128;
horzbar2 = uint8(horzbar2);
notes = '';
dmns=[8, barlengthP, barwidthP, 1];
savecx('./ctx/horzbar2.ctx', notes, dmns, horzbar2)


% gen orientation mapping movies
maskflag=1;
maskRad=4;
maskRad = maskRad*pixelsperdegree;
maskstd=0;
NumOri = 17;
refreshRate = 85; % Hz
frametime = 1/refreshRate; % seconds
movielength = 300; % ms
NumFrames = round(movielength*refreshRate/1000);
actualmovielength = NumFrames * frametime;
temporalfreq = 2; % cycles / second
deltaphase = temporalfreq * 2 * pi * frametime;
OrientationDiff = 2*pi/NumOri;
for thisOri = 1:NumOri
    orientation = (thisOri-1)* OrientationDiff;
    filename = ['./ctx/grat',num2str(round(orientation/2/pi*360)),'.ctx'];
    pattern = zeros(sizeX,sizeY,NumFrames);
    for thisFrame = 1:NumFrames
        phase = (thisFrame-1)*deltaphase;
        pattern(:,:,thisFrame) = cartesiangratings(sizeX,sizeY,meanLum,contrast,orientation,cyclesperimage,phase,supersampling,bgLum,maskflag,maskRad,maskstd);
    end
    dmns = [8, sizeX,sizeY,NumFrames];
    pattern = round(pattern*127)+128;
    pattern = uint8(pattern);
    savecx_movie_firstframe(filename, notes, dmns, pattern(:,:,1))
    for i = 2:NumFrames
        savecx_movie_succframe(filename, notes, dmns, pattern(:,:,i))
    end
end
