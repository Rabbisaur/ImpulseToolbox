function HistMat = WavDensity(Wav, Resolution, Fs, YLim, CStr)
% Resolution in ms; Fs in frequency/sec

NumWav = size(Wav, 1);
if Resolution<1/Fs*1000
    [InterpWav, InterpT] = Resampling(Wav, Fs, Resolution, 'cubic');
else
    InterpWav = Wav;
    InterpT = (0:size(Wav, 2)-1)*1/Fs*1000;
end
InterpMat = repmat(InterpT, [NumWav 1]);
YEdges = linspace(YLim(1), YLim(end), 100);
HistMat = Hist2(InterpMat(:,:), InterpWav(:,:), InterpT, YEdges);
% HistMat(HistMat==0) = nan;
MaxFreq = max(max(HistMat));
HPColor = imagesc(InterpT, YEdges, HistMat*5);
set(HPColor, 'CDataMapping', 'direct');
shading interp;
% CStr = 'k';
CData = Str2RGB(CStr);
CColor(:,1) = linspace(1, CData(1), min(256, MaxFreq));
CColor(:,2) = linspace(1, CData(2), min(256, MaxFreq));
CColor(:,3) = linspace(1, CData(3), min(256, MaxFreq));
colormap(CColor);
FreezeColors;
set(gca, 'YDir', 'normal');
axis square;
return;

function CData = Str2RGB(CStr)
switch CStr
    case 'k'
        CData = [0 0 0];
    case 'b'
        CData = [0 0 1];
    case 'r'
        CData = [1 0 0];
    case 'g'
        CData = [0 1 0];
    case 'c'
        CData = [0 1 1];
    case 'm'
        CData = [1 0 1];
end

function HistMat  = Hist2(x, y, xedges, yedges)
% function histmat  = Hist2(x, y, xedges, yedges)
% Extract 2D histogram data containing the number of events
% of [x , y] pairs that fall in each bin of the grid defined by 
% xedges and yedges. The edges are vectors with monotonically 
% non-decreasing values.
if nargin~=4
    error ('The four input arguments are required!');
end
if any(size(x) ~= size(y))
    error ('The size of the two first input vectors should be same!');
end
%# bin centers (integers)
xNumBins = numel(xedges);
yNumBins = numel(yedges);
%# map X/Y values to bin indices
Xi = round( interp1(xedges, 1:xNumBins, x, 'nereast') );
Yi = round( interp1(yedges, 1:yNumBins, y, 'nereast') );
%# limit indices to the range [1,numBins]
Xi = max( min(Xi,xNumBins), 1);
Yi = max( min(Yi,yNumBins), 1);
%# count number of elements in each bin
HistMat = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);