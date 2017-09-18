function HistMat = APDensity(Coeff, CStr)

% Contours
warning('off');
Gauss = gmdistribution.fit(Coeff(:,1:2), 1); clc;
[X, Y] = meshgrid(-80:60);
PDFVal = pdf(Gauss,[X(:) Y(:)]);
PDFVal = reshape(PDFVal, [size(X, 1) size(Y, 1)]);
ContData = contourc(X(1,:), Y(:,1), PDFVal, [.99 .5 .25]*max(max(PDFVal)));
CData = Str2RGB(CStr);
CColor(:,1) = linspace(1, CData(1), 100);
CColor(:,2) = linspace(1, CData(2), 100);
CColor(:,3) = linspace(1, CData(3), 100);
while 1
    if isempty(ContData)
        break;
    end
    PVal = ContData(1,1);
    NumPoint = ContData(2,1);
    XY = ContData(:, 1+(1:NumPoint));
    plot(XY(1,:), XY(2,:), 'color', CColor(round(PVal/max(max(PDFVal))*100),:));
    ContData(:,1:NumPoint+1) = [];
end

% %Imagesc histogram
% MaxWav = max(Coeff, [], 1);
% MinWav = min(Coeff, [], 1);
% XYLim = [MinWav-7*std(Coeff);  MaxWav+7*std(Coeff)];
% XEdges = linspace(XYLim(1,1), XYLim(2,1), 100);
% YEdges = linspace(XYLim(1,2), XYLim(2,2), 100);
% HistMat  = Hist2(Coeff(:,1), Coeff(:,2), XEdges, YEdges);
% HistMat(HistMat==0) = nan;
% HPColor = imagesc(XEdges, YEdges, HistMat*5);
% set(HPColor, 'CDataMapping', 'direct');
% shading interp;
% colormap(CColor);
% FreezeColors;
% set(gca, 'YDir', 'normal');

% %Scatter histogram
% [XEdges, YEdges] = meshgrid(XEdges, YEdges);
% scatter(XEdges(HistMat~=0), YEdges(HistMat~=0), 15, CColor(HistMat(HistMat~=0),:), 'fill', 'o');
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
    case 'y'
        CData = [1 1 0];
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
Xi = round( interp1(xedges, 1:xNumBins, x, 'cubic') );
Yi = round( interp1(yedges, 1:yNumBins, y, 'cubic') );
%# limit indices to the range [1,numBins]
Xi = max( min(Xi,xNumBins), 1);
Yi = max( min(Yi,yNumBins), 1);
%# count number of elements in each bin
HistMat = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);