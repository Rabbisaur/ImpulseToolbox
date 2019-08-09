function [locs, Pprominance]=peakseek(x,minpeakdist,peakwidth,sortflag)
baselength = 10;

% Find all maxima and ties
% tmp = diff(x);
% idx = tmp(1:end-1) > 0 & tmp(2:end) <=0;
% locs=find(idx)+1;
tmp = x(2:(end-1));
idx = tmp - x(1:end-2) > 0 & tmp - x(3:end) >= 0;
locs=find(idx)+1;

locs(locs < peakwidth + baselength) = [];
locs(locs + peakwidth + baselength > numel(x)) = [];

% find peak base
Pbase = zeros(numel(locs),1);
Pval = x(locs);
idx = 1:baselength;
for thisloc = 1:numel(locs)
    leftbaseIdx = round(locs(thisloc) - peakwidth/2 - idx);
    rightbaseIdx = round(locs(thisloc) + peakwidth/2 + idx);
    leftbase = x(leftbaseIdx);
    rightbase = x(rightbaseIdx);
    Pbase(thisloc) = sum([leftbase;rightbase])/baselength/2;
end
% if size(Pval,2) > size(Pval,1)
%     Pval = Pval';
% end
Pprominance = Pval - Pbase;
[~, idx] = sort(Pprominance,'descend');
locsSorted = locs(idx);

if minpeakdist>1
    counter = 0;
    while 1
        % starting from the highest and eliminate peaks that are too close
        counter = counter + 1;
        if counter > numel(locsSorted)
            break;
        end
        HiPeak = locsSorted(counter);
        eliminateIdx = abs(HiPeak - locsSorted) < minpeakdist & abs(HiPeak - locsSorted) ~=0;
        locsSorted(eliminateIdx) = [];
        idx(eliminateIdx) = [];
    end
end

% sort locs according to peak prominance
if sortflag
    Pprominance = Pprominance(idx);
    locs = locsSorted;
else
    locs = locs(sort(idx,'ascend'));
    Pprominance = Pprominance(sort(idx,'ascend'));
end
% if size(locs,2) > size(locs,1)
%     locs = locs';
% end
% if size(Pprominance,2) > size(Pprominance,1)
%     Pprominance = Pprominance';
% end

end