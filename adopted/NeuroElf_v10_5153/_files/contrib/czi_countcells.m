function count = czi_countcells(imfile, chn)

% requires neuroelf
n = neuroelf;

% input check
if nargin < 1 || ...
   ~ischar(imfile) || ...
    isempty(imfile)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2
    chn = 4;
end

% load image
try
    czi = xff(imfile);
    if numel(czi) ~= 1 || ...
       ~isxff(czi, 'czi')
        if numel(czi) == 1 && ...
            isxff(czi, true)
            czi.ClearObject;
        end
        error( ...
            'neuroelf:BadFileContent', ...
            'Bad CZI file content.' ...
        );
    end
    
    % make sure we have the channels we need
    ch = cell(1, 4);
    chc = 1;
    for cc = 1:numel(czi.RawContent)
        if ~isempty(czi.RawContent(cc).FileType) && ...
            strcmpi(czi.RawContent(cc).FileType, 'dvimage') && ...
            isstruct(czi.RawContent(cc).FileCooked) && ...
            isfield(czi.RawContent(cc).FileCooked, 'Image')
            ch{chc} = czi.RawContent(cc).FileCooked.Image;
            chc = chc + 1;
            if chc > 4
                break;
            end
        end
    end
    czi.ClearObject;
    if chc <= 4
        error( ...
            'neuroelf:BadFileContent', ...
            'Image must have at least 4 channels.' ...
        );
    end
catch ne_eo;
    rethrow(ne_eo);
end

% copy the channel we need
countc = double(ch{chn});
count = zeros(size(countc));

% minimal smoothing
ccs = [Inf, Inf; 1, 1; 1, 1; size(countc)];
smkmin = n.smoothkern(1.25);
smkmin = smkmin(smkmin>1e-4);
scountc = n.flexinterpn(countc, ccs, smkmin, 1, 0);

% remove bright spots (and their neighbors)
countc(scountc >= 240) = NaN;

% begin thresholding
tc = 1;
lsk = 0;
sks = 30;
skl = n.smoothkern(3.5);
skl(skl < 0.001) = [];
lcountc = n.flexinterpn(n.flexinterpn(countc, ...
    ccs, {skl, [0;1;0]}, {1, 1}, 0), ...
    ccs, {[0;1;0], skl}, {1, 1}, 0);
keeplooking = true;
while keeplooking

    % smooth data with current kernel
    if lsk ~= sks
        sk = n.smoothkern(sks);
        sk = sk(sk >= (0.001 * max(sk)));
        lsk = sks;
    end
    scountc = n.flexinterpn(n.flexinterpn(countc, ...
        ccs, {sk, [0;1;0]}, {1, 1}, 0), ...
        ccs, {[0;1;0], sk}, {1, 1}, 0);
    
    % find maximum value
    mval = max(scountc(:));
    
    % threshold image at 60%
    tcountc = (scountc >= (0.6 * mval) & count == 0);
    
    % cluster
    [cs, cv, ccl] = n.clustercoordsc(tcountc, 4);
    
    % for each cluster
    for cc = 1:numel(cs)
        
        % at least 250 pixels
        if cs(cc) >= (100 * sqrt(sks)) && ...
            cs(cc) <= 5000
            cvcc = find(cv(:) == cc);
            cvci = find(ccl(:, 4) == cc);
            cclc = ccl(cvci, 1:2);
            mcc = mean(cclc);
            mincrd = min(cclc, [], 1);
            maxcrd = max(cclc, [], 1);
            boxdata = lcountc(mincrd(1):maxcrd(1), mincrd(2):maxcrd(2));
            cvccdata = boxdata(sub2ind(size(boxdata), cclc(:, 1) - (mincrd(1) - 1), cclc(:, 2) - (mincrd(2) - 1)));
            fprintf('Cluster %03d: %5d pixels around (%f, %f)\n', tc, cs(cc), mcc(1), mcc(2));
            count(cvcc) = tc;
            countc(cvcc) = -abs(countc(cvcc));
            pause(0.01);
            tc = tc + 1;
        end
    end
    n.scaleimage(countc + abs(countc));
    
    % reduce size
    if mval < 24 || ...
        isempty(cs) || ...
        all(cs < (100 * sqrt(sks)) | cs > 5000)
        sks = 0.875 * sks;
        if sks < 3
            keeplooking = false;
        end
    end
end

% for each of the clusters, print out the mean of each of the channels
for cc = 1:(tc - 1)
    ccc = find(count(:) == cc);
    fprintf('Cluster %03d: Ch1: %5.1f / Ch2: %5.1f / Ch3: %5.1f / Ch4: %5.1f\n', ...
        cc, mean(ch{1}(ccc)), mean(ch{2}(ccc)), mean(ch{3}(ccc)), mean(ch{4}(ccc)));
end
