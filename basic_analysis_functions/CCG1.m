function [fRawCCG fShuffleCCG fShuffleCi] = CCG(iTrain1, iTrain2, iNormMethod, iShufMethod, nShuffle, iJitterWin)
% Compute CCG from spike trains of 2 neurons
%   Usage: [fRawCCG fShuffleCCG fShuffleCi] = CCG(iTrain1, iTrain2, iNormMethod, iShufMethod, nShuffle, iJitterWin)
%
% Inputs
%   iTrain1 - trial x SpikeCount train of neuron 1
%   iTrain2 - trial x SpikeCount train of neuron 2
%   iNormMethod - Normlization method
%     0/1/2/3/4--length/count1/count2/Geocount/number of spikes
%   iShufMethod - shuffle to correct stimulus-locked synchrony
%     0/1/2/3/4/5=>none/JPSTH/trial shuffle/jitter/pattern shift/re-assign
%   nShuffle - number of shuffle
%   iJitterWin - shuffle within a fixed number of points
%
% Outputs:
%   fRawCCG - raw CCG for each trial
%   fShuffleCCG - shuffle CCG
%   fShuffleCi - 95% two-side confidence interval
%
% External Includes:
%
% History:
%   create: 10/27/2011, MGChen@BNU
%   revised: 3/15/2012, MGChen&FWang@BNU

%% Check inputs
switch nargin
    case 2
        iNormMethod = 3;
        iShufMethod = 0;
        nShuffle = 50;
        iJitterWin = 25;
    case 3
        iShufMethod = 0;
        nShuffle = 0;
        iJitterWin = 25;
    case 4
        nShuffle = 0;
        iJitterWin = 25;
    case 5
        iJitterWin = 25;
    case 6
        %nothing
    otherwise
        disp('FAILURE: CCG accepts 2/3/4/5/6 input arguments.'); return;
end
assert(iJitterWin>=1, 'The iJitterWin must be larger than 0.');
%% Compute CCG & confidence interval
iTrain1 = single(iTrain1);
iTrain2 = single(iTrain2);
[nRow nCol1] = size(iTrain1);
[nRow nCol2] = size(iTrain2);
iJitterWin = min(round(iJitterWin), min(nCol1, nCol2));
switch iShufMethod
    case {0 1 2} %none/JPSTH/trial shuffle
        [fRawCCG fShuffleCCG fShuffleCi] = GetCCG(iTrain1, iTrain2, iNormMethod, iShufMethod, nShuffle); 
    case 3  %spike jitter
        fRawCCG = GetCCG(iTrain1, iTrain2, iNormMethod);
        iTrain = {iTrain1; iTrain2};
        fSampleCCG = zeros([nRow nCol1+nCol2-1 nShuffle], 'single');        
        for i = 1:nShuffle
            iShuffle = {rand(size(iTrain1), 'single'); rand(size(iTrain2), 'single')};
            for n = 1:2
                for j = 1:ceil(size(iTrain{n}, 2)/iJitterWin)
                    iTrainIdx = ((j-1)*iJitterWin+1) : min(j*iJitterWin, size(iTrain{n}, 2));
                    [B iShuffleIdx] = sort(iShuffle{n}(:, iTrainIdx), 2);
                    iShuffle{n}(:,iTrainIdx) = iShuffleIdx<=(sum(iTrain{n}(:,iTrainIdx), 2)*ones(1, numel(iTrainIdx)));
                end
            end
            fSampleCCG(:,:,i) = GetCCG(iShuffle{1}, iShuffle{2}, iNormMethod);
        end
        fSampleCCG = sort(fSampleCCG, 3);
        fShuffleCi = zeros([nRow nCol1+nCol2-1 2], 'single');
        fShuffleCi(:,:,1) = fSampleCCG(:,:,max(1, floor(nShuffle*0.025)), :);
        fShuffleCi(:,:,2) = fSampleCCG(:,:,max(1, floor(nShuffle*0.975)), :);
        fShuffleCCG = mean(fSampleCCG, 3);
    case 4  %pattern shift
        fRawCCG = GetCCG(iTrain1, iTrain2, iNormMethod);
        iTrain = {iTrain1; iTrain2};
        fSampleCCG = zeros([nRow nCol1+nCol2-1 nShuffle], 'single');  
        for i = 1:nShuffle
            iShuffle = {zeros(size(iTrain1), 'single'); zeros(size(iTrain2), 'single')};
            for n = 1:2
                for j = 1:ceil(size(iTrain{n}, 2)/iJitterWin)
                    iTrainIdx = ((j-1)*iJitterWin+1):min(j*iJitterWin, size(iTrain{n}, 2));
                    for k = 1:nRow
                        if sum(iTrain{n}(k, iTrainIdx))>0
                            iLoLim = find(iTrain{n}(k, iTrainIdx)==1, 1, 'first');
                            iHiLim = find(iTrain{n}(k, iTrainIdx)==1, 1, 'last');
                            iShift = (1-iLoLim):(numel(iTrainIdx)-iHiLim);
                            iShift = iShift(randi(numel(iShift), 1));
                            iShuffle{n}(k, iTrainIdx(1)-1+(iLoLim:iHiLim)+iShift) = iTrain{n}(k, iTrainIdx(1)-1+(iLoLim:iHiLim));
                        end
                    end
                end
            end               
            fSampleCCG(:,:,i) = GetCCG(iShuffle{1}, iShuffle{2}, iNormMethod);
        end
        fSampleCCG = sort(fSampleCCG, 3);
        fShuffleCi = zeros([nRow nCol1+nCol2-1 2], 'single');
        fShuffleCi(:,:,1) = fSampleCCG(:,:,max(1, floor(nShuffle*0.025)), :);
        fShuffleCi(:,:,2) = fSampleCCG(:,:,max(1, floor(nShuffle*0.975)), :);
        fShuffleCCG = mean(fSampleCCG, 3);
    case 5  %spike assignment, support maximum of 1 spike/bin
        fRawCCG = GetCCG(iTrain1, iTrain2, iNormMethod);
        iTrain = {iTrain1; iTrain2};
        fSampleCCG = zeros([nRow nCol1+nCol2-1 nShuffle], 'single');
        for i = 1:nShuffle
            iShuffle = {zeros(size(iTrain{1}), 'single'); zeros(size(iTrain{2}), 'single')};
            for n = 1:2
                for j = 1:ceil(size(iTrain{n}, 2)/iJitterWin)                
                    idx = ((j-1)*iJitterWin+1):min(j*iJitterWin, size(iTrain{n}, 2));                
                    nspkcol = sum(iTrain{n}(:, idx), 1);
                    nspkrow = sum(iTrain{n}(:, idx), 2);
                    while true
                        shuf = zeros([nRow numel(idx)]);
                        for k = 1:numel(idx)
                            spkinrow = mod( sum(shuf, 2), 10^6 );  %10^6==0
                            zeroidx = shuf==0;
                            shuf( logical( ((spkinrow +sum(zeroidx, 2))==nspkrow) *ones(1, numel(idx)) .*zeroidx) ) = 1;
                            shuf( logical( (spkinrow==nspkrow) *ones(1, numel(idx)) .*zeroidx) ) = 10^6;
                            spkincol = mod( sum(shuf, 1), 10^6 );
                            zeroidx = shuf==0;
                            shuf( logical(ones(nRow, 1) *((spkincol +sum(zeroidx, 1))==nspkcol) .*zeroidx) ) = 1;
                            shuf( logical(ones(nRow, 1) *(spkincol==nspkcol) .*zeroidx) ) = 10^6;
                            spkinrow = mod( sum(shuf, 2), 10^6 );  %10^6==0
                            zeroidx = shuf==0;
                            shuf( logical( ((spkinrow +sum(zeroidx, 2))==nspkrow) *ones(1, numel(idx)) .*zeroidx) ) = 1;
                            shuf( logical( (spkinrow==nspkrow) *ones(1, numel(idx)) .*zeroidx) ) = 10^6;
                            zeroidx = shuf==0;
                            if sum(sum(zeroidx))==0
                                break;
                            end
                            zeroidx = find(shuf(:,k)==0);
                            if ~isempty(zeroidx)  %at least 1 bin
                                zeroidx = zeroidx( randperm(numel(zeroidx)) );
                                shuf( zeroidx(1 : min(nspkcol(k)-mod(sum(shuf(:, k)), 10^6), numel(zeroidx))), k ) = 1;
                                shuf(shuf(:, k)==0, k) = 10^6;
                            end
                        end
                        shuf( shuf==10^6 ) = 0;
                        if isequal(nspkcol, sum(shuf, 1)) && isequal(nspkrow, sum(shuf, 2))
                            iShuffle{n}(:, idx) = shuf;
                            break;
                        end
                    end
                end
            end
            fSampleCCG(:,:,i) = GetCCG(iShuffle{1}, iShuffle{2}, iNormMethod);
        end
        fSampleCCG = sort(fSampleCCG, 3);
        fShuffleCi = zeros([nRow nCol1+nCol2-1 2], 'single');
        fShuffleCi(:,:,1) = fSampleCCG(:,:,max(1, floor(nShuffle*0.025)), :);
        fShuffleCi(:,:,2) = fSampleCCG(:,:,max(1, floor(nShuffle*0.975)), :);
        fShuffleCCG = mean(fSampleCCG, 3);
    otherwise
        error('The shuffle METHOD was not found.');
end

%=======================================================================
% Functions called, do CCG
%=======================================================================
function [fRawCCG fShiftCCG fShiftCi] = GetCCG(iTrain1, iTrain2, iNormMethod, ShifPredictor, nShift)
% Compute CCG(cross-correlogram) from spike trains of 2 neurons
%   Usage: [fCCG fShiftCCG fShiftCi] = get_ccg(iTrain1, iTrain2, iNormMethod, bShift, nShift)
%
% Inputs
%   iTrain1 - trial x SpikeCount train of neuron 1
%   iTain2 - trial x SpikeCount train of neuron 2
%   iNormMethod - Normlization method
%     0/1/2/3/4--length/count1/count2/Geocount/number of spikes
%   ShifPredictor - 0/1/2 - none/JPSTH/theretical shift predictor
%   nShift - number of shift
%
% Outputs:
%   fRawCCG - raw CCG for each trial
%   fShiftCCG - theretical shift predictor
%   fShiftCi - 95% two-side confidence interval for shift predictor
%
% External Includes:
%
% History:
%   create: 10/27/2011, MGChen@BNU
%   revised: 3/15/2012, MGChen&FWang@BNU
%% Check inputs
switch nargin
    case 2
        iNormMethod = 3;
        ShifPredictor = 0;
        nShift = 0;
    case 3
        ShifPredictor = 0;
        nShift = 0;
    case 4
        nShift = 0;
    case 5
        %nothing
    otherwise
        disp('FAILURE: Load_nev accepts 2/3/4/5 input arguments.'); return;
end
if ~isnumeric(iTrain1) || ~isnumeric(iTrain2)
    error('TRAIN1&TRAIN2 must be a numeric array.');
end
if size(iTrain1, 1)~=size(iTrain2, 1)
    error('TRAIN1&TRAIN2 must have same number of rows.');
end
if sum(sum(isnan(iTrain1)))~=0 || sum(sum(isnan(iTrain2)))~=0
    error('Some elements in TRAIN1/TRAIN2 are NAN.');
end
if sum(sum(iTrain1<0))~=0 || sum(sum(iTrain2<0))~=0
    error('Some elements in TRAIN1/TRAIN2 <0.');
end
if sum(sum(iTrain1>1))~=0 || sum(sum(iTrain2>1))~=0
%     error('Some elements in TRAIN1/TRAIN2 >1.');
end
%% Compute CCG, drag iTrain2
iTrain1 = single(iTrain1);
iTrain2 = single(iTrain2);
iTrain2lr = fliplr(iTrain2);
[nRow nCol1] = size(iTrain1);
[nRow nCol2] = size(iTrain2);
fRawCCG = nan(nRow, nCol2+nCol1-1, 'single');
for i = 1:nRow
    fRawCCG(i,:) = conv2(iTrain1(i,:), iTrain2lr(i,:));
end
if nCol1==nCol2
    iCount1 = [cumsum(iTrain1, 2), fliplr(cumsum(fliplr(iTrain1(:,2:end)), 2))];
    iCount2 = [cumsum(iTrain2lr, 2), fliplr(cumsum(iTrain2(:,1:end-1), 2))];
elseif nCol1>nCol2
    iCount1 = [cumsum(iTrain1, 2), fliplr(cumsum(fliplr(iTrain1(:,(end-nCol2+2):end)), 2))];
    iMinusCol = (nCol2+1):nCol1;
    iCount1(:,iMinusCol) = iCount1(:,iMinusCol) - iCount1(:,iMinusCol-nCol2);
    iCount2 = [cumsum(iTrain2lr, 2), sum(iTrain2, 2)*ones(1, nCol1-nCol2-1), fliplr(cumsum(iTrain2, 2))];    
elseif nCol1<nCol2
    iCount1 = [cumsum(iTrain1, 2), sum(iTrain1, 2)*ones(1, nCol2-nCol1-1), fliplr(cumsum(fliplr(iTrain1), 2))];
    iCount2 = [cumsum(iTrain2lr, 2), fliplr(cumsum(iTrain2(:,1:(nCol1-1)), 2))];
    iMinusCol = (nCol1+1):nCol2;
    iCount2(:,iMinusCol) = iCount2(:,iMinusCol) - iCount2(:,iMinusCol-nCol1);
end
switch iNormMethod
    case 0
        iMinCol = min(nCol1, nCol2);
        fNorm = ones([nRow 1])*[1:iMinCol ones([1 nCol1+nCol2-2*iMinCol])*iMinCol (iMinCol-1):-1:1];
    case 1
        fNorm = iCount1;
    case 2
        fNorm = iCount2;
    case 3
        fNorm = sqrt(iCount1.*iCount2);
    case 4
        fNorm = iCount1+iCount2-fRawCCG;
end
fRawCCG = fRawCCG./fNorm;
fRawCCG(isnan(fRawCCG)) = 0;
fShiftCCG = zeros(1, nCol1+nCol2-1, 'single');
fShiftCi = zeros(2, nCol1+nCol2-1, 'single');
switch ShifPredictor
    case 0  % None
    case 1  % JPSTH        
%         [fRawCCG fShiftCCG] = JPSTH(iTrain1, iTrain2);
        fRawCCG = mean(fRawCCG.*fNorm, 1);
        fShiftCCG = conv2(mean(iTrain1, 1), mean(iTrain2lr, 1));
        switch iNormMethod
            case 0
                iMinCol = min(nCol1, nCol2);
                fNorm = [1:iMinCol ones([1 nCol1+nCol2-2*iMinCol])*iMinCol (iMinCol-1):-1:1];
            case 1
                fNorm = mean(iCount1, 1);
            case 2
                fNorm = mean(iCount2, 1);
            case 3
                fNorm = sqrt(mean(iCount1, 1).*mean(iCount2, 1));
            case 4
                fNorm = mean(iCount1, 1)+mean(iCount2, 1)-fRawCCG;
        end
        fRawCCG = fRawCCG./fNorm;
        fRawCCG(isnan(fRawCCG)) = 0;
        fShiftCCG = fShiftCCG./fNorm;
        fShiftCCG(isnan(fShiftCCG)) = 0;
    case 2  % Theoretical shift predictor
        fTempCCG = zeros(nRow-1,nCol1+nCol2-1,'single');
        fShiftCCG = zeros(nRow-1,nCol1+nCol2-1,'single');
        for i = 2:nRow
            iTrain2lr = circshift(iTrain2lr,-1);
            for j = 1:nRow
                fTempCCG(j,:) = conv2(iTrain1(j,:), iTrain2lr(j,:));
            end
            iCount2 = circshift(iCount2,-1);
            switch iNormMethod
                case 0
                    iMinCol = min(nCol1, nCol2);
                    fNorm = [1:iMinCol ones([1 nCol1+nCol2-2*iMinCol])*iMinCol (iMinCol-1):-1:1];
                    fTempCCG = fTempCCG./fNorm;
                case 1
                    fTempCCG = fTempCCG./iCount1;
                case 2
                    fTempCCG = fTempCCG./iCount2;
                case 3
                    fTempCCG = fTempCCG./sqrt(iCount1.*iCount2);
                case 4
                    fTempCCG = fTempCCG./(iCount1+iCount2-fTempCCG);
            end
            fTempCCG(isnan(fTempCCG)) = 0;
            fShiftCCG(i,:) = sum(fTempCCG,1);
        end
        fShiftCCG = sum(fShiftCCG,1)/nRow/(nRow-1);
        
%         fFullCCG = zeros(nRow*(nRow-1),nCol1+nCol2-1,'single');
%         RowIdx = nchoosek(1:nRow, 2);
%         RowIdx = [RowIdx; fliplr(RowIdx)];
%         for i = 1:size(fFullCCG, 1)
%             fFullCCG(i,:) = conv2(iTrain1(RowIdx(i,1),:), iTrain2lr(RowIdx(i,2),:));
%         end
%         switch iNormMethod
%             case 0
%                 iMinCol = min(nCol1, nCol2);
%                 fNorm = repmat([1:iMinCol ones([1 nCol1+nCol2-2*iMinCol])*iMinCol (iMinCol-1):-1:1], [nRow*(nRow-1) 1]);
%             case 1
%                 fNorm = iCount1(RowIdx(:,1),:);
%             case 2
%                 fNorm = iCount2(RowIdx(:,2),:);
%             case 3
%                 fNorm = sqrt(iCount1(RowIdx(:,1),:).*iCount2(RowIdx(:,2),:));
%             case 4
%                 fNorm = iCount1(RowIdx(:,1),:)+iCount2(RowIdx(:,2),:)-fFullCCG;
%         end
%         fFullCCG = fFullCCG./fNorm;
%         fFullCCG(isnan(fFullCCG)) = 0;
%         fShiftCCG = mean(fFullCCG, 1);

%         if nShift>0 && nShift<nRow*(nRow-1)
%             iSample = randi(nRow-1, [nRow nShift])+(0:(nRow-1))'*(nRow-1)*ones(1, nShift);
%             fSampleCCG = zeros(nShift, nCol1+nCol2-1, 'single');
%             for i = 1:nShift
%                 fSampleCCG(i,:) = mean(fFullCCG(iSample(:,i),:), 1);
%             end
%             fSampleCCG = sort(fSampleCCG, 1);
%             fShiftCi(1, :) = fSampleCCG(max(1, floor(nShift*0.025)), :);
%             fShiftCi(2, :) = fSampleCCG(min(nShift, ceil(nShift*0.975)), :);
%         elseif nShift>0 && nShift>=nRow*(nRow-1)
%             fSampleCCG = zeros(nRow*(nRow-1), nCol1+nCol2-1, 'single');
%             for i = 1:nRow
%                 fSampleCCG = fSampleCCG + fFullCCG(randperm(nRow*(nRow-1)), :);
%             end
%             fSampleCCG = fSampleCCG/nRow;
%             fSampleCCG = sort(fSampleCCG, 1);
%             fShiftCi(1, :) = fSampleCCG(max(1, floor(nRow*(nRow-1)*0.025)), :);
%             fShiftCi(2, :) = fSampleCCG(min(nRow*(nRow-1), ceil(nRow*(nRow-1)*0.975)), :);
%         end
end