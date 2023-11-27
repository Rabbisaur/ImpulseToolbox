function [salpadata,OneInterpolationLength, AllinterpolationIdxLength]= run_salpa_paoloData(elecData,realpulseIdx, pulseIdx, interpolationLength, samplerate)
debug = 0;

% interpolationLength = 0.5/1000; %s
interpolationIdxLength = round(interpolationLength / (1/samplerate));


% interIdx = bsxfun(@plus,AllStimIdx,idx);
salpadata = zeros(size(elecData));
NumTrials = size(elecData,1);
for thisTrial = 1:NumTrials

    idx = (0:interpolationIdxLength)';
    interIdx = bsxfun(@plus,pulseIdx{thisTrial},idx);
    idx = interIdx < 1 | interIdx > size(salpadata,2);
    interIdx(idx) = [];
    %                 x = 1:size(polynomialfitARdata,1);
    %                 xx = x;
    %                 x(interIdx) = [];
    y = elecData(thisTrial,:);
    
    if debug
        figure, hold on, plot(y), plot(interIdx,ones(size(interIdx)),'.'),plot(pulseIdx{thisTrial},ones(size(pulseIdx{thisTrial})),'x')
    end

    y(interIdx) = inf;

    if debug
        plot(y)
    end

    %                 yy = interp1(x,y,xx);
    salparesult = salpa(y,'tau',60); % apply SALPA
    
    if debug
        plot(salparesult)
    end
    salpadata(thisTrial,:) = salparesult;
end
clearvars elecData

% get a fair interpolation parameter for all electrodes and all trials
AllinterpolationIdxLength = zeros(NumTrials,1);
for thisTrial = 1:NumTrials
    y = squeeze(salpadata(thisTrial,:));
    y(1) = find(~isnan(y),1,'first');
    y(end) = find(~isnan(y),1,'last');
    nanidx = isnan(y);
    nanidx2 = find(diff(nanidx)==1);
    nanidx3 = find(diff(nanidx)==-1);
    nanlength = abs(nanidx3-nanidx2);
    [~,idx1] = min(abs(pulseIdx{thisTrial} - realpulseIdx{thisTrial}(1)));
    [~,idx2] = min(abs(pulseIdx{thisTrial} - realpulseIdx{thisTrial}(end)));
    interpolationIdxLength = max(nanlength(idx1:idx2));
    if isempty(interpolationIdxLength)
        interpolationIdxLength = -1;
    end
    AllinterpolationIdxLength(thisTrial) = interpolationIdxLength;
end

% calculate one single interpolation length and use it for all trial
OneInterpolationLength = max(AllinterpolationIdxLength);
interpolationIdxLength = OneInterpolationLength;
% linear interpolation
for thisTrial = 1:NumTrials
    y = squeeze(salpadata(thisTrial,:));
    idx = (0:(interpolationIdxLength+1))';
    interIdx = bsxfun(@plus,pulseIdx{thisTrial}-1,idx);
    idx = interIdx < 1 | interIdx > size(salpadata,2);
    interIdx(idx) = [];

    x = 1:numel(y);

    xx = x;
    x(interIdx) = [];

    idx = isnan(y);
    y(idx) = 0;
    y(interIdx) = [];
    yy = interp1(x,y,xx);

    idx = isnan(yy);
    yy(idx) = 0;

    % fix the last value
%     yy(end) = 0;
    if debug
        plot(yy)
    end
    salpadata(thisTrial,:) = yy;
end
