function data = getMUA2(data,MUAfilters,samplerate,targetsamplerate,OneInterpolationLength,InterpolationIdx)

% transpose to time x trials
if size(data,1) < size(data,2)
    data = data';
end

% bandpass

data = filtfilt(MUAfilters.BMUAbp, MUAfilters.AMUAbp,data);

% full wave rectify
data = abs(data);

NumTrials = size(data,2);

% linear interpolation
if OneInterpolationLength > 0

    for thisTrial = 1:NumTrials
        y = squeeze(data(:,thisTrial));
        idx = (0:(OneInterpolationLength+1))';
        interIdx = bsxfun(@plus,InterpolationIdx{thisTrial}-1,idx);
        idx = interIdx < 1 | interIdx > size(data,1);
        interIdx(idx) = [];
        interIdx = round(interIdx);

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
        yy(end) = 0;

        data(:,thisTrial) = yy;
    end
end


% lowpass
data = filtfilt(MUAfilters.BMUAlow, MUAfilters.AMUAlow,data);

% downsample
downsamplerate = samplerate / targetsamplerate;
data = downsample(data,downsamplerate);