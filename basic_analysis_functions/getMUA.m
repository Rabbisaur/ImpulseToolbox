function data = getMUA(data,MUAfilters,samplerate,targetsamplerate,OneInterpolationLength,InterpolationIdx)

% % MUA bandpass filter
% Fn = RawDataSampleRate/2;
% Fbp=[500,5000];
% N  = 4;    % filter order
% [BMUAbp, AMUAbp] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn],'bandpass'); % BandPass
%
% % MUA lowpass filter
% Fn = RawDataSampleRate/2;
% Fbp=[200];
% N  = 4;    % filter order
% [BMUAlow, AMUAlow] = butter(N, min(Fbp)/Fn,'low'); % BandPass

% transpose to time x trials
data = data';

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

% transpose back to trial x time
data = data';