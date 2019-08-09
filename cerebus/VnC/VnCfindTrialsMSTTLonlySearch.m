function fileinfo = VnCfindTrialsMSTTLonlySearch(fileinfo,TTLchannel)

% define

disp('Searching for trials...')

elecfp = fopen(fileinfo.electrodeCachePath{TTLchannel},'r');
tmpdata = fread(elecfp,inf,'int16=>double');

Threshold = max(tmpdata)/2;
idx = tmpdata>Threshold;
idx = diff(idx) == 1;
idx = find(idx);
% 
% figure
% hold on
% plot(tmpdata)
% plot(idx,tmpdata(idx),'ro')

numTrials = numel(idx);
disp(['Got ', num2str(numTrials), ' trials.'])

fileinfo.trialInfo.trialAlignST = idx;