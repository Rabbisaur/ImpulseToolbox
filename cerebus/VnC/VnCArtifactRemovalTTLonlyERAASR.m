function instanceinfo = VnCArtifactRemovalTTLonlyERAASR(instanceinfo,trialParam,trialAlignSTfieldname)
disp('Removing stimulation artifacts...')

matversion = '-v6';
if ispc
    slash = '\';
else
    slash = '/';
end

% parameters
% SampleRate = instanceinfo.samplerate;
% low pass filter 150Hz
% Fl=150;
% Fn = SampleRate/2;
% N = 2;
% [BLFP, ALFP] = butter(N,Fl/Fn,'low'); % LFP low pass


numElec = instanceinfo.numElec;
SampleRate = double(instanceinfo.samplerate);

% remove artifact
fprintf('\n')
disp('Removing artifacts')

% StimStartTime = round(-trialParam.startT * SampleRate);
% remove artifact electrode by electrode
for thisElec = 1:min([numElec,96])
    fprintf('.')
    EID = instanceinfo.ElecOrder(thisElec);
    disp(['EID=',num2str(EID)])
    % load in data
    cmd = ['trialAlignST = instanceinfo.trialInfo.',trialAlignSTfieldname,';'];
    eval(cmd);
    electrodeTrialData = VnCGetTrialData2(instanceinfo, 1, EID,trialAlignST,trialParam);
    trialData = electrodeTrialData.trialData;
    ArrayTrialData(:,:,thisElec) = trialData;
end

% ERAASR algorithm
% data = data_trials_by_time_by_channels;
data = ArrayTrialData;

%% Setup ERAASR Parameters

opts = ERAASR.Parameters();
opts.Fs = SampleRate; % samples per second
Fms = opts.Fs / 1000; % multiply to convert ms to samples

opts.thresholdHPCornerHz = 250;
opts.thresholdChannel = 16;
opts.thresholdValue = -1000;

opts.alignChannel = 1;
opts.alignUpsampleBy = 10;
opts.alignWindowPre = Fms * 0.5;
opts.alignWindowDuration = Fms * 12;

% 60 ms stim, align using 20 ms pre start to 110 post
opts.extractWindowPre = Fms * 20;
opts.extractWindowDuration = Fms * 200;
opts.cleanStartSamplesPreThreshold = Fms * 0.5;
        
opts.cleanHPCornerHz = 10; % light high pass filtering at the start of cleaning
opts.cleanHPOrder = 4; % high pass filter order 
opts.cleanUpsampleBy = 1; % upsample by this ratio during cleaning
opts.samplesPerPulse = round(Fms * 1000/300.6942); % 3 ms pulses
opts.nPulses = 50;

opts.nPC_channels = 12;
opts.nPC_trials = 2;
opts.nPC_pulses = 6;

opts.omit_bandwidth_channels = 1;
opts.omit_bandwidth_trials = 1;
opts.omit_bandwidth_pulses = 1;

opts.alignPulsesOverTrain = false; % do a secondary alignment within each train, in case you think there is pulse to pulse jitter. Works best with upsampling
opts.pcaOnlyOmitted = true; % if true, build PCs only from non-omitted channels/trials/pulses. if false, build PCs from all but set coefficients to zero post hoc

opts.cleanOverChannelsIndividualTrials = false;
opts.cleanOverPulsesIndividualChannels = false;
opts.cleanOverTrialsIndividualChannels = false;

opts.cleanPostStim = true; % clean the post stim window using a single PCR over channels

% opts.showFigures = false; % useful for debugging and seeing well how the cleaning works
% opts.plotTrials = 1; % which trials to plot in figures, can be vector
% opts.plotPulses = 1; % which pulses to plot in figures, can be vector
% opts.figurePath = pwd; % folder to save the figures
% opts.saveFigures = false; % whether to save the figures
% opts.saveFigureCommand = @(filepath) print('-dpng', '-r300', [filepath '.png']); % specify a custom command to save the figure

%% Do alignment and cleaning procedure

dataCleaned = ERAASR.cleanTrials(data, opts);

figure
subplot(2,2,1), hold on
plot(data(:,:,33)','blue')
plot(dataCleaned(:,:,33)','red')
subplot(2,2,2), hold on
plot(data(:,:,43)','blue')
plot(dataCleaned(:,:,43)','red')
subplot(2,2,3), hold on
plot(data(:,:,53)','blue')
plot(dataCleaned(:,:,53)','red')
subplot(2,2,4), hold on
plot(data(:,:,63)','blue')
plot(dataCleaned(:,:,63)','red')

% save(''