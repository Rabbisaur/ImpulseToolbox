function instanceinfo = VnCfindTrialsNoSerialcode(instanceinfo,AlignBit)

disp('Searching for trials...')
% numInstances = numel(instanceinfo);
% for thisInstance = 1:numInstances
%     disp(['  Working on instance ',num2str(thisInstance),'/',num2str(numInstances),'...'])
    instanceinfo.nev.digitimestamps = double(instanceinfo.nev.digitimestamps);
    instanceinfo.nev.digidata = double(instanceinfo.nev.digidata);
    

        idx = instanceinfo.nev.digidata == 2^AlignBit;
        TrialTimeStamps = double(instanceinfo.nev.digitimestamps(idx));
        
    
    numTrials = numel(TrialTimeStamps);
    
    disp(['Got ', num2str(numTrials), ' trials.'])

    instanceinfo.trialInfo.trialAlignST = TrialTimeStamps;
% end