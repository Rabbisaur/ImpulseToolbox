function [EVENT, Trilist, TLelnames] = EXinf2(EVENT)
%[EVENT, St] = EXinf2(EVENT)
%usage : 
%called by Tdt2ml to retrieve info about events and 
%the timing data of strobe events from a TDT Tank
%
%If used in a batch file you must initialize these values:
%input: EVENT.Mytank = 'the tank you want to read from';
%       EVENT.Myblock = 'the block you want to read from';
%
%output: EVENT ;  a structure containing a lot of info 
%        Trilist ; an array containing timing info about all the trials
%             1st column contains stimulus onset times
%             2nd contains trial onset times, when the monkey starts to fixate
%             3rd saccade onset times
%             4th target onset times
%             5th correct(1) or not correct(0)
%             6th error (1) or no error(0)
%             7th micro stim times
%             8th 15 bit word value (0 - 2^15) for conditional stimulus data
%
% Chris van der Togt, 13/01/2006

matfile = [EVENT.Mytank EVENT.Myblock]; %name of file used to save event structure and trial list data

if exist([matfile '.mat'], 'file') > 0
  load(matfile)
  return
end


%E.UNKNOWN = hex2dec('0');  %"Unknown"
E.STRON = hex2dec('101');  % Strobe ON "Strobe+"
%E.STROFF = hex2dec('102');  % Strobe OFF "Strobe-"
%E.SCALAR = hex2dec('201');  % Scalar "Scalar"
E.STREAM = hex2dec('8101');  % Stream "Stream"
E.SNIP = hex2dec('8201');  % Snip "Snip"
%E.MARK = hex2dec('8801');  % "Mark"
%E.HASDATA = hex2dec('8000');  % has associated waveform data "HasData"

%event info indexes
I.SIZE   = 1;
%I.TYPE   = 2;
%I.EVCODE = 3;
I.CHAN   = 4;
%I.SORT   = 5;
I.TIME   = 6;
I.SCVAL  = 7;
%I.FORMAT  = 8;
I.HZ     = 9;
I.ALL    = 0;

F = figure('Visible', 'off');
H = actxcontrol('TTANK.X', [20 20 60 60], F);
H.ConnectServer('local','me');
H.OpenTank(EVENT.Mytank, 'R');
H.SelectBlock(EVENT.Myblock);


H.CreateEpocIndexing;
EVS = H.GetEventCodes(E.STREAM); %gets the long codes of event types
STRMS = size(EVS,2);
if ~isnan(EVS)
    for i = 1:STRMS;
        EVENT.strms{i} = H.CodeToString(EVS(i));
    end
    
    for j = 1:size(EVENT.strms,2)
        Epoch = char(EVENT.strms{j});
        Recnum = H.ReadEventsV(1000, Epoch, 0, 0, 0, 0, 'ALL'); %read in number of events      
    %call ReadEventsV before ParseEvInfoV !!!! I don't expct more than a
    %1000 channels per event
        EVENT.(Epoch).sampf = H.ParseEvInfoV(0, 1, I.HZ); %9 = sample frequency
        EVENT.(Epoch).size = H.ParseEvInfoV(0, 1, I.SIZE); %1 = number of samples * bytes (4??) 
        EVENT.(Epoch).channels = max(H.ParseEvInfoV(0, Recnum, I.CHAN)); %4 = number of channels
    end
end
%get snip events    
EVS = H.GetEventCodes(E.SNIP);
SNIPS = size(EVS,2);
if ~isnan(EVS)
    for i = 1:SNIPS;
            EVENT.snips{i} = H.CodeToString(EVS(i));
    end
    for j = 1:size(EVENT.snips,2)
        Epoch = char(EVENT.snips{j});
        Recnum = H.ReadEventsV(100000, Epoch, 0, 0, 0, 0, 'ALL'); %read in number of events      
        EVENT.(Epoch).sampf = H.ParseEvInfoV(0, 1, I.HZ); %9 = sample frequency
        EVENT.(Epoch).size = H.ParseEvInfoV(0, 1, I.SIZE); %1 = number of samples * bytes (4??) 
        %EVENT.(Epoch).channels = max(H.ParseEvInfoV(0, Recnum, I.CHAN)); %4 = number of channels
        
        Timestamps = H.ParseEvInfoV(0, Recnum, I.TIME); %6 = the time stamp
        Channel =    H.ParseEvInfoV(0, Recnum, I.CHAN);
        Chnm = max(Channel);
        EVENT.(Epoch).channels = Chnm;
        
        while Recnum == 100000
            Recnum = H.ReadEventsV(100000, Epoch, 0, 0, 0, 0, 'NEW'); %read in number of events 
            Timestamps = [Timestamps H.ParseEvInfoV(0, Recnum, I.TIME)];
            Channel = [Channel H.ParseEvInfoV(0, Recnum, I.CHAN)];
        end
        Times = cell(Chnm,1);
        for k = 1:Chnm
            Times(k) = {Timestamps(Channel == k)};
        end
        EVENT.(Epoch).times = Times; 
   end
end

EVS = H.GetEventCodes(E.STRON); %get long codes of STROBE on types
STRONS = size(EVS,2);
if ~isnan(EVS)
    for i = 1:STRONS;
        EVENT.strons{i} = H.CodeToString(EVS(i));           
    end
%get timing and other usefull data on the STREAM and STROBE ON events
%Epoch = {};
    for j = 1:size(EVENT.strons,2)
        Epoch = char(EVENT.strons{j});
        Recnum = H.ReadEventsV(10000, Epoch, 0, 0, 0, 0, 'ALL'); %read in number of events       
        TINFO = H.ParseEvInfoV(0, Recnum, I.TIME); %6 = the time stamp 
        if (strcmp(Epoch, 'word') || strcmp(Epoch, 'Word'))
         TINFO(2,:) = H.ParseEvInfoV(0, Recnum, I.SCVAL); %we also want the scalar value of the word epoch event
        end
      
        EVENT.info{j} = TINFO;
    end 
end

H.CloseTank;
H.ReleaseServer;
close(F)



Tz = [];

Word_INF = [];
Corr_INF = [];
Rewd_INF = [];
Trl_INF = [];
Sacc_INF = [];
Err_INF = [];
Targ_INF = [];
Stm_INF = [];
Micr_INF = [];
%obtain an index for each event code, these correspond with the arrays in
%the info matrix
for i = 1:length(EVENT.strons)
    switch EVENT.strons{i}
                    
        case 'word'
           WORD = i; 
            
            [Word_INF, Idx] = sort(EVENT.info{WORD}(1,:).');
            Word_INF(:,2) = EVENT.info{WORD}(2,Idx).';
                   
        case 'corr' 
            CORR = i;
            Corr_INF = sort(EVENT.info{CORR}.');
            
        case 'rewd'
            REWD = i;
            Rewd_INF = sort(EVENT.info{REWD}.');
            
        case 'tril'
            TRIL = i; 
            Trl_INF = sort(EVENT.info{TRIL}.');
            
        case 'stim'
            STIM = i;
            Stm_INF = sort(EVENT.info{STIM}.');           
            
        case 'targ'
            TARG = i;
            Targ_INF = sort(EVENT.info{TARG}.'); 
            
        case 'erro'
            ERROR = i;
            Err_INF = sort(EVENT.info{ERROR}.');
            
        case 'sacc'
            SACC = i;
            Sacc_INF = sort(EVENT.info{SACC}.');
            
        case 'Tzer'
            TZER = i;
            Tz =  sort(EVENT.info{TZER}.');  
            
        case 'micr'
            MICR = i;
            Micr_INF = sort(EVENT.info{MICR}.');
    end
end
            
            
          
%stim or catch
%these two time events are exclusive for each trial; 
%insert an extra column and add a token to characterize a stim or a catch
%then put them in one matrix and sort in time.

if isempty(Stm_INF), disp('Error no Stim')
end

if isempty(Corr_INF), disp('Warning no correct bit data')
end

if isempty(Sacc_INF), disp('Warning no sacc bit data')
end

if isempty(Word_INF), disp('Warning no word data')
end

if isempty(Err_INF), disp('Warning no error bit data')
end

if isempty(Rewd_INF), disp('Warning no reward bit data')
end

if isempty(Targ_INF), disp('Warning no targ bit data')
end

if isempty(Micr_INF), disp('Warning no micr bit data')
end

if isempty(Trl_INF), disp('Error no tril bit data')
else 

%add a dimension for the correct, saccade and word bit and whether this stimfix can be used (pass), 
%add the corresponding parameters at the right times to this 4 dimensional array.
TrlNm = length(Trl_INF);
Trilist = zeros(TrlNm,9); %stim onset, trial onset, saccade onset, target onset
                           %correct, error, micr, word
TLelnames = {'stim_onset', 'trial_onset', 'saccade_onset', 'target_onset', 'correct', 'reward', 'error', 'micro_stim_time', 'word'};

Trilist(:,2) = Trl_INF;
for i = 1:TrlNm   %go from trial to trial only for the selected indices
            if ~isempty(Stm_INF)
                if i < TrlNm
                    Ixj = find(Stm_INF > Trl_INF(i) & Stm_INF < Trl_INF(i+1));
                else
                    Ixj = find(Stm_INF > Trl_INF(i), 1, 'first');
                end
                if ~isempty(Ixj)
                    Trilist(i,1) = Stm_INF(Ixj(1));
                else
                    Trilist(i,1) = nan;
                end
            end
             if ~isempty(Sacc_INF)
                 if i < TrlNm
                    Ixm = find(Sacc_INF > Trl_INF(i) & Sacc_INF < Trl_INF(i+1));
                 else
                    Ixm = find(Sacc_INF > Trl_INF(i), 1, 'first');
                 end
                 if ~isempty(Ixm)
                     Trilist(i,3) = Sacc_INF(Ixm(1));  %saccade onset
                 else Trilist(i,3) = nan;
                 end    
             end
             
             if ~isempty(Targ_INF)
                 if i < TrlNm
                    Ixm = find(Targ_INF > Trl_INF(i) & Targ_INF < Trl_INF(i+1));
                 else
                    Ixm = find(Targ_INF > Trl_INF(i), 1, 'first');
                 end
                 if ~isempty(Ixm)
                     Trilist(i,4) = Targ_INF(Ixm(1));  %target onset
                 else Trilist(i,4) = nan;
                 end    
             end
             
            if ~isempty(Corr_INF)
                if i < TrlNm
                    Ixk = find(Corr_INF > Trl_INF(i) & Corr_INF < Trl_INF(i+1));
                else
                    Ixk = find(Corr_INF > Trl_INF(i), 1, 'first');
                end
                 if ~isempty(Ixk)
                     Trilist(i,5) = 1;          %correct trial
                 end                     
            end
        
             if ~isempty(Rewd_INF)
                if i < TrlNm
                    Ixk = find(Rewd_INF > Trl_INF(i) & Rewd_INF < Trl_INF(i+1));
                else
                    Ixk = find(Rewd_INF > Trl_INF(i));
                end
                 if (~isempty(Ixk) && length(Ixk) >= 2)
                     Trilist(i,6) = 1;         %manual reward in trial
                 end            
             end 
            
            if ~isempty(Err_INF)
                if i < TrlNm
                    Ixk = find(Err_INF > Trl_INF(i) & Err_INF < Trl_INF(i+1));
                else
                    Ixk = find(Err_INF > Trl_INF(i), 1, 'first');
                end
                 if ~isempty(Ixk)
                     Trilist(i,7) = 1;         %error in trial
                 end            
            end 
            
            if ~isempty(Micr_INF)
                if i < TrlNm
                    Ixk = find(Micr_INF > Trl_INF(i) & Micr_INF < Trl_INF(i+1));
                else
                    Ixk = find(Micr_INF > Trl_INF(i), 1, 'first');
                end
                 if ~isempty(Ixk)
                     Trilist(i,8) = Micr_INF(Ixk(1));  %micr stim onset
                 else Trilist(i,8) = nan;
                 end            
            end
            
            if ~isempty(Word_INF)
                 if i == 1
                    Ixk = find(Word_INF(:,1) < Trl_INF(i), 1, 'last');
                 else
                    Ixk = find(Word_INF(:,1) > Trl_INF(i-1) & Word_INF(:,1) < Trl_INF(i), 1, 'last');
                 end
                 if ~isempty(Ixk)
                     Trilist(i,9) = Word_INF(Ixk(1),2);  %conditional information
                 else Trilist(i,9) = nan;
                 end       
            end        
    
end


IxE = find(isnan(Trilist(:,9)) == 1);
if ~isempty(IxE)
    % IxW = find(isnan(Trilist(:,8)) == 0);  %select only those trials with a valid word info
    % Trilist = Trilist(IxW,:);
     errordlg(['Trials without wordinfo!!!! : ' num2str(IxE.')] )
end
%unpak the Word 
%apos = bitand(StFx(:,5), 31);
%astm = bitand(fix(StFx(:,5).*2^-5), 31);
%asel = bitand(fix(StFx(:,5).*2^-10), 31);

%StFx = [StFx(:,1:4), apos, astm, asel];

end

save(matfile, 'EVENT', 'Trilist', 'TLelnames')