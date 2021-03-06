function EVENT = Exinf3(EVENT)
%[EVENT, St] = Exinf3(EVENT)
%usage : 
%called by Tdt2ml to retrieve info about events and 
%the timing data of strobe events from a TDT Tank
%
%If used in a batch file you must initialize these values:
%input: EVENT.Mytank = 'the tank you want to read from';
%       EVENT.Myblock = 'the block you want to read from';
%
%output: EVENT ;  a structure containing a lot of info 

%
% Chris van der Togt, 29/05/2006
% 28-11-2006 :  changed to solve error saving to HDF, snip times only saved
%               for channels containing valid data

    matfile = [EVENT.Mytank EVENT.Myblock]; %name of file used to save event structure 
                                                            
    if exist([matfile '.mat'], 'file')
        load(matfile);
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
ALL = H.GetEventCodes(0);
AllCodes = cell(length(ALL),1);
for i = 1:length(ALL)
        AllCodes{i} = H.CodeToString(ALL(i));
end


EVS = H.GetEventCodes(E.STREAM); %gets the long codes of event types
STRMS = size(EVS,2);
strms = cell(STRMS,1);
if ~isnan(EVS)
    for i = 1:STRMS;
        strms{i} = H.CodeToString(EVS(i));
        IxC = find(strcmp(AllCodes, strms(i)));
        AllCodes(IxC) = [];
    end
    
    for j = 1:length(strms)
        Epoch = char(strms{j});
        Recnum = H.ReadEventsV(1000, Epoch, 0, 0, 0, 0, 'ALL'); %read in number of events      
    %call ReadEventsV before ParseEvInfoV !!!! I don't expct more than a
    %1000 channels per event

        T = H.ParseEvV(0, 1);
        EVENT.strms.(Epoch).size = size(T,1);    %number of samples in each event epoch
        EVENT.strms.(Epoch).sampf = H.ParseEvInfoV(0, 1, I.HZ); %9 = sample frequency  
        EVENT.strms.(Epoch).channels = max(H.ParseEvInfoV(0, Recnum, I.CHAN)); %4 = number of channels
        EVENT.strms.(Epoch).bytes = H.ParseEvInfoV(0, 1, I.SIZE); %1 = number of samples * bytes (4??) 

    end
end
%get snip events    

EVS = H.GetEventCodes(E.SNIP);
SNIPS = size(EVS,2);
snips = cell(SNIPS,1);
if ~isnan(EVS)
    for i = 1:SNIPS;
            snips{i} = H.CodeToString(EVS(i));
            %remove this item from allcodes
            IxC = find(strcmp(AllCodes, snips(i)));
            AllCodes(IxC) = [];
    end
    for j = 1:length(snips)
        Epoch = char(snips{j});
        Recnum = H.ReadEventsV(100000, Epoch, 0, 0, 0, 0, 'ALL'); %read in number of events
        T = H.ParseEvV(0, 1); 
        EVENT.snips.(Epoch).size = size(T,1); %number of samples per epoch event
        EVENT.snips.(Epoch).sampf = H.ParseEvInfoV(0, 1, I.HZ); %9 = sample frequency
        
        Timestamps = H.ParseEvInfoV(0, Recnum, I.TIME); %6 = the time stamp
        Channel =    H.ParseEvInfoV(0, Recnum, I.CHAN);
        EVENT.snips.(Epoch).bytes = H.ParseEvInfoV(0, 1, I.SIZE);
        Chnm = unique(Channel); %some channels may not contain data   
        [rw clm] = size(Chnm);
        if rw > clm        %make sure it is a row vector
            Chnm = Chnm.';
        end
        EVENT.snips.(Epoch).channels = {Chnm};

        
        while Recnum == 100000
            Recnum = H.ReadEventsV(100000, Epoch, 0, 0, 0, 0, 'NEW'); %read in number of events 
            Timestamps = [Timestamps H.ParseEvInfoV(0, Recnum, I.TIME)];
            Channel = [Channel H.ParseEvInfoV(0, Recnum, I.CHAN)];
        end
        Times = cell(1, length(Chnm));
        for k = 1:length(Chnm)
            Times(k) = {Timestamps(Channel == Chnm(k))};
        end
        EVENT.snips.(Epoch).times = Times; 
   end
end

%EVS = H.GetEventCodes(E.STRON); %get long codes of STROBE on types
%STRONS = size(EVS,2);
%strons = cell(STRONS,1);
%if ~isnan(EVS)
%    for i = 1:STRONS;
%        strons{i} = H.CodeToString(EVS(i));           
%    end
%get timing and other usefull data on the STREAM and STROBE ON events
%Epoch = {};
if ~isempty(AllCodes)
    strons = AllCodes;
    for j = 1:length(strons)
        Epoch = char(strons{j});
        Recnum = H.ReadEventsV(10000, Epoch, 0, 0, 0, 0, 'ALL'); %read in number of events       
        TINFO = H.ParseEvInfoV(0, Recnum, I.TIME); %6 = the time stamp 
        if (strcmp(Epoch, 'word') || strcmp(Epoch, 'Word'))
         TINFO(2,:) = H.ParseEvInfoV(0, Recnum, I.SCVAL); %we also want the scalar value of the word epoch event
        end
      
        EVENT.strons.(Epoch) = TINFO;
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

   Fnames = fieldnames(EVENT.strons);
   for i = 1:length(Fnames)
       switch Fnames{i}
           
            case 'word'     
            [Word_INF, Idx] = sort(EVENT.strons.word(1,:).');
            Word_INF(:,2) = EVENT.strons.word(2,Idx).';
                    
            case  'corr' 
            Corr_INF = sort(EVENT.strons.corr.');
        
            case 'rewd'
            Rewd_INF = sort(EVENT.strons.rewd.');
                    
            case  'tril'
            Trl_INF = sort(EVENT.strons.tril.');
          
            case  'stim'
            Stm_INF = sort(EVENT.strons.stim.');
                    
            case  'targ'
            Targ_INF = sort(EVENT.strons.targ.'); 
            
            case  'erro'
            Err_INF = sort(EVENT.strons.erro.');
            
            case  'Erro' %this should not be neccessary
            Err_INF = sort(EVENT.strons.Erro.');
          
            case 'sacc'
            Sacc_INF = sort(EVENT.strons.sacc.');
          
            case  'Tzer'
            Tz =  sort(EVENT.strons.Tzer.');  
          
            case 'micr'
            Micr_INF = sort(EVENT.strons.micr.');
        
       end
   end
   
    if isempty(Word_INF), disp('Warning no word events'),  end          
    if isempty(Corr_INF), disp('Warning no correct events'),  end
    if isempty(Rewd_INF), disp('Warning no reward events'),  end
    if isempty(Stm_INF), disp('Error no Stim on events'),  end
    if isempty(Targ_INF), disp('Warning no target on events'),  end
    if isempty(Err_INF), disp('Warning no error events'),  end
    if isempty(Sacc_INF), disp('Warning no saccade on events'),  end
    if isempty(Micr_INF), disp('Warning no micro stim events'),  end
    
if isempty(Trl_INF) && isempty(Tz), disp('Error no trial events (tril or Tz)')
else 

Names = {'stim_onset', 'trial_onset', 'saccade_onset', 'target_onset', 'correct', 'reward', 'error', 'micro_stim_time', 'word'};

%Tz is used only for compatibility with the first datatanks from Pieters
%projects
if ~isempty(Trl_INF)
    TrlNm = length(Trl_INF);
    Trilist = zeros(TrlNm,9);
    Trilist(:,2) = Trl_INF; %trial onsets(2)
else
    TrlNm = length(Tz);
    Trilist = zeros(TrlNm,9);
    Trilist(:,2) = Tz;
end


for i = 1:TrlNm   %go from trial to trial only for the selected indices
            if ~isempty(Stm_INF)
                %stimulus onsets (1)
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
                 %saccade onsets (3)
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
                 %target onsets (4)
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
                %corrects (5)
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
                 %rewards (6)
                if i < TrlNm
                    Ixk = find(Rewd_INF > Trl_INF(i) & Rewd_INF < Trl_INF(i+1));
                else
                    Ixk = find(Rewd_INF > Trl_INF(i));
                end
                 if ~isempty(Ixk)
                     Trilist(i,6) = length(Ixk);         %manual reward in trial
                 end                                     %if corr == 0 and rewd == 1
             end                                         %or corr == 1 and rewd > 1
            
            if ~isempty(Err_INF)
                %errors (7)
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
                %micr stim onsets (8)
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
                %words (9)
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

%save trial list in EVENT structure
EVENT.Trials.(Names{1}) = Trilist(:,1);
EVENT.Trials.(Names{2}) = Trilist(:,2);

if ~isempty(Sacc_INF)
    EVENT.Trials.(Names{3}) = Trilist(:,3);
end
if ~isempty(Targ_INF)
    EVENT.Trials.(Names{4}) = Trilist(:,4);
end
if ~isempty(Corr_INF)
    EVENT.Trials.(Names{5}) = Trilist(:,5);
end
if ~isempty(Rewd_INF)
    EVENT.Trials.(Names{6}) = Trilist(:,6);
end
if ~isempty(Err_INF)
    EVENT.Trials.(Names{7}) = Trilist(:,7);
end
if ~isempty(Micr_INF)
    EVENT.Trials.(Names{8}) = Trilist(:,8);
end
if ~isempty(Word_INF)
    EVENT.Trials.(Names{9}) = Trilist(:,9);
end

save(matfile, 'EVENT')
end

