function H5save(EVENT, Trilist)
%H5save(EVENT, Trilist)
%USAGE: hdf5 support for saving data from TDT Tanks
%IN : EVENT data structure, (see Exinfo3 )
%     Trilist array, first column should contain stimulus onset times
%     In the EVENT structure fields Triallngth and Start should be set.
%
%Chris van der Togt, 30-06-2006
%28-11-2006 : changed to save only snip channels containing values

%have all data members been set
if( ~isfield(EVENT, 'Triallngth') || ~isfield(EVENT, 'Start'))
    disp('Sorry first set Trial length and trial start')
    return
end


EVNT = hdf5.h5compound;

h5str = hdf5.h5string(EVENT.Mytank);
addMember(EVNT, 'Mytank')
setMember(EVNT, 'Mytank', h5str)
h5str = hdf5.h5string(EVENT.Myblock);
addMember(EVNT, 'Myblock')
setMember(EVNT, 'Myblock', h5str)
strm = [];
snip = [];

if isfield(EVENT, 'strms')
    STRMS = hdf5.h5compound;
    strm = fieldnames(EVENT.strms);
    for i = 1:length(strm)    
        SM = hdf5.h5compound;
        addMember(SM, 'size')
        setMember(SM, 'size', EVENT.strms.(strm{i}).size)
        addMember(SM, 'sampf')
        setMember(SM, 'sampf', EVENT.strms.(strm{i}).sampf)
        addMember(SM, 'channels')
        setMember(SM, 'channels', EVENT.strms.(strm{i}).channels)
        addMember(SM, 'bytes')
        setMember(SM, 'bytes', EVENT.strms.(strm{i}).bytes)
    
        addMember(STRMS, strm{i})
        setMember(STRMS, strm{i}, SM)
    end
    addMember(EVNT, 'strms')
    setMember(EVNT, 'strms', STRMS)
end %adding streams

if isfield(EVENT, 'snips')
    SNIPS = hdf5.h5compound;
    snip = fieldnames(EVENT.snips);
    for i = 1:length(snip)
        SP = hdf5.h5compound;
        addMember(SP, 'size')
        setMember(SP, 'size', EVENT.snips.(snip{i}).size)
        addMember(SP, 'sampf')
        setMember(SP, 'sampf', EVENT.snips.(snip{i}).sampf)
        addMember(SP, 'bytes')
        setMember(SP, 'bytes', EVENT.snips.(snip{i}).bytes)
        addMember(SP, 'channels')
        h5arr = hdf5.h5array(EVENT.snips.(snip{i}).channels{:});        
        setMember(SP, 'channels', h5arr)
        addMember(SP, 'times')
        SPT = hdf5.h5compound;
        for j=1:length(EVENT.snips.(snip{i}).channels{:})
            h5arr = hdf5.h5array(EVENT.snips.(snip{i}).times{j});
            addMember(SPT, num2str(j))
            setMember(SPT, num2str(j), h5arr)
        end
        setMember(SP, 'times', SPT)
        
        addMember(SNIPS, snip{i})
        setMember(SNIPS, snip{i}, SP)
    end
    addMember(EVNT, 'snips')
    setMember(EVNT, 'snips', SNIPS)       
end %add snips

if isfield(EVENT, 'strons')
    STRONS = hdf5.h5compound;
    stron = fieldnames(EVENT.strons);
    for i = 1:length(stron)
        h5arr = hdf5.h5array(EVENT.strons.(stron{i}));
        addMember(STRONS, stron{i})
        setMember(STRONS, stron{i}, h5arr)
    end
    addMember(EVNT, 'strons')
    setMember(EVNT, 'strons', STRONS)
end %add strobe on event info

if isfield(EVENT, 'Trials')
    TRIALS = hdf5.h5compound;
    trilnm = fieldnames(EVENT.Trials);
    for i = 1:length(trilnm)
        h5arr = hdf5.h5array(EVENT.Trials.(trilnm{i}));
        addMember(TRIALS, trilnm{i})
        setMember(TRIALS, trilnm{i}, h5arr)
    end
    addMember(EVNT, 'Trials')
    setMember(EVNT, 'Trials', TRIALS)
end

addMember(EVNT, 'Triallngth')
setMember(EVNT, 'Triallngth', EVENT.Triallngth)
addMember(EVNT, 'Start')
setMember(EVNT, 'Start', EVENT.Start)
if isfield(EVENT, 'CHAN')
    addMember(EVNT, 'CHAN')
    arr = hdf5.h5array(EVENT.CHAN);
    setName(arr, 'Selected channels');
    setMember(EVNT, 'CHAN', arr)
end
    
setName(EVNT, 'EVENT');

filename = [EVENT.Mytank '.h5'];
dset.Name = 'EVENT';
dset.Location = ['/' EVENT.Myblock];
if exist(filename, 'file')
    wmode = 'append';
    S = hdf5info(filename,'ReadAttributes', false);
    Topnm = size(S.GroupHierarchy.Groups,2);
    for i = 1:Topnm
        if strcmp(S.GroupHierarchy.Groups(i).Name, dset.Location)
            button = questdlg( 'This datset already exists. Choosing Overwrite completely deletes the HDF file', ...
                                'Overwrite or Cancel', 'Overwrite','Cancel','Cancel');
            if strcmp(button, 'Overwrite')
                wmode = 'overwrite';
            else
                return
            end
        end
    end
else
    wmode = 'overwrite';
end

hdf5write( filename, dset, EVNT, 'WriteMode', wmode)
disp('EVENT structure saved to HDF5');
%clear mex; %clears the file handel on the hdf file (bug in Matlab !!!!!!!!)
%retrieve all channels for envelop data
Chs = [strm; snip];
[Selection,ok] = listdlg('ListString', Chs,...
                         'ListSize', [120 length(strm)*15], ...
                         'Name', 'Save');
if(ok == 1)
    
    %only select trials with a stim onset
        Trials = Trilist(~isnan(Trilist(:,1)),:);
        TRLS = hdf5.h5array(Trials);
        setName(TRLS, 'Trials')    
        dset.Name = 'Trials';
        hdf5write( filename, dset, TRLS, 'WriteMode', 'append');
    
    
    for i = 1:length(Selection)
        if ismember(Chs(Selection(i)), strm)
            EVENT.Myevent = Chs{Selection(i)};
            SIG = Exd2(EVENT, Trials);
            SIGN = hdf5.h5array(SIG);
            setName(SIGN, EVENT.Myevent);
            dset.Name = EVENT.Myevent;
            hdf5write( filename, dset, SIGN, 'WriteMode', 'append');
            disp([ EVENT.Myevent ' saved to HDF5' ])
        
        elseif ismember(Chs(Selection(i)), snip)
            EVENT.Myevent = Chs{Selection(i)};
            SNP = Exsnip(EVENT, Trials);
            ni = size(SNP,1);           

            for i = 1:ni
                TMP = [SNP{i,:}];
                DATA = [TMP.data];   
                TIME = [TMP.time];

                Cmpd = hdf5.h5compound;
                addMember(Cmpd, 'time');
                setMember(Cmpd, 'time', hdf5.h5array(TIME) );
                addMember(Cmpd, 'data');
                setMember(Cmpd, 'data', hdf5.h5array(DATA) );
                        
                dset.Location = ['/' EVENT.Myblock '/' EVENT.Myevent ];
                dset.Name = num2str(i);
                hdf5write( filename, dset, Cmpd, 'WriteMode', 'append')
                

            end
            disp([ EVENT.Myevent ' saved to HDF5' ])
        end
    end
end
disp('finished saving to HDF5!')
