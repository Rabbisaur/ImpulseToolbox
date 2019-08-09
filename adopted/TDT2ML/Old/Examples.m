%Examples plotting psth's after retrieving data from HDF files

[EVENT, DATA] = H5read;

%to plot snip data
EVENT.Myevent = 'Snip'; 
EVENT.type = 'snips'; 
%EVENT.CHAN = (1:20);
snip = Exsniptimes(EVENT, DATA.Trials); 
Expsth(EVENT, snip)

%to plot Envl data
EVENT.Myevent = 'Envl'; 
EVENT.type = 'strms'; 
Expsth(EVENT, DATA.Envl)


%retrieving snip waveform data from a HDF file
import ttv.*
%first browse the hdf file visually
[FileName,PathName] = uigetfile('.h5');
f = fullfile(PathName, FileName);
S = hdf5info(f,'ReadAttributes', false);

T = tankview();
T.initmodel(strtok(FileName, '.'))
ParseObj(S.GroupHierarchy, T)
T.setmodel()
T.setVisible(true)
%Determine which datasets you would like to retrieve
%for example all snip trials for channel 1 would be
EVENT = eventparser(hdf5read(f, '/Block-4/EVENT'));
% values refer to         Block     Snip   channel  trials(All)
NumName = length(S.GroupHierarchy.Groups(3).Groups.Groups(9).Datasets);
loc = '/Block-4/Snip/c_2/t_';
Name = zeros(NumName,1);
for i = 1:NumName
    StrName = S.GroupHierarchy.Groups(3).Groups.Groups(9).Datasets(i).Name;
    Name(i) = str2num(char(regexp(StrName, '[^/t_]+$', 'match')));
end
[B,IX] = sort(Name);
snip = cell(NumName,1);
n = 1;
for i = IX'
pad = [loc num2str(Name(i))];
OBJ = hdf5read(f, pad);
snip{n} = eventparser(OBJ);
n = n+1;
end

%to plot the waveforms for a partcular trial
times = snip{1}.time;
data = snip{1}.data;
sampf = EVENT.snips.Snip.sampf;
x = ((0:29)/sampf)';
Len = length(times);
x = repmat(x, 1, Len);
for i = 1:Len
    x(:,i) = x(:,i) + times(i);
end
plot(x, data)

