function [SIG,T] = Exd_all(EVENT)

%
%Called by Tdt2ml to retrieve neurophysiological data (only stream data) after
%filtering the trigger array.
%Returns a matrix (SIG) :
%       1st dimension contains samples for the length of one trial
%       2nd dimension channels
%       3rd dimension number of trials
%
%usage in batch files:
%define the following variables in EVENT
%Input : EVENT.Myevent = (string) event  %must be a stream event
%        EVENT.type = 'strms'   %must be a stream event
%        EVENT.CHAN = (string) channel numbers

%Changed so that it extracts all the data.

%Chris van der Togt, 11/11/2005
%updated 08/04/2008
%

SIG = [];

EvCode = EVENT.Myevent;

F = figure('Visible', 'off');
H = actxcontrol('TTANK.X', [20 20 60 60], F);
H.ConnectServer('local', 'me');

if 0 == H.OpenTank(EVENT.Mytank, 'R')
    errordlg([EVENT.Mytank 'does not exist!!!'])
    H.CloseTank;
    H.ReleaseServer;
    close(F)
    return
end

H.SelectBlock(EVENT.Myblock);
H.CreateEpocIndexing;

Sampf = EVENT.strms(4).sampf; %sample frequency for this event
Evlngth = EVENT.strms(4).size; %number of samples in each epoch
Evtime = Evlngth/Sampf; %timespan of one event epoch plus one for safety
ChaNm = EVENT.strms(4).channels; %channels in block
if isfield(EVENT, 'CHAN') && length(EVENT.CHAN) <= ChaNm
    Chans = EVENT.CHAN;  %SELECTED CHANNELS
else
    Chans = 1:ChaNm;
end

Span = 100;
EVNUM = round((ChaNm *(Span+(2*Evtime))*Sampf)/Evlngth); %more event epochs than needed
Recnum = H.ReadEventsV(EVNUM, EvCode, 0, 0, 0.0, 0.0, 'ALL');

ChnIdx = H.ParseEvInfoV(0, Recnum, 4);      %channel number corresponding to event epoch
Times = H.ParseEvInfoV(0, Recnum, 6);
Data = H.ParseEvV(0, Recnum);   %event epoch data

t = Times(ChnIdx == 1);
i1 = 1;
for n = 1:length(t)-1
    st = t(n);
    ed = t(n+1);
    
    T(i1:i1+63) = linspace(st,ed,64);
    i1 = i1+Evlngth;
end
T = [T,linspace(t(end),t(end)+Evtime,64)];


for n = 1:ChaNm
    
    buf = Data(:,ChnIdx==n);
    SIG(:,n)  = reshape(buf,size(buf,1)*size(buf,2),1);
    
end





H.CloseTank;
H.ReleaseServer;
close(F)
