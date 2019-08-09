function [stimes,sdurs,ptime]=dchmmgetstatetimes(events,sseq,S)

%*************************************************************************
% Calculate state start and end times + durations
%*************************************************************************
% Input variables:
% events: the event train data
% sseq:   state sequence
% S:      # states
%*************************************************************************
% Output variables:
% stimes{i}: an array of start and end positions of state i's occurences
% sdurs{i}:  the lengths of these occurences
% ptime(i):  relative time of the neuron staying in state i, sum(ptime)=1
%*************************************************************************
% Copyright (C) 2005 Márton Danóczy & Richard Hahnloser
% Version: 1.1
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% any later version.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, get it from gnu.org/copyleft/gpl.html
%*************************************************************************

X=length(sseq); % #events
stimes=cell(1,S);
state=sseq(1);
sstart=events(1);

for t=2:X
    if sseq(t) ~= state
        send=events(t);
        stimes{state} = [ stimes{state} [sstart; send] ];
        sstart = send;
        state = sseq(t);
    end;
end;
stimes{state}=[stimes{state} [sstart; events(X+1)]];

for i=1:S
    sdurs{i} = diff(stimes{i},1);
end;

%calculate sdurs and ptime with omitting the first and the last state,
%since their durations are unknown
s1=sseq(1); se=sseq(end);
sdurs{s1}=sdurs{s1}(2:end);
sdurs{se}=sdurs{se}(1:end-1);

%calculate ptime with T=(beginning of last state)-(end of 1st state)
T=stimes{se}(1,end)-stimes{s1}(2,1);
for i=1:S
    ptime(i) = sum(sdurs{i})/T;
end;
