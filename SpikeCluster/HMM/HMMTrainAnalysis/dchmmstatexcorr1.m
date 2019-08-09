function sxcr=dchmmstatexcorr(f, t1st, tend, dt)

%*************************************************************************
% Calculate Pearson's correlation coefficient for two neurons
%*************************************************************************
% Input variables:
% f{n}:       state time array for one state, as provided by getstatetimes
%             The functions to be correlated are set to 1 inside the
%             investigated state described by f{n}, and 0 outside
% t1st, tend: start and end times
% dt:         time shift for each state sequence data
%*************************************************************************
% Output variables:
% sxcr:       correlation coefficient
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

if ~iscell(f) %autocorrelation
    f={f,f};
    nf=2;
else
    nf=numel(f);
    if nf==1 %autocorrelation
        f{2}=f{1};
        nf=2;
    end;
end;

if nargin==3
    dt=zeros(1,nf);
end;

if length(dt)==nf-1
    dt=[0 dt];
end;

T=tend-t1st;

mu=zeros(1,nf);
sigma=zeros(1,nf);
scomb=[]; index=[];

for i=1:nf
    
    %shift by dt
    s=f{i}+dt(i);    
    
    %delete parts before t1st
    ix=find(s(2,:)>t1st,1,'first');
    s=s(:,ix:end);
    if s(1,1)<t1st
        s(1,1)=t1st;
    end;
        
    %delete parts after tend
    ix=find(s(1,:)<tend,1,'last');
    s=s(:,1:ix);
    if s(2,end)>tend
        s(2,end)=tend;
    end;
        
    %we correlate sequences as defined by: 1 if in my state, 0 if outside
    tin=sum(diff(s,1)); %total state dwell time
    mu(i)=tin/T;        %mean of the sequence
    
    %the stdev of the sequence
    sigma(i)=sqrt( (1-mu(i))^2*tin + (0-mu(i))^2*(T-tin) );
    
    %flatten s
    s=s(:)';
    
    %cumulate values for sorting
    index=[index i*ones(1,length(s))]; %111111122222222 etc
    scomb=[scomb s];
end;

[scomb v]=sort(scomb); index=index(v);
clear v s tin ix;

numerator=0;
on=zeros(1,nf); %on: 1 if in my state, 0 if outside
on(index(1))=1; %first state is "on"
for i=2:length(scomb)
    numerator = numerator + (scomb(i)-scomb(i-1)) * prod(on-mu);
    on(index(i)) = ~on(index(i));
end;

numerator = numerator + (tend-scomb(end)) * prod(on-mu);

sxcr=numerator/prod(sigma);