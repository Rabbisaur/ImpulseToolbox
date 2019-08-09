function s=dchmm(infile,outfile,nstates,readchannels,nbins,jackknife,varargin)

%*************************************************************************
% Call all the DCHMM functions and save the output to a file
%*************************************************************************
% Input variables:
% infile:       input file name (.mat)
% outfile:      output file name, no output written if not specified
% nstates:      number of hidden states of the DCHMM to fit
% readchannels: which recording channel(s) to read in, default: all
% nbins:        how many bins to use for the CIF
% jackknife:    time array in seconds to cut out, default: none
% varargin:     passed on to dchmmtrain
%*************************************************************************
% Output variables:
% s:            structure containing all the data
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

tic;

nnames={'RA','HVc_{[i]}','HVC_{[RA]}','HVc_{[X]}','RA_{[i]}',...
        '','','','','','','','','','Nif_{[HVC]}','Nif_{[i]}',...
        'Nif_{[pseudo]}','field L','','Lman_{[RA]}','Lman_{[i]}'};

s.S=nstates;

inf = load(infile);
[fpath,s.fromfile,fext,fversn] = fileparts(infile);
if isfield(s,'scanrate')
    s.scanrate=inf.scanrate;
else
    s.scanrate=20000;
end;

if isempty(readchannels)
    readchannels=1:length(inf.ch);
end;

s.N=length(readchannels); % #neurons

if numel(nbins)==1
    nbins=ones(s.N,1)*nbins;
end;

%read & accumulate spike trains
s.spks=cell(1,s.N);       % the spike trains
s.t=[];                   % united spike trains
s.t1st=-Inf;              % start and end of time window ...
s.tend=-Inf;              % ... (where all the phases are known)
for n=1:s.N
    fromch=readchannels(n);
    s.nid(n)=inf.ch(fromch);
    s.nname{n}=nnames(s.nid(n));

    %read spike train of neuron n
    if iscell(inf.spks)
        spktemp=inf.spks{fromch}/s.scanrate;
    else
        spktemp=inf.spks(fromch,find(inf.spks(fromch,:)))/s.scanrate;
    end;
    
    %remove data between jacknife(1) and (2)
    if ~isempty(jackknife)
        spk1=spktemp(find(spktemp<jackknife(1)));
        spk2=spktemp(find(spktemp>jackknife(2)));
        if ~isempty(spk1) && ~isempty(spk2)
            spk2=spk2-(spk2(1)-spk1(end));
            spk2=spk2(2:end);
        end;
        spktemp=[ spk1 spk2 ];
    end;
   
    s.t1st=max(s.t1st,spktemp(2));
    s.tend=max(s.tend,spktemp(end));
    s.spks{n}=spktemp;
    s.t=[s.t spktemp];
end;
clear spktemp;

%sort & delete duplicate entries (where more than one neuron spiked)
s.t=unique(s.t);

%calculate delta
for n=1:s.N
    s.delta(n,:)=ismember(s.t,s.spks{n});
end;

%calculate phi
s.phi=zeros(s.N,length(s.t));
oldphi=zeros(s.N,1); %init all neurons with phi=0
for k=2:length(s.t)
    oldphi(find(s.delta(:,k-1)))=0;    %spike at k-1 => reset phase
    s.phi(:,k)=oldphi+s.t(k)-s.t(k-1); %phase=old phase + iei
    oldphi=s.phi(:,k);
end;
clear oldphi;

%strip data outside of borders
ix=find(s.t>=s.t1st & s.t<=s.tend);
s.t=s.t(ix);
s.delta=s.delta(:,ix);
s.phi=s.phi(:,ix);

%calculate interevent intervals
s.IEIs=diff(s.t); 

%set guess values for the HMM
s.r0=ones(s.S,1);
s.G0=(ones(s.S)-eye(s.S))./(s.S-1);

for n=1:s.N
    s.bins{n}=[0 logspace(-3, 1, nbins(n))];
    Q=length(s.bins{n})-1;
    %init every state of every neuron as a poisson process with unit rate
    s.lambda0{n}=ones(s.S,Q);
end;
%except for the last neuron, here we use gaussians
Q=length(s.bins{s.N})-1;
for i=1:s.S
    s.lambda0{s.N}(i,:)=1000*normpdf(1:Q,i/(s.S+1)*Q,0.05*Q);
end;

%train hmm
[s.r,s.lambda,s.G,s.loglik]=dchmmtrain(s.IEIs,s.delta,s.phi,s.r0,s.lambda0,s.G0,s.bins,varargin);

%get the posteriors
[s.post s.logpx]=dchmmdecode(s.IEIs,s.delta,s.phi,s.r,s.lambda,s.G,s.bins);

%get most likely state sequence
[s.stateseq, s.logpxs]=dchmmviterbi(s.IEIs,s.delta,s.phi,s.r,s.lambda,s.G,s.bins);

%histograms
for n=1:s.N
    s.ISIhist{n}=gethist(diff(s.spks{n}), s.bins{n});

    stateISIhist=zeros(s.S,length(s.ISIhist{n}));
    [hc,iphi]=histc(s.phi(n,:),s.bins{n});
    for k=1:length(s.IEIs)
        if s.delta(n,k+1)
            stateISIhist(s.stateseq(k),iphi(k+1)) = ...
            stateISIhist(s.stateseq(k),iphi(k+1)) + 1;
        end;
    end;
    s.stateISIhist{n}=stateISIhist;
end;

[s.stimes, s.sdurs, s.stimerel]=dchmmgetstatetimes(s.t,s.stateseq,s.S);

if length(outfile) > 0
    save(outfile,'-struct','s');
end;

toc;

function hc=gethist(data,bins)
hc=histc(data,bins);
hc=hc(1:end-1);
