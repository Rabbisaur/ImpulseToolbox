function [Pi, Ak, Bk] = dchmmprepare(x,delta,F,r,lambda,G,bins);

%*************************************************************************
% Calculate the hidden and observable state transition probabilities
%*************************************************************************
% Input variables:
% x:      interevent intervals (IEI)
% delta:  which channels emit a spike at the end of each IEI
% F:      the structure describing f^n_k as returned by dchmmquantizephase
% r(i):   transition rate constants of state i
% G(i,j): conditional transition probability from state i->j
% lambda: conditional intensity function (CIF)
% bins:   the borders of the bins for the CIFs of every neuron
%*************************************************************************
% Output variables:
% Pi(s):     probability to start in state s
% Ak(i,j,k): hidden state transition probability from state i at event k-1
%            to state j at event k
% Bk(i,k,n): observable state transition probability in neuron n from
%            the observable state right after event k-1 to the one at
%            event k, given hidden state S_k = i
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

S=length(r);     % #states
N=size(delta,1); % #neurons
K=length(x);     % #events

Ak=zeros(S,S,K);
for i=1:S
    Ak(i,i,:) = exp( -r(i)*x );
    for j=[1:i-1 i+1:S] %for all j \neq i
        Ak(i,j,:)=G(i,j) * (1-Ak(i,i,:));
    end;
end;

Bkn=zeros(S,K,N);
for n=1:N
    lambdan=lambda{n};
    bw=repmat(diff(bins{n}),S,1);

    for k=1:K
        
        i1=F.i1(n,k); i2=F.i2(n,k); 
        switch i2-i1
            case 0
                fnk=F.f2(n,k);
            case 1
                fnk=[F.f1(n,k) F.f2(n,k)];
            otherwise
                fnk=[F.f1(n,k) ones(1,i2-i1-1) F.f2(n,k)];
        end;
        fnk=repmat(fnk,S,1);
        %Since lambda is a firing rate of potentially >> 1 Hz, we have
        %to divide by an arbitrary frequency value in order to obtain
        %output probabilities in the intervsal [0;1]. Here we chose the
        %sampling rate of 20 kHz, i.e., the maximum possible firing rate 
        Bkn(:,k,n) = exp( ...
                       -sum(lambdan(:,i1:i2) .* fnk .* bw(:,i1:i2), 2) ...
                     ) .* (lambdan(:,i2)).^delta(n,k) / 20000;
    end;
end;

%The total probability is the joint prob of all neurons
Bk=prod(Bkn,3);

%T are the expected state durations
T=1./r; Pi=T/sum(T);