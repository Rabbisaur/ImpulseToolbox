function [sseq, lnpseq] = dchmmviterbi(x,delta,phi,r,lambda,G,bins)

%*************************************************************************
% Viterbi Algorithm for the Double Chain Hidden Markov Model
%*************************************************************************
% Input variables:
% x:      interevent intervals (IEI)
% delta:  which channels emit a spike at the end of each IEI
% phi:    phase of every neuron at the end of each IEI
% r(i):   Transition rate constants of state i
% G(i,j): Conditional transition probability from state i->j
% lambda: Conditional intensity function (CIF)
% bins:   the borders of the bins for the CIFs of every neuron
%*************************************************************************
% Output variables:
% sseq:   optimal state sequence
% lnpseq: log likelihood of this sequence
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


S=length(r);   % #states
N=size(phi,1); % #neurons
K=length(x);   % #events

%delta and phi relate to spkt, so strip the first entry to match the IEIs
delta0=delta(:,1); delta=delta(:,2:end);
phi0=phi(:,1); phi=phi(:,2:end);

%get the fractions of phi passed in every IEI
F = dchmmquantizephase(delta,phi,delta0,phi0,bins);

%get the hidden and observable state transition probs for every event
[Pi, Ak, Bk] = dchmmprepare(x,delta,F,r,lambda,G,bins);
clear F;

%we use logarithms to avoid numerical issues
warning off MATLAB:log:logOfZero;
lnPi=log(Pi); lnAk=log(Ak); lnBk=log(Bk); 
warning on MATLAB:log:logOfZero;
clear Pi Ak Bk;

%start in a state according to pi
lnV = lnPi + lnBk(:,1);

newlnV=zeros(size(lnV));
psi=zeros(S,K);

%find optimal sequences for every state at every event
for k=2:K
    for j=1:S
        [newlnV(j),psi(j,k)] = max( lnV + lnAk(:,j,k) + lnBk(j,k) );
    end;
    lnV=newlnV;
end;

%backtracking
sseq=zeros(1,K,'uint16');
[lnpseq,sseq(K)]=max(lnV);
for k=K:-1:2
    sseq(k-1)=psi(sseq(k),k);
end;
