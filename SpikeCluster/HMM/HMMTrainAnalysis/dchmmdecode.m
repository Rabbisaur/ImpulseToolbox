function [gamma, loglik] = chmmdecode(x,delta,phi,r,lambda,G,bins)

%*************************************************************************
% Get the posterior probabilities and the log lokelihood of a DCHMM
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
% gamma(s,k):  posterior probability for state s at event k
% loglik: log likelihood p(x|psi)
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


%delta and phi relate to spkt, so strip the first entry to match the IEIs
delta0=delta(:,1); delta=delta(:,2:end);
phi0=phi(:,1); phi=phi(:,2:end);

%the fractions of phi passed in every IEI
F = dchmmquantizephase(delta,phi,delta0,phi0,bins);

%the hidden and observable state transition probs for every event
[Pi, Ak, Bk] = dchmmprepare(x,delta,F,r,lambda,G,bins);

%scaled forward and backward probabilities
[fp, bp, scale] = dchmmfwbw(Pi, Ak, Bk);

%log likelihood ln p(x|psi)
loglik = sum(log(scale));

%posteriors p( S_k=s | x,psi )
gamma = fp .* bp;
