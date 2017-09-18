function [fp, bp, scale] = dchmmfwbw(Pi,Ak,Bk)

%*************************************************************************
% Baum-Welch forward-backward procedure
%*************************************************************************
% Input variables:
% Pi(s):     probability to start in state s
% Ak(i,j,k): hidden state transition probability from state i at event k-1
%            to state j at event k
% Bk(i,k,n): observable state transition probability in neuron n from
%            the observable state right after event k-1 to the one at
%            event k, given hidden state S_k = i
%*************************************************************************
% Output variables:
% fp, bp:    scaled forward and backward probabilities
% scale:     the scaling factors
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


S=size(Ak,1); % #states
K=size(Bk,2); % #events

%calculate scaled forward probabilities
fp=zeros(S,K);    %forward probabilities
scale=zeros(1,K); %scaling factor

%start in a state according to Pi
fp(:,1) = Pi .* Bk(:,1);
scale(1) = sum(fp(:,1));
fp(:,1) = fp(:,1) ./ scale(1);

for k=2:K
    fp(:,k) = ( fp(:,k-1)' * Ak(:,:,k) )' .* Bk(:,k);
    scale(k) = sum(fp(:,k));
    fp(:,k) = fp(:,k) ./ scale(k);
end;

%initialize diagonal matrices BkDiag(:,:,k) (needed for backward probs)
BkDiag=zeros(S,S,K);
for i=1:S
    BkDiag(i,i,:)=Bk(i,:);
end;

%calculate scaled backward probabilities
bp=zeros(S,K);     %backward probabilities
bp(:,K)=ones(S,1); %init with prob.=1 for all states

for k=K-1:-1:1
    bp(:,k) = ( Ak(:,:,k+1) * BkDiag(:,:,k+1) * bp(:,k+1) ) ./ scale(k+1);
end;
