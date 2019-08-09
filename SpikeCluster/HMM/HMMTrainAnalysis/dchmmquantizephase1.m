function F = dchmmquantizephase(delta,phi,delta0,phi0,bins)

%*************************************************************************
% Calculate the fractions of the phases passed in an IEI
%*************************************************************************
% Input variables:
% delta:  which channels emit a spike at the end of each IEI
% phi:    phase of every neuron at the end of each IEI
% delta0: value of delta before k=1
% phi0:   value of phi before k=1
% bins:   the borders of the bins for the CIFs of every neuron
%*************************************************************************
% Output variables:
% F.i1, F.i2: first and last nonzero index
% F.f1, F.f2: value of f^n_k at i1 and i2, respectively
% e.g., F.i1(n,k)=3 F.i2(n,k)=7 F.f1(n,k)=.2 F.f2(n,k)=.8
%       means f^n_k=[0, 0, .2, 1, 1, 1, .8, 0, 0, ...]
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

N=size(delta,1); % #neurons
K=size(delta,2); % #events

F.i1=zeros(N,K,'uint16');
F.i2=zeros(N,K,'uint16');
F.f1=zeros(N,K);
F.f2=zeros(N,K);

for n=1:N %for every neuron
    %get bin widths
    b=bins{n}; bw=diff(b); Q=length(bw);

    %get bin indices for the phases
    [dummy,pq0]=histc(phi0(n),b);
    [dummy,pq]=histc(phi(n,:),b);
    clear dummy;
    if any( pq > Q | pq <= 0 )
        error('dchmmquantize: Data outside of bins!')
    end;

    if delta0(n)
        oldphi=b(1); oldpq=1;
    else
        oldphi=phi0(n); oldpq=pq0;
    end;

    for k=1:K %for every event
        newphi=phi(n,k); newpq=pq(k);
        F.i1(n,k)=oldpq;
        F.i2(n,k)=newpq;            
        if newpq~=oldpq
            F.f1(n,k)=(b(oldpq+1)-oldphi)/bw(oldpq);
            F.f2(n,k)=(newphi-b(newpq))/bw(newpq);
        else
            F.f2(n,k)=(newphi-oldphi)/bw(newpq);
        end;
        if delta(n,k)
            oldphi=b(1); oldpq=1;
        else
            oldphi=newphi; oldpq=newpq;
        end;
    end;

end;
