function tau = dchmmtimerescale(x,sseq,delta,phi,lambda,bins)

%*************************************************************************
% Apply the Time-Rescaling Theorem to our DCHMM
%*************************************************************************
% See: E. Brown et al., The Time-Rescaling Theorem and Its Application to
% Neural Spike Train Data Analysis, in Neural Computation, Vol 14, Num 2
% http://neco.mitpress.org/cgi/content/abstract/14/2/325
%*************************************************************************
% Input variables:
% x:      interevent intervals (IEI)
% sseq:   state sequence
% delta:  which channels emit a spike at the end of each IEI
% phi:    phase of every neuron at the end of each IEI
% lambda: Conditional intensity function (CIF)
% bins:   the borders of the bins for the CIFs of every neuron
%*************************************************************************
% Output variables:
% tau:    IEIs of the time rescaled process, which should be exponentially
%         distributed with parameter 1 to verify the model
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


N=size(phi,1); % #neurons
K=length(x);   % #events

%delta and phi relate to spkt, so strip the first entry to match the IEIs
delta0=delta(:,1); delta=delta(:,2:end);
phi0=phi(:,1); phi=phi(:,2:end);

%get the fractions of phi passed in every IEI
F = dchmmquantizephase(delta,phi,delta0,phi0,bins);

for n=1:N
    temp=zeros(1,K);
    lambdan=lambda{n};
    bw=diff(bins{n});

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
        temp(k)=sum( lambdan(sseq(k),i1:i2) .* bw(i1:i2) .* fnk );
    end;
    
    temp=cumsum(temp);
    temp=temp(find(delta(n,:)));
    tau{n}=diff(temp);

end;