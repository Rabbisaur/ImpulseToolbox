function isi=CIF2ISI(lambda,bins)

%*************************************************************************
% Calculate the expected ISI pdf integrated between the bins
%*************************************************************************
% Input variables:
% lambda: Conditional intensity function (CIF)
% bins:   the borders of the bins for the CIFs of every neuron
%*************************************************************************
% Output variables:
% isi(k): ISIcdf(k)-ISIcdf(k-1)
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


S=size(lambda,1); % #states
Q=size(lambda,2); % #phases

bw=repmat(diff(bins),S,1); %widths of the bins

%isi(k)=exp( - \sum_{i=1}^{k-1} lambda(i)*bw(i)  ) * ( 1-exp(-lambda(k)*bw(k)) )
isi=exp( -cumsum(lambda.*bw, 2) );
isi=circshift(isi,[0 1]);
isi(:,1)=1;
isi=isi .* (-exp( -lambda.*bw ) + 1);

%if lambda=nan => isi=0
isi(find(isnan(isi)))=0;

%normalize isi
isi=isi./repmat(sum(isi,2),1,Q);