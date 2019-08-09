function [datafit,Amps,freqs,Fval,sig]=fitlinesc(data,params,p,plt,f0)
% fits significant sine waves to data (continuous data).
%
%  Inputs:  
%       data        (data in [N,C] i.e. time x channels/trials or a single
%       vector) - required.
%       params      
%           tapers :  [TW K] where TW is the time-bandwidth product and K is the number of tapers to be used 
%	        Fs 	        (sampling frequency) -- optional. Defaults to 1.
%           fpass       (frequency band to be used in the calculation in the form [fmin fmax])
%	    p		    (P-value to calculate error bars for) - optional. 
%                           Defaults to 0.05/N where N is the number of samples which
%	                 corresponds to a false detect probability of approximately 0.05.
%       plt         (y/n for plot and no plot respectively)
%       f0          frequencies at which you want to remove the
%                   lines - if unspecified the program
%                   will compute the significant lines
%
%
%  Outputs: 
%       datafit        (linear superposition of fitted sine waves)
%       Amps           (amplitudes at significant frequencies)
%       freqs          (significant frequencies)
%       Fval           (Fstatistic at all frequencies)
%       sig            (significance level for F distribution p value of p)
% adapted from Chronux on 5/24/2011

[N,C]=size(data);
if nargin < 2 || isempty(params); params=[]; end;
taper=params.tapers;
Fs=params.Fs;
if nargin < 3 || isempty(p);p=0.05/N;end;
if nargin < 4 || isempty(plt); plt='n'; end;
if nargin < 5; f0=[]; end;
[Fval,A,f,sig] = ftestc(data,params,p,plt);
if isempty(f0);
   fmax=findpeaks(Fval,sig);
   freqs=cell(1,C);
   Amps=cell(1,C);
   datafit=data;
   for ch=1:C;
       fsig=f(fmax(ch).loc);
       freqs{ch}=fsig;
       Amps{ch}=A(fmax(ch).loc,ch);
       Nf=length(fsig);
       datafit(:,ch)=exp(i*2*pi*(0:N-1)'*fsig/Fs)*A(fmax(ch).loc,ch)+exp(-i*2*pi*(0:N-1)'*fsig/Fs)*conj(A(fmax(ch).loc,ch));
   end;
else
   indx = zeros( length(f0) );
   for n=1:length(f0);
       [fsig,indx(n)]=min(abs(f-f0(n)));
   end;
   fsig=f(indx);
   for ch=1:C;
       freqs{ch}=fsig;
       Amps{ch}=A(indx,ch);
       Nf=length(fsig);
       datafit(:,ch)=exp(i*2*pi*(0:N-1)'*fsig/Fs)*A(indx,ch)+exp(-i*2*pi*(0:N-1)'*fsig/Fs)*conj(A(indx,ch));
   end;
end;

