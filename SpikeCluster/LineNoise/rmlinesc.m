function data=rmlinesc(data,params,p,plt,f0)
% removes significant sine waves from data (continuous data).
%
%  Inputs:  
% Note that units of Fs, fpass have to be consistent.
%       data        (data in [N,C] i.e. time x channels/trials or a single vector) - required.
%       params      structure containing parameters with fields
%           tapers :  [TW K] where TW is the time-bandwidth product and K is the number of tapers to be used 
%	        Fs 	        (sampling frequency) -- optional. Defaults to 1.
%           fpass       (frequency band to be used in the calculation in the form [fmin fmax])
%	    p		    (P-value for F-test) - optional. Defaults to 0.05/N
%	    where N is data length. This corresponds to a false detect
%	    probability of approximately 0.05
%
%       plt         (y/n for plot and no plot respectively)
%       f0          frequencies at which you want to remove the
%                   lines - if unspecified the program uses the f statistic
%                   to determine appropriate lines.
%
%  Outputs: 
%       data        (data with significant lines removed)
% adapted from Chronux on 5/24/2011

[N,C]=size(data);
if nargin < 2 || isempty(params); params=[]; end;
user_specified_pval=0;
if nargin < 3 || isempty(p);p=0.05/N; else; user_specified_pval=1; end;
if nargin < 4 || isempty(plt); plt='n'; end;
if nargin < 5; f0=[]; end;
if isempty(f0) && user_specified_pval==1; p=p/N; end;
[datafit,Amps,freqs,Fval,sig]=fitlinesc(data,params,p,'n',f0);
datan=data-datafit;
data=datan;