function [p,mu,stdev,n,x] = Undetected(w,threshes,criteria_func)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/201
%
% undetected - estimate fraction of events that did not reach threshold
%
% Usage:
%     [p,mu,stdev,n,x] = undetected(waveforms,threshes,criteria_func)
%
% Description:  
%    Estimates fraction of events that did not reach threshold by applying 
% the detection metric to each waveform and then fitting it with a Gaussian 
% that has a missing tail.
%
% The distribution of detection metric values is turned into a histogram.  
% A Gaussian is fit to the historgram to minimize the absolute error 
% between the Gaussian and the histogram for values above threshold.  
% The integral of this Gaussian that is below threshold is the estimate of 
% the fraction of missing events.
%
% Note that values are normalized so that the threshold is +/- 1.  The function
% attempts to preserve the sign of the original threshold, unless thresholds
% on different channels had different signs. In the case of multiple channels, 
% each channel is normalized so that the threshold has a magnitude of 1.  Then, 
% for each event, only the channel with the most extreme value of the detection 
% metric is used. 
%
% By default, this function assumes that a simple voltage crossing was used 
% for detection, but see "criteria_fun" below for alternatives. In the case 
% of a simple voltage threshold, note that the threshold is interpreted as 
% responding to crossings away from zero, i.e., negative thresholds 
% imply negative-going crossings and positive thresholds imply 
% positive-going crossings. 
%
%
% Input:
%   waveforms  - [Events X Samples X Channels] the waveforms of the cluster
%   threshes   - [1 X Channels] the threshold for each channel
%   criteria_func - Used to determine what the detection metric is on each
%                   waveform.  If this is the string "auto" or "manual" then
%                   it is assumed that a simple voltage threshold was used. 
%                   The detection criterion then is to divide each channel
%                   by its threhsold and use the maximum value.  Otherwise
%                   the criteria_func is assumed to be a function handle that
%                   takes in waveforms and threshes and returns the detection
%                   metric for each event [Events x 1].  The function will
%                   be called as
%                      criteria = criteria_func( waveforms, threshes)
%                   It is assumed that the values of criteria are normalized 
%                   to use a threshold value of + 1.
%
% Output:
%  p            - estimate of probability that a spike is missing because it didn't reach threshhold
%  mu           - mean estimated for gaussian fit
%  stdev        - standard deviation estimated for gaussian fit
%  n            - bin counts for histogram used to fit Gaussian
%  x            - bin centers for histogram used to fit Gaussian

% constant bin count
bins = 75;

% check for detection method
if isequal( criteria_func, 'auto') || isequal( criteria_func, 'manual' )
    % normalize all waveforms by threshold
    th(1,1,:) = threshes;
    w = w ./ repmat( th, [size(w,1) size(w,2) 1] );

    % get maximum value on each channel
    criteria = max( w(:,:), [], 2 );

else
   criteria = criteria_func( w, threshes);
end

% create the histogram values
global_max = max(criteria);
mylims = linspace( min(1,global_max),max(1,global_max),bins+1);
x = mylims +  (mylims(2) - mylims(1))/2;
n = histc( criteria, mylims );

% fit the histogram with a cutoff gaussian
m = mode_guesser(criteria, .05);    % use mode instead of mean, since tail might be cut off
[stdev,mu] = stdev_guesser(criteria, n, x, m); % fit the standard deviation as well

% Now make an estimate of how many spikes are missing, given the Gaussian and the cutoff
p = normcdf( 1,mu,stdev);

% attempt to keep values negative if all threshold values were negative
if all( threshes < 0 )
   mu = -mu;
   x = -x;
end

% fit the standard deviation to the histogram by looking for an accurate
% match over a range of possible values
function [stdev,m] = stdev_guesser( thresh_val, n, x, m)
% initial guess is juts the RMS of just the values below the mean
init = sqrt( mean( (m-thresh_val(thresh_val>=m)).^2  ) );
% try 20 values, within a factor of 2 of the initial guess
num = 20;
st_guesses = linspace( init/2, init*2, num );
m_guesses  = linspace( m-init,max(m+init,1),num);
for j = 1:length(m_guesses)
    for k = 1:length(st_guesses)
          b = normpdf(x,m_guesses(j),st_guesses(k));
          b = b *sum(n) / sum(b);
          error(j,k) = sum(abs(b(:)-n(:)));
    end        
end
% which one has the least error?
[val,pos] = min(error(:));
jpos = mod( pos, num ); if jpos == 0, jpos = num; end
kpos = ceil(pos/num);
stdev = st_guesses(kpos);
% refine mode estimate
m = m_guesses(jpos);

function m = mode_guesser(x,p)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
%
% mode_guesser - guess mode of the 
%
% Usage:
%    m = mode_guesser(x,p)
%
% Description:
%   Guesses mode by looking for location where the data is most tightly
% distributed.  This is accomplished by sorting the vector x and 
% looking for the p*100 percentile range of data with the least range.
%
% Input: 
%   x - [1 x M] vector of scalars
%
% Option input:
%   p - proportion of data to use in guessing the mode, defaults to 0.1
%
% Output:
%   m - guessed value of mode

%check for whether p is specified
if nargin < 2, p = .1; end
% determine how many samples is p proportion of the data
num_samples = length(x);
shift = round( num_samples * p );
% find the range of the most tightly distributed data
x = sort(x);
[val,m_spot] = min( x(shift+1:end) - x(1:end-shift) );
% use median of the tightest range as the guess of the mode
m = x( round(m_spot + (shift/2)) );    