function [f,y] = ez_powermeasure(x,Fs,HW)

%Measures the power spectrum of a vector (x) with sampling frequency (Fs) via FFT and returns the frequencies (f)
%and power (y). HW determnines whether a Hann window is applied to the
%data, default = 0; 
%The power returned is corrected for the application of teh Hann window by
%pre-multiplyign teh amplitudes by 2 (as the Hann window has a mean of 0.5) before converting to power.

%M.W.Self 2014

if nargin < 3
    HW = 0;
end

if HW
    x = x.*hann(length(x));
end
%Nyquist limit is half the samplign frequency
Nyq = Fs/2;
%Length of time-series (in secs)
L = length(x)./Fs;
%Number of samples
N = length(x);

%Make samples and time-base
smps = 0:1:(N-1);
t = smps./Fs;

%FFT ANALYSIS%%%%%%%%%%%%%%%%%%
%First do fft, remember to scale by N.
Y = fft(x)/N;
%the resulting frequencies will be in cycles/timeseries
%starting with a frequency of 0 (DC) going upto frequencies that are equal to the length of the time-series-1.
%So the 2nd pos in the vector contains a frequency of 1 cycle/timeseries which has a
%frequency of 1/L;  This will continue upto Nyquist limit.
%So we can now make our frequency x-axis by diving the number of samples by the length of the time-series;
f = smps/L;

if HW
    %This takes power as the absoulte value
    %If using a Hann window the amplituide will be decreased by half (as the
    %Hann window has a mean of 0.5), power has a more complex correction
    %(3/8) or something)
    y = abs(Y);
    y = y.*4; %2 for the negative frequencies and 2 for the Hann window
    y = y.^2;
else
    %If not using a Hann window the Complex conjugate gives an identical result to taking the squatre of the absoulte value
    %and is quicker
    %If Y = a+ib
    %Then CC = a-ib
    %If you work out the brackets this equals real^2+imag^2 = power
    y = Y.*conj(Y);
end

return