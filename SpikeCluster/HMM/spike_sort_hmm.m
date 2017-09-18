%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2008 Joshua Herbst
% 
% spike_sort_hmm is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or (at your option) any later version.
% 
% spike_sort_hmm is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
% PURPOSE.  See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with 
% this program. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This file (including hmm_Baum_Welch and hmm_Viterbi) implements a spike 
% sorting algorithm for a single neuron. The code is not optimized for 
% speed, but readability. To increase speed, the following should be 
% changed:
% - the transition matrix should be replaced by a hardcoded transition
% given by the ring structure of the hidden state space.
% - the normalization term in the Gaussian probabilities can be left away.
% - all the calculations can be done in the log-domain, i.e. normalizaion
% is no longer necessary.

close all;
clear;
load('observation.mat');
load('LFP.mat');
data=data(1:25000);
% data=data-mean(data);

% initialize parameters
p=1e-2;
sigma=std(data);
mu=zeros(60,1);
mu(10:20)=sin(0:pi/5:2*pi)*max(data)/5;

% run several iterations of the Baum-Welch algorithm to update model parameters
iterations=10;
for count=1:iterations
    [mu,sigma,p] = hmm_Baum_Welch(data,mu,sigma,p);
end 

% calculate the most likely hidden state sequence via the Viterbi algorithm
most_likely_state_sequence = hmm_Viterbi(data,mu,sigma,p);

% plot the original data and the data fit
most_likely_emission=mu(most_likely_state_sequence);

figure(1);clf;plot(data);hold on;plot(most_likely_emission,'r');