%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2008 Joshua Herbst
% 
% This file is part of spike_sort_hmm.
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


% 'hmm_Viterbi' calculates the most probable state sequence for the hidden
% variable, given the observation DATA and the model parameterd MU, SIGMA
% and P.

function most_likely_state_sequence = hmm_Viterbi(data,mu,sigma,p)

% the observation DATA must be a column vector, and of mean 0
s=size(data);
if s(1)>s(2)
    data=data';
    s=size(data);
end
if s(1)>1
    error('The observation DATA must be 1-dimensional.')
end

% the means MU of the Gaussians for each state must be a row vector;
s=size(mu);
if s(2)>s(1)
    mu=mu';
    s=size(mu);
end
if s(2)>1
    error('The means MU must be 1-dimensional.')
end

% check SIGMA
if prod(size(sigma))~=1
    error('The standard deviation SIGMA must be a scalar.')
end

% check P and initialize the transition matrix
if prod(size(p))~=1
    error('The spiking probability P must be a scalar.')
end
if or(p<=0,p>=1)
    error('The spiking probability P must be a probability. 0 <= P <= 1')
end

% determine dimensions and initialize state transition matrix
length_of_observation=size(data,2);
number_of_states=size(mu,1);
transition_matrix=zeros(number_of_states);
transition_matrix(1,1)=1-p;
transition_matrix(2,1)=p;
transition_matrix(1,end)=1;
transition_matrix(3:end,2:end-1)=eye(number_of_states-2);

% Initialize the variables for the Viterbi algorithm: the best sequence
% probabilities, the state history and the most probable state sequence.

best_sequence_prob=zeros(number_of_states,length_of_observation);
best_sequence_prob=zeros(number_of_states,length_of_observation);
best_sequence_prob(1,1)=1;
state_history=zeros(number_of_states,length_of_observation);
most_likely_state_sequence=zeros(1,length_of_observation);
most_likely_state_sequence(end)=1;

% calculate the best sequence probability and the history
for count=2:length_of_observation
    % Gaussian emission probability
    best_sequence_prob(:,count-1)=best_sequence_prob(:,count-1)...
        .*((1/(sqrt(2*pi)*sigma))*...
        exp(-(mu-data(count-1)).*(mu-data(count-1))/sigma^2));
    % ring structure transitions
    for state=1:number_of_states
        [best_sequence_prob(state,count) state_history(state,count)]...
            =max(transition_matrix(state,:)'...
            .*best_sequence_prob(:,count-1),[],1);        
    end
    % normalization
    best_sequence_prob(:,count)=best_sequence_prob(:,count)...
        /sum(best_sequence_prob(:,count));
end

% backtrack to find the most probable state sequence
for count=length_of_observation-1:-1:1
    most_likely_state_sequence(count)=...
        state_history(most_likely_state_sequence(count+1),count+1);
end