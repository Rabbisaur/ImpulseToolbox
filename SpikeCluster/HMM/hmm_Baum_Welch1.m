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

% 'hmm_Baum_Welch' calculates the posterior probabilities of the
% hidden states given the observation DATA and the model parameters MU,
% SIGMA and P. Based on the posterior probabilities, the model parameters
% are updated to MU_NEW, SIGMA_NEW and P_NEW.

function [mu_new,sigma_new,p_new] = hmm_Baum_Welch(data,mu,sigma,p)

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
    error('The spiking probability P must be a probability between 0 and 1.')
end

% determine dimensions and initialize state transition matrix
length_of_observation=size(data,2);
number_of_states=size(mu,1);

transition_matrix=zeros(number_of_states);
transition_matrix(1,1)=1-p;
transition_matrix(2,1)=p;
transition_matrix(1,end)=1;
transition_matrix(3:end,2:end-1)=eye(number_of_states-2);

% Initialize the variables for the Baum-Welch algorithm. The posterior
% probabilities, the forward variable alpha and the backward variablie
% beta.

posterior_probabilities=zeros(number_of_states,length_of_observation);
alpha=zeros(number_of_states,length_of_observation);
alpha(1,1)=1;
beta=zeros(number_of_states,length_of_observation);
beta(1,end)=1;

% calculate the forward probabilities alpha
for count=2:length_of_observation
    % Gaussian emission probability
    alpha(:,count-1)=alpha(:,count-1).*((1/(sqrt(2*pi)*sigma))*...
        exp(-(mu-data(count-1)).*(mu-data(count-1))/sigma^2));
    % ring structure transitions
    alpha(:,count)=transition_matrix*alpha(:,count-1);
    % normalization
    alpha(:,count,:)=alpha(:,count)/sum(alpha(:,count));
end

% calculate the backward probabilities beta
for count=length_of_observation-1:-1:1
    % Gaussian emission probability
    beta(:,count+1)=beta(:,count+1).*((1/(sqrt(2*pi)*sigma))*...
        exp(-(mu-data(count+1)).*(mu-data(count+1))/sigma^2));
    % ring structure transitions
    beta(:,count)=transition_matrix'*beta(:,count+1);
    % normalization
    beta(:,count)=beta(:,count)/sum(beta(:,count));
end

% calculate the posterior probabilities
posterior_probabilities=alpha.*beta;
posterior_probabilities=posterior_probabilities./...
    repmat(sum(posterior_probabilities,1),number_of_states,1);

% calculate the update of the parameters
mu_new=(posterior_probabilities*data')./sum(posterior_probabilities,2);
p_new=sum(posterior_probabilities(2,:))/sum(posterior_probabilities(1,:));
sigma_new=std(data-mu'*posterior_probabilities);