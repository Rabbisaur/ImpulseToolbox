function [r,lambda,G,loglik]=dchmmtrain(x,delta,phi,r0,lambda0,G0,bins,varargin)

%*************************************************************************
% Baum-Welch Training of the Double Chain Hidden Markov Model
%*************************************************************************
% Input variables:
% x:      interevent intervals (IEI)
% delta:  which channels emit a spike at the end of each IEI
% phi:    phase of every neuron at the end of each IEI
% ...0:   guess values for the parameters
% bins:   the borders of the bins for the CIFs of every neuron
%*************************************************************************
% Output variables:
% r(i):   Transition rate constants of state i
% G(i,j): Conditional transition probability from state i->j
% lambda: Conditional intensity function (CIF)
% loglik: log likelihood for each iteration
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

%defaults
tol=1e-6;
maxiter=500;
verbose=false;

%check optional parameters
if nargin > 6
    if iscell(varargin{1})
        varargin=varargin{1};
    end;
    okargs = {'tolerance','maxiterations','verbose'};
    for j=1:2:numel(varargin)
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k)
            error('mhmmtrain:BadParameter',...
                'Unknown parameter name:  %s.',pname);
        elseif length(k)>1
            error('mhmmtrain:BadParameter',...
                'Ambiguous parameter name:  %s.',pname);
        else
            switch(k)
                case 1  % tolerance
                    tol=pval;
                case 2  % max iterations
                    maxiter=pval;
                case 3  % verbose
                    if strmatch(upper(pval), strvcat('YES', 'ON', 'Y', 'TRUE'), 'exact')
                        verbose=true;
                    end;
            end;
        end;
    end;
end;

warning off MATLAB:divideByZero;

S=size(r0,1);    % #states
N=size(delta,1); % #neurons
K=length(x);     % #events
for n=1:N
    Q(n)=size(lambda0{n},2); % #bins
end;

%delta and phi relate to the spike train, so strip the first entries
%to match the IEIs
delta0=delta(:,1); delta=delta(:,2:end);
phi0=phi(:,1); phi=phi(:,2:end);

%get the fractions of phi passed in every IEI
F = dchmmquantizephase(delta,phi,delta0,phi0,bins);

%initialize estimates with the guess values
r=r0; G=G0; lambda=lambda0;

converged=false;
loglikold=-99999999;
loglik=NaN(1,maxiter);
r_est=NaN(S,maxiter);

for iter=1:maxiter;
  
    r_est(:,iter)=r;
    
    %get the hidden and observable state transition probs for every event
    [Pi, Ak, Bk] = dchmmprepare(x,delta,F,r,lambda,G,bins);

    %get scaled forward and backward probabilities
    [fp,bp,scale] = dchmmfwbw(Pi, Ak, Bk);
    gamma = fp .* bp;

    %reestimation of lambda
    for n=1:N %for all neurons
        sumgamma1=zeros(S,Q(n)); %sum(gamma) for all k where delta(k)=1
        sumfgamma=zeros(S,Q(n)); %sum(f*bw*gamma)
        bw=repmat(diff(bins{n}),S,1); %widths of the bins
       
        for k=1:K %for all events

            i1=F.i1(n,k); i2=F.i2(n,k);  
            switch i2-i1
                case 0
                    fnk=F.f2(n,k);
                case 1
                    fnk=[F.f1(n,k) F.f2(n,k)];
                otherwise
                    fnk=[F.f1(n,k) ones(1,i2-i1-1) F.f2(n,k)];
            end;

            %calculate the denominator of the reestimation formula ...
            sumfgamma(:,i1:i2) = sumfgamma(:,i1:i2) + gamma(:,k) * fnk;

            %...and now the numerator
            if delta(n,k)
                sumgamma1(:,i2)=sumgamma1(:,i2) + gamma(:,k);
            end;
            
        end;

        lambda{n} = sumgamma1 ./ (sumfgamma .* bw);
    end;
        
    %reestimate the new time constants r(i) and G(i,j)
    for i=1:S
        
        for j=1:S
            xi(i,j)  = sum( fp(i,1:K-1) .* ...
                       permute( Ak(i,j,2:K), [1 3 2] ) .* Bk(j,2:K) .* ... 
                       bp(j,2:K) ./ scale(2:K) );
            xix(i,j) = sum( fp(i,1:K-1) .* ...
                       permute( Ak(i,j,2:K), [1 3 2] ) .* Bk(j,2:K).* ...
                       bp(j,2:K) ./ scale(2:K) .* x(2:K) );
        end;
        
        J=[1:i-1 i+1:S];
        r(i) = sum(xi(i,J)) / ( xix(i,i) + 0.5 * sum(xix(i,J)) );
        for j=J
            G(i,j) = xi(i,j) / sum(xi(i,J));
        end;
    end;
    
    %get the log likelihood ln p(x|psi)
    loglik(iter)=sum(log(scale));
    
    %graphics
    if verbose
        figure(11);
        subplot(3,1,1);    
        plot(r_est');
        stitle=sprintf('Iter %d, ln p({\\bfx}|\\phi)=%.4f', iter, loglik(iter));
        for i=1:S
            stitle = [stitle sprintf(', r_%d=%.4fHz=(%.2fs)^{-1}',i,r(i),1./r(i))];
        end;
        title(stitle);
        for n=1:N
            subplot(3,N,N+n); newplot;
            hold on;
            plot(lambda0{n}',':');
            stairs(lambda{n}','LineWidth',2);
            hold off;
            
            subplot(3,N,2*N+n); newplot;
            hold on;
            plot(dchmmCIF2ISI(lambda0{n},bins{n})',':');
            stairs(dchmmCIF2ISI(lambda{n},bins{n})','LineWidth',2);
            hold off;
        end;    
        drawnow;
    end;
        
    llchange=1-loglik(iter)/loglikold;
    
    if llchange <= tol
        converged=true;
        break;
    end;
    loglikold=loglik(iter);
  
end; % end iterations

loglik=loglik(find(~isnan(loglik)));

if converged
    disp(sprintf('Algorithm converged after %d iterations.',iter))
else
    disp(sprintf('Algorithm failed to converge after %d iterations.',iter))
end;

warning on MATLAB:divideByZero;