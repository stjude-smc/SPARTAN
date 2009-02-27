function [LLtot,eA,eMU,eSIG,p0s,ps] = BWiterate( observations, A, mus, sigmas, p0s )
%
% Usage: [LLtot,eA,eMU,eSIG,p0s,ps] = BWiterate( Omat, A, mus, sigmas, p0s )
%
% Performs one iteration of adapted Baum-Welch algorithm to return log-likelihood
% and re-estimated model parameters [A,mus,sigmas,p0s] as well as overall fractional
% occupancy of each state ps (which is not a model term but we get for free)
%
% Note that data must be properly formatted in oberservations 'O' and repetion markers 'reps'
%
% DAB 2008.3.30
% DST 2008.6.29  Modified

[nTraces nFrames] = size(observations);
 
% if (length(reps) == 0) | (reps(1) == 0)
%   % Default is to use all reps
%   reps = 1:nTraces;
% end
% nTraces = nTraces;  %length(reps);

Nstates = length(mus);

LLtot = 0;
pstot = zeros(1,Nstates);
Osave = [];
Etot  = zeros(Nstates,Nstates);

for n = 1:nTraces
  trace = observations(n,:);

  % Trim trace to remove data after donor photobleaching
  lim = find(trace~=0, 1,'last');
  trace = trace(1:lim);
  
  NT = length(trace);
  Osave = [Osave trace(1:NT-1)];

  [E alphas LLc1] = BWtransition(trace,A,mus,sigmas,p0s);
  
  LLtot = LLtot + log(sum(alphas(NT,:)))+LLc1;%(NT);

  lambdas = zeros(NT-1,Nstates);
  Enorm = zeros(Nstates,Nstates);
  for t = 1:NT-1
    Enorm = Enorm + E{t};
	lambdas(t,:) = sum(E{t},2);
  end
  pstot = pstot + lambdas(1,:);
  
  % normalize
  for m = 1:Nstates
    Enorm(m,:) = Enorm(m,:)/sum(lambdas(:,m));
  end
  Etot = Etot+Enorm;
  
  %???
  if n == 1
	lamb_tot = lambdas;
  else
	lamb_tot(end+(1:NT-1),:) = lambdas;
  end
end

%% Re-estimate rate parameters (A matrix)
eA = Enorm;
% eA = Etot;


%% Re-estimate emmission parameters (stdev and mean)

eMU = zeros(1,Nstates);
eSIG = zeros(1,Nstates);

for n = 1:Nstates
  eMU(n) = Osave*lamb_tot(:,n) ./ sum(lamb_tot(:,n));
  eSIG(n) = sqrt((Osave-eMU(n)).^2*lamb_tot(:,n) ./ sum(lamb_tot(:,n)));
end



%% Calc....

p0s = pstot/nTraces;
LLtot = LLtot/nTraces;  % return LL per trial -- though not necessary...

% Calculate total occupancy of each state
ps = sum(lamb_tot)/sum(lamb_tot(:));

