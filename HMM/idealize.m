function [idealization,LL] = idealize(obs, model, start_p, trans_p)
% IDEALIZE   HMM discovery of underlying sequence of hidden states
%
%    [IDL,LL] = idealize( TRACES, MODEL, P0, A )
%    Produces an idealization (sequence of hidden states) that best explains
%    the observed data in TRACES, given a guassian emmission MODEL, 
%    starting probabilities (P0) and transition probabilities (A).
%    TRACES is a NxM matrix, where N is the number of traces and M is the 
%    length of each trace.  MODEL is a Sx2 matrix of state mean and stdev 
%    values. P0 is an Sx1 vector of starting probabilities for each state.
%    A is a matrix of probability of transiting from state i to state j 
%    at each time point (frame).
%
%    IDL is an NxM matrix of state assignments for each
%    datapoint in TRACES. LL is a 1xN vector of the log-likelihood 
%    ofeach trace, given the sequence of states (Viterbi path) and the model.

%   Copyright 2007-2016 Cornell University All Rights Reserved.



% Initialize emission probability distribution function
mu    = model(:,1);
sigma = model(:,2);
nStates = numel(mu);
start_p = reshape(start_p,[nStates,1]);

assert( ~all(trans_p(:)==0) && ~all(start_p==0) );

% Log transform version
% The small constant ensures values are less than unity.
start_p = log(start_p);
trans_p = log(trans_p);
prefactor = log( 0.001./(sqrt(2*pi).*sigma) );


% Generate Gaussian probability mass functions for each emission state.
% Using the PMF instead of the PDF is slower and less precise, but the
% numbers are properly normalized and can be interpreted as probabilities.
% nBins = 2000;
% minx = min( mu -5.*sigma );
% maxx = max( mu +5.*sigma );
% x = [-Inf linspace(minx,maxx,nBins-2) Inf];  %bin edges
% 
% pmf = zeros( numel(x)-1, nStates );
% for s=1:nStates,
%     cdf = normcdf(x,mu(s),sigma(s));
%     pmf(:,s) = cdf(2:end) - cdf(1:end-1);
% end


% Predict the sequence of hidden model states for each trace
nTraces = size(obs,1);
idealization = zeros( size(obs) );
LL = zeros(1,nTraces);

for i=1:nTraces,

    % Trim trace to remove data after donor photobleaching (where E=0).
    traceLen = find(obs(i,:)~=0, 1,'last');
    trace = obs(i,1:traceLen);
    
    if isempty(traceLen) || traceLen<2,
        continue;
    end
    
    % Precompute emmission probability matrix.
    Bx = zeros(nStates,traceLen);
    for s=1:nStates,
        % Using the PDF.
        %Bx(s,:) = 0.001*normpdf( trace, mu(s), sigma(s) );
        Bx(s,:) = prefactor(s) - 0.5*((trace-mu(s))./sigma(s)).^2;
        
        % Using the PMF version (see notes above):
        %[~,binIdx] = histc(trace,x);
        %Bx(s,:) = pmf(binIdx,s);
    end
    
    % Find the most likely viterbi path in model space
    [vPath, LL(i)] = forward_viterbi(start_p, trans_p, Bx);
    
    % Convert sequence of state assignments to dwell-times in each state
    % and add this new idealization to the output
    idealization(i,1:traceLen) = vPath;
end

end %function idealize


