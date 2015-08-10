function [dwt,idealization,offsets,LL] = idealize(obs, model, start_p, trans_p)
% IDEALIZE   HMM discovery of underlying sequence of hidden states
%
%    [DWT,IDL,OFFSETS,LL] = idealize( TRACES, MODEL, P0, A )
%    Produces an idealization (sequence of hidden states) that best explains
%    the observed data in TRACES, given a guassian emmission MODEL, 
%    starting probabilities (P0) and transition probabilities (A).
%    TRACES is a NxM matrix, where N is the number of traces and M is the 
%    length of each trace.  MODEL is a Sx2 matrix of state mean and stdev 
%    values. P0 is an Sx1 vector of starting probabilities for each state.
%    A is a matrix of probability of transiting from state i to state j 
%    at each time point (frame).
%
%    DWT is a 1xN cell array of dwell sequences (see loadDWT.m), each of
%    which is a 2xD matrix of state+dwell time pairs, where D is the number
%    of such dwells.  IDL is an NxM matrix of state assignments for each
%    datapoint in TRACES. OFFSETS is a 1xN vector of indexes into the raw 
%    data file (linearized TRACES).  LL is a 1xN vector of the log-likelihood 
%    ofeach trace, given the sequence of states (Viterbi path) and the model.
%



[nTraces,nFrames] = size(obs);


% Initialize emission probability distribution function
mu    = model(:,1);
sigma = model(:,2);
nStates = numel(mu);

start_p = reshape(start_p,[nStates,1]);

assert( ~all(trans_p(:)==0) && ~all(start_p==0) );


% Generate Gaussian probability mass functions for each emission state.
nBins = 10000;
minx = min( mu -5.*sigma );
maxx = max( mu +5.*sigma );
x = [-Inf linspace(minx,maxx,nBins-2) Inf];  %bin edges

pmf = zeros( numel(x)-1, nStates );
for s=1:nStates,
    pmf(:,s) = normcdf(x(2:end),mu(s),sigma(s)) - normcdf(x(1:end-1),mu(s),sigma(s));
end


% Predict the sequence of hidden model states for each trace
dwt = cell(1,nTraces);
idealization = zeros( size(obs) );
LL = zeros(1,nTraces);

for i=1:nTraces,

    % Trim trace to remove data after donor photobleaching (where E=0).
    trace = obs(i,:);
    traceLen = find(trace~=0, 1,'last');
    trace = trace(1:traceLen);
    
    if isempty(traceLen) || traceLen<2,
        dwt{i} = zeros(0,nStates);
        continue;
    end
    
    % Precompute emmission probability matrix.
    Bx = zeros(nStates,traceLen);
    for s=1:nStates,
        % Discretize data, assigning each observation to a bin in the pmf.
        [~,binIdx] = histc(trace,x);
        
        % Calculate the probability of each observation from the pmf.
        Bx(s,:) = pmf(binIdx,s);
    end
    
    % Find the most likely viterbi path in model space
    [vPath, LL(i)] = forward_viterbi(start_p, trans_p, Bx);
    
    % Convert sequence of state assignments to dwell-times in each state
    % and add this new idealization to the output
    idealization(i,1:traceLen) = vPath;
    dwt{i} = RLEncode(vPath);
end

% Add offsets to relate idealization back to raw data
offsets = (0:(nTraces-1))*nFrames;

end %function idealize


