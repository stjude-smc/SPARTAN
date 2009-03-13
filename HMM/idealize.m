function [idealization,offsets,LL] = idealize(obs, model, start_p, trans_p)
% IDEALIZE   HMM discovery of underlying sequence of hidden states
%
%    [IDEALIZATION,OFFSETS,LL] = idealize( TRACES, MODEL, P0, A )
%    Produces a sequence of hidden states (IDEALIZATION) that best explains
%    the observed data in TRACES, given a guassian emmission MODEL, 
%    starting probabilities (P0) and transition probabilities (A).  TRACES
%    is a NxM matrix, where N is the number of traces and M is the length
%    of each trace.  MODEL is a Sx2 matrix of state mean and stdev values.
%    P0 is an Sx1 vector of starting probabilities for each state.  A is a
%    matrix of probability of transiting from state i to state j at each
%    time point (frame).
%
%    Idealization is a 1xN cell array of dwell sequences, each of which is
%    a 2xD matrix of state+dwell time pairs, where D is the number of
%    such dwells.  OFFSETS is a 1xN vector of indexes into the raw data
%    file (linearized TRACES).  LL is a 1xN vector of the log-likelihood of
%    each trace, given the sequence of states (Viterbi path) and the model.

[nTraces,nFrames] = size(obs);
nStates = size(model,1);

% Initialize emission probability distribution function
mu    = model(:,1);
sigma = model(:,2);

% C = 1./(sqrt(2*pi)*sigma); %gaussian leading coefficient
C = 0.3989422804014327./sigma; %gaussian leading coefficient
D = 2*(sigma.^2);  %exponential denominator

% Predict the sequence of hidden model states for each trace
idealization = cell(1,nTraces);
LL = zeros(1,nTraces);

for i=1:nTraces,

    % Trim trace to remove data after donor photobleaching
    trace = obs(i,:);
    traceLen = length(trace);
%     traceLen = find(trace~=0, 1,'last');
%     trace = trace(1:traceLen);

    % Precompute emmission probability matrix
    Bx = zeros(nStates,traceLen);
    for s=1:nStates,
%         Bx(s,:) = C(s) .* exp(-  ((trace-mu(s)).^2) ./ D(s)  );
        Bx(s,:) = normpdf( trace, mu(s), sigma(s) )/6;
    end
    
    % Find the most likely viterbi path in model space
    % (totalLL is the probability of the observation sequence given the model)
    [totalLL, vPath, vLL] = forward_viterbi(start_p, trans_p, Bx);
    
    % Convert sequence of state assignments to dwell-times at each state
    % and add this new idealization to the output
    idl = RLEncode(vPath);
    LL(i) = vLL;
    
    % HACK: Remove last dwel in dark state.
    % Not doing so confuses MIL...
    if idl(end,1)==1,
        idl = idl(1:end-1,:);
    end
    idealization{i} = idl;
end

% Add offsets to relate idealization back to raw data
offsets = (0:(nTraces-1))*nFrames;

end %function idealize


