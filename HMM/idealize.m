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

% Strictly speaking, this code is incorrect:
%  The probability of observing any particular FRET value is zero
%  because there are infinitely more FRET values nearly the same.
%  In principle, one should calculate the probability of a FRET
%  value occuring in a specific range. In other words, one
%  should bin the distribution and distribute FRET values into
%  bins with set probabilities. The effect on the algorithm
%  here would be to multiply each datapoint by the precision
%  of computer representation of FRET values (eg, 1e-9).
%  The effect of the result would be to add M*log10(precision)
%  to the final log-likelihood and to slightly reduce the
%  precision of calculating probabilities. Since this realistically
%  adds nothing to the algorithm but extra work, I have not
%  implemented this strategy.
%
% Also note that multiplying any other normalization factor
%  to the emission probabilities also has no effect on
%  the final path.
%


[nTraces,nFrames] = size(obs);


% Initialize emission probability distribution function
mu    = model(:,1);
sigma = model(:,2);
nStates = numel(mu);

start_p = reshape(start_p,[nStates,1]);

assert( all(trans_p(:)>=0) && all(start_p>=0) );


% Predict the sequence of hidden model states for each trace
dwt = cell(1,nTraces);
idealization = zeros( size(obs) );
LL = zeros(1,nTraces);

for i=1:nTraces,

    % Trim trace to remove data after donor photobleaching (where E=0).
    trace = obs(i,:);
    traceLen = find(trace~=0, 1,'last');
    trace = trace(1:traceLen);
    
    if isempty(traceLen) || traceLen<1,
        dwt{i} = zeros(0,nStates);
        continue;
    end
    
    % Precompute emmission probability matrix.
    % This ist reated as if using descretized Gaussian distributions with
    % a bin size of 0.001. While the values are computed explicitly, the
    % probabilities should be very similar. While slightly slower, this
    % method is simple and easy in MATLAB.
    Bx = zeros(nStates,traceLen);
    for s=1:nStates,
        Bx(s,:) = 0.001*normpdf( trace, mu(s), sigma(s) );
    end
    
    % Find the most likely viterbi path in model space
    [vPath, vLL] = forward_viterbi(start_p, trans_p, Bx);
    
    % Convert sequence of state assignments to dwell-times at each state
    % and add this new idealization to the output
    idealization(i,1:traceLen) = vPath;
    dwt{i} = RLEncode(vPath);
    LL(i) = vLL;
end

% Add offsets to relate idealization back to raw data
offsets = (0:(nTraces-1))*nFrames;

end %function idealize


