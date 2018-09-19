function LL = mplIter(data, dt, model, params)
% Maximum Point Likelihood algorithm (MPL)
%
%   LL = mplIter(DATA, DT, p0, CLASSIDX, RATEMASK, PARAMS)
%   Runs one iteration of the Maximum Point Likelihood algorithm (MPL).
%   DATA is a matrix of FRET observations (one molecule trace per row).
%   DT is the experimental sampling interval of the data in seconds.
%   P0 is a vector of initial state probabilities.
%   CLASSIDX gives the class number of each state.
%   RATEMASK is true for rates that are to be varied in optimization, with
%     the same size as the full rate matrix (Q).
%   PARAMS is a vector of all optimizable parameters arranged as follows:
%     [mu1 mu2 ... stdev1 stdev2 ... k12 k13 .. k21 k23 .. k31 k32 .. ...].
%   The rate constants in this list are obtained as: Q(rateMask)'.
%
%   LL is the -log(likelihood) for current parameters: -log[ P(data|model) ]
%   This allows for use with minimizers (e.g., fmincon).
%  
%   See Qin et al (2000) Biophys J 79, pg. 1915-1927 for algorithm details.
%
%   See also: mplOptimize, bwOptimize, milOptimize, batchKinetics.

%   Copyright 2018 Cornell University All Rights Reserved.


narginchk(4,4);
nargoutchk(0,1);


%% Process input arguments

% Make sure vectors are the correct orientation.
p0 = reshape(model.p0, 1, model.nStates);
if isvector(data)
    data = reshape(data, 1, numel(data));
end


% Unpack parameters from fminunc input vector
nMu = sum(~model.fixMu);
nSigma = sum(~model.fixSigma);

mu = model.mu;
sigma = model.sigma;

mu(~model.fixMu)    = params( 1:nMu );
sigma(~model.fixSigma) = params( nMu + (1:nSigma) );

Q = model.rates;
I = logical(eye(model.nStates));
rateMask = ~I & Q~=0 & ~model.fixRates;
Q(rateMask) = params( nMu+nSigma + 1:end );

I = logical(eye(model.nStates));
Q(I) = -sum(Q,2);  %normalize so rows sum to zero

% Calculate discrete-time transition probabilities (A) from rate matrix
transitionProb = expm( Q*dt );

% Duplicate emission parameters for degenerate states
mu = mu(model.class);
sigma = sigma(model.class);


%% Calculate log likelihood and partial derivatives for each trace.
LL = 0;

for n=1:size(data,1)
    % Remove frames after donor photobleaching (which are precisely zero).
    trace = data(n,:);
    bleachFrame = find( trace==0, 1, 'first' );
    if ~isempty(bleachFrame)
        trace = trace(1:bleachFrame-1);
        if numel(trace)<5, continue; end  %skip extremely short traces
    end
    
    % Get partial probabilities using the forward-backward algorithm
    [LLtrace,~] = forwardBackward( p0, transitionProb, trace, mu, sigma );
    LL = LL+LLtrace;

end %for each trace

% Return opposite of LL and dLL since to convert maximizing LL into minimizing
% -LL for use with fminunc/fmincon.
LL = -LL;


end  %function mplIter


