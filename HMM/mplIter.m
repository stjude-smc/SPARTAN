function LL = mplIter(data, dt, p0, classidx, rateMask, params)
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


narginchk(6,6);
nargoutchk(0,1);


%% Process input arguments
nStates = numel(classidx);
nClass  = max(classidx);

% Make sure vectors are the correct orientation.
p0 = reshape(p0, 1, nStates);
if isvector(data)
    data = reshape(data, 1, numel(data));
end

% Unpack parameters from fminunc input vector
mu    = params( 1:nClass );
sigma = params( nClass + (1:nClass) );

Q = zeros(nStates);
Q(rateMask) = params( 2*nClass + 1:end );
I = logical(eye(nStates));
Q(I) = -sum(Q,2);

% Calculate discrete-time transition probabilities (A) from rate matrix
transitionProb = expm( Q*dt );

% Duplicate emission parameters for degenerate states
mu = mu(classidx);
sigma = sigma(classidx);


%% Calculate log likelihood and partial derivatives for each trace.
LL = 0;

for n=1:size(data,1)
    % Remove frames after donor photobleaching (which are precisely zero).
    trace = data(n,:);
    bleachFrame = find( trace==0, 1, 'first' );
    trace = trace(1:bleachFrame-1);
    nFrames = numel(trace);
    if nFrames<5, continue; end  %skip extremely short traces
    
    % Get partial probabilities using the forward-backward algorithm
    [LLtrace,~] = BWtransition( p0, transitionProb, trace, mu, sigma );
    LL = LL+LLtrace;

end %for each trace

% Return opposite of LL and dLL since to convert maximizing LL into minimizing
% -LL for use with fminunc/fmincon.
LL = -LL;


end  %function mplIter


