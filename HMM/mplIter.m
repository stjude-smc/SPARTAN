function [LL,dLL] = mplIter(data, dt, p0, classidx, rateMask, params)
% Maximum Point Likelihood algorithm (MPL)
%
%   [LL,dLL] = mplIter(DATA, DT, RATES, p0, CLASSIDX, PARAMS)
%   Runs one iteration of the Maximum Point Likelihood algorithm (MPL).
%   DATA are experimental observations (one molecule trace per row).
%   DT is the experimental sampling interval of the data in seconds.
%   P0 is a vector of initial state probabilities.
%   CLASSIDX gives the class number of each state.
%
%   PARAMS is a vector of all optimizable parameters arranged as follows:
%     [mu1 mu2 ... stdev1 stdev2 ... k12 k13 .. k21 k23 .. k31 k32 .. ...].
%
%   LL is the log likelihood for current parameters: log[ P(data|model) ].
%   dLL is the gradient of the LL function w/r/t each element in PARAMS.
%   Both values are negated to facilitate use with minimizers (e.g., fmincon).
%  
%   See Qin et al (2000) Biophys J 79, pg. 1915-1927 for algorithm details.
%
%   See also: mplOptimize, bwOptimize, milOptimize, batchKinetics.
%

%   Copyright 2018 Cornell University All Rights Reserved.

% disp(params);



%% Process input arguments
narginchk(6,6);
nargoutchk(0,2);

nStates = numel(p0);
p0 = reshape(p0, 1, nStates);
if isvector(data)
    data = reshape(data, 1, numel(data));
end
nTraces = size(data,1);

% Unpack parameters from fminunc input vector
mu     = params( 1:nStates );
sigma2 = params( nStates   + (1:nStates) );
sigma2(sigma2<1e-6) = 1e-7;  %steer optimizer from zero or negative sigma.
sigma = sqrt(sigma2);

Q = zeros(nStates);
Q(rateMask) = params( 2*nStates + 1:end );
I = logical(eye(nStates));
Q(I) = -sum(Q,2);


% Calculate transition probability matrix (A) using spectral matrices.
% See pg. 615 of 1995 book "Single Channel Recording". See eq. 6 in Qin 2000.
% expm(Q*dt) would be simpler, but we need the spectral matrices anyway.
[righteig,eigval] = eig(Q);
eigval = diag(eigval);
transitionProb = zeros(nStates);
spectralMatrix = cell(nStates,1);

for s=1:nStates
    lefteig = righteig^-1;
    spectralMatrix{s} = righteig(:,s) * lefteig(s,:);
    transitionProb = transitionProb  +  exp(eigval(s)*dt) * spectralMatrix{s};
end


%% Calculate log likelihood and partial derivatives for each trace.
LL = 0;

for n=1:nTraces
    % Remove frames after donor photobleaching (which are precisely zero).
    trace = data(n,:);
    bleachFrame = find( trace==0, 1, 'first' );
    trace = trace(1:bleachFrame-1);
    nFrames = numel(trace);
    if nFrames<10, continue; end  %skip extremely short traces
    
    % Calculate emmission probabilities at each timepoint
    observProb = zeros(nFrames, nStates);
    for i=1:nStates
        observProb(:,i) = exp(-0.5 * ((trace - mu(i))./sigma(i)).^2) ./ (sqrt(2*pi) .* sigma(i));
    end
    observProb = observProb+eps;  %prevent rows of all zeros
    
    % Get partial probabilities using the forward-backward algorithm
    [LLtrace,~] = BWtransition( p0, transitionProb, observProb );
    LL = LL+LLtrace;

end %for each trace

% Return opposite of LL and dLL since to convert maximizing LL into minimizing
% -LL for use with fminunc/fmincon.
LL = -LL;



end  %function mplIter


