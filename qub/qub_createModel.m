function model = qub_createModel(nStates)
% qub_createModel  Creates a simple starting point model
%     
%   [MODEL] = qub_createModel( NSTATES )
%   Creates a model object which can be used for optimization in either
%   Baum-Welch methods or QuB Methods (MIL).
%   See also qub_saveModel, qub_loadModel, qub_milOptimize, qub_skmIdealize
%
%   Member |  size | Description
%   --------------------------------------------------------------------------
%   .mu        1xN   Mean FRET value of each state
%   .sigma     1xN   Stdev of FRET values of each state
%   .rates     NxN   Matrix of rates (aka Q matrix); row=src, col=dest
%   .tp        NXN   Matrix of transition probabilities (aka A matrix)
%   .p0        1xN   Initial state probabilities
%   .pt        1xN   Time-averaged (steady-state) occupancy in each state
%   .fixMu     1xN   Set to 1 for to fix FRET value of each state 
%   .fixSigma  1xN   Set to 1 for to fix stdev of FRET for each state 
%   .fixRates  NxN   Set to 1 to fix each rate (NOT IMPLEMENTED!!!)
%   --------------------------------------------------------------------------
%   NOTE: pt is optional.
%   NOTE: fixRates does NOT work with the Baum-Welch algorithm!
%
%  http://www.qub.buffalo.edu

vectorSize = [1,nStates];

% MODEL PARAMETERS
model.mu = 0:1/(nStates+1):1;
model.mu = model.mu(1:end-1);

model.sigma = repmat( 0.061, vectorSize );

model.rates = repmat( 2.0, [nStates,nStates] );
model.rates( logical(eye(nStates,nStates)) ) = 0;

model.p0 = repmat( 1/nStates, vectorSize );

model.pt = model.p0;


% MODEL-SPECIFIC FITTING CONSTRATINGS
model.fixMu    = repmat( 0, vectorSize );
model.fixSigma = repmat( 0, vectorSize );
model.fixRates = repmat( 0, [nStates,nStates] );
