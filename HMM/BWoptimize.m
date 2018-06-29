function [idlTotal,optModel,LL] = BWoptimize( observations, sampling, model, paramInput )
% BWOPTIMIZE  Find HMM parameter values that best explain data
%
%    [LL,A,mu,sigma,p0] = BWoptimize( OBSERVATIONS, DT, MODEL, PARAMS )
%    Learn an optimal model 
%    Uses the Baum-Welch parameter optimization method to discovery HMM
%    parametrs which best fit the observed data (OBSERVATIONS).  The number
%    of states in the optimized model can be specified (nStates).
%    DT specifies the timestep (sampling interval) in sec.
%
%    See the definition of the MODEL structure in qub_createModel
%    for details of how to input initial conditions and constraints
%    for fitting. NOTE: .fixRates does not work with Baum-Welch.
%
%    OPTIONS can have any of the following specified:
%   Member   | size | Description
%   -----------------------------------------------------------------------
%   .maxItr     1x1   Max number of iterations before terminating
%   .convLL     1x1   Stop if LL converges within this threshold
%   .convGrad   1x1   Stop if all rates converge within this treshold
%   .showItr    1x1   Report on each iteration to the console (not impl.)
%   -----------------------------------------------------------------------
%   Other options, such as those specific to qub_milOptimize, are ignored.
%   Specifically, .seperately does NOT work.
%
%    Baum-Welch is sensitive to initial conditions of mu and sigma --
%    They may converge to different points depending on initial conditions
%    and convergence may be very slow.  A and pi, however, usually converge
%    very rapidly once mu and sigma converge.  Try several starting points
%    to be sure you have the best fit (highest LL).

%   Copyright 2007-2018 Cornell University All Rights Reserved.

% Format of the A-matrix:
%   1..N = lowest to highest FRET values (same order as mu, sigma, p0).
%   A(2,3) = Prob( i=2, j=3 )
%   A(3)   = A(3->1)  -- because MATLAB is column major
%
%   1->1  1->2  1->3  ..  1->N
%   2->1  2->2  2->3  ..  2->N
%   3->1  3->2  3->3  ..  3->N
%    ..    ..    ..   ..  .. 
%   N->1  N->2  N->3  ..  N->N
%


%% ----------- PARSE INPUT PARAMETER VALUES -----------
%FIXME: use inputParser

narginchk(3,4);
nargoutchk(0,2);

% Verify model input.
nStates = size(model.rates,1);
nClass  = numel(model.mu);
assert( nStates==nClass, 'SKM: aggregate states not supported.' );
assert(qub_verifyModel(model),'Invalid model');

% Set default values for any paramaters not specified.
params.maxItr   = 100;
params.convLL   = 1e-4;
% params.convGrad = 1e-3;  %not implemented yet.
params.quiet    = false;
params.fixRates = false; %FIXME: should ultimately be in the model.
params.zeroEnd  = false;
params.seperately = false;
params.exclude  = false( size(observations,1), 1 );

if nargin>=4
    params = mergestruct(params, paramInput);
end

% Check parameters
if isfield(model,'fixRates') && any(model.fixRates(:)),
    warning('BW:fixRates','Fixing specific rates not support!');
end
if params.seperately,
    warning('Individual trace analysis not yet supported!');
end
if isfield(params,'deadtime')
    warning('BW:deadtime','Deadtime correction not supported');
end
if isfield(params,'convGrad')
    warning('BW:convGrad','convGrad not yet implemented');
end



%% -----------------------  RUN BAUM-WELCH  ----------------------- %%
% Launch parallel pool for processig larger data sets in parallel.
wbh = waitbar(0,'Starting parallel pool...');
if numel(observations)>1e5 && cascadeConstants('enable_parfor'),
    pool = gcp;
    M = pool.NumWorkers;
else
    M = 0;
end

% Initialize parameter values
origSize = size(observations);  %before exlusions
observations = observations( ~params.exclude, : );
observations = max(-0.5, min(1.5,observations));  %remove outlier values
LL = zeros(0,1);
dL = Inf;
A  = model.calcA(sampling/1000);
p0 = to_row(model.p0);
[mu,mu_start]       = deal( to_row(model.mu) );
[sigma,sigma_start] = deal( to_row(model.sigma) );

% Run Baum-Welch optimization iterations until convergence
waitbar(0,wbh,'Running Baum-Welch...');

for n = 1:params.maxItr
    [LL(n),A,mu,sigma,p0] = BWiterate(observations,A,mu,sigma,p0,M);
   
    if any(isnan(A(:))) || any(isnan(mu)) || any(isnan(sigma)) || any(isnan(p0))
        disp('Warning: NaN parameter values in BWoptimize');
    end
  
    % Enforce crude re-estimation constraints.
    % FIXME: are there better ways apply constraints? 
    mu(model.fixMu) = mu_start(model.fixMu);
    mu(model.fixSigma) = sigma_start(model.fixSigma);
    
    % Update progress bar
    if n>1, dL = LL(n)-LL(n-1); end
    progress = max( n/params.maxItr, log10(dL)/log10(params.convLL) );
    waitbar(progress, wbh);
    
%     if ~params.quiet
        fprintf( '   Iter %d: %.2e %.2e\n', n, LL(n), dL);
%         disp( [mu' sigma' p0' A] );
%     end
    
    % Check for convergence
    if dL<0, disp('Warning: Baum-Welch is diverging!'); end
    if abs(dL)<=params.convLL, break; end
end
if n>=params.maxItr
    disp('Warning: Baum-Welch exceeded maximum number of iterations');
end
close(wbh);

% Save results
optModel = copy(model);
optModel.mu    = mu;
optModel.sigma = sigma;
optModel.p0    = p0;
A( logical(eye(size(A))) ) = 0;
optModel.rates = A/(sampling/1000);
LL = LL(end);

% Return idealized traces if requested
if nargin>=2
    idlTotal = zeros( origSize );
    fretModel = [to_col(optModel.mu) to_col(optModel.sigma)];
    idl = idealize( observations, fretModel, optModel.p0, model.calcA(sampling/1000) );
    idlTotal( ~params.exclude, :) = idl;
end


end %function BWoptimize




%% ----------- BAUM-WELCH PARAMETER ESTIMATION ROUTINE -----------

function [LLtot,eA,eMU,eSIG,p0s] = BWiterate( observations, A, mus, sigmas, p0s, M )
% Performs one iteration of adapted Baum-Welch algorithm to return log-likelihood
% and re-estimated model parameters [A,mus,sigmas,p0s] as well as overall fractional
% occupancy of each state ps (which is not a model term but we get for
% free)

assert( all(A(:)>=0), 'Invalid A matrix' );
assert( all(p0s(:)>=0), 'Invalid p0 vector' );

Nstates = length(mus);

% Get the photobleaching point of each trace
traceLengths = zeros( 1,size(observations,1) );

for i=1:size(observations,1),
    trace = observations(i,:);
    NT = find(trace~=0, 1,'last');
    if isempty(NT), NT=0; end
    traceLengths(i) = NT;
end

% Remove all traces that are very short
observations = observations( traceLengths>=5, : );
traceLengths = traceLengths( traceLengths>=5 );
nTraces = size(observations,1);

% intialize variables assigned in loop
LLtot = 0;
p0tot = zeros(1,Nstates);
Etot  = zeros(Nstates,Nstates);
lambda_cell = cell(nTraces,1);


parfor (n=1:nTraces, M)
  NT = traceLengths(n);
  obs = observations(n,:);

  % Calculate transition probabilities at each point in time using the
  % forward/backward algorithm. LLc1 is the final backward LL.
  [E,alphas,LLc1] = BWtransition( obs(1:NT), A, mus, sigmas, p0s );
  
%   if any( isnan(LLc1) )
%       fprintf('BW:NaN: NaN found at %d\n',n);
%       break;
%   end
  
  LLtot = LLtot + log(sum(alphas(NT,:)))+LLc1;

  % The re-estimated A matrix is simply the average of the transition
  % probability all every point in the dataset
  % (normalization at end of loop)
  lambdas = zeros(NT-1,Nstates);
  for t = 1:NT-1
    Etot  = Etot + E{t};
	lambdas(t,:) = sum(E{t},2); %probability of the state at time=t
  end
  p0tot = p0tot + lambdas(1,:); %unnormalized initial probability
  
  lambda_cell{n} = lambdas;
end


% lamb_tot = unnormalized probability of being in state i at time t
% with the concatinated observation sequence (Osave)
totalLength = sum(traceLengths-1);  %size of output array
starts = [1 cumsum(traceLengths-1)]-1;  %index into output array

lamb_tot = zeros(totalLength,Nstates);
Osave = zeros( totalLength, 1 );

for n=1:nTraces,
   NT = traceLengths(n);
   lamb_tot( starts(n)+(1:NT-1), : ) = lambda_cell{n};
   Osave( starts(n)+(1:NT-1) ) = observations( n, 1:NT-1 );
end
% lamb_tot = cat(1,lambda_cell{:});  %why doesn't this work?


% Calculate total occupancy of each state
% ps = sum(lamb_tot);
% ps = ps/sum(ps);


% Re-estimate rate parameters (A matrix)
% Enorm contains information only from the final trace and therefor
% would be a poor way to do this, but that's the way it came.

% This is calculated by summing the transition probability matrix
% (trellice) at every time point and normalizing.
for m = 1:Nstates
    Etot(m,:) = Etot(m,:)/sum(lamb_tot(:,m));
end
eA = Etot;


% Re-estimate emmission parameters (stdev and mean).
% Applies a weight (lambda_i = prob. of being in state i) to each
% datapoint to generate the mean and stdev for each state.
% Weights comes from forward/backward probabilities.
Osave = Osave';
eMU = zeros(1,Nstates);
eSIG = zeros(1,Nstates);

for n = 1:Nstates
  gamma = lamb_tot(:,n)/sum(lamb_tot(:,n)); %normalized weights

  eMU(n) = Osave * gamma;                      % =sum(Osave_i * z_i)
  eSIG(n) = sqrt( (Osave-eMU(n)).^2 * gamma );   % 1xN * Nx1 = 1x1
  
  if isnan(eSIG(n)) || isnan(eMU(n)),
    fprintf('BW:NaN: NaN found at %d\n',n);
    break;
  end
end

% Maintain an absolute minimum for sigma to prevent errors.
% This happens when only a single datapoint contributes to parameter
% estimates -- there is no variance.
% Ideally, a higher probability should be assigned to the solution
% with a state that is not occupied...
eSIG = max(0.02,eSIG);


p0s = p0tot/nTraces;
LLtot = LLtot/nTraces;  % LL per trial

assert( all(eA(:)>=0), 'Invalid estimated A matrix' );
assert( all(p0s>=0), 'Invalid estimated p0 vector' );


end %function BWiterate


