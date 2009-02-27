function [results,errorResults,meanResults] = BWoptimize( ...
          observations, mu_start, sigma_start, varargin )
% BWOPTIMIZE  Find HMM parameter values that best explain data
%
%    [LL,A,mu,sigma,p0] = BWoptimize( OBSERVATIONS, nStates )
%    Uses the Baum-Welch parameter optimization method to discovery HMM
%    parametrs which best fit the observed data (OBSERVATIONS).  The number
%    of states in the optimized model can be specified (nStates).  If it is
%    not, the number is automatically predicted (NOT IMPLEMENTED YET!)
%    
%    [LL,A,mu,sigma,p0] = BWoptimize( OBSERVATIONS, MU, SIGMA, ... )
%    Same as above, with initial values for state mean (MU) and stdev
%    (SIGMA) of emission.
%
%    OPTIONAL Parameter values:
%      'A_start'  ->  transition probability matrix (starting value)
%      'p0_start' ->  Initial probabilities (start value)
%
%    Initial conditions can be held constant during re-estimation using
%    the arguments: 'FixMu' and 'FixSigma'.  Each is a mask the same size
%    as the input parameter, with a value of zero allowing that value to
%    be re-estimated and 1 holds it constant.
% 
%    NOTE: this proceedure is very sensitive to initial conditions.  A and
%    P0 will generally converge to the same point, but the emission
%    parameters (mu,sigma) may not!  Try several starting points to be
%    sure you have the best fit (highest LL).

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

nTraces = size(observations,1);

%% ----------- PARSE INPUT PARAMETER VALUES -----------

% PARSE REQUIRED PARAMETER VALUES
if nargin<=2 && numel(mu_start)==1,
    nStates  = mu_start;
    mu_start = 0:(1/(nStates-1)):1;  % emmission means for each state   
else
    nStates = numel(mu_start);
end

if ~exist('sigma_start','var') || isempty(sigma_start),
    warning('BW:BWoptimize:sigma','Initial sigma values not provided, using default (0.061)');
    sigma_start = repmat(0.061, size(mu_start));
end

% PARSE OPTIONAL PARAMETER VALUES: initial kinetic parameter values
optargs = struct( varargin{:} );

if isfield(optargs,'A_start') % transition probability matrix
    A_start = optargs.A_start;
    assert( all(A_start(:)>0) );
else
    A_start = repmat(0.05, nStates,nStates);
end
A_start( find(eye(nStates)) ) = 1-sum(A_start,2);

if isfield(optargs,'p0_start') % state initial probabilities
    p0_start = optargs.p0_start;
else
    p0_start = A_start^10000;
    p0_start = p0_start(1,:);
    p0_start = p0_start/sum(p0_start);
end

% PARSE OPTIONAL PARAMETER VALUES: re-estimation constraints
if isfield(optargs,'FixMu')
    muMask = 1-optargs.FixMu;
else
    muMask = ones(1,nStates); %all are re-estimated (=1)
end

if isfield(optargs,'FixSigma')
    sigmaMask = 1-optargs.FixSigma;
else
    sigmaMask = ones(1,nStates); %all are re-estimated (=1)
end

% PARSE OPTIONAL PARAMETER VALUES: options
if isfield(optargs,'maxItr')
    maxItr = optargs.maxItr;
else
    maxItr = 100;
end

if isfield(optargs,'LLConv')
    LLConv = optargs.LLConv;
else
    LLConv = 1e-4;
end

if isfield(optargs,'bootstrapN')
    bootstrapN = optargs.bootstrapN;
else
    bootstrapN = 0;
end


%% ----------- RUN BAUM-WELCH, ESTIMATE ERROR W/ BOOTSTRAPPING -----------

if nTraces>1 && bootstrapN>1,
    h = waitbar(0,'Baum-Welch parameter estimation...');
end

params.maxItr = maxItr;
params.LLConv = LLConv;
params.muMask = muMask;
params.sigmaMask = sigmaMask;

% Obtain parameter estimates for complete sample.
% If bootstraping is disabled, this is all we do.
initialValues = {A_start mu_start sigma_start p0_start};

results = BWrun( observations, initialValues, params );


% Run boostrapping proceedure for error estimation -- 
% Initial conditions are set close to final results for speed.
if bootstrapN>1
    if nTraces>1,  waitbar( 1/(bootstrapN+1),h );  end
    
    initialValues = {results.A, results.mu, results.sigma, results.p0};

    for i=1:bootstrapN,
        % Generate bootstrap dataset (random set, with replacement)
        selection = floor(rand(nTraces,1)*nTraces)+1;
        bootstrapData = observations(selection,:);

        % Optimize parameter values, add results to output.
        bootstrapResults(i) = BWrun( bootstrapData, initialValues, params );
        
        if nTraces>1,  waitbar( (i+1)/(bootstrapN+1), h );  end
    end

    % Estimate errors
    fnames = fieldnames(bootstrapResults);

    for i=1:length(fnames)
        data = {bootstrapResults.(fnames{i})};
        [mu,sigma] = cellDist(data);

        meanResults.(fnames{i})  = mu;
        errorResults.(fnames{i}) = sigma;
    end
    
    if nTraces>1,  close(h);  end
else
    errorResults = struct([]);
    meanResults  = struct([]);
end

end %function BWoptimize



function [meanA,stdA] = cellDist( A )
% works with 1D or 2D input A{i}s.

assert( iscell(A) && ~isempty(A) );

N = numel(A);
B = zeros( N,numel(A{1}) );

for i=1:numel(A)
    assert( numel(A{i})==numel(A{1}) );
    B(i,:) = A{i}(:);
end

meanA = mean(B);
meanA = reshape(meanA,size(A{1}));

stdA  = std(B);
stdA = reshape(stdA,size(A{1}));

end






%% ----------- BAUM-WELCH ITERATOR -----------
function results = BWrun( observations, initialValues, params )

LL = zeros(0,1);

[A_start mu_start sigma_start p0_start] = initialValues{:};
[A mu sigma p0] = initialValues{:};


% Run Baum-Welch optimization iterations until convergence
for n = 1:params.maxItr

    [LL(n) A mu sigma p0 ps] = BWiterate2(observations,A,mu,sigma,p0);
   
    if any(isnan(A(:))) || any(isnan(mu)) || any(isnan(sigma)) || ...
        any(isnan(p0))  || any(isnan(ps))
        disp('NaN = BWoptimize');
    end
    
%     a = A(:)*25;
%     disp( sprintf( '%0.2f  ', a') );
  
    % Apply re-estimation constraints
    mu    = mu.*params.muMask       + mu_start.*(1-params.muMask);
    sigma = sigma.*params.sigmaMask + sigma_start.*(1-params.sigmaMask);
    
    disp( [sprintf( '   Iter %d: %f', n, LL(n)) sprintf('\t%.3f',mu)] );
    
    % Check for convergence of likelihood
    if n>1 && LL(n)-LL(n-1) <= params.LLConv,  break;  end
    
    drawnow;
end

% Save results
results.LL = LL(end);
results.A = A;
results.mu = mu;
results.sigma = sigma;
results.p0 = p0;
results.ps = ps;


end %function BWrun








%% ----------- BAUM-WELCH PARAMETER ESTIMATION ROUTINE -----------

function [LLtot,eA,eMU,eSIG,p0s,ps] = BWiterate2( observations, A, mus, sigmas, p0s )
% Performs one iteration of adapted Baum-Welch algorithm to return log-likelihood
% and re-estimated model parameters [A,mus,sigmas,p0s] as well as overall fractional
% occupancy of each state ps (which is not a model term but we get for
% free)
%
% DAB 2008.3.30
% DST 2008.6.29  Modified


Nstates = length(mus);

% Set maximal FRET to 1.0 to prevent overflows with zero probability
observations( observations>1.0 ) = 1.0;
observations( observations<-1.0 ) = -1.0;

% Get the photobleaching point of each trace
traceLengths = zeros( 1,size(observations,1) );

for i=1:size(observations,1),
    trace = observations(i,:);
    %NT = find(trace>=0.15, 1,'last');
    NT = find(trace~=0, 1,'last');
    if isempty(NT), NT=0; end
    traceLengths(i) = NT;
end

% Remove all traces that are very short
observations = observations( traceLengths>=5, : );
traceLengths = traceLengths( traceLengths>=5 );
nTraces = size(observations,1);

%
totalLength = sum(traceLengths-1);  %size of output array
starts = [1 cumsum(traceLengths-1)]-1;  %index into output array

% intialize variables assigned in loop
LLtot = 0;
p0tot = zeros(1,Nstates);
Etot  = zeros(Nstates,Nstates);
lamb_tot = zeros(totalLength,Nstates);
Osave = zeros( totalLength, 1 );

for n = 1:nTraces  %for each trace **with at least 5 datapoints**
  NT = traceLengths(n);
  obs = observations(n,1:NT);

  % Calculate 
  % LLc1 is the final backward LL.
  [E alphas LLc1] = BWtransition(obs,A,mus,sigmas,p0s);
  
  if any( isnan(LLc1) )
      disp( sprintf('BW:NaN: NaN found at %d',n) );
      break;
  end
  
  LLtot = LLtot + log(sum(alphas(NT,:)))+LLc1;

  % The re-estimated A matrix is simply the average of the transition
  % probability all every point in the dataset
  % (normalization at end of loop)
  lambdas = zeros(NT-1,Nstates);
%   Enorm = zeros(Nstates,Nstates);
  for t = 1:NT-1
%     Enorm = Enorm + E{t};
    Etot  = Etot + E{t};
	lambdas(t,:) = sum(E{t},2); %probability of the state at time=t
  end
  p0tot = p0tot + lambdas(1,:); %unnormalized initial probability
  
  % lamb_tot = unnormalized probability of being in state i at time t
  % with the concatinated observation sequence (Osave)
  lamb_tot( starts(n)+(1:NT-1), : ) = lambdas;  %really slow
  
  %
  Osave( starts(n)+(1:NT-1) ) = observations( n, 1:NT-1 );
end



% Re-estimate rate parameters (A matrix)
% Enorm contains information only from the final trace and therefor
% would be a poor way to do this, but that's the way it came.

% This is calculated by summing the transition probability matrix
% (trellice) at every time point and normalizing.
for m = 1:Nstates
    Etot(m,:) = Etot(m,:)/sum(lamb_tot(:,m));
end
eA = Etot;


% Re-estimate emmission parameters (stdev and mean)
Osave = Osave';

% Applies a weight (lambda_i = prob. of being in state i) to each
% datapoint to generate the mean and stdev for each state.
% Weights comes from forward/backward probabilities.
eMU = zeros(1,Nstates);
eSIG = zeros(1,Nstates);

for n = 1:Nstates
  gamma = lamb_tot(:,n)/sum(lamb_tot(:,n)); %normalized weights

  eMU(n) = Osave * gamma;                      % =sum(Osave_i * z_i)
  eSIG(n) = sqrt( (Osave-eMU(n)).^2 * gamma );   % 1xN * Nx1 = 1x1
  
  if isnan(eSIG(n)) || isnan(eMU(n)),
    disp( sprintf('BW:NaN: NaN found at %d',n) );
    break;
  end
end

% Maintain an absolute minimum for sigma to prevent errors.
% This happens when only a single datapoint contributes to parameter
% estimates -- there is no variance.
% Ideally, a higher probability should be assigned to the solution
% with a state that is not occupied...
eSIG(eSIG<0.02) = 0.02;

% Calc....
p0s = p0tot/nTraces;
LLtot = LLtot/nTraces;  % return LL per trial -- though not necessary...

% Calculate total occupancy of each state
ps = sum(lamb_tot)/sum(lamb_tot(:));

% disp( [ps' eMU' eSIG'] );
% disp('\n');

end
