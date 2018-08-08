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



%% ----------- PARSE INPUT PARAMETER VALUES -----------
%FIXME: use inputParser
tic;
narginchk(3,4);
nargoutchk(0,2);

% Verify model input.
nStates = size(model.rates,1);
nClass  = numel(model.mu);
assert( nStates==nClass, 'Aggregate states not supported.' );
[ok,msg] = model.verify();
if ~ok, error(msg);
elseif ~isempty(msg), warning(msg);
end

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

wbh = waitbar(0,'Running Baum-Welch...');

for n = 1:params.maxItr
    % Reestimate model parameters using the Baum-Welch algorithm
    [LL(n),A,mu,sigma,p0] = BWiterate(observations,A,mu,sigma,p0);
   
    if any(isnan(A(:))) || any( isnan(mu) | isnan(sigma) | isnan(p0))
        disp('Warning: NaN parameter values in BWoptimize');
    end
    if any(A(:)<0) || any(p0<0), error('Invalid A or p0'); end
  
    % Enforce crude re-estimation constraints.
    % FIXME: are there better ways apply constraints? 
    mu(model.fixMu) = mu_start(model.fixMu);
    mu(model.fixSigma) = sigma_start(model.fixSigma);
    
    % Update progress bar
    if n>1, dL = LL(n)-LL(n-1); end
    progress = max( n/params.maxItr, log10(dL)/log10(params.convLL) );
    if ~ishandle(wbh) || ~isvalid(wbh)
        error('spartan:op_cancelled','Operation cancelled by user');
    end
    waitbar(progress, wbh);
    
%     if ~params.quiet
        fprintf( '   Iter %d: %.5e %.2e\n', n, LL(n), dL);
        disp( [mu' sigma' p0' A] );
%     end
    
    % Check for convergence
    if dL<0, disp('Warning: Baum-Welch is diverging!'); end
    if abs(dL)<=params.convLL, break; end
end
if n>=params.maxItr
    disp('Warning: Baum-Welch exceeded maximum number of iterations');
end
close(wbh);


% Return idealized traces if requested
idlTotal = zeros( origSize );
idl = idealize( observations, [to_col(mu) to_col(sigma)], p0, A );
idlTotal( ~params.exclude, :) = idl;

% Save results
optModel = copy(model);
optModel.mu    = mu;
optModel.sigma = sigma;
optModel.p0    = p0;
A( logical(eye(size(A))) ) = 0;
optModel.rates = A/(sampling/1000);
LL = LL(end);

disp(toc);


end %function BWoptimize




%% ----------- BAUM-WELCH PARAMETER ESTIMATION ROUTINE -----------

function [LLtot,A,mu,sigma,p0] = BWiterate( observations, A, mu, sigma, p0 )
% Optimize model parameters using the Baum-Welch algorithm

nTraces = size(observations,1);
nStates = size(A,1);
[LLtot, p0tot, Etot] = deal(0);
[gamma_tot,obs_tot] = deal([]);

for n=1:nTraces
    obs = observations(n,:);

    % Remove donor bleaching at the end of each trace
    nFrames = find(obs~=0, 1,'last')-1;
    if ~isempty(nFrames)
        if nFrames<5, continue; end  %ignore very short traces
        obs = obs(1:nFrames);
    end

    % Calculate emmission probabilities at each timepoint
    B = zeros(nFrames, nStates);
    for i=1:nStates
        B(:,i) = exp(-0.5 * ((obs - mu(i))./sigma(i)).^2) ./ (sqrt(2*pi) .* sigma(i));
    end

    % Calculate transition probabilities at each point in time using the
    % forward/backward algorithm.
    [LL,~,~,gamma,E] = BWtransition( p0, A, B );

    LLtot = LLtot + LL;
    Etot = Etot+E;
    p0tot = p0tot + gamma(1,:);

    % Accumulate trace data and most likely state assignments for
    % emission parameter re-estimation below.
    gamma_tot = [gamma_tot gamma'];    %#ok<AGROW>
    obs_tot = [obs_tot obs];           %#ok<AGROW>
end
gamma_tot = gamma_tot';

% Normalize. FIXME: what about short traces that were skipped???
p0    = p0tot/nTraces;
LLtot = LLtot/nTraces;

% Reestimate transition probability matrix (A).
% SUM_t(E) is the expected number of each type of transition.
A = bsxfun(@rdivide, Etot, sum(Etot,2));   %normalized so rows sum to 1

% Reestimate emmission parameters (stdev and mean).
% Weighted average using gamma(t,i) = P(state i at time t) weights.
% FIXME: check whether mu and sigma should have simultaneous updates or not...
gamma_tot = bsxfun(@rdivide, gamma_tot, sum(gamma_tot)); 

for n = 1:size(gamma_tot,2)  %for each state
    mu(n) = obs_tot * gamma_tot(:,n);                      % =sum(data_i * gamma_i)
    sigma(n) = sqrt( (obs_tot-mu(n)).^2 * gamma_tot(:,n) );   % 1xN * Nx1 = 1x1
end

% Prevent stdev from converging to zero (e.g., single data point in state).
sigma = max(0.02,sigma);



end %function BWiterate


