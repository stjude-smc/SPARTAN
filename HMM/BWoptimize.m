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
[ok,msg] = model.verify();
if ~ok, error(msg);
elseif ~isempty(msg), warning(msg);
end

% Set default values for any paramaters not specified.
params.maxItr   = 100;
params.convLL   = 1e-5;
params.convGrad = 1e-5;
params.verbose  = true;
params.fixRates = false; %FIXME: should ultimately be in the model.
params.zeroEnd  = false;
params.seperately = false;
params.updateModel = false;
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



%% -----------------------  RUN BAUM-WELCH  ----------------------- %%

wbh = waitbar(0,'Running Baum-Welch...');

rateMask = model.rates~=0;

% Remove excluded traces from analysis.
origSize = size(observations);  %before exlusions
observations = observations( ~params.exclude, : );
observations = max(-0.5, min(1.5,observations));  %clip outlier values

LL = zeros(0,1);
dL = Inf;

% Monitor parameter values for debugging
%muIter    = zeros(model.nClasses, params.maxItr);
%sigmaIter = zeros(model.nClasses, params.maxItr);
%QIter     = zeros(sum(rateMask(:)),  params.maxItr);

% Initialize parameter values
A = model.calcA(sampling/1000);
p0 = to_row(model.p0);
mu = to_row(model.mu);
sigma = to_row(model.sigma);

for n = 1:params.maxItr
    % Reestimate model parameters using the Baum-Welch algorithm
    [LL(n),p0_new,A_new,mu_new,sigma_new] = BWiterate( p0, A, observations, mu, sigma, model.class );
   
    if any(isnan(A(:))) || any(isnan(p0)) || any( isnan(mu) | isnan(sigma) )
        disp('Warning: NaN parameter values in BWoptimize');
    end
    if any(A(:)<0) || any(p0<0), error('Invalid A or p0'); end
  
    % Enforce crude re-estimation constraints.
    % FIXME: are there better ways apply constraints? 
    mu(~model.fixMu)       = mu_new(~model.fixMu);
    sigma(~model.fixSigma) = sigma_new(~model.fixSigma);
    
    % Calculate step size (for convergence)
    step = norm( [p0_new(:)'-p0(:)' A_new(:)'-A(:)' mu_new(:)'-mu(:)' sigma_new(:)'-sigma(:)'] );
    
    % Update model parameters
    p0 = p0_new;
    A  = A_new;
    
    if params.updateModel
        % Modify input model for display in batchKinetics GUI.
        Q = logm(A) / (sampling/1000);
        model.rates(rateMask) = Q(rateMask);
        [model.mu, model.sigma, model.p0] = deal(mu,sigma,p0);
        
        %muIter(:,n) = mu;
        %sigmaIter(:,n) = sigma;
        %QIter(:,n) = Q(rateMask)';
    end
    
    % Update progress bar
    if n>1, dL = LL(n)-LL(n-1); end
    progress = max( [n/params.maxItr, log10(dL)/log10(params.convLL), log10(step)/log10(params.convGrad)] );
    if ~ishandle(wbh) || ~isvalid(wbh)
        disp('Baum-Welch manually terminated before convergence by user');
        break;
    end
    waitbar(progress, wbh);
    
    if params.verbose
        fprintf( '   Iter %d:\t%.5e\t%.2e\t%.2e\n', n, LL(n), dL, step);
        disp( [mu(model.class)' sigma(model.class)' p0' A] );
    end
    
    % Check for convergence
    if dL<0, disp('Warning: Baum-Welch is diverging!'); end
    if abs(dL)<=params.convLL, break; end
    if step<=params.convGrad, break; end
end
if n>=params.maxItr
    disp('Warning: Baum-Welch exceeded maximum number of iterations');
end

if ishandle(wbh) && isvalid(wbh)
    close(wbh);
end

% Idealize traces using optimized model
imu    = mu( model.class );
isigma = sigma( model.class );
idl = idealize( observations, [to_col(imu) to_col(isigma)], p0, A );
classes = [0; to_col(model.class)];
idl = reshape( classes(idl+1), size(idl) );

% Give zeros to excluded traces to indicate they were not idealized.
idlTotal = zeros( origSize );
idlTotal( ~params.exclude, :) = idl;


% Save results
optModel = copy(model);
[optModel.mu, optModel.sigma, optModel.p0] = deal(mu, sigma, p0);
LL = LL(end);

% Convert transition probabilities to rates, setting any rates to zero
% that are not connected in the original model. There doesn't seem to be
% a way to enforce this constraint during optimization with Baum-Welch.
% NOTE: A/dt gives approximate rates, but logm(A)/dt seems more correct.
optModel.rates = logm(A) / (sampling/1000);
if any(optModel.rates<0)
    optModel.rates = (A-eye(size(A))) / (sampling/1000);
    disp('logm(A) return negative rates');
end
optModel.rates( model.rates==0 ) = 0;

disp(toc);


end %function BWoptimize




%% ----------- BAUM-WELCH PARAMETER ESTIMATION ROUTINE -----------

function [LLtot,p0,A,mu,sigma] = BWiterate( p0, A, observations, mu, sigma, classidx )
% Optimize model parameters using the Baum-Welch algorithm

nTraces = size(observations,1);
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

    % Calculate transition probabilities at each point in time using the
    % forward/backward algorithm.
    [LL,~,~,gamma,E] = forwardBackward( p0, A, obs, mu(classidx), sigma(classidx) );

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

% Combine probabilities from each class of degenerate states
nClass = max(classidx);
gamma_combined = zeros( size(gamma_tot,1), nClass );

for c=1:max(classidx)
    gamma_combined(:,c) = sum( gamma_tot(:,classidx==c), 2 );
end

% Reestimate emmission parameters (stdev and mean).
% Weighted average using gamma(t,i) = P(state i at time t) weights.
% FIXME: check whether mu and sigma should have simultaneous updates or not...
gamma_combined = bsxfun(@rdivide, gamma_combined, sum(gamma_combined)); 

for c = 1:nClass
    mu(c) = obs_tot * gamma_combined(:,c);                      % =sum(data_i * gamma_i)
    sigma(c) = sqrt( (obs_tot-mu(c)).^2 * gamma_combined(:,c) );   % 1xN * Nx1 = 1x1
end

% Prevent stdev from converging to zero (e.g., single data point in state).
sigma = max(0.01,sigma);



end %function BWiterate


