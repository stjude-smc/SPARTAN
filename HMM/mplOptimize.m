function [idl,optModel,LL] = mplOptimize(fret, dt, model, optionsInput)
% mplOptimize  Maximum Point Likelihood model optimization
%
%   [optModel,LL] = mplOptimize(DATA, MODEL, OPTIONS)
%   Finds a model (optModel) that maximizes the probability of the experimental
%   data (fret values) given the model parameters, expressed as the log 
%   likelihood (LL). For algorithm details, see the mplIter function, which 
%   implements the likelihood function and its partial derivatives for 
%   optimization with a standard optimizer (fmincon here).
%
%   This algorithm is similar to maximum interval likelihood (MIL) in that it
%   directly optimizes the log likelihood using analytical partial derivatives,
%   but that function only considers dwell-time information.
%
%   DATA is a TracesFret object, providing the input data in the 'fret' field.
%   MODEL is a QubModel object.
%   OPTIONS is a struct with additional settings (all fields optional):
%     .maxIter:  maximum number of iterations (200)
%     .convLL:   termination condition (tolerance) in LL values.
%     .convGrad: termination condition for parameter values.
%     .verbose:  display information about each iteration for debugging.
%     .quiet:    do not display waitbar or any diagnostic text.
%
%   See also: mplIter, milOptimize, bwOptimize, batchKinetics.

%   Copyright 2018 Cornell University All Rights Reserved.


narginchk(3,4);
nargoutchk(1,3);

nStates = model.nStates;
I = logical(eye(nStates));
rateMask = ~I & model.rates~=0;
nRates = sum(rateMask(:));
dt = dt/1000;  %convert to seconds/frame


% Define default optional arguments, mostly taken from fmincon
options.maxIter  = 200;
options.convLL   = 10^-6;   %OptimalityTolerance in fmincon (sort of)
options.convGrad = 10^-8;  %StepTolerance in fmincon
options.verbose  = true;
if nargin>=4
    options = mergestruct(options, optionsInput);
end
if isfield(options,'exclude') && any(options.exclude)
    error('Excluding traces not supported by MPL yet');
end

% Construct options for fmincon.
fminopt = optimoptions('fmincon', 'Algorithm','trust-region-reflective');  %'SQP');  %,'interior-point');
if options.verbose
    fminopt.Display='iter';
    fminopt.OutputFcn = @outfun;
end

try
    % Legacy version for MATLAB R2015a and before
    fminopt.MaxIter = options.maxIter;
    fminopt.TolX    = options.convGrad;
    fminopt.TolFun  = options.convLL;
    fminopt.GradObj = 'on';
catch
    fminopt.MaxIterations            = options.maxIter;
    fminopt.StepTolerance            = options.convGrad;
    fminopt.OptimalityTolerance      = options.convLL;
    fminopt.SpecifyObjectiveGradient = 'true';
end

% Define optimization function and initial parameter values.
% NOTE: mplIter optimizes the variance (sigma^2).
optFun = @(x)mplIter(fret(1:100,:), dt, model.p0, model.class, rateMask, x);
x0 = [ model.mu(:)'  model.sigma(:)'.^2  model.rates(rateMask)'  ];

% Constrained version.
% NOTE: fmincon doesn't like x0=lb or x0=ub in any parameter and will 'fix'
% things in ways that can break the bounds, so I added -eps to the minimum rate
% in hopes that the actual minimum is 0.
lb = [ -0.5*ones(1,nStates)  0.01^2*ones(1,nStates)   -eps*ones(1,nRates) ];
ub = [  1.5*ones(1,nStates)  0.15^2*ones(1,nStates)  10/dt*ones(1,nRates) ];
[optParam,LL] = fmincon( optFun, x0, [],[],[],[],lb,ub,[],fminopt );

% Save results
optModel = copy(model);  %do not modify model in place
% optModel.p0    = beta(:,1);   %FIXME: get p0 from forward-backward directly
optModel.mu    = optParam(1:nStates);
optModel.sigma = sqrt( optParam(nStates + (1:nStates)) );

% rates = exp( optParam(2*nStates + 1:end) );
rates = optParam(2*nStates + 1:end);
optModel.rates(rateMask) = rates;
optModel.rates(I) = 0;

% Idealize traces if requested.
idl = idealize( fret, [to_col(optModel.mu) to_col(optModel.sigma)], optModel.p0, optModel.calcA(dt) );


end %function mplOptimize




%%
function stop = outfun(x,optimValues,state)
% Called in each iteration of fmincon optimizer to track the progress of
% optimization for debugging

persistent X;
persistent dX;
stop=false;  %if true, optimizer terminates early.

switch state
    case 'init'
        X  = zeros( 1000, numel(x) );
        dX = zeros( 1000, numel(x) );
          
    case 'iter'
        % Keep track of parameter values at each iteration
        X(optimValues.iteration+1,:)  = x;
        dX(optimValues.iteration+1,:) = optimValues.gradient;
          
    case 'done'
        % Once optimizer completes, plot how parameters change over iterations
        X  = X(1:optimValues.iteration,:);
        dX = dX(1:optimValues.iteration,:);
        
        figure;
        for i=1:numel(x)
            ax(1,i) = subplot(2,numel(x),i);
            plot(1:optimValues.iteration, X(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
            if i==1
                ylabel('Param. Value');
            end
            %title('mu2');
            
            ax(2,i) = subplot(2,numel(x),numel(x)+i);
            plot(1:optimValues.iteration, dX(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
            if i==1
                xlabel('Iteration');
                ylabel('Gradient');
            end
        end
        linkaxes(ax(:), 'x');
        xlim(ax(1), [1,optimValues.iteration]);
end


end %function outfun
