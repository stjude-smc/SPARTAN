function [optModel,LL] = milOptimize(dwt, dt, model, optionsInput)
% mplOptimize  Maximum Point Likelihood model optimization
%
%   [optModel,LL] = milOptimize(DWT, DT, MODEL, OPTIONS)
%   Finds a model (optModel) that maximizes the probability of the experimental
%   dwell times (dwt) given the model parameters, expressed as the log 
%   likelihood (LL). For algorithm details, see the milIter function, which 
%   implements the likelihood function and its partial derivatives for 
%   optimization with a standard optimizer (fmincon here).
%
%   DWT is a cell array, one element per trace, each of which contains class
%   numbers and dwell times (seconds) in the first and second column, resp.
%   DT is the measurement dead time in seconds (not used!).
%   MODEL is a QubModel object.
%   OPTIONS is a struct with additional settings (all fields optional):
%     .maxIter:  maximum number of iterations (200)
%     .convLL:   termination condition (tolerance) in LL values.
%     .convGrad: termination condition for parameter values.
%     .verbose:  display information about each iteration for debugging.
%
%   See also: milIter, mplIter, mplOptimize, bwOptimize, batchKinetics.

%   Copyright 2018 Cornell University All Rights Reserved.


narginchk(3,4);
nargoutchk(1,2);


% Construct a mask to select only rates from connected states.
% FIXME: this excludes connections where ONE rate is zero...
nStates = model.nStates;
I = logical(eye(nStates));
rateMask = ~I & model.rates~=0;
nRates = sum(rateMask(:));

% Define default optional arguments, mostly taken from fmincon
options.maxIter  = 100;
options.convLL   = 10^-8;   %OptimalityTolerance in fmincon (sort of)
options.convGrad = 10^-8;  %StepTolerance in fmincon (sort of)
options.verbose  = true;
if nargin>=4
    options = mergestruct(options, optionsInput);
end
if isfield(options,'exclude') && any(options.exclude)
    error('Excluding traces not supported by MIL yet');
end

% Construct options for fmincon.
fminopt = optimoptions('fmincon');
fminopt.DiffMaxChange = 1;
if options.verbose
    fminopt.Display='iter';
    fminopt.OutputFcn = @outfun;
end

try
    % Legacy version for MATLAB R2015a and before
    fminopt.MaxIter = options.maxIter;
    fminopt.TolX    = options.convGrad;
    fminopt.TolFun  = options.convLL;
%     fminopt.GradObj = 'on';
catch
    fminopt.MaxIterations            = options.maxIter;
    fminopt.StepTolerance            = options.convGrad;
    fminopt.OptimalityTolerance      = options.convLL;
%     fminopt.SpecifyObjectiveGradient = true;
end

% If requested, remove dark state dwells at end of each trace.
% This can make fitting converge much faster.
% FIXME: assumes dark state is lowest class value.
for i=1:numel(dwt)
    if dwt{i}(end,1)==1
        dwt{i}(end,:) = [];
    end
end
dwt( cellfun(@isempty,dwt) ) = [];  %remove any now-empty traces

% Define optimization function and initial parameter values.
% NOTE: mplIter optimizes the variance (sigma^2).
optFun = @(x)milIter(dwt, dt, model.p0, model.class, rateMask, x);
x0 = model.rates(rateMask)';

% Constrained version.
% NOTE: fmincon doesn't like x0=lb or x0=ub in any parameter and will 'fix'
% things in ways that can break the bounds, so I added -eps to the minimum rate
% in hopes that the actual minimum is 0.
lb =  -eps*ones(1,nRates);
ub = 10/dt*ones(1,nRates);
[optParam,LL] = fmincon( optFun, x0, [],[],[],[],lb,ub,[],fminopt );

% Save results
optModel = copy(model);  %do not modify model in place
optModel.rates(:) = 0;
optModel.rates(rateMask) = optParam;


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
        %disp(x);
          
    case 'done'
        % Once optimizer completes, plot how parameters change over iterations
        X  = X(1:optimValues.iteration,:);
        dX = dX(1:optimValues.iteration,:);
        %[idxin,idxout] = find(rateMask);  %state numbers for each transition type.
        
        figure;
        for i=1:numel(x)
            ax(1,i) = subplot(2,numel(x),i);
            plot(1:optimValues.iteration, X(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
            if i==1
                ylabel('Param. Value');
            end
            %title( sprintf('k%d->%d',idxin(i),idxout(i)) );
            
            ax(2,i) = subplot(2,numel(x),numel(x)+i);
            plot(1:optimValues.iteration, dX(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
            if i==1
                xlabel('Iteration');
                ylabel('Gradient');
            end
        end
        linkaxes(ax(:), 'x');
        try
            xlim(ax(1), [1,optimValues.iteration]);
        catch
        end
end


end %function outfun
