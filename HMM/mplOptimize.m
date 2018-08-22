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
nClass  = model.nClasses;

I = logical(eye(nStates));
rateMask = ~I & model.rates~=0;
nRates = sum(rateMask(:));
dt = dt/1000;  %convert to seconds/frame


% Define default optional arguments, mostly taken from fmincon
options = struct('maxIter',200,  'convLL',10^-6, 'convGrad',10^-6, ...
                 'verbose',true, 'updateModel',false);
if nargin>=4
    options = mergestruct(options, optionsInput);
end
if isfield(options,'exclude') && any(options.exclude)
    warning('Excluding traces not supported by MPL yet');
end

% Construct options for fmincon.
fminopt = optimoptions('fmincon', 'UseParallel',cascadeConstants('enable_parfor') );
if options.verbose
    fminopt.Display='iter';
    fminopt.OutputFcn = @outfun;
end
fminopt.MaxIter = options.maxIter;
fminopt.TolX    = options.convGrad;
fminopt.TolFun  = options.convLL;

% Run fmincon optimizer, with loose contraints to aid convergence.
optFun = @(x)mplIter(fret, dt, model.p0, model.class, rateMask, x);
x0 = [ model.mu(:)' model.sigma(:)' model.rates(rateMask)' ];

lb = [ -0.3*ones(1,nClass)  0.01*ones(1,nClass)     zeros(1,nRates) ];
ub = [  1.2*ones(1,nClass)  0.12*ones(1,nClass) 3/dt*ones(1,nRates) ];
[optParam,LL] = fmincon( optFun, x0, [],[],[],[],lb,ub,[],fminopt );

% Save results
optModel = copy(model);  %do not modify model in place
optModel.mu    = optParam(1:nClass);
optModel.sigma = optParam(nClass + (1:nClass));
optModel.rates(:) = 0;
optModel.rates(rateMask) = optParam(2*nClass + 1:end);

% Idealize traces if requested.
% FIXME: idealize should properly handle degenerate states
imu    = optModel.mu( optModel.class );
isigma = optModel.sigma( optModel.class );
idl = idealize( fret, [to_col(imu) to_col(isigma)], optModel.p0, optModel.calcA(dt) );
classes = [0; to_col(model.class)];
idl = reshape( classes(idl+1), size(idl) );



%%
function stop = outfun(x,optimValues,state)
% Called in each iteration of fmincon optimizer to track the progress of
% optimization for debugging

persistent X;
persistent dX;
persistent wbh;
stop=false;  %if true, optimizer terminates early.
itr = optimValues.iteration;

switch state
    case 'init'
        X  = zeros( 1000, numel(x) );
        dX = zeros( 1000, numel(x) );
        wbh = waitbar(0,'Running MPL...');
          
    case 'iter'
        % Keep track of parameter values at each iteration
        X(itr+1,:)  = x;
        dX(itr+1,:) = optimValues.gradient;
        
        if options.updateModel
            model.mu    = x( 1:nClass );
            model.sigma = x( nClass + (1:nClass) );
            model.rates(rateMask) = x( 2*nClass + 1:end );
            drawnow;
        end
        
        progress = max( itr/options.maxIter, log10(optimValues.stepsize)/log10(options.convLL) );
        if ~ishandle(wbh) || ~isvalid(wbh)
            stop=true; %user closed waitbar
            return;
        end
        if itr>1,  waitbar( max(0,min(1,progress)), wbh );  end
          
    case 'done'
        if ishandle(wbh), close(wbh); end
        if ~options.verbose, return; end
        
        % Once optimizer completes, plot how parameters change over iterations
        X  = X(1:itr,:);
        dX = dX(1:itr,:);
        
        figure;
        for i=1:numel(x)
            ax(1,i) = subplot(2,numel(x),i);
            plot(1:itr, X(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
            if i==1
                ylabel('Param. Value');
            end
            %title('mu2');
            
            ax(2,i) = subplot(2,numel(x),numel(x)+i);
            plot(1:itr, dX(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
            if i==1
                xlabel('Iteration');
                ylabel('Gradient');
            end
        end
        linkaxes(ax(:), 'x');
        xlim(ax(1), [1,itr+1]);
        disp(X);
end

end %function outfun




end %function mplOptimize
