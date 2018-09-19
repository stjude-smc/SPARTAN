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

rateMask = ~eye(model.nStates) & model.rates & ~model.fixRates;
dt = dt/1000;  %convert to seconds/frame


% Define default optional arguments, mostly taken from fmincon
options = struct('maxItr',200,   'convLL',10^-6, 'convGrad',10^-6, ...
                 'verbose',true, 'updateModel',false);
if nargin>=4
    options = mergestruct(options, optionsInput);
end
% if isfield(options,'exclude') && any(options.exclude)
%     warning('Excluding traces not supported by MPL yet');
% end

% Construct options for fmincon.
fminopt = optimoptions('fmincon', 'UseParallel',cascadeConstants('enable_parfor') );
if options.verbose
    fminopt.Display='iter';
    fminopt.OutputFcn = @outfun;
end
fminopt.MaxIter = options.maxItr;
fminopt.TolX    = options.convGrad;
fminopt.TolFun  = options.convLL;

% Run fmincon optimizer, with loose contraints to aid convergence.
idl = zeros(size(fret));
fret = fret(~options.exclude,:);
optFun = @(x)mplIter(fret, dt, model, x);

nMu = sum(~model.fixMu);
nSigma = sum(~model.fixSigma);
nRates = sum(rateMask(:));

mu = model.mu(~model.fixMu);
sigma = model.sigma(~model.fixSigma);
rates = model.rates(rateMask);
x0 = [mu(:)' sigma(:)' rates(:)'];

lb = [ -0.3*ones(1,nMu)  0.01*ones(1,nSigma)      zeros(1,nRates) ];
ub = [  1.2*ones(1,nMu)  0.12*ones(1,nSigma) 10/dt*ones(1,nRates) ];
[optParam,LL] = fmincon( optFun, x0, [],[],[],[],lb,ub,[],fminopt );

% Save results
optModel = copy(model);  %do not modify model in place
optModel.mu(~model.fixMu)       = optParam(1:nMu);
optModel.sigma(~model.fixSigma) = optParam(nMu + (1:nSigma));
optModel.rates(rateMask)        = optParam(nMu+nSigma + 1:end);

% Idealize traces if requested.
% FIXME: idealize should properly handle degenerate states
imu    = optModel.mu( optModel.class );
isigma = optModel.sigma( optModel.class );
idl_sel = idealize( fret, [to_col(imu) to_col(isigma)], optModel.p0, optModel.calcA(dt) );
classes = [0; to_col(model.class)];
idl_sel = reshape( classes(idl_sel+1), size(idl_sel) );
idl(~options.exclude,:) = idl_sel;



%%
function stop = outfun(x,optimValues,state)
% Called in each iteration of fmincon optimizer to track the progress of
% optimization for debugging

% persistent X;
% persistent dX;
persistent wbh;
stop=false;  %if true, optimizer terminates early.
itr = optimValues.iteration;

switch state
    case 'init'
%         X  = zeros( 1000, numel(x) );
%         dX = zeros( 1000, numel(x) );
        wbh = waitbar(0,'Running MPL...');
          
    case 'iter'
        % Keep track of parameter values at each iteration
%         X(itr+1,:)  = x;
%         dX(itr+1,:) = optimValues.gradient;
        
        if options.updateModel
            model.mu(~model.fixMu)       = x(1:nMu);
            model.sigma(~model.fixSigma) = x(nMu + (1:nSigma));
            model.rates(rateMask)        = x(nMu+nSigma + 1:end);
            drawnow;
        end
        
        progress = max( itr/options.maxItr, log10(optimValues.stepsize)/log10(options.convLL) );
        if ~ishandle(wbh) || ~isvalid(wbh)
            stop=true; %user closed waitbar
            return;
        end
        if itr>1,  waitbar( max(0,min(1,progress)), wbh );  end
          
    case 'done'
        if ishandle(wbh), close(wbh); end
        if ~options.verbose, return; end
        
        % For debugging: plot change in parameters during optimization
%         X  = X(1:itr,:);
%         dX = dX(1:itr,:);
%         
%         figure;
%         for i=1:numel(x)
%             ax(1,i) = subplot(2,numel(x),i);
%             plot(1:itr, X(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
%             if i==1
%                 ylabel('Param. Value');
%             end
%             %title('mu2');
%             
%             ax(2,i) = subplot(2,numel(x),numel(x)+i);
%             plot(1:itr, dX(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
%             if i==1
%                 xlabel('Iteration');
%                 ylabel('Gradient');
%             end
%         end
%         linkaxes(ax(:), 'x');
%         xlim(ax(1), [1,itr+1]);
%         disp(X);
end

end %function outfun




end %function mplOptimize
