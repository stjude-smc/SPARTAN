function [optModel,LL,exitflag] = milOptimize(dwt, dt, model, optionsInput)
% milOptimize  Maximum Interval Likelihood (MIL) model optimization
%
%   [optModel,LL] = milOptimize(DWT, DT, MODEL, OPTIONS)
%   Finds a model (optModel) that maximizes the probability of the experimental
%   dwell times (dwt) given the model parameters, expressed as the log 
%   likelihood (LL). For algorithm details, see the milIter function, which 
%   implements the likelihood function for optimization with fmincon.
%
%   DWT is a cell array, one element per trace, each of which contains class
%     numbers and dwell times (seconds) in the first and second column, resp.
%   DT is the measurement sampling interval in seconds.
%   MODEL is a QubModel object, providing initial parameter values.
%
%   OPTIONS is a struct with any additional settings (all fields optional):
%     .maxIter:  maximum number of iterations (150)
%     .convLL:   termination tolerance in LL values.
%     .convGrad: termination tolerance for parameter step size.
%     .verbose:  print information about each iteration for debugging.
%     .updateModel: modify input model object in each iteration (false).
%                Enables rate constants to be viewed during optimization.
%     .removeBleaching: remove final dwell in dark state (true).
%
%   See also: milIter, fmincon, batchKinetics, mplOptimize, bwOptimize.

%   Copyright 2018 Cornell University All Rights Reserved.


narginchk(3,4);
nargoutchk(1,3);


% Remove traces that were excluded from analysis
dwt = dwt( ~cellfun(@isempty,dwt) );
    
% Define default optional arguments, mostly taken from fmincon
options = hmmopt(mfilename);
if nargin>=4
    options = mergestruct(options, optionsInput);
end

% Construct options for fmincon.
fminopt = optimoptions('fmincon', 'UseParallel',options.UseParallel );
if options.verbose
    fminopt.Display='iter';
    fminopt.OutputFcn = @outfun;
end
fminopt.MaxIter = options.maxItr;
fminopt.TolX    = options.convGrad;
fminopt.TolFun  = options.convLL;

% Construct a mask to select only rates from connected states.
% FIXME: this excludes connections where ONE rate is zero...
nStates = model.nStates;
I = logical(eye(nStates));
rateMask = ~I & model.rates~=0 & ~model.fixRates;
nRates = sum(rateMask(:));

% Sanity checks
temp = vertcat(dwt{:});
assert( all(temp(:,1)<=model.nClasses), 'Idealization does not match current model!' );
assert( nRates>0, 'MIL requires model with at least two states' );

% Run fmincon optimizing with weak constraints to avoid negative rates.
% (fminunc works just as well; negative rates just cause harmless restarts).
optFun = @(x)milIter(dwt, dt, model, rateMask, x);
x0 = model.rates(rateMask)';
lb =     zeros(1,nRates);
% ub = 10/dt*ones(1,nRates);
[optParam,LL,exitflag] = fmincon( optFun, x0, [],[],[],[],lb,[],[],fminopt );

% Save results
optModel = copy(model);  %do not modify model in place
optModel.rates(rateMask) = optParam;



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
        wbh = waitbar(0,'Running MIL...');
          
    case 'iter'
%         % Keep track of parameter values at each iteration
%         X(itr+1,:)  = x;
%         dX(itr+1,:) = optimValues.gradient;
        
        if options.updateModel
            model.rates(rateMask) = x;
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
        
        % For debugging, plot how parameters change over iterations
%         X  = X(1:itr,:);
%         dX = dX(1:itr,:);
%         [idxin,idxout] = find(rateMask);  %state numbers for each transition type.
%         
%         figure;
%         for i=1:numel(x)
%             ax(1,i) = subplot(2,numel(x),i);
%             plot( 1:itr, X(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
%             if i==1
%                 ylabel('Param. Value');
%             end
%             title( sprintf('k%d->%d',idxin(i),idxout(i)) );
%             
%             ax(2,i) = subplot(2,numel(x),numel(x)+i);
%             plot( 1:itr, dX(:,i)', 'k.-', 'MarkerFaceColor',[1 0 1]);
%             if i==1
%                 xlabel('Iteration');
%                 ylabel('Gradient');
%             end
%         end
%         linkaxes(ax(:), 'x');
%         xlim(ax(1), [1,itr+1]);
end

end %function outfun



end %function mplOptimize

