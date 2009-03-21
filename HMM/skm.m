function [dwt,model,LL] = skm( data, sampling, initialModel, params )
% SKM  Crude model re-estimation using iterative idealization
% 
%   [DWT,NEW_MODEL] = SKM( DATA, SAMPLING, MODEL, params )
%   Optimizes the given FRET/kinetic MODEL soas to maximize the 
%   likelihood data given the model. DATA is a NxM matrix of
%   N FRET traces of M datapoints in length each. MODEL is a typical
%   model specification, as defined in createModel.m.  SKM
%   returns the optimal model (NEW_MODEL) and the idealization
%   with maximum likelihood (DWT).
%
%   NOTE that constraints on FRET values and stdev specified in
%   the model file WILL be enforced in the fitting with SKM.
%   Constraints on kinetics will be ignored!
%   
%   The following params may be specified: FIXME...
%    - maxItr (100):  maximum number of iterations before terminating
%    - convLL (1e-2): stop iterating when LL converges within this limit
%    - 
%   

% DEPENDS: idealize, forward_viterbi, countEvents??

% NOTE: aggregated states not yet supported!!!!

% NOTE for now we use the standard proceedure:
% 1. optimize for whole file and idealize
% 2. truncate each trace to last instance of non-zero FRET
% 3. optimize a model for each trace individually,
%    fixing mean FRET and stdev values.
% 4. idealize each trace with its optimal model...


if nargin<3,
    error('SKM: not enough input arguments');
end

if nargin<4,
    params = struct([]);
end



%% Initialize algorithm

% Verify model correctness...
assert(qub_verifyModel(initialModel),'SKM: Invalid initial model');
nStates = numel(initialModel.mu);

if ~isfield(initialModel,'fixMu')
    initialModel.fixMu = zeros(nStates,1);
end
if ~isfield(initialModel,'fixSigma')
    initialModel.fixSigma = zeros(nStates,1);
end

% Parse model constraints...
if isfield(initialModel,'params') && any(initialModel.fixRates(:)),
    warning('SKM:fixRates','Fixing specific rates not support!');
end

% Parse convergence criteria
if ~isfield(params,'maxItr')
    params(1).maxItr = 100;
end

if ~isfield(params,'convLL')
    params.convLL = 1e-4;
end

if isfield(params,'convGrad')
    warning('SKM:convGrad','convGrad not yet implemented');
end

if ~isfield(params,'quiet'),
    params.quiet = 0;
end

if ~isfield(params,'convLL')
    params.convLL = 1e-4;
end

% Modify data to fall within specified ranges...
data(data>1) = 1;
data(data<-0.3) = -0.3;



%% Run the SKM algorithm

% Here we have several choices.
% If params.seperately = 
% NO:  Optimize all the data together and return a single model.
% YES: Optimize each trace individually, returning a model array
%      and a single idealization combining all results.

if isfield(params,'seperately') && params.seperately==1,
    [nTraces,nFrames] = size(data);
    dwt = cell(nTraces,1);
    LL  = zeros(nTraces,1);
    
    % Optimize each trace seperately, using the 
    for n=1:nTraces,
        [newDWT,model(n),newLL] = runSKM( data(n,:), ...
                                     sampling, initialModel, params );
        dwt{n} = newDWT{1};
        LL(n) = newLL(end);
    end
else
    % Optimize a single model and idealize all data using this model.
    [dwt,model,LL] = runSKM(data, sampling, initialModel, params);
end


end %function skm







%% SKM CORE METHOD
function [dwt,model,LL] = runSKM(data, sampling, initialModel, params)

[nTraces,nFrames] = size(data);
nStates = numel(initialModel.mu);

% Setup initial conditions
itr = 1; %number of iterations so far
LL = [];

model = initialModel;
mu = reshape(model.mu,nStates,1);
sigma = reshape(model.sigma,nStates,1);
A = model.rates*(sampling/1000);
A( logical(eye(nStates)) ) = 1-sum(A,2);
%A = model.A;
p0 = reshape(model.p0,nStates,1);

while( itr < params.maxItr ),

    % Idealize the data using the Viterbi algorithm
    % NOTE: this is the slowest step and should be optimized!
    [dwt,idl,offsets,vLL] = idealize( data, [mu sigma], p0, A );
    LL(itr) = sum(vLL)/nTraces;
    
    % Display intermediate progress...
    if ~params.quiet
        if itr==1,
            disp( sprintf('%d: %f', itr, LL(itr) ));
        else
            disp( sprintf('%d: %f (%f)', itr, LL(itr), LL(itr)-LL(itr-1) ));
            if (LL(end)-LL(end-1))<0,
                warning('SKM: LL is decreasing...');
            end
        end
        
        disp( [mu sigma p0 A] );
    end
    
    
    % Re-estimate FRET model parameters using idealization
    for state=1:nStates,
        edata = data(idl==state);
        if length(edata)<1, continue; end
        
        if ~model.fixMu(state)
            mu(state) = mean( edata );
        end
        if ~model.fixSigma(state)
            sigma(state) = std( edata );
            sigma(state) = max(sigma(state),0.01);
        end
    end
    
    
    % Re-estimate kinetic model parameters using idealization by:
    % counting number of each type of transition occuring in idl
    % and dividing by the total time spent in the source state.
    A = estimateA(idl,nStates);
    
    
    % Re-estimate initial probabilities.
    % Did the paper specify how this should be done?
    for state=1:nStates,
        % count number of times the first datapoint is in this state
        % and normalize to number of traces...
        p0(state) = sum( idl(:,1)==state )/nTraces;
    end
    p0 = p0/sum(p0);

    
    % Check for convergence
    if length(LL)>1 && abs(LL(end)-LL(end-1))<params.convLL,
        break;
    end
    
    itr = itr+1;

end %for each iteration...

% disp(vLL);

% Save final parameter values back into model
model.mu = mu;
model.sigma = sigma;
model.A = A;  %model.Q = ...
model.p0 = p0;

if ~params.quiet,
    disp( sprintf('SKM: Finished after %d iterations with LL=%f',itr,LL(end)) );
end



end %skm core method..


%%
function idlFinal = dwtToIdl( dwt, traceLen )
% NOTE how truncated data is handled?

nTraces = numel(dwt);
idlFinal = zeros(traceLen,nTraces);
    
for dwtID=1:nTraces,
    states = dwt{dwtID}(:,1);
    times  = double(dwt{dwtID}(:,2));
    assert(all(states>0));

    ends = [0; cumsum(times)];
    for j=1:numel(states),
        idlFinal( (ends(j)+1):ends(j+1),dwtID ) = states(j);
    end
end

idlFinal = idlFinal';

end


%%
function [A2] = estimateA( idl,nStates )
% Counts each transition type and saves in EVENTS,
% where EVENTS(i,j) is the number of transitions from state i to state j.
% NOTE: if idealizing very little data, some states may not be occupied.
% This will lead to a row in A that is all zeros and cause division by 0.
% TO avoid this, 

[nTraces,nFrames] = size(idl);

A = zeros( nStates, nStates );

for n=1:nTraces
    for t=2:nFrames,
        prev = idl(n,t-1);
        curr = idl(n,t);
        if curr==0, break; end  %end of sequence.

        A(prev,curr) = A(prev,curr)+1;
    end
end

normFact = repmat(sum(A,2),1,nStates);
normFact(normFact==0) = 1; %prevent division by zero

A2 = A./normFact;
assert(all(A2(:)>=0),'SKM: invalid A-matrix');

end
