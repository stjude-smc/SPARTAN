function [dwt,model,LL,offsets] = skm( data, sampling, initialModel, params )
% SKM  Crude model re-estimation using iterative idealization
% 
%   [DWT,NEW_MODEL,LL,OFFSETS] = SKM( DATA, SAMPLING, MODEL, params )
%   Optimizes the given FRET/kinetic MODEL soas to maximize the 
%   likelihood data given the model. DATA is a NxM matrix of
%   N FRET traces of M datapoints in length each. MODEL is a typical
%   model specification, as defined in createModel.m.  SKM
%   returns the optimal model (NEW_MODEL) and the idealization
%   with maximum likelihood (DWT). SAMPLING is in ms.
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

% DEPENDS: idealize, forward_viterbi


if nargin<3,
    error('SKM: not enough input arguments');
end

if nargin<4,
    params = struct([]);
end



%% Initialize algorithm

% Verify model correctness...
[goodModel,msg] = qub_verifyModel(initialModel);
assert(goodModel,['SKM: Invalid initial model: ' msg]);
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

if ~isfield(params,'fixRates')
    params.fixRates = 0;
end

% Modify data to fall within specified ranges...
data(data>10) = 10;
data(data<-1) = -1;


%% Run the SKM algorithm

[nTraces,nFrames] = size(data);

% Here we have several choices.
% If params.seperately = 
% NO:  Optimize all the data together and return a single model.
% YES: Optimize each trace individually, returning a model array
%      and a single idealization combining all results.

if isfield(params,'seperately') && params.seperately==1,
    % Start the matlab thread pool if not already running. perfor below will
    % run the calculations of the available processors. The speed up is not
    % that significant (2.5-fold for 4 cores), but useful.
    if matlabpool('size')==0,  matlabpool;  end
    
    dwt = cell(nTraces,1);
    LL  = zeros(nTraces,1);
    
    % Optimize each trace seperately, using the 
    wbh = waitbar(0,'Idealizing traces seperately,..');
    
    parfor n=1:nTraces,
    %for n=1:nTraces,   %use this line instead of single-thread.
        [newDWT,model(n),newLL] = runSKM( data(n,:), ...
                                     sampling, initialModel, params );
        dwt{n} = newDWT{1};
        LL(n) = newLL(end);
        
        %if mod(n,25)==0,
        %    waitbar(n/nTraces,wbh);
        %end
    end
    close(wbh);
    
    offsets = nFrames*((1:nTraces)-1);
else
    % Optimize a single model and idealize all data using this model.
    [dwt,model,LL,offsets] = runSKM(data, sampling, initialModel, params);
end


% Add dwell in zero-state until end of trace, if requested
if isfield(params,'zeroEnd') && params.zeroEnd==1,
    for i=1:nTraces,
        states = dwt{i}(:,1);
        times  = dwt{i}(:,2);
        if numel(states)<1, continue; end
        
        remainder = nFrames-sum(times)-1;
        if remainder<=0, continue; end
        
        if states(end)==1
            times(end) = times(end)+remainder;
        else
            states = [states ; 1];
            times  = [times ; remainder];
        end
        
        dwt{i} = [states times];
    end
end


end %function skm







%% SKM CORE METHOD
function [dwt,model,LL,offsets] = runSKM(data, sampling, initialModel, params)

[nTraces,nFrames] = size(data);
nStates = size(initialModel.rates,1);
nClass  = numel(initialModel.mu);


% Setup initial conditions
itr = 1; %number of iterations so far
LL = [];

model = initialModel;
mu = reshape(model.mu,nClass,1);
sigma = reshape(model.sigma,nClass,1);
p0 = reshape(model.p0,nStates,1);
classes = [0; model.class];


% Calculate correct Q matrix from rates...
Q = model.rates;
Q( logical(eye(size(Q))) ) = 0;
Q( logical(eye(size(Q))) ) = -sum( Q,2 );
dt = (sampling/1000); %convert from ms to seconds.
A = expm( Q.*dt );

while( itr < params.maxItr ),

    % Idealize the data using the Viterbi algorithm (slow step)
    imu    = mu( model.class );
    isigma = sigma( model.class );
    
    [dwt,idl,offsets,vLL] = idealize( data, [imu isigma], p0, A );
    LL(itr) = sum(vLL)/nTraces;
    
    
    % Display intermediate progress...
    if ~params.quiet
        if itr==1,
            disp( sprintf('%d: %f', itr, LL(itr) ));
        else
            disp( sprintf('%d: %f (%f)', itr, LL(itr), LL(itr)-LL(itr-1) ));
            if (LL(end)-LL(end-1))<0,
                disp('SKM Warning: LL is decreasing...');
            end
        end
        
        disp( [imu isigma p0 A] );
    end
    
    
    % Re-estimate kinetic model parameters using idealization by: 
    % counting number of each type of transition occuring in idl
    % and dividing by the total time spent in the source state.
    % (this is a slow step - should be optimized).
    if params.fixRates<1,
        A = estimateA(idl,nStates);
    end
    
    
    % Re-estimate initial probabilities.
    % Did the paper specify how this should be done?
    for state=1:nStates,
        % count number of times the first datapoint is in this state
        % and normalize to number of traces...
        p0(state) = sum( idl(:,1)==state )/nTraces;
    end
    p0 = p0/sum(p0);
    
    
    % Convert state-list idealization to class-list.
    idlClasses = classes( idl+1 );
    idlClasses = reshape( idlClasses, size(idl) );
    for i=1:size(idl,1)
        trace = idlClasses(i,:);
        if all(trace<=0), continue; end  %nothing here...
        dwt{i} = RLEncode(  trace(1:find(trace>0,1,'last'))  );
    end
    
    
    % Re-estimate FRET model parameters using idealization
    for class=1:nClass,
        % Get all datapoints assigned to <class>
        edata = data( idlClasses==class );
        
        if length(edata)<1, continue; end
        
        % Resestimate parameters as their expectation from observations.
        if ~model.fixMu(class)
            mu(class) = mean( edata );
        end
        if ~model.fixSigma(class)
            sigma(class) = std( edata );
            sigma(class) = max(sigma(class),0.01); %prevent convergace to 0
        end
    end
    
    
    % Check for convergence
    if length(LL)>1 && abs(LL(end)-LL(end-1))<params.convLL,
        break;
    end
    
    itr = itr+1;

end %for each iteration...

% Convert state list into class list for idealization
idl = classes( idl+1 );
reshape( idlClasses, size(idl) );


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
function [A2] = estimateA( idl,nStates )
% Counts each transition type and saves in EVENTS,
% where EVENTS(i,j) is the number of transitions from state i to state j.
% NOTE: if idealizing very little data, some states may not be occupied.
% This will lead to a row in A that is all zeros and cause division by 0.
% TO avoid this, use a small value on the diagonal instead of zero, which
% will get normalized to 1 (can transition to itself) upon normalization.
%

[nTraces,nFrames] = size(idl);

A = 0.01*eye( nStates );

for n=1:nTraces
    
    for t=2:nFrames,
        prev = idl(n,t-1);
        curr = idl(n,t);
        if curr==0, break; end  %end of sequence.

        A(prev,curr) = A(prev,curr)+1;
    end
end

normFact = repmat(sum(A,2),1,nStates); %sum over rows.

A2 = A./normFact;
assert( all(A2(:)>=0) & ~all(A2(:)==0), 'SKM: invalid A-matrix' );

end
