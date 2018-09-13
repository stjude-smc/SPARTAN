function [idlTotal,model,LL] = skm( data, sampling, initialModel, params )
% SKM  Crude model re-estimation using iterative idealization
% 
%   [DWT,NEW_MODEL,LL,OFFSETS] = SKM( DATA, SAMPLING, MODEL, params )
%   Idealizes the FRET traces in the NxM matrix DATA, assigns the state of the
%   system at each point in time. The fitting process optimizes the initial
%   QubModel object MODEL as a starting point,  it using the segmental K-means
%   algorithm.
%
%   Optimizes the given FRET/kinetic MODEL soas to maximize the 
%   likelihood data given the model. DATA is a NxM matrix of
%   N FRET traces of M datapoints in length each. MODEL is a typical
%   model specification, as defined in createModel.m.  SKM
%   returns the optimal model (NEW_MODEL) and the idealization
%   with maximum likelihood (DWT). SAMPLING is in ms.
%
%   The following model parameters will affect the fitting procedure:
%    - fixMu:     Cx1 logical vector; fix class model mean FRET value.
%    - fixSigma:  Cx1 logical vector; fix class model stdev FRET value.
%    - fixRates:  SxS logical array -- NOT IMPLEMENTED!
%   
%   The following params may be specified:
%    - maxItr (100):       maximum number of iterations before terminating
%    - convLL (1e-2):      stop iterating when LL converges within this limit
%    - quiet  (false):     if true, avoid displaying console output or waitbar.
%    - fixRates (false):   if true, fix all rates at their initial values.
%    - zeroEnd  (false):   explicitly add a zero-state dwell at end of traces
%    - seperately (true):  if true, fit a model for each trace individually.
%
%   See also: QubModel, idealize, forward_viterbi.

%   Copyright 2007-2017 Cornell University All Rights Reserved.



%% Process input arguments

narginchk(3,4);
nargoutchk(0,4);


% Set default values for any paramaters not specified.
defaultParams.maxItr   = 100;
defaultParams.convLL   = 1e-4;
defaultParams.quiet    = false;
defaultParams.fixRates = false; %FIXME: should ultimately be in the model.
defaultParams.zeroEnd  = false;
defaultParams.seperately = true;
defaultParams.exclude  = false( size(data,1), 1);

params = mergestruct(defaultParams, params);


% Convert struct model input into a QubModel object
if isstruct(initialModel),
    initialModel = QubModel(initialModel);
elseif ~isa(initialModel,'QubModel')
    error('Invalid model input');
end


% Warnings for any features not yet implemented.
if any(initialModel.fixRates(:)),
    warning('SKM:fixRates','Fixing specific rates not supported!');
end

if isfield(params,'convGrad')
    warning('SKM:convGrad','convGrad not yet implemented.');
end


origSize = size(data);  %size before removing outliers
data = max(-1, min(10,data));    %remove outliers
data = data(~params.exclude,:);  %remove manually excluded traces
[nTraces,nFrames] = size(data);


%% Run the SKM algorithm
% Optimize each trace individually, returning a model array
if params.seperately,
    if ~params.quiet,
        wbh = parfor_progressbar(nTraces,'Idealizing traces separately...');
    else
        wbh = [];
    end
    
    % Use multi-process execution only for large datasets.
    if nTraces*nFrames > 1e4 && cascadeConstants('enable_parfor'),
        pool = gcp;
        M = pool.NumWorkers;
    else
        M = 0;  % Single-thread execution.
    end
    
    LL  = zeros(nTraces,1);
    idl = zeros(size(data));
    
    % Optimize each trace seperately.
    parfor (n=1:nTraces,M)
        [idl(n,:),model(n),newLL] = runSKM( data(n,:), ...
                                     sampling, initialModel, params );
        LL(n) = newLL(end);
        
        if ~isempty(wbh) && mod(n,10)==0,
            wbh.iterate(10);
        end
    end
    if ~isempty(wbh), close(wbh); end
    
% Optimize one model with all traces together.
else
    [idl,model,LL] = runSKM(data, sampling, initialModel, params);
end


idlTotal = zeros( origSize );
idlTotal( ~params.exclude, :) = idl;


% % Add dwell in zero-state until end of trace, if requested
if params.zeroEnd
    idlTotal(idlTotal==0)=1;
    idlTotal(params.exclude) = 0;
end


end %function skm







%% SKM CORE METHOD
function [idlClasses,model,LL] = runSKM(data, sampling, initialModel, params)

nTraces = size(data,1);
nStates = initialModel.nStates;
nClass  = initialModel.nClasses;


% Setup initial conditions
itr = 1; %number of iterations so far
LL = [];

model = copy(initialModel);
mu    = to_col(model.mu);
sigma = to_col(model.sigma);
p0    = to_col(model.p0);
classes = [0; to_col(model.class)];
A  = model.calcA(sampling/1000);  %transition probability matrix.

while( itr < params.maxItr ),

    % Idealize the data using the Viterbi algorithm (slow step).
    % idl is the state (NOT class) assignment at each point.
    imu    = mu( model.class );
    isigma = sigma( model.class );
    
    [idl,vLL] = idealize( data, [imu isigma], p0, A );
    LL(itr) = sum(vLL)/nTraces;
    
    
    % Display intermediate progress...
%     if ~params.quiet
%         if itr==1,
%             fprintf('%d: %f\n', itr, LL(itr) );
%         else
%             fprintf('%d: %f (%f)\n', itr, LL(itr), LL(itr)-LL(itr-1) );
%             if (LL(end)-LL(end-1))<0,
%                 disp('SKM Warning: LL is decreasing...');
%             end
%         end
%         
%         disp( [imu isigma p0 A] );
%     end
    
    
    % Re-estimate transition probability matrix (kinetic parameters).
    if params.fixRates<1,
        A = estimateA( idlToDwt(idl), nStates );
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
    
    
    % Re-estimate FRET model parameters using idealization
    for class=1:nClass,
        % Get all datapoints assigned to <class>
        edata = data( idlClasses==class );
        
        if numel(edata)<1, continue; end
        
        % Resestimate parameters as their expectation from observations.
        if ~model.fixMu(class)
            mu(class) = sum(edata)/numel(edata);
        end
        if ~model.fixSigma(class)
            sigma(class) = sqrt(var(edata));
            sigma(class) = max(sigma(class),0.01); %prevent convergace to 0
        end
    end
    
    
    % Check for convergence
    if numel(LL)>1 && abs(LL(end)-LL(end-1))<params.convLL,
        break;
    end
    
    itr = itr+1;

end %for each iteration...


% Save final parameter values back into model
model.mu = mu;
model.sigma = sigma;
model.p0 = p0;

% Convert transition probabilities to rates, setting any rates to zero
% that are not connected in the original model.
model.rates = logm(A) / (sampling/1000);
model.rates( initialModel.rates==0 ) = 0;

% if ~params.quiet,
%     fprintf('SKM: Finished after %d iterations with LL=%f\n',itr,LL(end));
% end



end %skm core method..



%%
function [A2] = estimateA( dwt,nStates )
%ESTIMATEA: Transition probability matrix from observed state sequence.
%
%   A = estimateA(DWT,nStates) estimates a transition probability matrix
%   (A) from observed dwell-times (DWT), where each element A(i,k) is the
%   probability at each frame of transition from state i to state j. This
%   is estimated simply as the number of transitions observed from i to j
%   divided by the total number of transitions out of i, including self
%   transitions (ie, the total number of frames in state i).

% Add a small constant for cases where some states are not occupied at all.
% The normalization would otherwise give a divide by zero error.
A = 0.01*eye(nStates);

for i=1:numel(dwt),
    states = dwt{i}(:,1);
    times  = dwt{i}(:,2);
    
    for d=1:numel(times),
        % Add self transitions (staying in the same state).
        cur = states(d);
        A(cur,cur) = A(cur,cur)+times(d);
        
        % Add transitions to another state (not self).
        if d<numel(times),
            next = states(d+1);
            A(cur,next) = A(cur,next)+1;
        end
    end
end

% Normalize to total time in each state.
normFact = repmat(sum(A,2),1,nStates); %sum over rows.
A2 = A./normFact;

assert( all(A2(:)>=0) & ~all(A2(:)==0), 'SKM: invalid A-matrix' );

end









