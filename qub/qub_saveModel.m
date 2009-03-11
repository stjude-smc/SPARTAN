function outputTree = qub_saveModel(model,modelFilename)
% qub_loadModel  Loads a model file created by QuB
%     
%   [DMODEL] = qub_saveModel( filename )
%   Saves a qub model file (.qmf)
%   See qub_createModel, qub_loadModel, qub_milOptimize, qub_skmIdealize
%
%  http://www.qub.buffalo.edu


% MEX CODE
% if ~exist(modelFilename,'file')
%     error('Model file doesn''t exist');
% end

if numel(model.mu)<2,
    error('Model must have at least 2 states');
end

% Verify model integrity
[status,msg] = qub_verifyModel(model);
if ~status,
    error( ['Model is not valid: ' msg] );
end


%% Load a simple default model as a starting point
outputTree = qub_loadTree('../default.qmf');

% Update FRET parametes
nStates = numel(model.mu);

outputTree.Amps.data(1:nStates) = model.mu;
outputTree.Stds.data(1:nStates) = model.sigma;


% Generate states and save initial probabilities
model.p0 = zeros(1,nStates);
outputTree.States.State = struct([]); %clear out existing stuff

clear s;
for i=1:nStates,
    s(i).y = 46.6;
    s(i).x = i*15;
    s(i).Gr = 0;
    s(i).Class = i;
    s(i).Pr = model.p0(i);
end
outputTree.States.State = s;


% Generate rate connections
ratePrototype = outputTree.Rates.Rate(1);
outputTree.Rates.Rate = struct([]);
nPairs = 0;
clear r;

for i=1:nStates,
    for j=1:nStates,
        if i==j || i>j, continue; end
        
        nPairs = nPairs+1;
        
        r(nPairs) = ratePrototype;
        r(nPairs).States.data = [i,j];
        r(nPairs).k0.data = [model.rates(i,j) model.rates(j,i)];
    end
end
outputTree.Rates.Rate = r;


% Generate rate constraints
clear f;

if isfield(model,'fixRates')
    outputTree.Constraints.FixRate = struct([]);
    
    for i=1:nStates,
        for j=1:nStates,
            if i==j || i>j, continue; end
            if ~model.fixRates(i,j), continue; end

            nPairs = nPairs+1;
            f(nPairs).data = [i,j];
            f(nPairs).HasValue.data = 0;
            f(nPairs).Value = 0.0;
        end
    end
    
    outputTree.Constraints.FixRate = f;
end


% Save the resulting to QUB_Tree .qmf file
if nargin>1,
    qub_saveTree(outputTree,modelFilename);
end






