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
%outputTree = qub_loadTree('default.qmf');
if isfield(model,'qubTree')
    outputTree = model.qubTree;
else
    outputTree = qub_loadTree('/home/dsterry/code/default_models/default.qmf');
end
if isfield( outputTree,'VRevs' ),
    outputTree = rmfield(outputTree,'VRevs');
end

% Update FRET parametes
nStates = size(model.rates,1);
nClass = numel(model.mu);
outputTree.Amps.data(1:nClass) = model.mu;
outputTree.Stds.data(1:nClass) = model.sigma;


% Generate states and save initial probabilities
s = outputTree.States.State;
for i=1:nStates,
    assert( s(i).Class.data+1 == model.class(i), 'Class-state mismatch' );
    s(i).Pr.data = model.p0(i);
end
outputTree.States.State = s;


% Generate rate connections
r = outputTree.Rates.Rate;
for i=1:numel(r),
    st = r(i).States.data+1;
    src = st(1);  dst = st(2);
    r(i).k0.data = [model.rates(dst,src) model.rates(src,dst)];
end
outputTree.Rates.Rate = r;


% Generate rate constraints
% if isfield(model,'fixRates')
%     outputTree.Constraints.FixRate = struct( ...
%                                 'data',[],'HasValue',[],'Value',[] );
%     f = struct([]);
%     
%     for i=1:nStates,
%         for j=1:nStates,
%             if i==j || i>j, continue; end
%             if ~model.fixRates(i,j), continue; end
% 
%             nPairs = nPairs+1;
%             f(nPairs).data = [i,j];
%             f(nPairs).HasValue.data = 0;
%             f(nPairs).Value = 0.0;
%         end
%     end
%     
%     outputTree.Constraints.FixRate = f;
% end


% Save the resulting to QUB_Tree .qmf file
if nargin>1,
    qub_saveTree(outputTree,modelFilename,'ModelFile');
end






