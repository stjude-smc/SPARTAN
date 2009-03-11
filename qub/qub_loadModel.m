function model = qub_loadModel(modelFilename)
% qub_loadModel  Loads a model file created by QuB
%     
%   [DMODEL] = qub_loadModel( filename )
%   Loads a qub model file (.qmf) -- specifically, this converts the
%   QUB_Tree object imbedded in a .qmf file into a structure, as
%   defined in qub_createModel.
%
%  http://www.qub.buffalo.edu


% MEX CODE
% if exist('qub_loadModel_MEX')~=3
%     error('Can''t find QuB compatibility layer files');
% end

if ~exist(modelFilename,'file')
    error('Model file doesn''t exist');
end

%[model.mu,model.sigma,model.p0s,model.rates,model.fixRates] = ...
%                                        qub_loadModel_MEX(modelFilename);


treeModel = qub_loadTree( modelFilename );

% Load initial probabilities
nStates = numel(treeModel.States.State);
mode.nStates = nStates;
model.p0 = zeros(1,nStates);

for i=1:nStates,
    model.p0(i) = treeModel.States.State(i).Pr.data;
end

if abs(sum(model.p0)-1)>0.02,
    warning('qub_loadModel:p0norm','Initial probabilities (p0) not normalized?');
end
    
% Load FRET parameters
model.mu = treeModel.Amps.data(1:nStates);
model.sigma = treeModel.Stds.data(1:nStates);

% Load rates into Q matrix
nRatePairs = numel( treeModel.Rates.Rate );
model.rates = zeros(nStates,nStates);

for i=1:nRatePairs,
    rate = treeModel.Rates.Rate(i);
    
    src = rate.States.data(1)+1;
    dst = rate.States.data(2)+1;
    
    model.rates(src,dst) = rate.k0.data(1);
    model.rates(dst,src) = rate.k0.data(2);
    
    %in principle, we could use dk0 for errors here
end

% NOTE: rate constraints not yet implemented!!!!
% Files that have constraints seem to crash the engine....
if isfield(treeModel,'Constraints') && isfield(treeModel.Constraints,'FixRate')
    nFixedRates = numel( treeModel.Constraints.FixRate );
    model.fixRates = zeros( size(model.rates) );
    
    for i=1:nFixedRates
        cons = treeModel.Constraints.FixRate(i);
        src = cons.data(1)+1;
        dst = cons.data(2)+1;
        
        model.fixRates(src,dst) = 1;
        
        if cons.HasValue.data~=0,
            warning('qub_loadModel:HasValue', ...
                       'Fixing rates to particular value not supported');
        end
    end    
end








