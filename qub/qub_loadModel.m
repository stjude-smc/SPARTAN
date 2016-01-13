function model = qub_loadModel(modelFilenameInput)
% qub_loadModel  Loads a model file created by QuB
%     
%   [DMODEL] = qub_loadModel( filename )
%   Loads a qub model file (.qmf) -- specifically, this converts the
%   QUB_Tree object imbedded in a .qmf file into a structure, as
%   defined in qub_createModel.
%
%  http://www.qub.buffalo.edu

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Disable these pointless warnings in qubtree/treestruct.cpp.
warning off qubtree:PointerFieldsNotSupported
warning off qubtree:MatrixFieldsNotSupported



persistent modelFilename;

% If no model filename is given, prompt the user for it.
% The selection will be remembered because modelFilename is persistent.
if nargin<1,
    [f,p] = uigetfile('*.qmf','Select a model file...',modelFilename);
    if f==0,
        model = [];
        return;
    end
    modelFilename = fullfile(p,f);
else
    modelFilename = modelFilenameInput;
end

if isstruct(modelFilename),
    disp('qub_loadModel: Recieved struct as input, assuming it is a ModelFile tree');
    treeModel = modelFilename;
    %verify model here
else
    if ~exist(modelFilename,'file')
        error('Model file doesn''t exist');
    end

    % Load QuBTree object saved to disk representing a model.
    treeModel = qub_loadTree( modelFilename );
end
    

% Load initial probabilities
nStates = numel(treeModel.States.State);
model.nStates = nStates;
model.p0 = zeros(nStates,1);
model.class = zeros(nStates,1);

for i=1:nStates,
    model.p0(i) = treeModel.States.State(i).Pr.data;
    model.class(i) = treeModel.States.State(i).Class.data +1;
end

if abs(sum(model.p0)-1)>0.02,
    warning('qub_loadModel:p0norm','Initial probabilities (p0) not normalized?');
end
    
% Load FRET parameters
nClass = max(model.class);
model.nClasses = nClass;
model.mu = treeModel.Amps.data(1:nClass);
model.sigma = treeModel.Stds.data(1:nClass);

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
        
%         if cons.HasValue.data~=0,
%             warning('qub_loadModel:HasValue', ...
%                        'Fixing rates to particular value not supported');
%         end
    end    
end


model.qubTree = treeModel;
if exist('modelFilename','var')
    model.filename = modelFilename;
end





