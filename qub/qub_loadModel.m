function model = qub_loadModel(input)
% qub_loadModel  Loads a model file created by QuB
%     
%   MODEL = qub_loadModel( FILENAME ) loads the QuB-compatible .qmf file
%   specified in the string FILENAME, constructing a new QubModel modelect
%   MODEL.
%
%   NOTE: this function is provided for compatibility with QuB, but it is
%   no longer recommended since not all features in the .qmf format are
%   supported.
%
%   See also: QubModel, qub_saveModel.

%   Copyright 2007-2018 Cornell University All Rights Reserved.


narginchk(1,1);
nargoutchk(0,1);


% Load QUB_Tree modelect from .qmf file and convert to a struct
[~,~,e] = fileparts(input);
assert( strcmpi(e,'.qmf'),   'Invalid model file type. Should be .qmf');
assert( exist(input,'file')==2, 'File does not exist or is not accessible');

warning off qubtree:PointerFieldsNotSupported
warning off qubtree:MatrixFieldsNotSupported

model = QubModel();
treeModel = qub_loadTree( input );


% Construct model from QUB_Tree values
ns = numel(treeModel.States.State);
[model.p0,model.class,model.x,model.y] = deal( zeros(ns,1) );
states = treeModel.States.State;

for i=1:ns
    model.p0(i)    = states(i).Pr.data;
    model.class(i) = states(i).Class.data +1;
    model.x(i)     = states(i).x.data;
    model.y(i)     = states(i).y.data;
end

nc = max(model.class);
model.mu    = treeModel.Amps.data(1:nc);
model.sigma = treeModel.Stds.data(1:nc);
model.rates = zeros(ns);

for i=1:numel( treeModel.Rates.Rate )
    conn = treeModel.Rates.Rate(i);
    src  = conn.States.data(1)+1;
    dst  = conn.States.data(2)+1;

    model.rates(src,dst) = conn.k0.data(1);
    model.rates(dst,src) = conn.k0.data(2);
end


% Extract supported rate constraints. NOTE that many constraint
% types are not supported and silently ignored!
[model.fixMu,model.fixSigma] = deal( false(1,nc) );
model.fixRates = false(ns);

if isfield(treeModel,'Constraints') && isfield(treeModel.Constraints,'FixRate')
    for i=1:numel( treeModel.Constraints.FixRate )
        pair = treeModel.Constraints.FixRate(i).data+1;
        model.fixRates( pair(1), pair(2) ) = true;
    end    
end


end %function qub_loadModel












