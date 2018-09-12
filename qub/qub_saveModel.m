function qub_saveModel(model, filename)
% qub_loadModel  Loads a model from a .qmf format file from QuB
%     
%   qub_saveModel( MODEL, FILENAME ) saves the QubModel object MODEL to
%   a file with the path given in the string FILENAME (ending in .qmf).
%   The QUB_Tree object constructed to save to the .qmf file is saved in
%   the 'qubTree' property of the input model object. The 'filename'
%   property is also updated.
%
%   NOTE: this function is provided for compatibility with QuB, but it is
%   no longer recommended since not all features in the .qmf format are
%   supported.
%   
%   See QubModel, qub_loadModel.

%   Copyright 2007-2018 Cornell University All Rights Reserved.


narginchk(2,2);
model.verify();


if ~isempty(model.qubTree)
    % Use the existing QUB_Tree object if available
    outputTree = model.qubTree;
    try
        outputTree = rmfield(outputTree,'VRevs');
    catch
    end
else
    % If no template is available, construct a new QUB_Tree struct.
    s = struct();
    outputTree = struct('States',s,'Rates',s,'Constraints',s,'ConstraintsAmpVar',s, ...
              'ExtraKineticPars',s, 'ChannelCount',int32(1), 'Amps',0:9, 'Stds',repmat(0.06,1,10), ...
              'Conds',zeros(1,10), 'NAr',zeros(1,10), 'Ars',s, 'VRev',int32(0));

    outputTree.Properties = struct('ColorBack',16777215,'ColorLine',0, 'ColorRate',0, ...
        'ColorSelected',255, 'ColorFrame',16777215, 'ColorPanel',-2147483633, ...
        'StateSize',5, 'LineWidth',0, 'AlignToGrid',1, 'UseGlobalCond',0, ...
        'ScrollRates',1, 'DiagonalRates',0, 'ShowK1',0, 'EnforceConstraints',1, ...
        'MarginH',5, 'MarginV',5);

    outputTree = datanode(outputTree); %Add .data and .dataType fields
end


% Update FRET parametes from current model
outputTree.Amps.data(1:model.nClasses) = model.mu;
outputTree.Stds.data(1:model.nClasses) = model.sigma;

% Generate states and save initial probabilities.
% Class numbers are zero-based.
s = struct('x',num2cell(model.x), 'y',num2cell(model.y), ...
           'Class',num2cell(int32(model.class-1)), 'Pr',num2cell(model.p0), ...
           'Gr',num2cell(zeros(size(model.p0))) );
outputTree.States.State = datanode(s);

% Generate rate connections
rateTemplate = struct('States',0:1, 'k0',[0 0], 'k1',[0 0], 'dk0',[0 0], ...
     'dk1',[0 0], 'P',[0 0], 'Q',[0 0], 'PValue',struct(),'QValue',struct() );
rateTemplate.PNames.PName = {'Ligand','Ligand'};
rateTemplate.QNames.QName = {'Voltage','Voltage'};
rateTemplate.RateFormats.RateFormat = {'',''};
rateTemplate = datanode(rateTemplate);

outputTree.Rates = [];
conn = model.connections;
for i=1:size(conn,1)
    st = conn(i,:); %src and dst state numbers.
    outputTree.Rates.Rate(i) = rateTemplate;
    outputTree.Rates.Rate(i).States.data = int32(st-1); %zero-based
    outputTree.Rates.Rate(i).k0.data = ...
               [ model.rates(st(1),st(2)) model.rates(st(2),st(1)) ];
end

% Generate rate constraints
% FIXME: this does not seem to be compatible with QuB.
% if any(model.fixRates(:))
%     [src,dst] = find(model.fixRates);
%     pairs = int32([src dst]-1);
%     outputTree.Constraints = datanode( struct('FixRate',pairs) );
% end

% Save the resulting to QUB_Tree .qmf file
qub_saveTree(outputTree, filename, 'ModelFile');
model.qubTree = outputTree;


end







function node = datanode(input)
% Convert struct to a format acceptable to qub_saveTree().
% Converts leaf nodes to structs with .data and .dataType elements.

    if ischar(input)
        node.data = input;
        node.dataType = 3;
        
    elseif isfloat(input)
        node.data = input;
        node.dataType = 13;
        
    elseif isinteger(input) || islogical(input)
        node.data = input;
        node.dataType = 9;
        
    elseif iscell(input)
        for i=1:numel(input)
            node(i) = datanode(input{i});
        end
        
    elseif isstruct(input)
        node = struct();
        fn = fieldnames(input);
        for i=1:numel(input)
            for f=1:numel(fn)
                node(i).(fn{f}) = datanode(input(i).(fn{f}));
            end
        end
        
    else
        error('Invalid type');
    end
    
end


