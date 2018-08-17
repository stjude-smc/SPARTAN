classdef QubModel < matlab.mixin.Copyable
% HMM model object designed for representing a .qmf model file generated in
% QuB. Even if QuB isn't involved, we still use this format for all model
% definitions. The filename cannot be changed once the object is created;
% a new object should be created to load another file. The exception is
% when saving a file -- the model is now tied to the newly save file.
%
%   QubModel(N)     creates a model with N states/classes.
%   QubModel(FILE)  loads the given .qmf file.
%   QubModel(MODEL) creates a copy of a QubModel object or similar struct.
%
%   See also: qub_loadTree, qub_saveTree.

%   Copyright 2007-2016 Cornell University All Rights Reserved.

% TODO: support basic constraint types from .qmf format.

 
properties (SetAccess=public, GetAccess=public, SetObservable)
    % Model parameters
    class;     %class number for each state (1-based)
    p0;        %initial probabilities
    mu;        %mean FRET value
    sigma;     %stdev of FRET (noise)
    rates;     %rate constants (i,j) is i->j
    
    % Fitting constraints.
    % If 1, value is not re-estimated/changed during model optimization.
    % QuB supports many types of constraints and none are really supported.
    % This variables are only (currently) set in MATLAB code. FIXME
    fixRates;
    fixMu;
    fixSigma;
    
    %Internal display data (must be public for QubModelViewer)
    x = [];
    y = [];
end

% Model properties derived from model parameters (above).
properties (SetAccess=immutable, GetAccess=public, Dependent)
    nStates;
    nClasses;
    connections;
end

properties (SetAccess=protected, GetAccess=public)
    % Full path and name of the model file that was loaded. This cannot be
    % changed by the user.
    filename = [];
    
    % Structure containing the .qmf format tree of all model information.
    % This includes many parameters we don't use but QuB expects.
    qubTree;
end

properties (SetAccess=protected, GetAccess=protected, Transient, Hidden)
    % UpdateModel event listener
    % FIXME: loading (incl. w/ parfor) will need to recreate listener!
    updateListener;
    rateUpdateListener;
    muteListeners = false;
end

events
    % Event triggered after model parameters have been changed, and the model
    % is verified to be valid.
    UpdateModel;
    UpdateRates;
end


methods
    %%%%%%%%%%%%%%%%%%%  CONSTRUCTOR & SERIALIZATION  %%%%%%%%%%%%%%%%%%%%
    
    function obj = QubModel( input )
    % Create a new QubModel object. Inputs can be a .qmf file, a struct with the
    % same fields as a QubModel, or a QubModel object to copy.
        
        % parfor seems to create empty objects it doesn't use. Silently ignore.
        if nargin<1, return;  end
        
        % Copy another QubModel object
        if nargin==1 && isa(input,'QubModel')
            obj = copy(input);
            obj.verify();
            return;
        
        % Load from a model struct (see qub_loadModel.m)
        elseif isstruct(input),
            m = input;
            
        % Load from file
        elseif ischar(input),
            m = qub_loadModel( input );
            obj.filename = input;
            
        % New model with a given number of states
        elseif isscalar(input) && numel(input)==1,
            m = qub_createModel(input);
            
        else
            error('Unexpected input for QubModel constructor');
        end
        
        % Set default values for optional fields
        N = length(m.rates);
        obj.fixRates = false( size(m.rates) );
        obj.fixMu    = false( size(m.mu)    );
        obj.fixSigma = false( size(m.sigma) );
        if ~isfield(m,'class'),
            obj.class = (1:length(m.rates))';  %assumes no degenerate states!
        end
        
        % Copy all relevant fields from input to object properties.
        mco   = ?QubModel;
        props = {mco.PropertyList.Name};
        
        for i=1:numel(props),
            if isfield(m,props{i}) && ~mco.PropertyList(i).Dependent,
                obj.(props{i}) = m.(props{i});
            end
        end
        
        % Extract model display settings from qubTree.
        if isempty(obj.x) || isempty(obj.y),
            if ~isempty(obj.qubTree)
                obj.x = zeros( obj.nStates,1 );
                obj.y = zeros( obj.nStates,1 );
                for i=1:obj.nStates,
                    obj.x(i) = obj.qubTree.States.State(i).x.data;
                    obj.y(i) = obj.qubTree.States.State(i).y.data;
                end
            else
                % Default locations, scaled below.
                obj.x = (1:N)';
                obj.y = obj.class;
            end

            % Rescale coordinates to properly fit in the box.
            obj.x = obj.x-min(obj.x);
            obj.y = obj.y-min(obj.y);
            obj.x = 75*obj.x/max(obj.x) +10;
            obj.y = 75*obj.y/max(obj.y) +12;
        end
        
        % Verify the model parameters make sense.
        obj.verify();
        obj.updateListener = addlistener(obj, {'mu','fixMu','sigma','fixSigma','rates'}, ...
                                        'PostSet',@obj.UpdateModel_Callback);
        obj.rateUpdateListener = addlistener(obj, {'rates'}, 'PostSet',@obj.UpdateRates_Callback);
    end
    
    function UpdateModel_Callback(obj,varargin)
        if ~obj.muteListeners
            notify(obj,'UpdateModel');
        end
    end
    function UpdateRates_Callback(obj,varargin)
        if ~obj.muteListeners
            notify(obj,'UpdateRates');
        end
    end
    
    
    %% ----------------------   SERIALIZATION   ---------------------- %%
    function save( model, fname )
        % Save the model to file. The code is essentially the same as
        % qub_saveModel, but that code can't be used directly since it
        % assumes the model is a struct, not an class object.
        narginchk(2,2);
        
        % If no filename is given, update the original model file.
        if nargin>=2,
            model.filename = fname;
        end
        model.verify();

        % Use loaded QuB tree or create a new one with default values.
        % NOTE: QuB expects index variables to have integer type!
        if ~isempty(model.qubTree)
            outputTree = model.qubTree;
            
            if isfield(outputTree,'VRevs'),
                outputTree = rmfield(outputTree,'VRevs');
            end
        else
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
        nClass = numel(model.mu);
        outputTree.Amps.data(1:nClass) = model.mu;
        outputTree.Stds.data(1:nClass) = model.sigma;
        
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
        for i=1:size(conn,1),
            st = conn(i,:); %src and dst state numbers.
            outputTree.Rates.Rate(i) = rateTemplate;
            outputTree.Rates.Rate(i).States.data = int32(st-1); %zero-based
            outputTree.Rates.Rate(i).k0.data = ...
                       [ model.rates(st(1),st(2)) model.rates(st(2),st(1)) ];
        end

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
        qub_saveTree(outputTree, model.filename, 'ModelFile');
        model.qubTree = outputTree;
    end
    
    
    function revert(model,varargin)
    % Revert to the state of the file before any modifications.
    
        % Wait until model is finalized to notify listeners
        model.muteListeners = true;
        
        newmodel = QubModel(model.filename);
        mco   = ?QubModel;
        props = {mco.PropertyList.Name};
        props = props(~[mco.PropertyList.Dependent] & ~[mco.PropertyList.Constant] & ~[mco.PropertyList.Hidden]);

        for i=1:numel(props),
            model.(props{i}) = newmodel.(props{i});
        end

        model.muteListeners = false;
        notify(model,'UpdateModel');  %inform listeners model has changed.
        notify(model,'UpdateRates');
    end
    
    
    
        
    %% ---------------------   GET/SET METHODS   --------------------- %%
    function n = get.nStates( model )
        n = numel(model.class);
    end
    
    function n = get.nClasses( model )
        n = numel(model.mu);
    end
    
    function C = get.connections(model)
        % Matrix with each row being state pairs with non-zero rates
        [row,col] = find(model.rates>0);
        C = unique( sort([row col],2), 'rows' );
    end
    
    function tf = isfield(model, fname)
        tf = ismember(fname, properties(model));
    end
    
    function tf = isempty(this)
        tf = isempty(this.class);
    end
    
    
    function addState(model, newClass, newP0, newX,newY)
    % Add a new state to the model
        if nargin<2, newClass=1; end
        if nargin<3,
            newP0 = double(model.nStates==0);  %set to 1.0 if empty model
        end
        if nargin<5, newX=50; newY=50; end
        
        enableListener(model.updateListener, false);
        N = model.nStates;
        
        % Add a state with default settings
        model.class(N+1) = newClass;
        model.p0(N+1) = newP0;
        model.rates(N+1,N+1) = 0;
        model.fixRates(N+1,N+1) = false;
        model.x(N+1) = newX;
        model.y(N+1) = newY;
        
        % If this is a new class, use reasonable defaults.
        if newClass>numel(model.mu)
            model.addClass(newClass);
        else
            notify(model,'UpdateModel');  %inform listeners model has changed.
        end
        model.verify();
        enableListener(model.updateListener, true);
    end
    
    function addClass(model, newClass, newMu, newSigma)
        % Update properties to include a new class, using reasonable defaults.
        % FIXME: what if newClass is not contiguous with existing states?
        % FIXME: what about empty models?
        narginchk(2,4);
        
        model.muteListeners = true;
        
        if nargin<3,
            newMu = max(model.mu) + 0.1;
        end
        if nargin<4,
            newSigma = model.sigma(end);
        end
        
        model.mu(newClass) = newMu;
        model.sigma(newClass) = newSigma;
        model.fixMu(newClass) = false;
        model.fixSigma(newClass) = false;
        
        model.verify();
        model.muteListeners = false;
        notify(model,'UpdateModel');  %inform listeners model has changed.
    end
    
    function removeState(model,id)
    % Add a new state to the model
        model.muteListeners = true;
        
        % Remove state from variables
        fields = {'class','p0','x','y'};
        for i=1:numel(fields)
            if ~isempty(model.(fields{i})),
                model.(fields{i})(id) = [];
            end
        end
        model.p0 = model.p0/sum(model.p0);  %re-normalize
        
        model.rates(id,:) = [];
        model.rates(:,id) = [];
        model.fixRates(id,:) = [];
        model.fixRates(:,id) = [];
        
        model.verify();
        model.muteListeners = false;
        notify(model,'UpdateModel');  %inform listeners model has changed.
    end
    
    % Verify model is self-consistent and valid (see qub_verifyModel).
    % isValid is true if ok, or false if there is a problem. str is an error
    % message. models can have non-fatal problems (b=1, but a str is given).
    function [isValid,str] = verify( model, throwErrors )
        str = [];
        
        % Ensure all publicly-accessible fields are set.
        % An instance without a qubTree is valid, but cannot be displayed/saved.
        mco   = ?QubModel;
        props = {mco.PropertyList.Name};
        sa    = {mco.PropertyList.SetAccess};
        
        if any( cellfun(@isempty,props) & strcmpi(sa,'public') ),
            str = 'Some properties are empty or not defined!';
        end
        
        if model.nClasses<2,
            % Is this still a problem??
            str = 'Model must have at least 2 states';
        end
        
        c = cellfun(@numel,{model.mu,model.sigma,model.fixMu,model.fixSigma});
        if ~(  numel(unique(c))==1 && numel(model.p0)==numel(model.class) && ...
               all(size(model.rates)==size(model.fixRates))  && ...
               all(numel(model.p0)==size(model.rates))  ),
            str = 'Model parameter size mismatch';
        end

        r = model.rates( ~logical(eye(model.nStates)) );  %ignore diagonal elements.
        if any( r<0 | isinf(r) ),
            str = 'Negative or infinite rates not allowed.';
        end
        
        isValid = isempty(str);
        
        % These checks only produce warnings
        if ~isempty(model.p0) && abs(sum(model.p0)-1) > 0.01
            str = 'p0 values not normalized';
        end        
        
        % If throwErrors is set, go ahead and produce the error or warnings
        % instead of letting the caller decide.
        if nargin==1 || throwErrors,
            if ~isValid, error(['Invalid model: ' str]); end
            if isValid && ~isempty(str),  warning(str);  end
        end
        
        % Force correct orientation for all properties.
        fields = {'p0','mu','sigma','fixMu','fixSigma'};
        for i=1:numel(fields),
            model.(fields{i}) = to_col(model.(fields{i}));
        end
    end
    
    % Matrix of rate constants normalized so that rows sum to zero.
    function Q = calcQ(model)
        Q = model.rates;
        I = logical(eye(size(Q)));
        Q(I) = 0;
        Q(I) = -sum( Q,2 );
    end
    
    % Probability of transitioning within the time step dt (in seconds).
    % Same size and indexing as the rate matrix Q.
    function A = calcA( model, dt )
        A = expm( model.calcQ().*dt );
    end
    
    % Eqiulibrium state probabilities (p0) assuming the system is ergodic.
    function p = calcEquilibriumP0( model )
        
        if any( sum(model.rates)==0 & sum(model.rates,2)'==0 ),
            % If there are isolated parts of the model, how do we know which 
            % one the system starts in? A more sophisticated test is needed 
            % to find isolated segments of a model (multiple states).
            warning('Isolated states may produce unexpected p0 estimates.');
        end
        
        % See Single-Channel Recording (1997), ISBN 978-1-4419-1230-5,
        % pg. 597 (eq. 17, section 3.2.2, chapter 20).
        U = ones(1, model.nStates);
        S = [ model.calcQ() U' ];
        p = U * (S * S')^-1;
    end
    
    
end %public methods

end  %classdef





function node = datanode(input)
% Convert struct to a format acceptable to qub_saveTree().
% Converts leaf nodes to structs with .data and .dataType elements.

    if ischar(input),
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
        for i=1:numel(input),
            for f=1:numel(fn)
                node(i).(fn{f}) = datanode(input(i).(fn{f}));
            end
        end
        
    else
        error('Invalid type');
    end
    
end





