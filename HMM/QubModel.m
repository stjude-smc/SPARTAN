classdef QubModel < matlab.mixin.Copyable
%   FRET observation and kinetic model for HMM analysis and optimization.
%
%   QubModel(N)     creates a model with N states/classes.
%   QubModel(FILE)  loads a model from file (.model or .qmf)
%   QubModel(MODEL) creates a copy of a QubModel object.
%
%   See also: qub_loadModel, qub_saveModel, QubModelViewer.

%   Copyright 2007-2018 Cornell University All Rights Reserved.

 
properties (SetAccess=public, GetAccess=public, SetObservable)
    % Model parameters
    class;     %class number for each state (1-based)
    p0;        %initial probabilities
    mu;        %mean FRET value
    sigma;     %stdev of FRET (noise)
    rates;     %rate constants (i,j) is i->j
    
    % If true, do not optimize the specified parameter.
    % These logical arrays are the same size as rates, mu, sigma arrays.
    fixRates;
    fixMu;
    fixSigma;
    
    % State locations for display (0..100; see QubModelViewer)
    x = [];
    y = [];
    
    % Structure containing the .qmf format tree of all model information.
    % This includes many parameters we don't use but QuB expects.
    % This is only updated when a .qmf file is saved to disk.
    qubTree;
    
    % If true, prevents events from being triggered.
    % Use this to avoid sending events while updating multiple parameters.
    muteListeners = false;
end

properties (GetAccess=public, SetAccess=protected, Transient)
    % Full path and name of the model file that was loaded.
    filename;
end

% Model properties derived from model parameters (above).
properties (SetAccess=immutable, GetAccess=public, Dependent)
    nStates;
    nClasses;
    connections;
end

properties (SetAccess=protected, GetAccess=protected, Transient, Hidden)
    % UpdateModel event listener
    % FIXME: loading (incl. w/ parfor) will need to recreate listener!
    updateListener;
    rateUpdateListener;
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
        
        if nargin<1, return; end  %Empty object constructor
        
        % Copy another QubModel object.
        if nargin==1 && isa(input,'QubModel')
            obj = copy(input);
            obj.verify();
            return;
            
        % Create a new model only specifying the number of states
        elseif isscalar(input)
            obj.initModel(input);
            
        % Load from file
        elseif ischar(input)
            [~,~,e] = fileparts(input);
            if strcmpi(e,'.qmf')
                obj = qub_loadModel(input);
            elseif strcmpi(e,'.model') || strcmpi(e,'.mat')
                temp = load(input,'-mat');
                if temp.version>=2, error('Model file version format is not supported'); end
                obj = temp.model;
            else
                error('Invalid model file format');
            end
            obj.filename = input;
            
        else
            error('Unexpected input for QubModel constructor');
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
    
    function initModel(obj,ns)
        % Initialize a valid model given only the number of states
        obj.class = (1:ns)';
        obj.mu    = (0:1/ns:1-(1/ns))';
        obj.sigma = ones(ns,1) * 0.06;
        obj.p0    = ones(ns,1) / ns;
        obj.rates = ones(ns);
        obj.rates( eye(ns)==1 ) = 0;

        obj.fixMu    = false(ns,1);
        obj.fixSigma = false(ns,1);
        obj.fixRates = false(ns);

        [obj.x, obj.y] = deal( linspace(11,86,ns)' );
    end
    
    
    function save( model, filename )
        % Save the model to file
        [~,~,e] = fileparts(filename);
        if strcmpi(e,'.qmf')
            qub_saveModel(model, filename);
        elseif strcmpi(e,'.model') || strcmpi(e,'.mat')
            version = 1.0; %#ok<NASGU>
            save(filename, 'model','version');
        else
            error('Unrecognized model output format');
        end
        model.filename = filename;
    end
    
    
    function revert(model,varargin)
    % Revert to the state of the file before any modifications.
    
        if isempty(model.filename)
            warning('Cannot revert temporary models');
            return;
        end
    
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




