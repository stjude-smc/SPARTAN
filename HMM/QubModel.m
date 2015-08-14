classdef QubModel < handle
% HMM model object designed for representing a .qmf model file generated in
% QuB. Even if QuB isn't involved, we still use this format for all model
% definitions. The filename cannot be changed once the object is created;
% a new object should be created to load another file. The exception is
% when saving a file -- the model is now tied to the newly save file.
%

% TODO: support basic constraint types from .qmf format.

 
properties (SetAccess=public, GetAccess=public)
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
end


% Model properties derived from model parameters (above).
properties (SetAccess=private, GetAccess=public, Dependent)
    nStates;
    nClasses;
end % Dependent properties

properties (SetAccess=protected, GetAccess=public)
    % Full path and name of the model file that was loaded. This cannot be
    % changed by the user.
    filename = [];
    
    % Display data
    x = [];
    y = [];
    
    % Structure containing the .qmf format tree of all model information.
    % This includes many parameters we don't use but QuB expects.
    qubTree;
end



methods
    %%%%%%%%%%%%%%%%%%%  CONSTRUCTOR & SERIALIZATION  %%%%%%%%%%%%%%%%%%%%
    
    function obj = QubModel( fname )
        % Create a model by loading a qmf (QuB model) file.        
        assert( nargin>=1, 'Empty models are not allowed' );
        m = qub_loadModel( fname );
        
        obj.filename = fname;
        obj.qubTree  = m.qubTree;
        obj.p0       = m.p0;
        obj.class    = m.class;
        obj.mu       = m.mu;
        obj.sigma    = m.sigma;
        obj.rates    = m.rates;
        
        % Set default values for constraint matrices if not given.
        if isfield(m,'fixRates'),
            obj.fixRates = m.fixRates;
        else
            obj.fixRates = zeros( size(m.rates) );
        end
        
        if isfield(m,'fixMu'),
            obj.fixMu = m.fixMu;
        else
            obj.fixMu = zeros( size(m.mu) );
        end
        
        if isfield(m,'fixSigma'),
            obj.fixSigma = m.fixSigma;
        else
            obj.fixSigma = zeros( size(m.sigma) );
        end
        
        % Get display settings from the tree.        
        obj.x = zeros( obj.nStates,1 );
        obj.y = zeros( obj.nStates,1 );
        for i=1:obj.nStates,
            obj.x(i) = obj.qubTree.States.State(i).x.data;
            obj.y(i) = obj.qubTree.States.State(i).y.data;
        end

        % Rescale coordinates to properly fit in the box.
        obj.x = obj.x-min(obj.x);
        obj.y = obj.y-min(obj.y);
        obj.x = 75*obj.x/max(obj.x) +10;
        obj.y = 75*obj.y/max(obj.y) +12;
        
        % Verify the model parameters make sense.
        obj.verify();
    end
    
    
    function save( model, fname )
        % Save the model to file. The code is essentially the same as
        % qub_saveModel, but that code can't be used directly since it
        % assumes the model is a struct, not an class object.
        
        % If no filename is given, update the original model file.
        if nargin<2,
            fname = model.filename;
        end        
        
        % Verify model integrity before saving to prevent later errors.
        model.verify();

        % The tree has many elements we don't use, but QuB needs, so either
        % update the original tree for this model file or create a new one
        % from a template model (default.qmf). Not guaranteed to work
        % correctly if the template has a different number of states?
        if ~isempty(model.qubTree)
            outputTree = model.qubTree;
        else
            outputTree = qub_loadTree('default.qmf');
        end
        
        if isfield(outputTree,'VRevs'),
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
            s(i).x.data = model.x(i);
            s(i).y.data = model.y(i);       
        end
        outputTree.States.State = s;


        % Generate rate connections
        r = outputTree.Rates.Rate;
        for i=1:numel(r),
            st = r(i).States.data+1;
            r(i).k0.data = [model.rates(st(1),st(2)) model.rates(st(2),st(1))];
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
            qub_saveTree(outputTree,fname,'ModelFile');
        end
        
        % This file now really represents the new file and not the old one.
        model.filename = fname;
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%  GET/SET METHODS  %%%%%%%%%%%%%%%%%%%%%%%%%%

    function n = get.nStates( model )
        n = numel(model.class);
    end
    
    function n = get.nClasses( model )
        n = numel(model.mu);
    end
        
    % Verify model is self-consistent and valid (see qub_verifyModel).
    % b is true if ok, or false if there is a problem. str is an error
    % message. models can have non-fatal problems (b=1, but a str is given).
    function [isValid,str] = verify( model, throwErrors )
        str = [];
        
        if any( cellfun(@isempty,properties(model)) ),
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

        if any( model.rates<0 | isinf(model.rates) ),
            str = 'Negative or infinite rates not allowed.';
        end
        
        isValid = isempty(str);
        
        % These checks only produce are only warnings
        if abs(sum(model.p0)-1) > 0.01
            str = 'p0 values not normalized';
        end        
        
        % If throwErrors is set, go ahead and produce the error or warnings
        % instead of letting the caller decide.
        if nargin==1 || throwErrors,
            if ~isValid, error(['Invalid model: ' str]); end
            if isValid && ~isempty(str),  warning(str);  end
        end
    end
    
    % Calculate the transition probability matrix (A), which is a discrete
    % version of the rate matrix (Q). dt is the timestep in seconds.
    function A = calcA( model, dt )        
        Q = model.rates;
        I = logical(eye(size(Q)));
        Q(I) = 0;
        Q(I) = -sum( Q,2 );
        A = expm( Q.*dt );
    end
    
    % Calculate the initial state probabilities (p0) expected if the system
    % is in equilibrium. This is calculated from the rate matrix.
    % Since this approximation uses the (discrete) probability matrix A,
    % the time step (dt, in seconds) is also needed.
    % FIXME: there is probably a way to directly calculate from rates.
    function p = calcEquilibriumP0( model, dt )
        A = model.calcA( dt );
        p = A^10000;
        p = p(1,:);
        p = p/sum(p);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   GUI METHODS   %%%%%%%%%%%%%%%%%%%%%%%%%
    function parent = showModel(model,parent)
    % Displays a representation of the model similar to the one used in
    % QuB. The user can click on the states or rate constants to change
    % them and this updates the object's properties automatically.
    % PARENT is the GUI object to draw in. If not specified, a new figure
    % is created.
    % FIXME: this displays the model /as it was originally/, but the
    % parameters can be modified after loading!
        
    % Verify the model makes sense before trying to draw it.
    model.verify();
    tree = model.qubTree;
    
    % Formatting options:
    boxsize = 6;
    linewidth = 2;
    % boxsize = tree.Properties.StateSize*6/20;  % state square size in pixels
    % linewidth = tree.Properties.LineWidth+2;
    %                k       r       b     dark g      y       m
    colors     = {[0 0 0],[1 0 0],[0 0 1],[0 0.7 0],[1 1 0],[1 0 1]};  % class colors
    textColors = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],  [0 0 0],[0 0 0]};  % class colors

    textFormat = {'Parent',parent,'FontSize',8, 'FontWeight','bold', ...
              'HorizontalAlignment','center', 'VerticalAlignment','middle' };
    lineFormat = {'Parent',parent, 'Color','k', 'LineWidth',linewidth };
    
    % Create a window for the model, or use one if given.
    if nargin<2,
        figure;
        parent = axes;
    end
    set(gcf,'WindowButtonMotionFcn',@figButtonMotion,'WindowButtonUpFcn',@dropObject); %FIXME
    draggedBox = [];  %no box is being dragged right now.
    cla; hold on;

    % Get coordinates of states
    nStates = model.nStates;
    c = model.class;

    % Draw all specified transitions with rates. Each element specifies both
    % the forward and reverse rates for a particular transition type.
    nRates = numel( tree.Rates.Rate );
    hRate = zeros( nRates,2 );
    hLine = zeros( nRates,1 ); %main line and two direction arrows

    for i=1:nRates,
        states = tree.Rates.Rate(i).States.data+1; %state numbers are zero-based.
        k0  = [ model.rates(states(1),states(2)) model.rates(states(2),states(1)) ];

        % Some new models made by QuB have a second column, unknown purpose!
        if all( size(states)>1 ),
            assert( all(size(states)==2) );
            states = states(:,1);
        end
        
        hLine(i) = line( 0, 0, lineFormat{:} );

        % Display rate numbers
        hRate(i,1) = text( 0,0, num2str(k0(1)), textFormat{:}, 'Color',[0 0 0] );
        hRate(i,2) = text( 0,0, num2str(k0(2)), textFormat{:}, 'Color',[0 0 0] );

        % Add callbacks for changing the rate constants.
        set( hRate(i,1), 'ButtonDownFcn', {@editRate,states       } );
        set( hRate(i,2), 'ButtonDownFcn', {@editRate,states([2 1])} );
    end

    % Create the lines connecting states and rate number text next to them.
    hBox  = zeros(nStates,1);  %handles to the graphics objects
    hText = zeros(nStates,1);
    
    xx = [-0.5 0.5  0.5 -0.5]*boxsize;  %coordinates defining a generic box.
    yy = [ 0.5 0.5 -0.5 -0.5]*boxsize;       

    for i=1:nStates,
        % Draw box for each state. 
        hBox(i) = patch( xx+model.x(i), yy+model.y(i), colors{c(i)},'Parent',parent );

        % State number (which may be different from class number)
        hText(i) = text( model.x(i),model.y(i), num2str(i), textFormat{:}, ...
                                'Color',textColors{c(i)},'Parent',parent );

        % Add a callback to the box so the properties can be seen and changed.
        set( hBox(i), 'ButtonDownFcn', {@editState,i} );
        set( hText(i),'ButtonDownFcn', {@editState,i} );
    end
    
    % Position the lines and text in the correct places.
    moveLines();

    axis equal;
    set(parent,'ydir','reverse');  %mimic orientation in QuB
    set(parent,'YTick',[]);
    set(parent,'XTick',[]);
    xlim([5 90]);
    ylim([7 90]);
    
    
    
    %%%%------  DRAWING FUNCTIONS   ------%%%%

    function moveLines(varargin)
    % Positions the lines connecting states and rate numbers next to them. This
    % is called both to initially display the model and whenever a state box is
    % moved by the user.
    r = 0.7*boxsize;  %distance from the line to draw text
    
    for i=1:nRates,
        states = tree.Rates.Rate(i).States.data+1; %state numbers are zero-based.

        % Some new models made by QuB have a second column, unknown purpose!
        if all( size(states)>1 ),
            assert( all(size(states)==2) );
            states = states(:,1);
        end

        state_x = model.x(states);
        state_y = model.y(states);

        % Display the rate numbers
        % Use some fancy geometry to position the numbers above and below the
        % line regardless of their orientation.
        t = atan2( diff(state_y), diff(state_x) );
        mx = mean(state_x);  my = mean(state_y);  %center between states

        % Draw the connecting lines for each rate. The fancy math is shorten it
        % to fit between the boxes. The other lines are arrowheads that
        % indicate which rate is described by the numbers.
        len = sqrt( diff(state_x)^2 + diff(state_y)^2 )/2 - 0.75*boxsize;
        line_x = [mx-len*cos(t) mx+len*cos(t)];
        line_y = [my-len*sin(t) my+len*sin(t)]; %these define the main line connecting states

        line_xx = [line_x(1)+(boxsize/2)*cos(t+(40*pi/180)) line_x(1) line_x(2) line_x(2)-(boxsize/2)*cos(t+(40*pi/180))];
        line_yy = [line_y(1)+(boxsize/2)*sin(t+(40*pi/180)) line_y(1) line_y(2) line_y(2)-(boxsize/2)*sin(t+(40*pi/180))];
        set( hLine(i), 'XData', line_xx );
        set( hLine(i), 'YData', line_yy );

        % Display rate numbers
        set( hRate(i,1), 'Position',[mx-r*cos(t+pi/2),my-r*sin(t+pi/2)] );
        set( hRate(i,2), 'Position',[mx+r*cos(t+pi/2),my+r*sin(t+pi/2)] );
    end
    
    end %function moveLines
    
    
    
    %%%%------  GUI CALLBACK FUNCTIONS   ------%%%%

    function editRate( hObject, ~, rateID )
    % Called whenever one of the rate labels is clicked.
    % NOTE: this does not update the qubTree object!!

    % Ask the user for the new value for the rate constant.
    initRate = model.rates( rateID(1), rateID(2) );
    text = sprintf('Enter a new rate constant (%d->%d)',rateID(1),rateID(2));
    a = inputdlg( text, 'Edit model parameters', 1, {num2str(initRate)} );

    % Save the value. It is saved in the showModel() function's scope because
    % this is a nested function (right?)
    if ~isempty(a), 
        a = str2double(a);
        model.rates( rateID(1), rateID(2) ) = a;
        set(hObject,'String', num2str(a) );
    end

    end %function editRate


    
    function editState( ~, ~, stateID )
    % Called whenever one of state boxes is clicked.
    % NOTE: this does not update the qubTree object!!

    % If the user left clicks, allow the box to be moved. If it is a right-click
    % or something else, let the user edit the settings for the state.
    if strcmpi( get(gcf,'SelectionType'), 'normal' ),
        draggedBox = stateID;
        return;
    end
    
    classID = model.class(stateID);
    initVal = { model.p0(stateID), model.mu(classID), model.sigma(classID) };

    prompt = {'Enter start prob:','Enter mean value:','Enter stdev:'};
    defaults = cellfun( @num2str, initVal, 'UniformOutput',false );

    a = inputdlg( prompt, 'Edit model parameters', 1, defaults );

    % Save the value. It is saved in the showModel() function's scope because
    % this is a nested function (right?)
    if ~isempty(a),
        a = cellfun(@str2double,a);

        % Verify and update class number. The mu/sigma vectors would have
        % to be updated too, so this is disabled for now.
        %if( floor(a(1))~=a(1) ),
        %    warning('Invalid class number');
        %else
        %    model.class(stateID) = a(1);
        %
        %    % If the class number is changed, update the box color as well.
        %    set( hBox(stateID), 'FaceColor', colors{a(1)} );
        %    set( hText(stateID), 'Color', textColors{a(1)} );
        %end

        model.p0(stateID)    = a(1);
        model.mu(classID)    = a(2);
        model.sigma(classID) = a(3);
    end

    end %function editState

    
    % Called when the user is dragging one of the state boxes
    function figButtonMotion(varargin)
        if ~isempty(draggedBox)
            % Get the Mouse Location
            curr_pt = get(parent,'CurrentPoint');
            curr_pt = max(5, min(90,curr_pt) );
            
            model.x(draggedBox) = curr_pt(1,1);
            model.y(draggedBox) = curr_pt(1,2);
            
            % Change the Position of the Patch
            set( hBox(draggedBox), 'XData',xx+curr_pt(1,1), 'YData',yy+curr_pt(1,2) );
            set( hText(draggedBox), 'Position',[curr_pt(1,1),curr_pt(1,2)] );
            
            % Update lines and rates
            moveLines();
        end
    end

    % Called when the state being dragged is released.
    function dropObject(varargin)
        draggedBox = [];
    end

    end %function showModel
    
    
end %methods



end  %classdef


