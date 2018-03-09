classdef QubModelViewer < handle
% QubModelViewer   Display a graphical, interactive QuB model interface.
%
%   Displays a representation of the model similar to the one used in
%   QuB. The user can click on the states or rate constants to change
%   them and this updates the object's properties automatically.
%   PARENT is the GUI object to draw in. If not specified, a new figure
%   is created.
%
%      QubModelViewer(model) displays the given QubModel object.
%      QubModelViewer(filename) loads the given .qmf file as a new QubModel.
%      QubModelViewer(..., axes) specifies the target axes to draw in. If none
%          given, a new figure will be created.
%   
%   See also: QubModel, batchKinetics, simulate.

%   Copyright 2007-2016 Cornell University All Rights Reserved.

% FIXME: destroy viewer if target axes destroyed?



%% ---------------------  CONSTRUCT HANDLES OBJECT  --------------------- %%

properties (SetAccess=protected, GetAccess=public)
    ax;               %target axes in which model is drawn
    model;            %QubModel object represented by GUI
    
    hBox  = [];       %handles to state box GUI objects
    hText = [];       %handles to text in state box
    hLine = [];       %handles to lines between states
    hRate = [];       %handles to text above rates
    
    draggedBox = [];  %state box currently being dragged (if any)
end

properties (GetAccess=public, Constant)
    % Display settings
    %                k       r       b     dark g      y       m
    colors     = {[0 0 0],[1 0 0],[0 0 1],[0 0.7 0],[1 1 0],[1 0 1]};  % class colors
    textColors = {[1 1 1],[1 1 1],[1 1 1],[1 1 1],  [0 0 0],[0 0 0]};  % class colors
    boxsize = 6;
    boxx = [-0.5 0.5  0.5 -0.5]*6;  %coordinates defining a generic box.
    boxy = [ 0.5 0.5 -0.5 -0.5]*6;
end




methods
    %% ------------------------- CONSTRUCTOR ------------------------- %%
    function this = QubModelViewer(input, target)
    % Create a new viewer.
    
    narginchk(1,2); nargoutchk(0,1);

    % Initialize the internal QubModel
    if ischar(input)
        if exist(input,'file'),
            this.model = QubModel(input);
        else
            error('Input file does not exist');
        end

    elseif isa(input,'QubModel')
        this.model = input;

    else
        error('Invalid QubModelViewer initialization arguments');
    end

    % Verify the model makes sense before trying to draw it.
    this.model.verify();
    
    % Look for initial axes argument. Create a new figure if none given.
    if nargin>=2
        assert( isscalar(target) && ishghandle(target) );  %isgraphics(varargin{1},'axes')
        this.ax = target;
        hf = ancestor(target,'figure');
    else
        hf = figure;
        this.ax = axes('Parent',hf);
    end

    % Setup internal variables.
    set(hf,'WindowButtonMotionFcn',@this.figButtonMotion,  ...
           'WindowButtonUpFcn',    @this.dropObject  );
    
    this.redraw();
    
    if ~isempty(this.model.filename)
        [~,f] = fileparts(this.model.filename);
    else
        f = 'New unsaved model';
    end
    title(this.ax, f, 'interpreter', 'none');
    
    end %constructor
    
    
    
    %% ----------------------- DRAWING FUNCTIONS ----------------------- %%
    
    function redraw(this)
    % Completely redraw all GUI elements

    % Formatting options:
    textFormat = {'Parent',this.ax,'FontSize',8, 'FontWeight','bold', ...
              'HorizontalAlignment','center', 'VerticalAlignment','middle' };
    lineFormat = {'Parent',this.ax, 'Color','k', 'LineWidth',2 };

    % Create a window for the model, or use one if given.
    this.draggedBox = [];  %no box is being dragged right now.
    cla(this.ax); hold(this.ax, 'on');

    % Draw all specified transitions with rates. Each element specifies both
    % the forward and reverse rates for a particular transition type.
    conn = this.model.connections;
    nRates = size(conn,1);

    this.hRate = zeros( nRates, 2 );
    this.hLine = zeros( nRates, 1 ); %main line and two direction arrows

    for i=1:nRates,
        states = conn(i,:);
        k0  = [ this.model.rates(states(1),states(2)) this.model.rates(states(2),states(1)) ];

        this.hLine(i) = line( 0,0, lineFormat{:} );

        % Display rate numbers
        this.hRate(i,1) = text( 0,0, num2str(k0(1),2), textFormat{:}, 'Color',[0 0 0] );
        this.hRate(i,2) = text( 0,0, num2str(k0(2),2), textFormat{:}, 'Color',[0 0 0] );

        % Add callbacks for changing the rate constants.
        set( this.hRate(i,1), 'ButtonDownFcn', {@this.editRate_callback,states}         );
        set( this.hRate(i,2), 'ButtonDownFcn', {@this.editRate_callback,states([2 1])}  );
    end

    % Create the lines connecting states and rate number text next to them.
    this.hBox  = zeros(this.model.nStates,1);  %handles to the graphics objects
    this.hText = zeros(this.model.nStates,1);

    for i=1:this.model.nStates,
        c = this.model.class(i);

        % Draw box for each state. 
        this.hBox(i) = patch( this.boxx+this.model.x(i), this.boxy+this.model.y(i), ...
                                 this.colors{c}, 'Parent',this.ax, 'UserData',i );

        % State number (which may be different from class number)
        this.hText(i) = text( this.model.x(i),this.model.y(i), num2str(i), textFormat{:}, ...
                       'Color',this.textColors{c}, 'Parent',this.ax, 'UserData',i );
    end

    % Add a callback to the box so the properties can be seen and changed.
    set([this.hBox; this.hText], 'ButtonDownFcn', @this.stateClicked_callback );

    hFig = ancestor(this.ax,'figure');
    menu = uicontextmenu(hFig);
    uimenu( menu, 'Label','State properties...', 'Callback', @this.editState_callback   );
    uimenu( menu, 'Label','Class properties...', 'Callback', @this.editClass_callback  );
    uimenu( menu, 'Label','Connect to...',       'Callback', @this.connect_callback     );
    uimenu( menu, 'Label','Remove state',        'Callback', @this.removeState_callback );
    set([this.hBox; this.hText], 'UIContextMenu', menu);


    % Position the lines and text in the correct places.
    this.moveLines();

    axis(this.ax,'equal');
    set(this.ax,'ydir','reverse');  %mimic orientation in QuB
    set(this.ax,'YTick',[]);
    set(this.ax,'XTick',[]);
    xlim(this.ax,[5 90]);
    ylim(this.ax,[7 90]);

    %---- Add context menu to for additional options
    menu = uicontextmenu(hFig);
    uimenu( menu, 'Label','Add state', 'Callback',@this.addState_callback );
    uimenu( menu, 'Label','Calculate and apply equilibrium p0', ...
                             'Callback',@this.calcEqP0_callback );
    %uimenu( menu, 'Label','Enforce loop balance' );  %see Colquhoun 2004, Biophys J 86, p. 3510. Minimum spanning tree method.
    uimenu( menu, 'Label','Revert to saved', 'Callback', @this.revert_callback, 'Separator','on' );
    uimenu( menu, 'Label','Save model as...', 'Callback', @this.save_callback );
    set(this.ax, 'UIContextMenu', menu);

    end %function showModel
    
    
    
    function moveLines(this)
    % Positions the lines connecting states and rate numbers next to them. This
    % is called both to initially display the model and whenever a state box is
    % moved by the user.

    r = 0.6*this.boxsize;  %distance from the line to draw text

    conn = this.model.connections;
    nRates = size(conn,1);

    for i=1:nRates,
        states = conn(i,:);
        state_x = this.model.x(states);
        state_y = this.model.y(states);

        % Display the rate numbers
        % Use some fancy geometry to position the numbers above and below the
        % line regardless of their orientation.
        t = atan2d( diff(state_y), diff(state_x) );
        mx = mean(state_x);  my = mean(state_y);  %center between states

        % Draw the connecting lines for each rate. The fancy math is shorten it
        % to fit between the boxes. The other lines are arrowheads that
        % indicate which rate is described by the numbers.
        len = sqrt( diff(state_x)^2 + diff(state_y)^2 )/2 - 0.75*this.boxsize;
        line_x = [mx-len*cosd(t) mx+len*cosd(t)];
        line_y = [my-len*sind(t) my+len*sind(t)]; %these define the main line connecting states

        line_xx = [line_x(1)+(this.boxsize/2)*cosd(t+40) line_x(1) line_x(2) line_x(2)-(this.boxsize/2)*cosd(t+40)];
        line_yy = [line_y(1)+(this.boxsize/2)*sind(t+40) line_y(1) line_y(2) line_y(2)-(this.boxsize/2)*sind(t+40)];
        set( this.hLine(i), 'XData', line_xx );
        set( this.hLine(i), 'YData', line_yy );

        % Display rate numbers
        set( this.hRate(i,1), 'Position',[mx-r*cosd(t+90),my-r*sind(t+90)], ...
                              'Rotation', rflct(180-t) );
        set( this.hRate(i,2), 'Position',[mx+r*cosd(t+90),my+r*sind(t+90)], ...
                              'Rotation', rflct(180-t) );
    end

    end %function moveLines

    



    %% -----------------------   CALLBACK FUNCTIONS   ----------------------- %%

    function editRate_callback(this, hObject, ~, rateID)
    % Called whenever one of the rate labels is clicked.
    
    % Ask the user for the new value for the rate constant.
    initRate = this.model.rates( rateID(1), rateID(2) );
    text = sprintf('Enter a new rate constant (%d->%d)',rateID(1),rateID(2));
    a = inputdlg( text, 'Edit model parameters', 1, {num2str(initRate)} );
    if isempty(a), return; end  %user hit cancel

    % Save the value.
    a = str2double(a);
    this.model.rates( rateID(1), rateID(2) ) = a;
    set(hObject,'String', num2str(a,2) );

    % Redraw if a connection was removed.
    if this.model.rates( rateID(1), rateID(2) )==0 && this.model.rates( rateID(2), rateID(1) )==0
        this.redraw();
    end

    end %function editRate


    function editState_callback(this, ~,~, stateID)
    % Called whenever one of state boxes is clicked.
        if nargin<4,  stateID = get(gco,'UserData');  end

        prompt = {'Start probability:','Class number:'};
        initVal = { this.model.p0(stateID), this.model.class(stateID) };
        defaults = cellfun( @num2str, initVal, 'UniformOutput',false );
        a = inputdlg( prompt, sprintf('Edit state %d',stateID), 1, defaults );
        if isempty(a), return; end

        % Check validity of inputs
        a = cellfun(@str2double,a);
        if a(1)<0 || a(1)>1 || ~ismember(a(2),1:10) || a(2)>this.model.nClasses+1,
            warndlg('Invalid parameter values');
            return;
        end

        % If this is a new class, use reasonable defaults.
        if a(2) > this.model.nClasses,
            this.model.addClass( a(2) );
        end
        this.model.p0(stateID)    = a(1);
        this.model.class(stateID) = a(2);

        % If the class number is changed, update the box color as well.
        set( this.hBox(stateID),  'FaceColor',this.colors{a(2)} );
        set( this.hText(stateID), 'Color',this.textColors{a(2)} );
       
    end %function editState


    function editClass_callback(this,varargin)
    % Change the class number for selected state (context menu callback)
        classID = this.model.class( get(gco,'UserData') );

        % Prompt user for new values
        prompt = {'Mean:','Stdev:'};
        initVal = { this.model.mu(classID), this.model.sigma(classID) };
        defaults = cellfun( @num2str, initVal, 'UniformOutput',false );
        a = inputdlg( prompt, sprintf('Edit class %d',classID), 1, defaults );

        % Check for valid values and save
        if ~isempty(a),
            a = cellfun(@str2double,a);
            if any(isnan(a)),
                warndlg('Invalid values');
                return;
            end

            this.model.mu(classID)    = a(1);
            this.model.sigma(classID) = a(2);
        end
    end


    function addState_callback(this,varargin)
    %
        this.model.addState();
        newState = numel(this.model.class);
        this.redraw();

        % Get the Mouse Location
        curr_pt = get(this.ax,'CurrentPoint');
        curr_pt = max(5, min(90,curr_pt) );

        this.model.x(newState) = curr_pt(1,1);
        this.model.y(newState) = curr_pt(1,2);
        xx = this.boxx + curr_pt(1,1);
        yy = this.boxy + curr_pt(1,2);

        % Change the Position of the Patch
        set( this.hBox(newState), 'XData',xx, 'YData',yy );
        set( this.hText(newState), 'Position',[curr_pt(1,1),curr_pt(1,2)] );
    end


    function removeState_callback(this,varargin)
    % State box context menu option to remove selected state.
        stateID = get(gco,'UserData');
        this.hBox(stateID)  = [];
        this.hText(stateID) = [];
        this.model.removeState(stateID);
        this.redraw();
    end
    
    
%     function removeConnction_callback(this,varargin)
%     % Remove a specified connection, setting rates to zero
%         
%     end


    function connect_callback(this,varargin)
    % Ask user which states to connect
        dst = str2double( inputdlg('Target state:','Connect states',1) );
        if isempty(dst)||isnan(dst), return; end  %user hit cancel

        src = get(gco,'UserData');

        % Add connection by setting non-zero values in rate matrix,
        % but do not modify if the connection already exists
        if this.model.rates(src,dst)==0 && this.model.rates(dst,src)==0
            this.model.rates(src,dst) = 1;
            this.model.rates(dst,src) = 1;
        end
        
        this.redraw();
    end


    function calcEqP0_callback(this,varargin)
    % Calculate state probabilities once the model reaches equilibrium.
    % This parameter is not displayed, so no need to redraw GUI.
        this.model.p0 = this.model.calcEquilibriumP0();
    end


    function save_callback(this,varargin)
    % Save the current model to file.
        if ~isempty(this.model.filename)
            [f,p] = uiputfile(this.model.filename);
        else
            [f,p] = uiputfile([pwd filesep '*.qmf']);
        end
        if ischar(f),
            fname = fullfile(p,f);
            this.model.save(fname);
            
            [~,f] = fileparts(fname);
            title(this.ax, f, 'interpreter', 'none');
        end
    end
    
    
    function revert_callback(this,varargin)
    % Save the current model to file.
        if ~isempty(this.model.filename)
            this.model.revert();
            this.redraw();
        end
    end




    %% ------------------   STATE DRAG-n-DROP CALLBACKS   ------------------ %%

    function figButtonMotion(this,varargin)
    % Called when the user is dragging one of the state boxes
        if ~isempty(this.draggedBox)
            % Get the Mouse Location
            curr_pt = get(this.ax,'CurrentPoint');
            curr_pt = max(5, min(90,curr_pt) );

            this.model.x(this.draggedBox) = curr_pt(1,1);
            this.model.y(this.draggedBox) = curr_pt(1,2);
            xx = this.boxx+curr_pt(1,1);
            yy = this.boxy+curr_pt(1,2);

            % Change the Position of the Patch
            set( this.hBox(this.draggedBox), 'XData',xx, 'YData',yy );
            set( this.hText(this.draggedBox), 'Position',[curr_pt(1,1),curr_pt(1,2)] );

            % Update lines and rates
            this.moveLines();
        end
    end

    function stateClicked_callback(this,varargin)
    % State box was clicked. Start dragging if single left click.
        if strcmpi( get(gcf,'SelectionType'), 'normal' ),
            this.draggedBox = get(gco,'UserData');
        end
    end

    function dropObject(this,varargin)
    % Called when the state being dragged is released.
        this.draggedBox = [];
    end


end %methods

end %classdef



function t = rflct(t)
    if t>90&&t<270, t=t-180; end
end





