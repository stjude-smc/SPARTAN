classdef TraceListViewer < handle
% TraceListViewer   Display an interactive list of fret traces.
%
%   viewer = TraceListViewer(Ax, sldY, sldX) creates a new trace list viewer
%   that will draw traces in the axes AX. sldY and sldX are handles to slider
%   objects to scroll through traces and zoom in on time, respectively.
%
%   To display trace data, the following data fields can be set:
%
%    - data:  TracesFret object containing traces to display (required)
%    - dataFilename: path to associated .traces file (optional)
%    - idl:   state assignment matrix / idealized traces (optional)
%    - model: QubModel object for showing model fret value lines (optional)
%
%   The following optional arguments control what is displayed:
%
%    - nTracesToShow: number of traces to display at once
%    - showStateMarkers: show model fret value markers as colored, dotted lines
%    - dataField: which field in the TracesFret object data to plot
%    - exclude: logical array, true if trace should be excluded from analysis
%
%   NOTE: The viewer does not update automatically when properties are changed.
%   To draw the viewer with new data, call the redraw() method. For minor
%   updates like adding an idealization or excluding traces, call the
%   showTraces() method. To update fret value markers, call showModelLines().
%
%   See also: batchKinetics, QubModelViewer.

%   Copyright 2017 Cornell University All Rights Reserved.



%% ---------------------------  PROPERTIES  --------------------------- %%

properties (SetAccess=public, GetAccess=public)
    data;             % Traces object represented by this GUI
    dataFilename;     % Full path to .traces file associated with data
    model;            % QubModel object for drawing model fret value markers
    idl;              % State assignment matrix (idealized traces)
    
    % Display settings that can be modified by user.
    nTracesToShow = 6;        % Number of traces to display in viewer
    showStateMarkers = true;  % Draw dotted lines to mark model fret values
    dataField = 'fret';       % Which field of target Traces object to show
    exclude = [];             % Logical array marking traces to exlude from analysis
end

properties (SetAccess=protected, GetAccess=public)
    ax;                      % Target axes in which traces are drawn
    sldTraces  = [];         % Handles to Y axis slider to scroll through traces
    sldTracesX = [];         % Handles to X axis slider zoom in on traces
    sldTracesListener = [];  % Listener that detects when scrollbars are altered.
    
    hFretLine   = [];        % Array of Line objects showing FRET traces,
    hIdlLine    = [];        %    idealization traces,
    hTraceLabel = [];        %    trace numbers.
    contextMenu = [];
    
    traceColor  = [0 0 1];   % RGB color of traces,
    idlColor    = [1 0 0];   %     idealization
end




methods
    %% -------------------------  CONSTRUCTOR  ------------------------- %%
    function this = TraceListViewer(target, hSld, hSldX)
    % Create a new viewer.
    
    narginchk(3,3); nargoutchk(0,1);
    assert( isscalar(target) && ishghandle(target) );  %isgraphics(varargin{1},'axes')
    
    this.ax = target;
    this.sldTraces = hSld;
    this.sldTracesX = hSldX;
    
    % Setup slider listeners for scrolling through data
    sl(1) = addlistener( hSld,  'Value', 'PostSet',@(h,e)this.showTraces() );
    sl(2) = addlistener( hSldX, 'Value', 'PostSet',@this.sldTracesX_Callback );
    this.sldTracesListener = sl;
    enableListener(sl, false);
    set( get(this.ax,'Parent'), 'WindowScrollWheelFcn', @this.wheelScroll_callback );
    
    % Add context menu for trace options.
    menu = uicontextmenu( ancestor(this.ax,'figure'), 'Visible','off' );
    uimenu( menu, 'Label','Display settings...', 'Callback',@this.mnuDisplaySettings_Callback );
    uimenu( menu, 'Label','Include all traces', 'Callback',@(h,e)this.mnuIncludeAll_Callback(false), 'Separator','on' );
    uimenu( menu, 'Label','Exclude all traces', 'Callback',@(h,e)this.mnuIncludeAll_Callback(true) );
    uimenu( menu, 'Label','Load selection list...', 'Callback', @this.mnuLoadSelList_Callback );
    uimenu( menu, 'Label','Save selection list...', 'Callback', @this.mnuSaveSelList_Callback );
    set( allchild(menu), 'Enable','off' );
    set(this.ax, 'UIContextMenu', menu);
    this.contextMenu = menu;
    
    xlabel(this.ax,'Time (s)');
    
    end %constructor
    
    
    
    
    %% ----------------------  DRAWING FUNCTIONS  ---------------------- %%
    
    function redraw(this, varargin)
    % (re-)initialize trace viewer, creating plot objects, etc.
    % Run this method when a new data file is loaded or nTracesToShow or
    % dataField properties are changed.
    %
    % Draws traces from current file in the trace viewer panel.
    % Creates line objects for data plotting, which are updated when the user
    % scrolls or makes other changes by calling showTraces() below.
    % 
    % sldTraces: scrolls vertical list of traces. Value is trace shown at top.
    %            Max means the Value showing the first traces (starting with 1).

    % Reset figure 
    cla(this.ax);
    [this.hFretLine, this.hIdlLine, this.hTraceLabel] = deal([]);
    set( allchild(this.contextMenu), 'Enable','off' );
    
    if isempty(this.data), return; end  %no data loaded.

    this.exclude = false(this.data.nTraces,1);
    
    % Update scroll bar limits, disabling listeners to prevent errors.
    enableListener(this.sldTracesListener, false);
    set( this.sldTraces, 'Enable',onoff(this.data.nTraces>this.nTracesToShow), ...
                         'Min',min(this.data.nTraces,this.nTracesToShow),  ...
                         'Max',this.data.nTraces, 'Value',this.data.nTraces );
    set(this.sldTracesX, 'Min',10, 'Max',this.data.nFrames, 'Value',this.data.nFrames);
    enableListener(this.sldTracesListener, true);

    % Determine colors for data display
    if strcmpi(this.dataField,'fret'),
        this.traceColor = [0 0 1];
    elseif strcmpi(this.dataField,'fret2')
        this.traceColor = [1 0 1];
    else
        idxch = find(strcmpi(this.dataField,this.data.channelNames) );
        if ~isempty(idxch) && any(idxch==this.data.idxFluor) && ...
                              isfield(this.data.fileMetadata,'wavelengths')
            this.traceColor = Wavelength_to_RGB( this.data.fileMetadata.wavelengths(idxch) );
        else
            this.traceColor = [0 0 0];
        end
    end

    % Setup axes for plotting traces.
    time = this.data.time/1000;
    xlim( this.ax, [0 time(end)] );
    dt = this.data.sampling/2/1000; %put idl transitions between FRET datapoints.

    for i=1:this.nTracesToShow,
        y_offset = 1.18*(this.nTracesToShow-i) +0.2;

        plot( this.ax, time([1,end]), y_offset+[0 0], 'k:', 'HitTest','off' );  %baseline marker

        % Draw traces, idealizations, and trace number labeles.
        z = y_offset + zeros(1,this.data.nFrames);
        this.hFretLine(i) = plot( this.ax, time, z, 'HitTest','off' );
        this.hIdlLine(i)  = stairs( this.ax, time-dt, z, 'r-', 'HitTest','off' );

        this.hTraceLabel(i) = text( 0.98*time(end),y_offset+0.1, '', ...
                   'Parent',this.ax, 'BackgroundColor','w', ...
                   'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
                   'ButtonDownFcn',@this.traceLabel_Callback, ...
                   'UIContextMenu',this.contextMenu);
    end

    ylim(this.ax, [0 1.2*this.nTracesToShow]);

    this.showTraces();
    this.showModelLines();
    set( allchild(this.contextMenu), 'Enable','on' );

    end %function redraw


    
    function showTraces(this, varargin)
    % Update ax to show the current subset -- called by sldTraces.
    % i is the index into the lines in the viewer.
    % idx is the index into the traces in the whole file.

    if isempty(this.data), return; end
    
    idxStart = get(this.sldTraces,'Max')-floor(get(this.sldTraces,'Value'));
    nToShow = min(this.nTracesToShow, this.data.nTraces);
    exText = {'',' Excluded'};

    for i=1:nToShow
        idx = i+idxStart;
        y_offset = 1.18*(this.nTracesToShow-i) +0.2;
        ex = this.exclude(idx);

        % Redraw FRET and idealization traces
        ydata = this.data.(this.dataField)(idx,:);
        ydata = y_offset + min(1.15, max(-0.15,ydata) );  %clip outliers, position w/i viewer
        set( this.hFretLine(i), 'YData',ydata, 'Color',min(1,this.traceColor+0.65*ex) );

        if ~isempty(this.idl)
            ydata = y_offset + this.idl(idx,:);
            set( this.hIdlLine(i), 'YData',ydata, 'Color',min(1,this.idlColor+0.65*ex) );
        end

        traceLabel = sprintf('%d%s',idx, exText{ex+1});
        set( this.hTraceLabel(i), 'String',traceLabel, 'UserData',idx );
    end

    set( this.hIdlLine, 'Visible',onoff(~isempty(this.idl)) );
    
    % Clear final lines if there isn't enough data to fill them.
    set( this.hFretLine(nToShow+1:end), 'YData',nan(1,this.data.nFrames) );
    set( this.hIdlLine(nToShow+1:end),  'YData',nan(1,this.data.nFrames) );
    set( this.hTraceLabel(nToShow+1:end), 'String','' );

    end %function showTraces
    
    
    
    function showModelLines(this)
    % Draw dotted lines to indicate model FRET values in trace viewer panel.
    % Should only be called when new data are loaded or model is modified.

    delete( findall(this.ax,'Tag','ModelMarker') );  %clear existing markers
    
    if ~this.showStateMarkers || isempty(this.data) || isempty(this.model)
        return;
    end

    nToShow = min(this.nTracesToShow, this.data.nTraces);
    time = this.data.time([1,end])/1000;
    mu = this.model.mu;
    h = [];

    % Redraw model FRET value markers.
    for i=1:nToShow
        y_offset = 1.18*(this.nTracesToShow-i) +0.2;

        for k=2:this.model.nClasses
            h(end+1) = plot( this.ax, time, y_offset+mu([k k]), ':', ...
                             'Color',QubModelViewer.colors{k}, 'Tag','ModelMarker', 'HitTest','off' ); %#ok<AGROW>
        end
    end

    uistack(h,'bottom');  %draw markers below data.

    end %function showModelLines





    %% ----------------------   CALLBACK FUNCTIONS   ---------------------- %%

    function mnuDisplaySettings_Callback(this, varargin)
    % Change display settings (e.g., number of traces displayed).
    if isempty(this.data), return; end
    
    opt = struct('dataField',this.dataField,'nTracesToShow',this.nTracesToShow, 'showStateMarkers',this.showStateMarkers);
    prompt = {'Data field', 'Number of traces to show', 'Show model FRET values over traces'};
    types = { this.data.channelNames(this.data.idxFret), @(x)(x==round(x)&x>0) };  %'total' field inaccessable...
    opt = settingdlg(opt, fieldnames(opt), prompt, types);

    if ~isempty(opt),
        this.dataField = opt.dataField;
        this.nTracesToShow = opt.nTracesToShow;
        this.showStateMarkers = opt.showStateMarkers;
        this.redraw();
    end

    end %function mnuDisplaySettings_Callback

    

    function sldTracesX_Callback(this, varargin)
    % Called when X axis slider is adjusted. Zoom in on early times.
    if isempty(this.data), return; end
    xlimit = floor(get(this.sldTracesX,'Value'));
    set( this.ax, 'XLim',[0 this.data.time(xlimit)/1000] );

    % Ensure trace number labels stay in the same spot
    for i=1:this.nTracesToShow,
        p = get(this.hTraceLabel(i), 'Position');
        p(1) = 0.98*this.data.time(xlimit)/1000;
        set( this.hTraceLabel(i), 'Position',p );
    end
    end %function sldTracesX_Callback


    
    function wheelScroll_callback(this, ~, eventData)
    % Mouse wheel scrolling moves the trace viewer pane up and down.
    % The event is triggered at the figure level.
    loc = get(this.sldTraces, 'Value')-3*eventData.VerticalScrollCount;
    loc = min( loc, get(this.sldTraces,'Max') );
    loc = max( loc, get(this.sldTraces,'Min') );
    set(this.sldTraces, 'Value', loc);  %triggers listener, updating viewer.
    end %function wheelScroll_callback



    function traceLabel_Callback(this, hObject, ~)
    % Executes when user clicks on trace number text in trace viewer panel.
    % Togger whether the exclude/include the trace in analysis.
    if strcmpi( get(gcbf,'SelectionType'), 'alt' ), return; end  %right-click passes through

    idxTrace = get(hObject,'UserData');
    this.exclude(idxTrace) = ~this.exclude(idxTrace);
    this.showTraces();
    end %function traceLabel_Callback



    function mnuIncludeAll_Callback(this, value)
    % Include or exclude all traces in trace viewer.
    this.exclude(:) = value;
    this.showTraces();
    end %function mnuIncludeAll_Callback



    function mnuLoadSelList_Callback(this, varargin)
    % Load text file listing which traces to include for analysis.
    % FIXME: check for out of bound traces?

    if ~isempty(this.dataFilename)
        [p,f] = fileparts( this.dataFilename );
    else
        p=pwd; f='';
    end
    fname = getFile( fullfile(p,[f '_sel.txt']), 'Load selection list' );

    if ~isempty(fname)
        fid = fopen(fname,'r');
        idx = fscanf(fid, '%d');
        fclose(fid);
        this.exclude(:) = true;
        this.exclude(idx) = false;
        this.showTraces();
    end

    end %function mnuLoadSelList_Callback

    

    function mnuSaveSelList_Callback(this, varargin)
    % Save text file listing which traces to include for analysis.

    if ~isempty(this.dataFilename)
        [p,f] = fileparts( this.dataFilename );
    else
        p=pwd; f='';
    end
    [f,p] = uiputfile( fullfile(p,[f '_sel.txt']), 'Save selection list' );
    if isequal(f,0), return; end  %user hit cancel

    fid = fopen( fullfile(p,f), 'w');
    fprintf( fid, '%d ', find(~this.exclude) );
    fclose(fid);

    end %function mnuSaveSelList_Callback

    

end %methods

end %classdef








