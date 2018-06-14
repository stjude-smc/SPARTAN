classdef MovieViewer < handle
% MovieViewer   Show and play TIFF movies with separated fluorescence fields
%
%   viewer = MovieViewer( FNAME, PARAMS );
%   AX = viewer.show();
%
%   Display a TIFF movie, given by the path in the string FNAME, with controls
%   for adjusting intensity levels and playing the movie.
%   PARAMS struct must contain: geometry, wavelengths, chNames, chDesc, 
%     nAvgFrames, bgBlurSize (see @MovieParser/openStk.m for details).
%   AX is an array of the axes in which the movie images are displayed.
%
%   See also: gettraces, MovieMontage, MovieParser, Movie_TIFF, showMovie.

%   Copyright 2018 Cornell University All Rights Reserved.



%% ---------------------------  PROPERTIES  --------------------------- %%

properties (SetAccess=public, GetAccess=public)
end

properties (SetAccess=protected, GetAccess=public)
    movie;         % Movie object to be viewed
    params;        % Gettraces parameters (for geometry field). FIXME
    background;    % Movie background estimated from MovieParser/openStk.m
    %stk_top;       % Average of first 10 frames used for molecule picking.
    
    ax;            % Handles to subplot axes in montage
    hImg;          % Handle to image viewers
    btnPlay;       % Handle to 'Play' button
    sldScrub;      % Handle to scroll bar for scrubbing through time
    sldIntensity;  % Handle to scroll bar for adjusting intensity levels
    edTime;        % Text box showing current time
end


properties (Dependent)
    curFrame;  % Frame number (1..nFrames) being displayed
    curTime;   % Time (in seconds) of frame currently being displayed.
end



methods
    %% ---- Constructor ---- %%
    function this = MovieViewer(movieInput, paramsInput)
    % Create a new movie viewer window.
    
    narginchk(0,2); nargoutchk(0,1);
    
    % Load movie from file, prompting user if no arguments given.
    if nargin<1
        movieInput = getFile('*.tif');
        if isempty(movieInput), return; end  %user hit cancel
    end
    
    if ischar(movieInput)
        this.movie = Movie.load(movieInput);
    elseif isa(movieInput,'Movie')
        this.movie = movieInput;
    else
        error('Invalid input. Must be filename or Movie object');
    end
    
    % If no parameters given, show raw movie (no subfields)
    if nargin<2
        constants = cascadeConstants;
        paramsInput = constants.gettraces_profiles(1);
    end
    
    if ~isfield(paramsInput,'idxField')
        [val,idx] = sort( paramsInput.geometry(:) );
        paramsInput.idxFields = idx(val>0);
    end
    this.params = paramsInput;
    
    end %function load
    
    
    function load(this, movieInput)
    % Load a different movie using the same settings.
    % We assume this is the exact same type of movie (even the same number of
    % frames). Usually run from sorttraces. FIXME.
        this.movie = Movie.load(movieInput);
        sldScrub_Callback(this, this.sldScrub);  % Redraw field images
    end %function load
    
    
    function close(this)
        try
            close( get(this.ax,'Parent') );
            delete(this);
        catch
        end
    end
    

    function delete(this)
        try
            delete( get(this.ax(1),'Parent') );
        catch
        end
    end %function delete
    
    
    
    %% ------------------ Display viewer window ------------------ %%
    function axOut = show(this, hPanel, varargin)
    % Create a new figure displaying movie frames.
    
    hFig = figure;
    colormap(hFig, gettraces_colormap);
    
    % Use a movie parser to extract fields and approximate background image.
    parser = MovieParser(this.movie, this.params);
    this.background = parser.background;
    
    % Get intensity scale
    sort_px = sort( parser.stk_top{1}(:) );
    val = sort_px( floor(0.99*numel(sort_px)) );
    high = min( ceil(val*10), 32000 );  %uint16 maxmimum value
    
    % Create a uipanel for axes to sit in.
    % FIXME: handle to panel as argument (for gettraces integration)
    if nargin<2
        hPanel = uipanel( hFig, 'Position',[0.05  0.15  0.9  0.85], 'BorderType','none' );
    end

    % Create axes for sub-fields (listed in column-major order, like stk_top)
    axopt = {'Visible','off'};
    
    switch numel(this.params.idxFields)
    case 1
        this.ax    = axes( hPanel, 'Position',[0.05  0.15 0.9  0.85], axopt{:} );

    case 2
        this.ax    = axes( hPanel, 'Position',[0     0    0.45 0.95], axopt{:} );  %L
        this.ax(2) = axes( hPanel, 'Position',[0.5   0    0.45 0.95], axopt{:} );  %R

    case {3,4}
        this.ax    = axes( hPanel, 'Position',[0.0   0.5  0.325 0.47], axopt{:} );  %TL
        this.ax(2) = axes( hPanel, 'Position',[0     0    0.325 0.47], axopt{:} );  %BL
        this.ax(3) = axes( hPanel, 'Position',[0.335 0.5  0.325 0.47], axopt{:} );  %TR
        this.ax(4) = axes( hPanel, 'Position',[0.335 0    0.325 0.47], axopt{:} );  %BR

    otherwise
        error('Invalid field geometry');
    end

    % Show fluorescence fields for all channels
    % NOTE: all properties are listed in wavelength order. Use idxFields to
    % translate to the order of physical channels, which may be different.
    axopt = {'YDir','reverse', 'Color',get(hFig,'Color'), 'Visible','off'};

    for i=1:numel(parser.stk_top),
        target = this.ax( this.params.idxFields(i) );
        this.hImg(i) = image( target, parser.stk_top{i}, 'CDataMapping','scaled' );
        set( target, 'UserData',i, 'CLim',[0 val], axopt{:} );
        %colormap(target, gettraces_colormap);
    end
    setAxTitles(this);
    linkaxes( this.ax );
    axOut = this.ax;
    
%     if abs(log2( size(total,2)/size(total,1) )) < 1
        % For roughly symmetric movies, allows for better window resizing.
        axis( this.ax, 'image' );  %adjust axes size to match image
%     else
%         % Allows for better zooming of narrow (high time resolution) movie.
%         axis( [ax handles.axTotal], 'equal' );  %keep the axes size fixed.
%     end

    % Add intensity and time scroll bars and play button.
    sldStyle = {'style','slider','units','normalized'};
    this.sldIntensity = uicontrol(sldStyle{:}, 'position',[0.01 0.1 0.025 0.8], 'Callback',@this.sldIntensity_Callback);
    set(this.sldIntensity,'min',0, 'max',high, 'value',val);
    
    this.sldScrub = uicontrol(sldStyle{:}, 'position',[0.185 0.05 0.65 .05], 'callback',@this.sldScrub_Callback);
    set( this.sldScrub, 'Min',1, 'Max',this.movie.nFrames, 'Value',1, ...
           'SliderStep',[1/this.movie.nFrames,0.02] );
    
    this.edTime = uicontrol('Style','Edit', 'Units','normalized', 'Enable','off', ...
            'Position',[0.85 0.05 0.1 0.05], 'String','0 s');
       
    this.btnPlay = uicontrol('style','pushbutton', 'units','normalized', 'String','Play', ...
            'Position',[0.08 0.05 0.08 .05], 'Callback',@this.btnPlay_Callback );
    
    end %function show
    
    
    function setAxTitles(this)
    % Set axes titles from imaging profile settings, including colors.

    % Create new titles
    p = this.params;
    for i=1:numel(p.idxFields)  %i is channel index
        chColor = Wavelength_to_RGB( p.wavelengths(i) );

        title( this.ax( p.idxFields(i) ), ...
               sprintf('%s (%s) #%d',p.chNames{i},p.chDesc{i},i), ...
               'BackgroundColor',chColor, 'FontSize',10, 'Visible','on', ...
               'Color',(sum(chColor)<1)*[1 1 1] ); % White text for dark backgrounds.
    end
    
    end %FUNCTION setTitles
    
    
    
    function highlightPeaks(this, coords)
    % Draw a circle to highlight a molecule's PSF in each channel.
    % Coords is a cell array, one per channels, with x in first col and y in
    % second, with molecules listed in rows.
    
    if ~all(ishandle(this.ax)), return; end  %window closed?
    delete( findall(this.ax,'type','Line') );
    
    [ny,nx] = size( this.background{1} );

    for i=1:numel(coords)
        x = rem( coords{i}(:,1), nx);
        y = rem( coords{i}(:,2), ny);
            
        % Translate coordinates from stitched movie to subfield
        if numel(this.params.idxFields)>1
            viscircles( this.ax(i), [x y], 3, 'EdgeColor','w' );
        else
            viscircles( this.ax, [x y], 3, 'EdgeColor','w' );
        end
    end
    
    % If zoomed in, recenter on around new peak.
    dx = diff(get(this.ax(end), 'XLim'))/2;
    dy = diff(get(this.ax(end), 'YLim'))/2;
    if dx<nx/2 || dy<ny/2
        set( this.ax(end), 'XLim',[x-dx x+dx], 'YLim',[y-dy y+dy] );
    end
    
    end %function highlightPeaks

    
    
    
    %% ----------------------- Callback Functions ----------------------- %%
    
    % --- Executes on slider movement.
    function sldIntensity_Callback(this, hObject, varargin)
    % Update axes color limits from new slider value

    minimum = get(hObject,'min');
    maximum = get(hObject,'max');

    val = get(hObject,'value');
    val = max(val,minimum+1); %prevent errors in GUI
    maximum = max(val,maximum);

    set( hObject, 'Value',val );
    set( hObject, 'max',maximum );
    set( this.ax, 'CLim',[minimum val] );

    end %FUNCTION sldIntensity_Callback



    % --- Executes on slider movement.
    function sldScrub_Callback(this, varargin)
    % Allows user to scroll through the movie

    % Stop any currently playing movies.
    set(this.btnPlay,'String','Play');

    % Read frame of the movie
    idxFrame = this.curFrame;
    allgeo = true( size(this.params.geometry) );
    fields = subfield( this.movie, allgeo, idxFrame );

    for f=1:numel(fields)
        field = single(fields{f}) - this.background{f};
        set( this.hImg(f), 'CData',field );
    end
    set( this.edTime, 'String',[num2str(this.curTime) ' s'] );

    end %FUNCTION sldScrub_Callback



    % --- Executes on button press in btnPlay.
    function btnPlay_Callback(this, hObject, varargin)
    % Play movie
    % FIXME: why not call sldScrub_Callback instead?

    % Clicking when the button is labeled 'stop' causes the loop below to terminate.
    if strcmpi( get(hObject,'String') ,'Play' )
        set(hObject,'String','Stop');
    else
        set(hObject,'String','Play');
        return;
    end

    startFrame = this.curFrame;
    allgeo = true( size(this.params.geometry) );

    for i=startFrame:this.movie.nFrames
        fields = subfield( this.movie, allgeo, i );
        
        for f=1:numel(fields)
            field = single(fields{f}) - this.background{f};
            set( this.hImg(f), 'CData',field );
        end

        set( this.edTime, 'String',sprintf('%.2f s',this.curTime) );
        set(this.sldScrub,'Value',i);
        drawnow;

        % Terminate early if the user clicks the 'Stop' button or closes window
        if ~ishandle(hObject) || strcmpi( get(hObject,'String'), 'Play' )
            return;
        end
    end

    set(hObject,'String','Play');
    
    end %FUNCTION btnPlay_Callback
    
    
    
    %% ---------------------- Dependent Properties ---------------------- %%
    
    function t = get.curTime(this)
        t = this.movie.timeAxis(this.curFrame)/1000;
    end
    
    function f = get.curFrame(this)
        f = round( get(this.sldScrub,'Value') );
    end
    
    
end %methods



end %classdef


