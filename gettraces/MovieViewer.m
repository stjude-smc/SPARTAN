classdef MovieViewer < handle
% MovieViewer   Display a montage of movies
%
%   
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
%   See also: gettraces, MovieParser, Movie_TIFF, showMovie.

%   Copyright 2017 Cornell University All Rights Reserved.


% size of 'movies' property dictates how it is displayed

%TODO: function suggestions:
% - centerQuad
% - layout options (set number of rows/columns)


%% ---------------------------  PROPERTIES  --------------------------- %%

properties (SetAccess=public, GetAccess=public)
end

properties (SetAccess=protected, GetAccess=public)
    filenames = {};  % Path to each .tif movie file
    movies    = {};  % MovieParser objects
    
    ax;            % Handles to subplot axes in montage
    hImg;          % Handle to image viewers
    btnPlay;       % Handle to 'Play' button
    sldScrub;      % Handle to scroll bar for scrubbing through time
    sldIntensity;  % Handle to scroll bar for adjusting intensity levels
end




methods
    %% ---- Constructor ---- %%
    function this = MovieViewer(varargin)
    % Create a new movie viewer window.
    
    narginchk(0,1); nargoutchk(0,1);
    
    this.load(varargin{:});
    %this.show();
    
    end %constructor
    
    
    
    %% ---- Load movies from file ---- %%
    function load(this, movieInput)
    % Load movies from input
    % FIXME: could be a cell array of Movie objects.
    % FIXME: use Movie as superclass, handle all subtypes generically.
    
    % If not specified, ask the user for input.
    if nargin<2
        movieInput = getFiles('*.tif');
        this.filenames = movieInput;
    end
    
    % Process list of file name from user.
    if iscell(movieInput)
        constants = cascadeConstants;
        this.movies = cell( size(movieInput) );
    
        for i=1:numel(movieInput)
            this.movies{i} = MovieParser( movieInput{i}, constants.gettraces_profiles(1) );
        end

    elseif all(isa(movieInput,'MovieParser'))
        this.movies = movieInput;

%     else
%         error('Invalid movie input');
    end
        
    end %function load
    
    
    
    
    %% ------------------ Display montage viewer window ------------------ %%
    function show(this, varargin)
    % Create figure
    
    % Set arrangement of panels
    % FIXME: gracefully handle empty panels
    if numel(varargin)>0
        this.movies = reshape(this.movies, varargin{:});
    end
    
    
    hFig = figure;
    
    if ~isempty(this.filenames)
        titles = trimtitles(this.filenames);
    end
    
    
    % Get intensity scale
    sort_px = sort( this.movies{1}.stk_top{1}(:) );
    val = sort_px( floor(0.99*numel(sort_px)) );
    high = min( ceil(val*10), 32000 );  %uint16 maxmimum value
    
    
    % Draw movie viewer axes
    [nr,nc] = size(this.movies);
    axopt = {'YDir','reverse', 'Color',get(hFig,'Color'), 'Visible','off'};
    
    for i=1:numel(this.movies)
        if isempty(this.movies{i}), continue; end
        
        [row,col] = ind2sub( size(this.movies), i );
        this.ax(i) = axes( hFig, 'Position',[0.05+0.95*(col-1)/nc  0.15+0.85*(row-1)/nr  (0.9/nc)  (0.85/nr)] );

        frame = single(this.movies{i}.movie.readFrames(1)) - this.movies{i}.background{1};
        this.hImg(i) = image( frame, 'CDataMapping','scaled' );
        set( this.ax(i), 'CLim',[0 val], axopt{:} );
        colormap(this.ax(i), gettraces_colormap);
        
        if numel(this.filenames)>=i
            title( this.ax(i), titles{i}, 'Visible','on' );
        end
    end
    
    linkaxes( to_row(this.ax) );
    
%     if abs(log2( size(total,2)/size(total,1) )) < 1
        % For roughly symmetric movies, allows for better window resizing.
        axis( this.ax, 'image' );  %adjust axes size to match image
%     else
%         % Allows for better zooming of narrow (high time resolution) movie.
%         axis( [ax handles.axTotal], 'equal' );  %keep the axes size fixed.
%     end


    % Add titles



    % Add intensity and time scroll bars and play button.
    sldStyle = {'style','slider','units','normalized'};
    this.sldIntensity = uicontrol(sldStyle{:}, 'position',[0.01 0.1 0.025 0.8], 'Callback',@this.sldIntensity_Callback);
    set(this.sldIntensity,'min',0, 'max',high, 'value',val);
    
    this.sldScrub = uicontrol(sldStyle{:}, 'position',[0.2 0.05 0.7 .05], 'callback',@this.sldScrub_Callback);
    set( this.sldScrub, 'Min',1, 'Max',this.movies{1}.nFrames, 'Value',1, ...
           'SliderStep',[1/this.movies{1}.nFrames,0.02] );
 
    this.btnPlay = uicontrol('style','pushbutton', 'units','normalized', 'String','Play', ...
            'Position',[0.08 0.05 0.08 .05], 'Callback',@this.btnPlay_Callback );
    
    end %function show
    
    
    
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
    function sldScrub_Callback(this, hObject, varargin)
    % Allows user to scroll through the movie

    % Stop any currently playing movies.
    set(this.btnPlay,'String','Play');

    idxFrame = round( get(hObject,'Value') );

    for i=1:numel(this.movies)
        frame = single(this.movies{i}.movie.readFrames(idxFrame)) - this.movies{i}.background{1};
        set( this.hImg(i), 'CData',frame );
    end

    end %FUNCTION sldScrub_Callback



    % --- Executes on button press in btnPlay.
    function btnPlay_Callback(this, hObject, varargin)
    % Play movie

    % Clicking when the button is labeled 'stop' causes the loop below to terminate.
    if strcmpi( get(hObject,'String') ,'Play' )
        set(hObject,'String','Stop');
    else
        set(hObject,'String','Play');
        return;
    end

    startFrame = round( get(this.sldScrub,'Value') );

    for i=startFrame:this.movies{1}.nFrames

        for j=1:numel(this.movies)
            frame = single(this.movies{j}.movie.readFrames(i)) - this.movies{j}.background{1};
            set( this.hImg(j), 'CData',frame );
        end

        set(this.sldScrub,'Value',i);
        drawnow;

        % Terminate early if the user clicks the 'Stop' button.
        if strcmpi( get(hObject,'String'), 'Play' ),  return;  end
    end

    set(hObject,'String','Play');

    end %FUNCTION btnPlay_Callback
    
    
    
    
end %methods




end %classdef


