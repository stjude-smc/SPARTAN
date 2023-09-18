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
    chExtractor;   % Encapsulates movie data; splits into spectral channels
    subtractBGImage = true;  % Automatically subtract background image
end

properties (SetAccess=protected, GetAccess=public)
    hFig;          % Handle to figure
    ax;            % Handles to subplot axes in montage
    hImg;          % Handle to image viewers
    btnPlay;       % Handle to 'Play' button
    sldScrub;      % Handle to scroll bar for scrubbing through time
    sldIntensity;  % Handle to scroll bar for adjusting intensity scaling
    txtMaxIntensity;  % Handle to text box for adjusting intensity scaling
    edTime;        % Text box showing current time
end


properties (Dependent)
    curFrame;  % Frame number (1..nFrames) being displayed
    curTime;   % Time (in seconds) of frame currently being displayed.
end



methods
    %% ---- Constructor ---- %%
    function this = MovieViewer(input)
    % Create a new movie viewer window.
    
    % Load movie from file. ChannelExtractor will use the movie file's
    % metadata to determine how to split fields into channels. If no
    % metadata available, it will default to single-channel.
    if nargin<1, input=[]; end
    
    if ischar(input) || isa(input,'Movie')
        this.chExtractor = ChannelExtractor(input);
    elseif isa(input,'ChannelExtractor')
        this.chExtractor = input;
    else
        error('Input must be Movie, path to movie file, or ChannelExtractor');
    end
    
    end %constructor
    
    
    % If MovieViewer object is deleted, close the associated window.
    % This is for command line use, where it makes sense, but causes
    % problems for other uses.
%     function close(this)
%         try
%             close(this.hFig);
%             delete(this);
%         catch
%         end
%     end
%     
%     % Destructor
%     function delete(this)
%         try
%             close(this.hFig);
%         catch
%         end
%     end %function delete
    
    
    
    %% ------------------ Display viewer window ------------------ %%
    function varargout = show(this, hPanel, varargin)
    % Create a new figure displaying movie frames.
    
    [varargout{1:nargout}] = deal([]);
    
    % Create a uipanel for axes to sit in.
    % FIXME: handle to panel as argument (for gettraces integration)
    if nargin<2
        this.hFig = figure;
        hPanel = uipanel( this.hFig, 'Position',[0 0 1 1], 'BorderType','none' );
    else
        this.hFig = ancestor(hPanel,'figure');
    end
    colormap(this.hFig, gettraces_colormap);
    delete( findall(hPanel,'type','axes') );
    
    % Use an average of beginning frames for initial image and thresholds
    stk_top = this.chExtractor.stk_top;
    px_max = zeros(numel(stk_top),1);
    
    for i=1:numel(stk_top)
        sort_px = sort( stk_top{i}(:) );
        px_max(i) = sort_px( floor(0.99*numel(sort_px)) );
    end
    
    val = max(px_max);
    high = min( ceil(val*10), 32000 );  %uint16 maxmimum value

    % Create axes for sub-fields (listed in column-major order, like stk_top)
    % Position is: L,B,W,H
    axopt = {'Visible','off', 'Parent',hPanel};
    nCh = numel(stk_top);
    
    switch nCh
    case 1
        this.ax    = axes( 'Position',[0.07 0.06 0.9 0.9], axopt{:} );

    case 2
        this.ax    = axes( 'Position',[0.07  0.06  0.3 0.9], axopt{:} );  %L
        this.ax(2) = axes( 'Position',[0.38  0.06  0.3 0.9], axopt{:} );  %R
        this.ax(3) = axes( 'Position',[0.69  0.06  0.3 0.9], axopt{:} );

    case {3,4}
        this.ax    = axes( 'Position',[0.07  0.51  0.3 0.45], axopt{:} );  %TL
        this.ax(2) = axes( 'Position',[0.07  0.06  0.3 0.45], axopt{:} );  %BL
        this.ax(3) = axes( 'Position',[0.38  0.51  0.3 0.45], axopt{:} );  %TR
        this.ax(4) = axes( 'Position',[0.38  0.06  0.3 0.45], axopt{:} );  %BR
        
        geo = this.chExtractor.fieldArrangement;
        if isempty(geo) || size(geo,3)>1
            % If channels have no spatial coding, use Quad-View convention
            this.ax = this.ax( [2 3; 1 4] );
            if nCh==3
                if this.chExtractor.channels(1).wavelength>500
                    this.ax(1) = [];
                else
                    this.ax(end) = [];
                end
            end
        else
            % If channels are tiled in space, try to keep the same physical
            % layout here as in the actual frames.
            temp = zeros(nCh,1);  %new axes list
            for i=1:numel(geo)
                if geo(i)==0, continue; end
                temp( geo(i) ) = this.ax(i);
            end
            this.ax = temp(temp~=0);  %ax(i) now directly corresponds to stk_top(i)
        end
        
        this.ax(end+1) = axes( 'Position',[0.69  0.25 0.3 0.45], axopt{:} );
        
    otherwise
        error('Invalid field geometry');
    end
    this.ax = to_row(this.ax);

    % Show fluorescence fields for all channels
    % NOTE: all properties are listed in wavelength order.
    axopt = {'YDir','reverse', 'Color',get(this.hFig,'Color'), 'Visible','off'};

    this.hImg = [];
    for i=1:nCh
        this.hImg(i) = image( stk_top{i}, 'CDataMapping','scaled', 'Parent',this.ax(i) );
        set( this.ax(i), 'UserData',i, 'CLim',[0 val], axopt{:} );
    end
    
    % Show total fluorescence channel
    if nCh>1
        total = sum( cat(3,stk_top{:}), 3);
        this.hImg(end+1) = image( total, 'CDataMapping','scaled', 'Parent',this.ax(end) );
        set(this.ax(end), 'CLim',[0 val*2], axopt{:} );
    end
    
%     if abs(log2( size(total,2)/size(total,1) )) < 1
        % For roughly symmetric movies, allows for better window resizing.
        axis( this.ax, 'image' );  %adjust axes size to match image
%     else
%         % Allows for better zooming of narrow (high time resolution) movie.
%         axis( ax, 'equal' );  %keep the axes size fixed.
%     end
    
    setAxTitles(this);
    linkaxes( this.ax );
    if nargout>0, varargout{1}=this.ax; end
    zoom on;

    % Add intensity and time scroll bars and play button. (L,B,W,H)
    % Don't re-create if left from a previous instance of MovieViewer.
    if isempty(this.sldIntensity)
        style = {'Parent',hPanel, 'units','normalized', 'style'};
        this.sldIntensity = uicontrol(style{:}, 'slider', 'Position',...
               [0.03 0.15 0.035 0.78], 'Callback',@this.sldIntensity_Callback);

        this.sldScrub = uicontrol(style{:}, 'slider', 'position', ...
                [0.185 0.05 0.65 .05], 'callback',@this.sldScrub_Callback);

        this.edTime = uicontrol(style{:},'edit', 'Enable','inactive', ...
                'Position',[0.85 0.05 0.1 0.05], 'String','0 s');

        this.btnPlay = uicontrol(style{:},'pushbutton', 'String','Play', ...
                'Position',[0.08 0.05 0.08 .05], 'Callback',@this.btnPlay_Callback );

        this.txtMaxIntensity = uicontrol(style{:},'edit', ...
                'Position',[0.01 0.93 0.08 0.05], 'Callback',@this.txtMaxIntensity_Callback );
    end
    set( this.sldIntensity,'min',0, 'max',high, 'value',val);
    set( this.txtMaxIntensity, 'String',sprintf('%.0f',val) );
    set( this.sldScrub, 'Min',1, 'Max',this.chExtractor.nFrames, 'Value',1 );
    set( this.sldScrub, 'SliderStep', [1 10]./this.chExtractor.nFrames );
        
    zoom(this.hFig, 'on');
    
    end %function show
    
    
    
    function setAxTitles(this)
    % Set axes titles from imaging profile settings, including colors.

    for i=1:this.chExtractor.nChannels  %i is channel index
        ch = this.chExtractor.channels(i);

        if isempty(ch.name)
            title( this.ax(i), '' );
        else
            if isfield(ch,'role') && ~isempty(ch.role)
                text = sprintf('%s (%s)',ch.name,ch.role);
            else
                text = ch.name;
            end
            
            chColor = Wavelength_to_RGB( double(ch.wavelength) );
            title( this.ax(i), text, 'BackgroundColor',chColor, 'FontSize',10, 'Visible','on', ...
                   'Color',(sum(chColor)<1)*[1 1 1] ); % White text for dark backgrounds.
        end
    end
    title(this.ax(end),'Total Intensity', 'FontSize',10, 'Visible','on');
    
    end %FUNCTION setTitles
    
    
    
    function highlightMolecule(this, coords)
    % Used by sorttraces (showMovie).
    % Draw a circle to highlight a molecule's PSF in each channel.
    % Coords is a cell array, one per channels, with x in first col and y in
    % second, with molecules listed in rows.
    
    if ~all(ishandle(this.ax)), return; end  %window closed?
    delete( findall(this.ax,'type','Line') );
    
    for i=1:numel(coords)
        x = rem( coords{i}(:,1), this.chExtractor.nX);
        y = rem( coords{i}(:,2), this.chExtractor.nY);
        
        % Translate coordinates from stitched movie to subfield
        viscircles( this.ax(i), [x y], 3, 'EdgeColor','w' );
    end
    
    % If zoomed in, recenter on around new peak.
    dx = diff(get(this.ax(end), 'XLim'))/2;
    dy = diff(get(this.ax(end), 'YLim'))/2;
    if dx<this.chExtractor.nX/2 || dy<this.chExtractor.nY/2
        set( this.ax(end), 'XLim',[x-dx x+dx], 'YLim',[y-dy y+dy] );
    end
    
    end %function highlightMolecule
    
    
    
    function highlightPeaks(this, selected, rejected, selectedTotal, rejectedTotal)
    % Used by gettraces (MovieParser).
    % Draw a circle to highlight a molecule's PSF in each channel.
    % Coords is an array of molecules, dimension (x/y), fields.
    
    if ~all(ishandle(this.ax)), return; end  %window closed?
    delete( findall(this.ax,'type','Line') );
    style = {'LineStyle','none','marker','o'};
    
    if size(selected,3)>1
        for i=1:size(selected,3)
            line( selected(:,1,i), selected(:,2,i), ...
                    style{:}, 'color','w', 'Parent',this.ax(i) );
            line( rejected(:,1,i), rejected(:,2,i), ...
                    style{:}, 'color',[0.4,0.4,0.4], 'Parent',this.ax(i) );
        end
    end
    
    % Draw markers on selection points (total intensity composite image).
    if nargin>2
        line( selectedTotal(:,1), selectedTotal(:,2), ...
                style{:}, 'color','y',  'Parent',this.ax(end) );
        line( rejectedTotal(:,1), rejectedTotal(:,2), ...
                style{:}, 'color',[0.4,0.4,0.0], 'Parent',this.ax(end) );
    end
    
    end %function highlightMolecule
    

    
    function hidePeaks(this)
        delete(findobj(this.hFig,'type','line'));
    end %function hidePeaks
    
    
    
    %% ----------------------- Callback Functions ----------------------- %%
    
    % --- Executes on slider movement.
    function sldIntensity_Callback(this, hObject, varargin)
    % Update axes color limits from new slider value

    val = get(hObject,'value');
    set(this.txtMaxIntensity,'String', sprintf('%.0f',val));
    this.txtMaxIntensity_Callback(this.txtMaxIntensity);

    end %FUNCTION sldIntensity_Callback
    
    
    
    % --- Executes on intensity text box change
    function txtMaxIntensity_Callback(this, hObject, varargin)
    % Update axes color limits from new slider value

    val = str2double( get(hObject,'String') );
    minimum = get(this.sldIntensity,'min');
    maximum = get(this.sldIntensity,'max');

    val = max(val,minimum+1); %prevent errors in GUI
    maximum = max(val,maximum);

    set( this.sldIntensity, 'Value',val );
    set( this.sldIntensity, 'max',maximum );
    set( this.ax, 'CLim',[minimum val] );
    
    if numel(this.ax)>1
        set( this.ax(end), 'CLim',[minimum*2 val*2] );
    end
    
    set(this.txtMaxIntensity,'String', sprintf('%.0f',val));

    end %FUNCTION sldIntensity_Callback



    % --- Executes on slider movement.
    function sldScrub_Callback(this, hObject, varargin)
    % Allows user to scroll through the movie

    % Stop any currently playing movies.
    set(this.btnPlay,'String','Play');

    % Read frame of the movie
    idx = floor(get(hObject,'Value'));
    fields = this.chExtractor.read(idx);

    for f=1:numel(fields)
        field = single(fields{f});
        if this.subtractBGImage, field = field - this.chExtractor.background{f}; end
        set( this.hImg(f), 'CData',field );
    end
    set( this.edTime, 'String',sprintf('%.2f s',this.chExtractor.timeAxis(idx)/1000) );
    
    % Color the time text by laser wavelength (if metadata available)
    wavelength = this.chExtractor.lasersActive(idx);
    if isempty(wavelength), wavelength=0; end
    set( this.edTime, 'ForegroundColor', 0.8*Wavelength_to_RGB(wavelength(1)) );
    
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
    set( this.edTime, 'ForegroundColor', [0 0 0] );

    startFrame = floor(get(this.sldScrub,'Value'));
    if startFrame==this.chExtractor.nFrames, startFrame=1; end

    for i=startFrame:this.chExtractor.nFrames
        fields = this.chExtractor.read(i);
        
        for f=1:numel(fields)
            field = single(fields{f});
            if this.subtractBGImage, field = field - this.chExtractor.background{f}; end
            set( this.hImg(f), 'CData',field );
        end

        set( this.edTime, 'String',sprintf('%.2f s',this.chExtractor.timeAxis(i)/1000) );
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
        t = this.chExtractor.timeAxis(this.curFrame)/1000;
    end
    
    function f = get.curFrame(this)
        f = round( get(this.sldScrub,'Value') );
    end
    
    function set.curFrame(this,val)
        set(this.sldScrub,'Value',floor(val));
    end
    
    
end %methods



end %classdef


