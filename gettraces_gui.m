function varargout = gettraces_gui(varargin)
% GETTRACES M-file for gettraces.fig
%      
%      Converts Metamorph .stk files into single-molecule fluorescence and
%      FRET traces. The picking algorithm was taken from gui3.m, written by
%      Harold Kim. I fixed problems with the batch mode, and reconfigured
%      the startup to open .stks instead of .pmas. I also reconfigured the
%      output, so that now the result is a new .traces file with a unique 
%      id for each molecule donor fluorescence, acceptor fluorescence, and
%      fret for each molecule.
%
%      To run in batch mode, first load a single .stk file and adjust the
%      threshold and txtOverlap rejection. Then select batch mode and choose a
%      directory. All the .stk files in that directory will be converted to
%      separate .traces files using the established threshold and txtOverlap
%      rejection. Batch mode is not recursive.
% 

% Last Modified by GUIDE v2.5 10-Jul-2014 15:44:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gettraces_OpeningFcn, ...
                   'gui_OutputFcn',  @gettraces_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --------------------- GUI INITIALIZATION --------------------- %
% --- Executes just before gettraces is made visible.
function gettraces_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gettraces (see VARARGIN)



% Load colormap for image viewer
fid=fopen('colortable.txt','r');
colortable = fscanf(fid,'%d',[3 256]);
handles.colortable = colortable'/255;
fclose(fid);


% Initialize GUI if gettraces is being launched for the first time.
if ~isfield(handles,'params')
    % Choose default command line output for gettraces
    handles.output = hObject;
    
    % Setup dropdown list of configuration profiles.
    constants = cascadeConstants();
    set( handles.cboGeometry, 'String',{constants.gettraces_profiles.name}, ...
                              'Value',constants.gettraces_defaultProfile );
    
    % Setup default values for parameter values -- 2-color FRET.
    params = constants.gettraces_profiles(constants.gettraces_defaultProfile);
    handles.alignment = [];  %current alignment parameters (status)
    handles.params = params;
    
    % Set up GUI elements to reflect the internal parameter values.
    handles = cboGeometry_Callback(handles.cboGeometry, [], handles);
end

% Update handles structure
guidata(hObject, handles);

% gettraces may be called from sorttraces to load the movie associated with
% a particular trace. The first argument is then the filename of the movie
% file and the second argument is the x-y coordinate of the trace.
% FIXME: this still doesn't work very well. And we don't have enough
% information to determine which profile to load, so we default to
% something simple that will work (1 big field). We may need more from
% fileMetadata to figure it out.
if numel(varargin) > 0,
    handles.stkfile = varargin{1};
    traceMetadata = varargin{2};
    
    % Get the coordinates for all of the fluorescence channels.
    fields = fieldnames(traceMetadata);
    xs = find(  ~cellfun( @isempty, strfind(fields,'_x') )  );
    ys = find(  ~cellfun( @isempty, strfind(fields,'_y') )  );
    
    x = zeros(0,2);  y = zeros(0,2);
    
    for i=1:numel(xs),
        x = [ x  traceMetadata.(fields{xs(i)}) ];
        y = [ y  traceMetadata.(fields{ys(i)}) ];
    end
    
    handles.total_x = x;  handles.total_y = y;
    handles.rtotal_x = zeros(0,2);
    handles.rtotal_y = zeros(0,2);
    
    % Choose the appropriate profile and configure the GUI.
    handles.params.geometry = 1;
    set( handles.cboGeometry, 'Value', 1 );
    handles = cboGeometry_Callback(handles.cboGeometry,[],handles);
    
    handles.num = length(x)/numel(xs);
    set(handles.nummoles,'String',num2str(handles.num));
    
    % Load stk file
    handles = OpenStk( handles.stkfile, handles, hObject );
    set(handles.getTraces,'Enable','on');
    
    % Highlight the selected trace.
    highlightPeaks(handles);
end


% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = gettraces_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------- OPEN SINGLE MOVIE (CALLBACK) ---------------- %

% --- Executes on button press in openstk.
function openstk_Callback(hObject, ~, handles)
% hObject    handle to openstk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get filename of input data from user. If multiple files are selected,
% they are considered sections (groups of frames) from a larger movie.
% This happens with very large (sCMOS) movies > 2GB in size.
[datafile,datapath] = uigetfile( '*.stk;*.tif*', 'Choose a movie file', ...
                                 'MultiSelect','on' );

if ~iscell(datafile),
    if datafile==0, return; end  %user hit cancel
end

handles.stkfile = strcat(datapath,datafile);

% Load the movie
handles = OpenStk( handles.stkfile, handles, hObject );

% Update GUI now that data has been loaded.
guidata(hObject,handles);



% --------------------- OPEN SINGLE MOVIE --------------------- %
function handles = OpenStk(filename, handles, hObject)
% Load movie data from file. filename is a cell array of files.

if isempty(filename), return; end
if ~iscell(filename), filename = {filename}; end

[p,f,e] = fileparts( filename{1} );
fnameText = fullfile(p, [f e]);

% Trancate name if too long ot fit into window without wrapping.
if numel(fnameText)>90,
    fnameText = ['...' fnameText(end-90:end)];
end

if numel(filename)>1,
    fnameText = [fnameText ' (multiple files)'];
end

set( handles.txtFilename,          'String',fnameText);
set( handles.txtOverlapStatus,     'String', '' );
set( handles.txtIntegrationStatus, 'String', '' );
set( handles.txtPSFWidth,          'String', '' );
set( handles.txtAlignStatus,      'String', '' );
set( handles.txtAlignWarning, 'Visible','off' );
set(handles.nummoles,'String','');

set(gcf,'pointer','watch'); drawnow;

% Clear the original stack to save memory
if isappdata(handles.figure1,'stkData')
    rmappdata(handles.figure1,'stkData');
end

% Load movie data
[stkData] = gettraces( filename );
handles.stk_top = stkData.stk_top;

% Since the image stack is very large, it is stored in ApplicationData
% instead of GUIData for memory efficiency
setappdata(handles.figure1,'stkData', stkData);


% Setup slider bar (adjusting maximum value in image, initially 2x max)
% low = min(min(handles.stk_top));
low=0;
sort_px = sort(handles.stk_top);
high = sort_px( floor(0.99*numel(sort_px)) );
high = min( ceil(high*2), 32000 );
val = (low+high)/3;

set(handles.scaleSlider,'min',low);
set(handles.scaleSlider,'max',high); %uint16 maxmimum value
set(handles.scaleSlider,'value', val);
set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));

%
image_t = handles.stk_top-stkData.background;
[nrow,ncol] = size(image_t);
fields = {};  ax = [];

chNames = upperFirst( handles.params.chNames );
chNames = strrep(chNames, 'Factor','Factor Binding'); %for display only!


%---- Show fields for Single-Channel (full-chip) recordings.
if handles.params.geometry==1,
    handles.total_t = image_t;
    
    % Show full field of view.
    handles.axTotal = subplot( 1,1,1, 'Parent',handles.panView, 'Position',[0.2 0 0.6 0.95] );
    
    
%---- Show fields for Dual-Channel (half-chip, L/R) recordings.
elseif handles.params.geometry==2,    
    fields{1} = image_t(:,1:ncol/2);  %left, donor
    fields{2} = image_t(:,(ncol/2)+1:end);  %right, acceptor
    handles.total_t = fields{1}+fields{2};

    % Setup axes
    handles.axDonor    = subplot( 1,3,1, 'Parent',handles.panView, 'Position',[0.025 0 0.3 0.95] );
    handles.axAcceptor = subplot( 1,3,2, 'Parent',handles.panView, 'Position',[0.350 0 0.3 0.95] );
    ax = [handles.axDonor,handles.axAcceptor];
        
    % Show total intensity image
    handles.axTotal = subplot( 1,3,3, 'Parent',handles.panView, 'Position',[0.675 0 0.3 0.95] );
    
    
%---- Show fields for Quad-Channel (quarter-chip, L/R/U/D) recordings.
elseif handles.params.geometry>2,
    % Split up channels and combine, assuming perfect alignment.
    fields{1} = image_t( 1:nrow/2, 1:ncol/2 );              %upperLeft
    fields{2} = image_t( 1:nrow/2, (ncol/2)+1:end );        %upperRight
    fields{3} = image_t( (nrow/2)+1:end, 1:ncol/2 );        %lowerLeft
    fields{4} = image_t( (nrow/2)+1:end, (ncol/2)+1:end );  %lowerRight
    handles.total_t = fields{1} + fields{2} + fields{3} + fields{4};
    
    % Setup axes
    handles.axUL = subplot( 2,3,1, 'Parent',handles.panView, 'Position',[0.025 0.5   0.25 0.5] );
    handles.axUR = subplot( 2,3,2, 'Parent',handles.panView, 'Position',[0.35  0.5   0.25 0.5] );
    handles.axLL = subplot( 2,3,4, 'Parent',handles.panView, 'Position',[0.025 0.025 0.25 0.5] );
    handles.axLR = subplot( 2,3,5, 'Parent',handles.panView, 'Position',[0.35  0.025 0.25 0.5] );
    ax = [handles.axUL handles.axUR handles.axLL handles.axLR];

    % Also plot the total intensity image.
    handles.axTotal = subplot( 2,3,3, 'Parent',handles.panView, 'Position',[0.65 0.25 0.25 0.5] );
end


% Show fluorescence fields for all channels
chColors = Wavelength_to_RGB(handles.params.wavelengths);

for i=1:numel(fields),
    % i is the index of physical CCD chip locations.
    % idxCh is the corresponding index into list of channels (there may be none).
    idxCh = find( handles.params.idxFields==i ); 
        
    imshow( fields{i}, [low*2 (high+low)], 'Parent',ax(i) );
    colormap(handles.colortable);  zoom on;

    if ~isempty(idxCh) && ~isempty(chNames{idxCh}),
        % Give each field a title with the background color matching the
        % wavelength of that channel.        
        h = title( ax(i), [chNames{idxCh} ' (' handles.params.chDesc{idxCh} ') #' num2str(i)], ...
                                 'BackgroundColor',chColors(idxCh,:) );

        % Use white text for very dark background colors.
        if sum(chColors(idxCh,:)) < 1,
            set(h,'Color',[1 1 1]);
        end
    end
end

% Link axes so zooming one zooms all.
linkaxes( [ax handles.axTotal] );

% Show total fluorescence channel
imshow( handles.total_t, [low*2 (high+low)], 'Parent',handles.axTotal );
colormap(handles.colortable);  zoom on;
title(handles.axTotal,'Total Intensity');

% Finish up
set(gcf,'pointer','arrow');
set(handles.getTraces,'Enable','on');
set(handles.btnMetadata,'Enable','on');
guidata(hObject,handles);


%end function OpenStk


function names = upperFirst(names)
% Make the first letter of a string uppercase

if ~iscell(names)
    names = {names};
end

for i=1:numel(names),
    if ~isempty( names{i} ),
        names{i} = [upper(names{i}(1)) names{i}(2:end)];
    end
end
    
% end


% --------------- OPEN ALL STK'S IN DIRECTORY (BATCH) --------------- %

% --- Executes on button press in batchmode.
function batchmode_Callback(hObject, ~, handles)
% hObject    handle to batchmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get input parameter values
skipExisting = get(handles.chkOverwrite,'Value');
recursive = get(handles.chkRecursive,'Value');


% Get location of files for gettraces to process
direct=uigetdir('','Choose directory:');
if direct==0, return; end
disp(direct);

% Get list of files in current directory (option: and all subdirectories)
if recursive
    movieFilenames  = rdir( [direct filesep '**' filesep '*.stk'] );
    movieFilenames  = [ movieFilenames ; rdir([direct filesep '**' filesep '*.tif*']) ];
else
    movieFilenames  = rdir( [direct filesep '*.stk'] );
    movieFilenames  = [ movieFilenames ; rdir([direct filesep '*.tif*']) ];
end

nFiles = length(movieFilenames);



% ---- For each file in the user-selected directory

nTraces  = zeros(nFiles,1); % number of peaks found in file X
existing = zeros(nFiles,1); % true if file was skipped (traces file exists)

% Show progress information
% h = waitbar(0,'Extracting traces from movies...');
set(handles.txtProgress,'String','Creating traces, please wait...');

% For each file...
for i=1:nFiles
    stk_fname = movieFilenames(i).name;
    handles.stkfile = stk_fname;
    
    % Skip if previously processed (.traces file exists)
    [p,name] = fileparts(stk_fname);
    traceFname = fullfile(p, [name '.rawtraces']);
    
    if skipExisting && exist(traceFname,'file'),
        disp( ['Skipping (already processed): ' stk_fname] );
        existing(i) = 1;
        continue;
    end
    
    % Load STK file
    handles = OpenStk(handles.stkfile,handles, hObject);
    
    % Pick molecules using default parameter values
    handles = getTraces_Callback(hObject, [], handles);
    
    % Save the traces to file
    saveTraces_Callback(hObject, [], handles);
    
    % Update progress information
    text = sprintf('Creating traces: %.0f%%', 100*(i/nFiles) );
    set(handles.txtProgress,'String',text);
    nTraces(i) = handles.num;
    
    guidata(hObject,handles);
%     waitbar(i/nFiles, h);
end
% close(h);



% ----- Create log file with results
log_fid = fopen( fullfile(direct,'gettraces.log'), 'w' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');

output = evalc('disp(handles.params)');
fprintf(log_fid, '%s', output);
%FIXME: structure parameters are not displayed here (alignment!)

% Log list of files processed by gettraces
fprintf(log_fid,'\n%s\n\n%s\n%s\n\n%s\n',date,'DIRECTORY',direct,'FILES');

for i=1:nFiles
    if existing(i),
        fprintf(log_fid, 'SKIP %s\n', movieFilenames(i).name);
    else
        fprintf(log_fid, '%.0f %s\n', nTraces(i), movieFilenames(i).name);
    end
end

fclose(log_fid);



% ----- Update GUI
set(handles.txtProgress,'String','Finished.');
guidata(hObject,handles);







% ------------------------ PICK INTENSITY PEAKS ------------------------ %

function handles = getTraces_Callback(hObject, ~, handles)
%

%----- Find peak locations from total intensity

% Do nothing if no file has been loaded. This function may be triggered by
% changing the settings fields before a file is open.
if ~isfield(handles, 'stkfile')
    return;
end

set(gcf,'pointer','watch'); drawnow;

% Locate single molecules
stkData = getappdata(handles.figure1,'stkData');
% params = handles.params;
% params.alignment = handles.alignment; %apply loaded alignment if any
[stkData,peaks] = gettraces( stkData, handles.params );
% axes( handles.axTotal );

% The alignment may involve shifting (or distorting) the fields to get a
% registered donor+acceptor field. Show this distorted imaged so the user
% can see what the algorithm is doing.
% axes(handles.axTotal);
val = get(handles.scaleSlider,'value');
minimum = get(handles.scaleSlider,'min');
val = max(val,minimum+1);

imshow( stkData.total_t, [minimum*2 val*2], 'Parent',handles.axTotal );
colormap( handles.axTotal, handles.colortable );
title( handles.axTotal, 'Total Intensity' );


% If no alignment data given (for example in single-channel recordings),
% don't display any status messages.
set( handles.txtAlignWarning, 'Visible','off' );

if ~isfield(stkData,'alignStatus') || isempty(stkData.alignStatus),
    set( handles.txtAlignStatus, 'String', '' );
    handles.alignment = [];
    
% Display alignment status to inform user if realignment may be needed.
% Format: translation deviation (x, y), absolute deviation (x, y)
else
    a = stkData.alignStatus;
    handles.alignment = a;
    
    if isfield(a,'quality'),
        text = 'Alignment applied:';
    else
        text = 'Alignment deviation:';
    end
    
    for i=1:numel(a),
        if isempty(a(i).theta),
           continue;   %ignore donor field alignment to itself
        end
        
        text = [text sprintf('\n%0.1f (x), %0.1f (y), %0.2f (rot), %.1f%% (s), %0.1f (dev)', ...
                       [a(i).dx a(i).dy a(i).theta 100*(a(i).sx-1) a(i).abs_dev] )  ];
    end
    
    % If the alignment quality (confidence) is low, warn the user.
    if isfield(a,'quality'),
        if any( [a.quality]<1.1 & [a.quality]>0 ),
            text = [text sprintf('\nLow quality alignment!')];
        end
    end
    
    set( handles.txtAlignStatus, 'String', text );

    % Color the text to draw attention to it if the alignment is bad.
    if any( [a.abs_dev] > 0.25 ),
        d = max( [a.abs_dev] );
        set( handles.txtAlignStatus, 'ForegroundColor', [(3/2)*min(2/3,d) 0 0] );
    else
        set( handles.txtAlignStatus, 'ForegroundColor', [0 0 0] );
    end
    
    % Total misalignment (no corresponding peaks) can give a relatively low
    % alignment since the algorithm can only search a 1px neighborhood. So
    % when things get close to 1, we need a big warning. 
    if any( [a.abs_dev] >=0.7 ),
        % Put up some big red text that's hard to ignore in the total
        % fluorescence image. FIXME: the warning should be different if the
        % software alignment was applied and the residual is low.
        set( handles.txtAlignWarning, 'Visible','on' );
    end
end


% Get locations also without overlap rejection to estimate the number of
% molecules that are overlapping. This can be used to give the user a
% warning if the density is too high (here by showing it in red).
percentOverlap = stkData.fractionOverlapped*100;

set(  handles.txtOverlapStatus, 'String', ...
      sprintf('%0.0f%% molecules overlapped', percentOverlap)  );

if percentOverlap>30,
    set( handles.txtOverlapStatus, 'ForegroundColor', [0.9 0 0] );
else
    set( handles.txtOverlapStatus, 'ForegroundColor', [0 0 0] );
end


% Get (approximate) average fraction of fluorescence collected within the
% integration window of each molecule. Set the text color to red where the
% intensity is not well collected at the current integration window size.
eff = 100*stkData.integrationEfficiency(:,handles.params.nPixelsToSum);
eff = mean(eff);
set(  handles.txtIntegrationStatus, 'String', ...
      sprintf('%0.0f%% intensity collected', eff)  );

if eff<70,
    set( handles.txtIntegrationStatus, 'ForegroundColor', [0.9 0 0] );
else
    set( handles.txtIntegrationStatus, 'ForegroundColor', [0 0 0] );
end


% Estimate the peak width from pixel intensity distribution.
eff = stkData.integrationEfficiency;
decay = zeros( size(eff,1), 1 ); %number pixels to integrate to get 70% intensity integrated.

for i=1:size(eff,1),
    decay(i) = find( eff(i,:)>=0.7, 1, 'first' );%default 3-color channel assignments.
end

set(  handles.txtPSFWidth, 'String', ...
                         sprintf('PSF size: %0.1f px', mean(decay))  );
                     
if mean(decay) > handles.params.nPixelsToSum,
    set( handles.txtPSFWidth, 'ForegroundColor', [0.9 0 0] );
else
    set( handles.txtPSFWidth, 'ForegroundColor', [0 0 0] );
end


%----- Graphically show peak centers
handles.x = peaks(:,1);
handles.y = peaks(:,2);
handles.total_x = stkData.total_peaks(:,1);
handles.total_y = stkData.total_peaks(:,2);

handles.rx = stkData.rejectedPicks(:,1);
handles.ry = stkData.rejectedPicks(:,2);
handles.rtotal_x = stkData.rejectedTotalPicks(:,1);
handles.rtotal_y = stkData.rejectedTotalPicks(:,2);

highlightPeaks( handles );

% Update GUI controls
handles.num = numel(handles.total_x);
set( handles.nummoles, 'String', sprintf('%d (of %d)',handles.num, ...
                                    numel(handles.rtotal_x)+handles.num) );

set(handles.saveTraces,'Enable','on');

set(gcf,'pointer','arrow');
guidata(hObject,handles);

% end function



function highlightPeaks(handles)
% Draw circles around each selected fluorescence spot, defined in handles.x
% and handles.y. FIXME: using line() to draw circles may be inefficient.
% Try other methods. The big thing is EraseMode.
%

    
% Parameters for drawing circles around each detected molecule.
style1 = {'LineStyle','none','marker','o','color','w'};
style2 = {'LineStyle','none','marker','o','color','y'};

% Parameters for overlapped molecules that won't actually be recorded.
% These are darker to de-emphasize them, but yet still be visible.
style1b = {'LineStyle','none','marker','o','color',[0.4,0.4,0.4]};
style2b = {'LineStyle','none','marker','o','color',[0.4,0.4,0.0]};


[nrow,ncol] = size(handles.stk_top);

% Clear any existing selection markers from previous calls.
delete(findobj(gcf,'type','line'));

    
if handles.params.geometry==2, %dual-channel
    % Assign each spot to the corresponding field image.
    indL  = find( handles.x<=(ncol/2) );
    indR  = find( handles.x> (ncol/2) );
    rindL = find( handles.rx<=(ncol/2) );
    rindR = find( handles.rx> (ncol/2) );
    
    % Draw markers on selection points (donor side)
    axes(handles.axDonor);
    line( handles.x(indL),   handles.y(indL),   style1{:}  );
    line( handles.rx(rindL), handles.ry(rindL), style1b{:} );

    % Draw markers on selection points (acceptor side)
    axes(handles.axAcceptor);
    line( handles.x(indR)-(ncol/2),   handles.y(indR),  style1{:}   );
    line( handles.rx(rindR)-(ncol/2), handles.ry(rindR), style1b{:} );
    
elseif handles.params.geometry>2, %quad-channel
    % Assign each spot to a quadrant. FIXME: can this be a loop?
    indUL = find( handles.x<=(ncol/2) & handles.y<=(nrow/2) );
    indLL = find( handles.x<=(ncol/2) & handles.y> (nrow/2) );
    indLR = find( handles.x> (ncol/2) & handles.y> (nrow/2) );
    indUR = find( handles.x> (ncol/2) & handles.y<=(nrow/2) );
    
    rindUL = find( handles.rx<=(ncol/2) & handles.ry<=(nrow/2) );
    rindLL = find( handles.rx<=(ncol/2) & handles.ry> (nrow/2) );
    rindLR = find( handles.rx> (ncol/2) & handles.ry> (nrow/2) );
    rindUR = find( handles.rx> (ncol/2) & handles.ry<=(nrow/2) );
    
    axes(handles.axUL);
    line( handles.x(indUL),   handles.y(indUL),   style1{:}  );
    line( handles.rx(rindUL), handles.ry(rindUL), style1b{:} );
    
    axes(handles.axLL);
    line( handles.x(indLL),   handles.y(indLL)-(nrow/2),   style1{:}  );
    line( handles.rx(rindLL), handles.ry(rindLL)-(nrow/2), style1b{:} );
    
    axes(handles.axLR);
    line( handles.x(indLR)-(ncol/2),   handles.y(indLR)-(nrow/2),   style1{:}  );
    line( handles.rx(rindLR)-(ncol/2), handles.ry(rindLR)-(nrow/2), style1b{:} );
    
    axes(handles.axUR);
    line( handles.x(indUR)-(ncol/2),   handles.y(indUR),   style1{:}  );
    line( handles.rx(rindUR)-(ncol/2), handles.ry(rindUR), style1b{:} );
end

% Draw markers on selection points (total intensity composite image).
axes(handles.axTotal);
line( handles.total_x,  handles.total_y,  style2{:}  );
line( handles.rtotal_x, handles.rtotal_y, style2b{:} );


% end function highlightPeaks




% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

% --- Executes on button press in saveTraces.
function saveTraces_Callback(~, ~, handles)

if iscell(handles.stkfile),
    filename = handles.stkfile{1};
else
    filename = handles.stkfile;
end

% Integrate fluorophore point-spread functions, generate fluorescence
% traces, and save to file.
stkData = getappdata(handles.figure1,'stkData');
gettraces( stkData, handles.params, filename );
clear stkData;






% --------------------- MISC. GUI CALLBACK FUNCTIONS --------------------- %

% --- Executes on slider movement.
function scaleSlider_Callback(hObject, ~, handles)
% Update axes color limits from new slider value
val = get(hObject,'value');
minimum = get(hObject,'min');

val = max(val,minimum+1); %prevent errors in GUI

if handles.params.geometry==1, %Single-channel recordings
    set( handles.axTotal,    'CLim',[minimum val] );
    
elseif handles.params.geometry==2, %Dual-channel recordings
    set( handles.axDonor,    'CLim',[minimum val] );
    set( handles.axAcceptor, 'CLim',[minimum val] );
    set( handles.axTotal,    'CLim',[minimum*2 val*2] );
elseif handles.params.geometry>2,
    set( handles.axUL, 'CLim',[minimum val] );
    set( handles.axUR, 'CLim',[minimum val] );
    set( handles.axLL, 'CLim',[minimum val] );
    set( handles.axLR, 'CLim',[minimum val] );
    set( handles.axTotal, 'CLim',[minimum*2 val*2] );
end

set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));

guidata(hObject,handles);


% --- Peak selection total intensity threshold specification
function txtIntensityThreshold_Callback(hObject, ~, handles)
% Update gettraces parameters using specified values
text = get(hObject,'String');
if ~isempty( text )
    handles.params.don_thresh = str2double(text);
elseif isfield(handles.params,'don_thresh');
    handles.params = rmfield( handles.params,'don_thresh' );
end

% Re-pick molecules with new settings.
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);


% --- Overlap rejection threshold specification
function txtOverlap_Callback(hObject, ~, handles)
% Update gettraces parameters using specified values
handles.params.overlap_thresh = 0;

text = get(hObject,'String');
if ~isempty( text )
    handles.params.overlap_thresh = str2double(text);
end

% Re-pick molecules with new settings.
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);


% --- Integration window size specification
function txtIntegrationWindow_Callback(hObject, ~, handles)
% Update gettraces parameters using specified values
text = get(hObject,'String');
if ~isempty( text )
    handles.params.nPixelsToSum = floor( str2double(text) );
elseif isfield(handles.params,'nPixelsToSum');
    handles.params = rmfield( handles.params,'nPixelsToSum' );
end

% Re-pick molecules with new settings.
% This is only really necessary to update the GUI status (% intensity
% collected).
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);




% --- Executes on button press in chkSaveLocations.
function chkSaveLocations_Callback(hObject, ~, handles)
% Update gettraces parameters using specified values
handles.params.saveLocations = get(handles.chkSaveLocations,'Value');
guidata(hObject,handles);




function txtMaxIntensity_Callback(hObject, ~, handles)
% Update axes color limits from new slider value
val = str2double( get(hObject,'String') );
minimum = get(handles.scaleSlider,'min');
maximum = get(handles.scaleSlider,'max');

val = max(val,minimum+1); %prevent errors in GUI
maximum = max(val,maximum);

set( handles.scaleSlider, 'Value',val );
set( handles.scaleSlider, 'max',maximum );

if handles.params.geometry==1, %Single-channel recordings
    set( handles.axTotal,    'CLim',[minimum val] );
    
elseif handles.params.geometry==2, %Dual-channel recordings
    set( handles.axDonor,    'CLim',[minimum val] );
    set( handles.axAcceptor, 'CLim',[minimum val] );
    set( handles.axTotal,    'CLim',[minimum*2 val*2] );
    
elseif handles.params.geometry>2,
    set( handles.axUL, 'CLim',[minimum val] );
    set( handles.axUR, 'CLim',[minimum val] );
    set( handles.axLL, 'CLim',[minimum val] );
    set( handles.axLR, 'CLim',[minimum val] );
    set( handles.axTotal, 'CLim',[minimum*2 val*2] );
end

set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));

guidata(hObject,handles);




% --- Executes on selection change in cboGeometry.
function handles = cboGeometry_Callback(hObject, ~, handles)
%

% Get parameter values associated with the selected profile.
% Warning: if cascadeConstants is changed to add a new profile or rearrange
% profiles, this can have unpredictable effects...
constants = cascadeConstants;
sel = get(hObject,'Value');
params = constants.gettraces_profiles(sel);
handles.params = params;


% Reset all GUI elements to reflect the default settings in the currently
% selected profile.
if ~isfield(params,'don_thresh') || params.don_thresh==0,
    set( handles.txtIntensityThreshold,'String','' );
else
    set( handles.txtIntensityThreshold,'String',num2str(params.don_thresh) );
end

set( handles.txtOverlap,           'String', num2str(params.overlap_thresh)   );
set( handles.txtIntegrationWindow, 'String', num2str(params.nPixelsToSum)     );
set( handles.txtPhotonConversion,  'String', num2str(params.photonConversion) );
set( handles.chkRecursive, 'Value', params.recursive    );
set( handles.chkOverwrite, 'Value', params.skipExisting );
set( handles.cboAlignMethod, 'Value',1 ); %params.alignment.method

if handles.params.geometry==1, %Single-channel recordings
    set( handles.txtDACrosstalk,    'Enable','off', 'String','', 'Visible','on' );
    set( handles.btnCrosstalk,      'Visible','off'  );
    set( handles.btnSaveAlignment,  'Enable','off' );
    set( handles.btnLoadAlignment,  'Enable','off', 'Value',0 );
else  %Multi-channel recordings
    set( handles.btnSaveAlignment,  'Enable','on' );
    set( handles.btnLoadAlignment,  'Enable','on', 'Value',0 );
end

% For multi-color FRET, the matrix of possible crosstalk values can't
% be handled in a little text box, so use a button for a dialog with
% all of the values listed instead (see btnCrosstalk_Callback).
if handles.params.geometry>1 && numel(params.crosstalk)==1,
    set( handles.txtDACrosstalk, 'Enable','on', 'Visible','on', ...
                                    'String',num2str(params.crosstalk) );
    set( handles.btnCrosstalk,   'Visible','off'  );
elseif handles.params.geometry>2,
    set( handles.txtDACrosstalk, 'Visible','off' );
    set( handles.btnCrosstalk,   'Visible','on'  );
end
    

% If a movie has already been loaded, reload movie with new setup.
if isfield(handles,'stkfile'),
    handles = OpenStk( handles.stkfile, handles, hObject );
end

guidata(hObject,handles);




function txtDACrosstalk_Callback(hObject, ~, handles)
% 
handles.params.crosstalk = str2double( get(hObject,'String') );
guidata(hObject,handles);





function txtPhotonConversion_Callback(hObject, ~, handles)
%
handles.params.photonConversion = str2double( get(hObject,'String') );
guidata(hObject,handles);



% --- Executes on button press in btnMetadata.
function btnMetadata_Callback(hObject, ~, handles)
% Display a simple diaglog with the MetaMorph metadata for the first frame.
% 

stkData = getappdata(handles.figure1,'stkData');

if isempty(stkData) || ~isfield(stkData.movie.header,'MM')
    % This movie doesn't have MetaMorph metadata (should never happen).
    set( hObject, 'Enable','off' );
    return;
end

% Grab extra MetaMorph fields outside the MM struct.
metadata = stkData.movie.stkHeader.MM(1);

if isfield( stkData.movie.stkHeader, 'MM_wavelength' ),
    wv = stkData.movie.stkHeader.MM_wavelength;
    metadata.MM_wavelength = wv(wv>100);
end

% Convert structure data into text to display.
fields = fieldnames(metadata);
flen = max( cellfun(@numel,fields) );
nFields = numel(fields);

output = cell( nFields, 1 );

for i=1:nFields, %for each field
    fname = fields{i};
    data  = metadata.(fname);
    
    if isnumeric(data),
        data = num2str(data);
    end
    
    output{i} = sprintf( ['%' num2str(flen) 's:  %s'], fname, data );
end


% Display the dialog.
msgbox( output, 'MetaMorph metadata' );






% --- Executes on button press in btnLoadAlignment.
function btnLoadAlignment_Callback(hObject, ~, handles)
% Load software alignment settings previously saved to file. The file
% the "align" structure defined in gettraces, including dx, dy, theta, etc.
%

% This should only be used in multi-color experiments.
assert( handles.params.geometry>1 );


pressed = get(hObject,'Value');

% Load an alignment file
if pressed == get(hObject,'Max')  %toggle is pressed: load alignment.
    [f,p] = uigetfile('*.mat','Select an alignment settings file');
    alignFilename = [p f];
    if f==0, return; end
    
    try
        % Overwrite alignment settings with those in the file.
        input = load(alignFilename);
        handles.params.alignment = input.alignment;
        handles.alignment = input.alignment; %???
    catch e,
        % If the file is invalid, give a warning and reset the button so
        % that is as if nothing happened.
        disp( ['Invalid alignment file: ' e.message] );
        set( hObject, 'Value',get(hObject,'Min') );
        return;
    end
    
    % 4) Disable alignment controls and set to checked.
    handles.params.alignMethod = 2;
    set( handles.cboAlignMethod, 'Value',2 );
    
    
% Unload the current alignment and reset to the normal state.
elseif pressed == get(hObject,'Min')    
    % Reset alignment parameters back to defaults.
    constants = cascadeConstants;
    sel = get(handles.cboGeometry,'Value');
    p = constants.gettraces_profiles(sel);
    
    handles.params.alignMethod = p.alignMethod;
    set( handles.cboAlignMethod, 'Value',p.alignMethod );
    
    if isfield(handles.params,'alignment'),
        handles.params = rmfield(handles.params,'alignment');
    end
    handles.alignment = [];
end


% Re-pick molecules with new settings.
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);


%end function btnLoadAlignment_Callback




% --- Executes on button press in btnSaveAlignment.
function btnSaveAlignment_Callback(~, ~, handles)
% Save current software alignment settings (which may be set to do nothing
% at all) to file so they can be reloaded later.
%

assert( isfield(handles,'alignment') && ~isempty(handles.alignment) && handles.params.geometry>1 );

% Verify there is a valid software alignment 
% if ~isfield(handles,'alignment') || isempty(handles.alignment),
%     set(handles.btnSaveAlignment,'Enable','off');
%     return;
% end

[f,p] = uiputfile('*.mat','Save software alignment settings','align.mat');
alignment = rmfield( handles.alignment, {'quality'} );
save( [p f], 'alignment' );


%end function btnSaveAlignment_Callback




% --- Executes on button press in btnCrosstalk.
function btnCrosstalk_Callback(hObject, ~, handles)
% When there are more than 2 channels, the crosstalk is more than just a
% scalar and can't be represented in the text box easily, so this button
% will launch a dialog to show all the possible parameter values and allow
% the user to change them.
%

params = handles.params;

assert( params.geometry>1 && numel(params.crosstalk)>1 );


% Only show crosstalk parameters for channels that are adjacent in
% wavelength space. Crosstalk in other channels is generally negligable.
% This simplifies the input. The channels are listed in the order of
% wavelength.
[wl,idx] = sort(params.wavelengths);

prompts = strcat( params.chNames(idx(1:end-1)'), ' (', params.chDesc(idx(1:end-1)'), ') -> ', ...
                 params.chNames(idx(2:end)'),   ' (', params.chDesc(idx(2:end)'),   '):' );

defaults = cell( numel(wl)-1, 1 );
for i=1:(numel(wl)-1),
    defaults{i} = num2str( params.crosstalk(idx(i),idx(i+1)) );
end

result = inputdlg( prompts, 'Enter crosstalk values', 1, defaults );

if isempty(result),
    return;
end

% Process the user input. We assume they are numbers.
for i=1:(numel(wl)-1),
    c = str2double( result{i} );

    if isnan(c) || c>1 || c<0,
        fprintf( 'Error: invalid crosstalk value %d ignored (%s)', ...
                                                  i, result{i} );
    else
        handles.params.crosstalk( idx(i), idx(i+1) ) = c;
    end
end %for each channel pair.

disp(handles.params.crosstalk);
guidata(hObject,handles);

%end function btnCrosstalk_Callback




% --- Executes on selection change in cboAlignMethod.
function cboAlignMethod_Callback(hObject, ~, handles)
% 

% methods = contents = cellstr(get(hObject,'String'));
sel = get(hObject,'Value');
% method = contents{sel};

handles.params.alignMethod = sel;

% FIXME: if the users selects "load alignment", this should ask for a file.

% Re-pick molecules with new settings.
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);


% END FUNCTION cboAlignMethod_Callback










