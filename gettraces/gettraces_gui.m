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

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% Last Modified by GUIDE v2.5 02-Sep-2015 14:55:40


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
handles.colortable = gettraces_colormap();

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
    
    set( handles.figure1, 'Name', ['gettraces (version ' constants.version ')'] );
end

% Update handles structure
guidata(hObject, handles);

% gettraces may be called from sorttraces to load the movie associated with
% a particular trace. The first argument is then the filename of the movie
% file and the second argument is the x-y coordinate of the trace.
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
function openstk_Callback(hObject, ~, handles)  %#ok<DEFNU>
% hObject    handle to openstk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get filename of input data from user. If multiple files are selected,
% they are considered sections (groups of frames) from a larger movie.
% This happens with very large (sCMOS) movies > 2GB in size.
[datafile,datapath] = uigetfile( '*.stk;*.tif;*.tiff', 'Choose a movie file', ...
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

% Remove "-file00x" extension for multi-file TIFF stacks.
[p,f,e] = fileparts( filename{1} );
f = regexprep(f,'-file[0-9]*$','');


% If a single file is selected, look for the others in multi-file TIFFs.
% If the user selected multiple files, we assume that they got all of them.
if numel(filename)==1,
    d = dir( [f '*.tif*'] );
    if numel(d)>1,
        d = regexpi( {d.name}, [f '(-file[0-9]*)?\.tiff?$'], 'match' );
        filename = [d{:}];
    end
end

% Trancate name if too long ot fit into window without wrapping.
fnameText = fullfile(p, [f e]);

if numel(fnameText)>90,
    fnameText = ['...' fnameText(end-90:end)];
end

if numel(filename)>1,
    fnameText = [fnameText ' (multiple files)'];
end

set( handles.txtFilename,          'String',fnameText);
set( handles.txtOverlapStatus,     'String', '' );
set( handles.txtIntegrationStatus, 'String', '' );
set( handles.txtWindowOverlap,     'String', '' );
set( handles.txtPSFWidth,          'String', '' );
set( handles.txtAlignStatus,      'String', '' );
set( handles.txtAlignWarning, 'Visible','off' );
set(handles.nummoles,'String','');

set(handles.figure1,'pointer','watch'); drawnow;

% Clear the original stack to save memory
if isappdata(handles.figure1,'stkData')
    rmappdata(handles.figure1,'stkData');
end

% Load movie data
[stkData] = gettraces( filename, handles.params );
handles.stk_top = stkData.stk_top;

% Since the image stack is very large, it is stored in ApplicationData
% instead of GUIData for memory efficiency
setappdata(handles.figure1,'stkData', stkData);


% Setup slider bar (adjusting maximum value in image, initially 2x max)
% low = min(min(handles.stk_top));
sort_px = sort(handles.stk_top(:));
low=0;
val = sort_px( floor(0.98*numel(sort_px)) );
high = min( ceil(val*10), 32000 );  %uint16 maxmimum value

set(handles.scaleSlider,'min',low);
set(handles.scaleSlider,'max',high);
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
        
    imshow( fields{i}, [low val], 'Parent',ax(i) );

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
imshow( handles.total_t, [low*2 val*2], 'Parent',handles.axTotal );
colormap(handles.axTotal,handles.colortable);
zoom(handles.axTotal,'on');
title(handles.axTotal,'Total Intensity');

% Finish up
set(handles.figure1,'pointer','arrow');
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
function batchmode_Callback(hObject, ~, handles,direct)
% hObject    handle to batchmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% direct     target directory location to look for new files

% Get input parameter values
skipExisting = get(handles.chkOverwrite,'Value');
recursive = get(handles.chkRecursive,'Value');

% Get location of files for gettraces to process
if nargin>=4 && exist(direct,'dir'),
    % Get a fresh copy of handles. The one passed in arguments is an old
    % copy made when the timer was created. Kind of ugly...
    handles = guidata(handles.chkAutoBatch);
else
    direct=uigetdir('','Choose directory:');
    if direct==0, return; end
    disp(direct);
end

% Get list of files in current directory (option: and all subdirectories)
movieFiles = regexpdir(direct,'^.*\.(tiff?|stk)$',recursive);
nFiles = length(movieFiles);

% Wait for 100ms to give sufficient time for polling file sizes in the
% main loop below.
pause(0.1);


% ---- For each file in the user-selected directory
nTraces  = zeros(nFiles,1); % number of peaks found in file X
existing = zeros(nFiles,1); % true if file was skipped (traces file exists)

% Show progress information
% h = waitbar(0,'Extracting traces from movies...');
set(handles.txtProgress,'String','Running...');

% For each file...
for i=1:nFiles
    stk_fname = movieFiles(i).name;
    handles.stkfile = stk_fname;
    
    % Skip if previously processed (.traces file exists)
    [p,name] = fileparts(stk_fname);
    traceFname = fullfile(p, [name '.rawtraces']);
    
    if skipExisting && exist(traceFname,'file'),
        %disp( ['Skipping (already processed): ' stk_fname] );
        existing(i) = 1;
        continue;
    end
    
    % Poll the file to make sure it isn't changing.
    % This could happen when a file is being saved during acquisition.
    d = dir(stk_fname);
    if movieFiles(i).datenum ~= d(1).datenum,
        disp( ['Skipping (save in process?): ' stk_fname] );
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
    text = sprintf('Running: %.0f%%', 100*(i/nFiles) );
    set(handles.txtProgress,'String',text);
    nTraces(i) = handles.num;
    
    guidata(hObject,handles);
%     waitbar(i/nFiles, h);
end
% close(h);



% ----- Create log file with results
log_fid = fopen( fullfile(direct,'gettraces.log'), 'wt' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');

output = evalc('disp(handles.params)');
fprintf(log_fid, '%s', output);
%FIXME: structure parameters are not displayed here (alignment!)

% Log list of files processed by gettraces
fprintf(log_fid,'\n%s\n\n%s\n%s\n\n%s\n',date,'DIRECTORY',direct,'FILES');

for i=1:nFiles
    if existing(i),
        fprintf(log_fid, 'SKIP %s\n', movieFiles(i).name);
    else
        fprintf(log_fid, '%.0f %s\n', nTraces(i), movieFiles(i).name);
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

set(handles.figure1,'pointer','watch'); drawnow;

% Locate single molecules
stkData = getappdata(handles.figure1,'stkData');
% params = handles.params;
% params.alignment = handles.alignment; %apply loaded alignment if any
[stkData,peaks] = gettraces( stkData, handles.params );

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
    
    if handles.params.alignMethod==1,
        text = 'Alignment deviation:';
    else
        text = 'Alignment applied:';
    end
    
    tableData = get(handles.tblAlignment,'Data');
    fmt = {'% 0.2f','% 0.2f','% 0.2f','% 0.2f %%','%0.2f','%0.2f'};  %sprintf formats for each field
    
    for i=1:numel(a),
        if ~isempty(a(i).theta),
            %FIXME: quality is not defined for all alignments (!)
            %FIXME: indexing here assumes donor is first; it may not be!
            row = [a(i).dx a(i).dy a(i).theta 100*(a(i).sx-1) a(i).quality a(i).abs_dev];
            tableData(i-1,:) = arrayfun( @(i) sprintf(fmt{i},row(i)), 1:numel(fmt), 'Unif',false );
            
            if a(i).quality==0,
                tableData{i-1,5}='';
            end
        end
    end
    
    set( handles.tblAlignment, 'Data',tableData(1:numel(a)-1,:) );
    set( handles.tblAlignment, 'RowName',handles.params.chDesc(2:end) );
    
    % If the alignment quality (confidence) is low, warn the user.
    % FIXME: 
    if isfield(a,'quality'),
        if any( [a.quality]<1.1 & [a.quality]>0 ),
            text = [text sprintf(' (low confidence!)')];
        end
    end
    
    set( handles.txtAlignStatus, 'String', text );

    % Color the text to draw attention to it if the alignment is bad.
    % FIXME: this should depend on the nhood/window size. 1 px may be small.
    % FIXME: try to color individual rows according to degree of misalignment,
    % for example using HTML tags in text (rg, <html><b><font color='red',
    % or color="#FF00FF"). <center> tag might also be useful.
    if any( [a.abs_dev] > 0.25 ),
        d = max( [a.abs_dev] );
        set( handles.tblAlignment,   'ForegroundColor', [(3/2)*min(2/3,d) 0 0] );
        set( handles.txtAlignStatus, 'ForegroundColor', [(3/2)*min(2/3,d) 0 0] );
    else
        set( handles.tblAlignment,   'ForegroundColor', [0 0 0] );
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
      sprintf('Molecules rejected: %0.0f%%', percentOverlap)  );

if percentOverlap>30,
    set( handles.txtOverlapStatus, 'ForegroundColor', [0.9 0 0] );
else
    set( handles.txtOverlapStatus, 'ForegroundColor', [0 0 0] );
end


% Determine the fractin of integration windows that overlap, as a measure of the
% density of molecules on the surface.
percentWinOverlap = mean(stkData.fractionWinOverlap*100);
% percentTracesWinOverlap = 100*sum(stkData.fractionWinOverlap>0)/numel(stkData.fractionWinOverlap);

set(  handles.txtWindowOverlap, 'String', ...
      sprintf('Residual win. overlap: %0.1f%%', percentWinOverlap)  );

if percentWinOverlap>10,
    set( handles.txtWindowOverlap, 'ForegroundColor', [0.9 0 0] );
else
    set( handles.txtWindowOverlap, 'ForegroundColor', [0 0 0] );
end




% Get (approximate) average fraction of fluorescence collected within the
% integration window of each molecule. Set the text color to red where the
% intensity is not well collected at the current integration window size.
eff = 100*stkData.integrationEfficiency(:,handles.params.nPixelsToSum);
eff = nanmean(eff);
set(  handles.txtIntegrationStatus, 'String', ...
      sprintf('Intensity collected: %0.0f%% ', eff)  );

if eff<70,
    set( handles.txtIntegrationStatus, 'ForegroundColor', [0.9 0 0] );
else
    set( handles.txtIntegrationStatus, 'ForegroundColor', [0 0 0] );
end


% Estimate the peak width from pixel intensity distribution.
eff = stkData.integrationEfficiency;
eff = eff( ~any(isnan(eff')), : );  %ignore NaN values, which can happen in with empty fields.
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
set(handles.btnHidePicks,'Enable','on');

set(handles.figure1,'pointer','arrow');
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
delete(findobj(handles.figure1,'type','line'));

    
if handles.params.geometry==2, %dual-channel
    % Assign each spot to the corresponding field image.
    indL  = find( handles.x<=(ncol/2) );
    indR  = find( handles.x> (ncol/2) );
    rindL = find( handles.rx<=(ncol/2) );
    rindR = find( handles.rx> (ncol/2) );
    
    % Draw markers on selection points (donor side)
    axes(handles.axDonor);
    line( handles.x(indL),   handles.y(indL),   style1{:},  'Parent',handles.axDonor );
    line( handles.rx(rindL), handles.ry(rindL), style1b{:}, 'Parent',handles.axDonor );

    % Draw markers on selection points (acceptor side)
    axes(handles.axAcceptor);
    line( handles.x(indR)-(ncol/2),   handles.y(indR),   style1{:},  'Parent',handles.axAcceptor );
    line( handles.rx(rindR)-(ncol/2), handles.ry(rindR), style1b{:}, 'Parent',handles.axAcceptor );
    
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
    line( handles.x(indUL),   handles.y(indUL),   style1{:},  'Parent',handles.axUL  );
    line( handles.rx(rindUL), handles.ry(rindUL), style1b{:}, 'Parent',handles.axUL );
    
    axes(handles.axLL);
    line( handles.x(indLL),   handles.y(indLL)-(nrow/2),   style1{:},  'Parent',handles.axLL  );
    line( handles.rx(rindLL), handles.ry(rindLL)-(nrow/2), style1b{:}, 'Parent',handles.axLL );
    
    axes(handles.axLR);
    line( handles.x(indLR)-(ncol/2),   handles.y(indLR)-(nrow/2),   style1{:},  'Parent',handles.axLR  );
    line( handles.rx(rindLR)-(ncol/2), handles.ry(rindLR)-(nrow/2), style1b{:}, 'Parent',handles.axLR );
    
    axes(handles.axUR);
    line( handles.x(indUR)-(ncol/2),   handles.y(indUR),   style1{:},  'Parent',handles.axUR  );
    line( handles.rx(rindUR)-(ncol/2), handles.ry(rindUR), style1b{:}, 'Parent',handles.axUR );
end

% Draw markers on selection points (total intensity composite image).
axes(handles.axTotal);
line( handles.total_x,  handles.total_y,  style2{:},  'Parent',handles.axTotal  );
line( handles.rtotal_x, handles.rtotal_y, style2b{:}, 'Parent',handles.axTotal );


% end function highlightPeaks




% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

% --- Executes on button press in saveTraces.
function saveTraces_Callback(~, ~, handles)

set(handles.figure1,'pointer','watch'); drawnow;

if iscell(handles.stkfile),
    filename = handles.stkfile{1};
else
    filename = handles.stkfile;
end

% Remove file extension for multi-part movies.
filename = regexprep(filename,'-file[0-9]*.[A-Za-z0-9]*$','');

% Integrate fluorophore point-spread functions, generate fluorescence
% traces, and save to file.
stkData = getappdata(handles.figure1,'stkData');
gettraces( stkData, handles.params, filename );

clear stkData;
set(handles.figure1,'pointer','arrow'); drawnow;






% --------------------- MISC. GUI CALLBACK FUNCTIONS --------------------- %

% --- Executes on slider movement.
function scaleSlider_Callback(hObject, ~, handles)  %#ok<DEFNU>
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
function txtIntensityThreshold_Callback(hObject, ~, handles)  %#ok<DEFNU>
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
function txtOverlap_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Update gettraces parameters using specified values

input = str2double(get(hObject,'String'));
if isnan(input),
    % If an invalid number is entered, reset the value to what it was.
    set(hObject,'String',handles.params.overlap_thresh);
    return;
else
    handles.params.overlap_thresh = input;
end

% Re-pick molecules with new settings.
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);


% --- Integration window size specification
function txtIntegrationWindow_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Update gettraces parameters using specified values
nPx = str2double(get(hObject,'String'));
if isnan(nPx),
    % If an invalid number is entered, reset the value to what it was.
    set(hObject,'String',handles.params.nPixelsToSum);
    return;
else
    handles.params.nPixelsToSum = floor(nPx);
end

% Re-pick molecules with new settings.
% This is only really necessary to update the GUI status (% intensity
% collected).
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);




function txtSettings_Callback(hObject, ~, handles, paramName)  %#ok<DEFNU>
% User changed one of the values int he "Analysis Settings" panel.
% paramName is passed to identify which one and the matching parameter.
% Only for buttons that have no side effects or special features.

inputstr = get(hObject,'String');
%if inputstr is empty, set the parameter to empty for automatic. TODO

input = str2double( inputstr );
if isnan(input),
    % Reset field for invalid numbers, presumably to a valid value.
    set( hObject, 'String', num2str(handles.params.(paramName)) );
else
    handles.params.(paramName) = input;
    guidata(hObject,handles);
end

% END FUNCTION txtSettings_Callback





% --- Executes on button press in chkSaveLocations.
function chkSaveLocations_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Update gettraces parameters using specified values
handles.params.saveLocations = get(handles.chkSaveLocations,'Value');
guidata(hObject,handles);




function txtMaxIntensity_Callback(hObject, ~, handles)  %#ok<DEFNU>
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

% If running, stop the "auto detect" timer. Otherwise, it may be triggered by
% the change in settings.
fileTimer = timerfind('Name','gettraces_fileTimer');
if ~isempty(fileTimer),
    stop(fileTimer);
    delete(fileTimer);
    set(handles.chkAutoBatch,'Value',0);
end

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
set( handles.cboAlignMethod,   'Value',1 ); %params.alignment.method
set( handles.edScaleAcceptor, 'String',num2str(params.scaleAcceptor) );

if handles.params.geometry==1, %Single-channel recordings
    set( handles.txtDACrosstalk,    'Enable','off', 'String','', 'Visible','on' );
    set( handles.btnCrosstalk,      'Visible','off'  );
    set( handles.btnSaveAlignment,  'Enable','off' );
    set( handles.btnLoadAlignment,  'Enable','off' );
else  %Multi-channel recordings
    set( handles.btnSaveAlignment,  'Enable','on' );
    set( handles.btnLoadAlignment,  'Enable','on' );
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



% --- Executes on button press in btnMetadata.
function btnMetadata_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Display a simple diaglog with the MetaMorph metadata for the first frame.
% 

stkData = getappdata(handles.figure1,'stkData');

if isempty(stkData) || ~isfield(stkData.movie.header,'MM')
    if isfield( stkData.movie.header,'ImageDescription' ),
        msgbox( stkData.movie.header.ImageDescription, 'ImageDescription' );
    else
        % This movie doesn't have any metadata
        set( hObject, 'Enable','off' );
    end
    return;
end

% Grab extra MetaMorph fields outside the MM struct.
metadata = stkData.movie.header.MM(1);

if isfield( stkData.movie.header, 'MM_wavelength' ),
    wv = stkData.movie.header.MM_wavelength;
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

% Reset the GUI to its previous state until alignment file is correctly loaded.
set( handles.cboAlignMethod, 'Value',handles.params.alignMethod );

% This should only be used in multi-color experiments.
assert( handles.params.geometry>1 );

% Load an alignment file
[f,p] = uigetfile('*.mat','Select an alignment settings file');
alignFilename = [p f];
if f==0, return; end

try
    input = load(alignFilename);
    
    % If the alignment was made with version 2.8 or earlier, transformations are
    % not about the center of the image (as they are in 2.9 and greater).
    % This is a big difference. Old alignments will not work correctly.
    if isfield(input.alignment(2).tform,'tdata'),
        error('gettraces:oldAlignment','Alignment files from version 2.8 and earlier are not supported.');
    elseif ~isa(input.alignment(2).tform,'affine2d'),
        error('gettraces:badAlignTform','Unrecognized tform class');
    end
    
    % Overwrite alignment settings with those in the file.
    handles.params.alignment = input.alignment;
    handles.alignment = input.alignment; %???
catch e,
    % If the file is invalid, give a warning and reset the button so
    % that is as if nothing happened.
    warning('gettraces:invalidAlignFile', ['Invalid alignment file: ' e.message]);
    return;
end

% Update GUI (alignment was successfully loaded) and pick molecules again.
handles.params.alignMethod = 2;
set( handles.cboAlignMethod, 'Value',2 );

handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);


%end function btnLoadAlignment_Callback




% --- Executes on button press in btnSaveAlignment.
function btnSaveAlignment_Callback(~, ~, handles)  %#ok<DEFNU>
% Save current software alignment settings (which may be set to do nothing
% at all) to file so they can be reloaded later.
%

assert( isfield(handles,'alignment') && ~isempty(handles.alignment) && handles.params.geometry>1 );

% Verify there is a valid software alignment 
% if ~isfield(handles,'alignment') || isempty(handles.alignment),
%     set(handles.btnSaveAlignment,'Enable','off');
%     return;
% end

[p,f] = fileparts(handles.stkfile);
alignfile = fullfile( p, [f '_align.mat'] );

[f,p] = uiputfile('*.mat','Save software alignment settings',alignfile);

if f,
    % FIXME: should also remove abs_dev.
    alignment = rmfield( handles.alignment, {'quality'} );   %#ok<NASGU>
    save( [p f], 'alignment' );
end


%end function btnSaveAlignment_Callback




% --- Executes on button press in btnCrosstalk.
function btnCrosstalk_Callback(hObject, ~, handles)  %#ok<DEFNU>
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
function cboAlignMethod_Callback(hObject, ~, handles)  %#ok<DEFNU>
% 

% methods = contents = cellstr(get(hObject,'String'));
sel = get(hObject,'Value');
% method = contents{sel};

if sel==2
    % Load alignment from file.
    % FIXME: if this fails or the user hits cancel, the dropdown has "load
    % alignment" selected, but is still in the previous state.
    btnLoadAlignment_Callback(hObject, [], handles);
else
    % Re-pick molecules with new settings.
    handles.params.alignMethod = sel;
    handles = getTraces_Callback(hObject, [], handles);
    guidata(hObject,handles);
end


% END FUNCTION cboAlignMethod_Callback



% --- Executes on button press in btnHidePicks.
function btnHidePicks_Callback(~, ~, handles)  %#ok<DEFNU>
% Hide the circles drawn to indicate molecule locations so the field of
% view image is more visible. They will show up again if the "Pick Peaks"
% button is clicked.
delete(findobj(handles.figure1,'type','line'));



% --- Executes on button press in chkAutoBatch.
function chkAutoBatch_Callback(hObject, ~, handles)  %#ok<DEFNU>
% 
    
% If another timer is running, stop it.
fileTimer = timerfind('Name','gettraces_fileTimer');

if ~isempty(fileTimer),
    stop(fileTimer);
    delete(fileTimer);
    %disp('Timer deleted.');
end

% Start a new timer if requested
if get(hObject,'Value') == get(hObject,'Max'),
    % Ask the user for a directory location
    targetDir = uigetdir('','Choose directory:');
    if targetDir==0,
        set( hObject, 'Value', get(hObject,'Min') );
        return;
    end
    disp(targetDir);
    
    % Force "skip processed data" setting so we don't end up trying to
    % reprocess every movie in every iteration.
    set(handles.chkOverwrite,'Value',1);
    
    % Start a thread that will periodically check for new movies every 5
    % seconds and process them automatically.
    %disp('Timer started.');
    fileTimer = timer('ExecutionMode','fixedSpacing','StartDelay',1, ...
                              'Name','gettraces_fileTimer', 'TimerFcn', ...
                              {@updateFileTimer,hObject,targetDir}, ...
                              'StopFcn',{@stopFileTimer,hObject}, ...
                              'Period',5.0,'BusyMode','drop');
    start(fileTimer);
    %FIXME: add an error/stop function to clear the checkbox.
end

% guidata(hObject,handles);

% END FUNCTION chkAutoBatch_Callback


function stopFileTimer(~,~,hObject)
% This function is called when there is an error during the timer callback
% or when the timer is stopped.
handles = guidata(hObject);
set(handles.chkAutoBatch,'Value',0);

% END FUNCTION stopFileTimer


function updateFileTimer(~,~,hObject,targetDir)
% This function runs each time the timer is fired, looking for any new
% movies that may have appeared on the path.
% disp('Timer fired');
batchmode_Callback( hObject, [], guidata(hObject), targetDir );


% END FUNCTION updateFileTimer
