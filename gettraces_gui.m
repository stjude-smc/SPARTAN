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
%               -JBM 12/06
%
% Changes:
% 
% 0. Much faster run time, less memory usage (AppData), comments
% 1. txtoverlap distances now calculated for centroid of fluor peak.
% 2. Full bg correction used for peak picking
% 3. Peak search uses plus-shaped window instead of 3x3
% 
%           -DT 10/22

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help gettraces

% Last Modified by GUIDE v2.5 06-Mar-2013 15:03:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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
function gettraces_OpeningFcn(hObject, eventdata, handles, varargin)
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

    % Setup default values for parameter values -- 2-color FRET.
    constants = cascadeConstants();
    params = constants.gettracesDefaultParams;

    if ~isfield(params,'don_thresh') || params.don_thresh==0,
        set( handles.txtIntensityThreshold,'String','' );
    else
        set( handles.txtIntensityThreshold,'String',num2str(params.don_thresh) );
    end
    set( handles.txtOverlap,           'String', num2str(params.overlap_thresh)   );
    set( handles.txtIntegrationWindow, 'String', num2str(params.nPixelsToSum)     );
    set( handles.txtDACrosstalk,       'String', num2str(params.crosstalk)        );
    set( handles.txtPhotonConversion,  'String', num2str(params.photonConversion) );
    
    set( handles.chkAlignTranslate, 'Value', params.alignTranslate );
    set( handles.chkAlignRotate,    'Value', params.alignRotate    );
    
    set( handles.chkRecursive, 'Value', params.recursive    );
    set( handles.chkOverwrite, 'Value', params.skipExisting );
    
    set( handles.cboGeometry, 'Value', params.geometry );

    handles.params = params;
end

% Update handles structure
guidata(hObject, handles);

% gettraces may be called from sorttraces to load the movie associated with
% a particular trace. The first argument is then the filename of the movie
% file and the second argument is the x-y coordinate of the trace.
% FIXME: this way of doing it generates warnings in sorttraces, so there is
% something not quite right!!
% FIXME: there is some missing information for 3-color implementation.
% Channel names (and the order) are really important, but they are not in
% traceMetadata!
if numel(varargin) > 0,
    handles.stkfile = varargin{1};
    traceMetadata = varargin{2};
    
    % Determine imaging geometry and build list of coordinates from metadata.
    geometry=1;
    handles.x = traceMetadata.donor_x;
    handles.y = traceMetadata.donor_y;
    
    if isfield(traceMetadata,'acceptor_x'),
        geometry=2;
        handles.x = [handles.x traceMetadata.acceptor_x];
        handles.y = [handles.y traceMetadata.acceptor_y];
    end
    
    if isfield(traceMetadata,'factor_x'),
        geometry=3;
        handles.x = [handles.x traceMetadata.factor_x];
        handles.y = [handles.y traceMetadata.factor_y];
    end
    
    handles.params.geometry = geometry;
    set( handles.cboGeometry, 'Value', geometry );
    
    handles.num = numel(handles.x)/handles.params.geometry;
    set(handles.nummoles,'String',num2str(handles.num));
    
    % Load file
    handles = OpenStk( handles.stkfile, handles, hObject );
    set(handles.getTraces,'Enable','on');
    
    % Highlight the selected trace.
    highlightPeaks(handles);
end


% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = gettraces_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------- OPEN SINGLE MOVIE (CALLBACK) ---------------- %

% --- Executes on button press in openstk.
function openstk_Callback(hObject, eventdata, handles)
% hObject    handle to openstk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get filename of input data from user
[datafile,datapath]=uigetfile( ...
    '*.stk;*.stk.bz2;*.tif*','Choose a movie file');
if datafile==0, return; end

handles.stkfile = strcat(datapath,datafile);

% Load the movie
handles = OpenStk( handles.stkfile, handles, hObject );

% Update GUI now that data has been loaded.
guidata(hObject,handles);



% --------------------- OPEN SINGLE MOVIE --------------------- %
function handles = OpenStk(filename, handles, hObject)

constants = cascadeConstants;

[p,f,e] = fileparts(filename);
fnameText = [p filesep f e];

% Trancate name if too long ot fit into window without wrapping.
if numel(fnameText)>90,
    fnameText = ['...' fnameText(end-90:end)];
end

set( handles.txtFilename,          'String',fnameText);
set( handles.txtOverlapStatus,     'String', '' );
set( handles.txtIntegrationStatus, 'String', '' );
set( handles.txtPSFWidth,          'String', '' );
set( handles.txtAlignStatus,      'String', '' );
set( handles.txtAlignWarning, 'Visible','off' );

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
high = max(max(handles.stk_top));
high = min( ceil(high*1.5), 32000 );
val = (low+high)/2;

set(handles.scaleSlider,'min',low);
set(handles.scaleSlider,'max',32000); %uint32 maxmimum value
set(handles.scaleSlider,'value', val);
set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));

%
image_t    = handles.stk_top-stkData.background;
[nrow,ncol] = size(image_t);


%---- Show fields for Single-Channel (full-chip) recordings.
if handles.params.geometry==1,
    handles.params.chNames = constants.gettraces_chNames2;
    handles.params.chDesc  = constants.gettraces_chDesc2;
    
    handles.total_t = image_t;
    
    % Show full field of view.
    handles.axTotal = subplot( 1,1,1, 'Parent',handles.panView, 'Position',[0.2 0 0.6 0.95] );
    imshow( image_t, [low (high+low)] );
    colormap(handles.colortable);  zoom on;
    title('Single-Channel');
    
%---- Show fields for Dual-Channel (half-chip, L/R) recordings.
elseif handles.params.geometry==2,
    handles.params.chNames = constants.gettraces_chNames2;
    handles.params.chDesc  = constants.gettraces_chDesc2;
    
    donor_t    = image_t(:,1:ncol/2);
    acceptor_t = image_t(:,(ncol/2)+1:end);
    total_t    = donor_t+acceptor_t;
    handles.total_t = total_t;
    
    chNames = upperFirst( handles.params.chNames );

    % Show donor image
    handles.axDonor = subplot( 1,3,1, 'Parent',handles.panView, 'Position',[0.025 0 0.3 0.95] );
    imshow( donor_t, [low (high+low)/2] );
    colormap(handles.colortable);  zoom on;
    title( [chNames{1} ' (' handles.params.chDesc{1} ')'] );

    % Show acceptor image
    handles.axAcceptor = subplot( 1,3,2, 'Parent',handles.panView, 'Position',[0.350 0 0.3 0.95] );
    imshow( acceptor_t, [low (high+low)/2] );
    colormap(handles.colortable);  zoom on;
    title( [chNames{2} ' (' handles.params.chDesc{2} ')'] );

    % Show total intensity image
    handles.axTotal = subplot( 1,3,3, 'Parent',handles.panView, 'Position',[0.675 0 0.3 0.95] );
    imshow( total_t, [low*2 (high+low)] );
    colormap(handles.colortable);  zoom on;
    title('Total');

    linkaxes( [handles.axDonor handles.axAcceptor handles.axTotal] );
    
%---- Show fields for Quad-Channel (quarter-chip, L/R/U/D) recordings.
elseif handles.params.geometry>2,
    % Split up channels and combine, assuming perfect alignment.
    upperLeft  = image_t( 1:nrow/2, 1:ncol/2 );
    upperRight = image_t( 1:nrow/2, (ncol/2)+1:end );
    lowerLeft  = image_t( (nrow/2)+1:end, 1:ncol/2 );
    lowerRight = image_t( (nrow/2)+1:end, (ncol/2)+1:end );
    total_t = upperLeft + upperRight + lowerLeft + lowerRight;
    
    handles.total_t = total_t;
    
    % Get channel names and make them easier to read.
    handles.params.chNames = constants.gettraces_chNames4;
    handles.params.chDesc  = constants.gettraces_chDesc4;
    chNames = upperFirst( handles.params.chNames );
    %chNames = strrep(chNames, 'Factor','Factor Binding');
    
    % Create axes for each of the fluorescence channels.
    handles.axUL    = subplot( 2,3,1, 'Parent',handles.panView, 'Position',[0.025 0.5 0.25 0.5] );
    imshow( upperLeft, [low*2 (high+low)] );
    colormap(handles.colortable);  zoom on;
    if ~isempty(chNames{1})
        title( [chNames{1} ' (' handles.params.chDesc{1} ')'] );
    end
    
    handles.axUR    = subplot( 2,3,2, 'Parent',handles.panView, 'Position',[0.35 0.5 0.25 0.5] );
    imshow( upperRight, [low*2 (high+low)] );
    colormap(handles.colortable);  zoom on;
    if ~isempty(chNames{4})
        title( [chNames{4} ' (' handles.params.chDesc{4} ')'] );
    end
    
    handles.axLR    = subplot( 2,3,5, 'Parent',handles.panView, 'Position',[0.35 0.025 0.25 0.5] );
    imshow( lowerRight, [low*2 (high+low)] );
    colormap(handles.colortable);  zoom on;
    if ~isempty(chNames{3})
        title( [chNames{3} ' (' handles.params.chDesc{3} ')'] );
    end
    
    handles.axLL    = subplot( 2,3,4, 'Parent',handles.panView, 'Position',[0.025 0.025 0.25 0.5] );
    imshow( lowerLeft, [low*2 (high+low)] );
    colormap(handles.colortable);  zoom on;
    if ~isempty(chNames{2})
        title( [chNames{2} ' (' handles.params.chDesc{2} ')'] );
    end
    
    handles.axTotal = subplot( 2,3,3, 'Parent',handles.panView, 'Position',[0.65 0.25 0.25 0.5] );
    imshow( total_t, [low*2 (high+low)] );
    colormap(handles.colortable);  zoom on;
    title('Total');
    
    linkaxes( [handles.axUL handles.axUR handles.axLL handles.axLR handles.axTotal] );
end


% Finish up
set(handles.getTraces,'Enable','on');
set(handles.btnMetadata,'Enable','on');
guidata(hObject,handles);



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
function batchmode_Callback(hObject, eventdata, handles)
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
    movieFilenames  = rdir( [direct filesep '**' filesep '*.stk*'] );
    movieFilenames  = [ movieFilenames ; rdir([direct filesep '**' filesep '*.tif*']) ];
else
    movieFilenames  = rdir( [direct filesep '*.stk*'] );
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
    stk_fname = strrep(stk_fname,'.bz2','');
    [p,name] = fileparts(stk_fname);
    traceFname = [p filesep name '.rawtraces'];
    
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
log_fid = fopen( [direct filesep 'gettraces.log'], 'w' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');

names = fieldnames(  handles.params );
vals  = struct2cell( handles.params );

for i=1:numel(names),
    if iscell( vals{i} )
        f = repmat( '%s, ', 1,numel(vals{i}));
        f = f(1:end-2);
        fprintf(log_fid, ['  %15s:  ' f '\n'], names{i}, vals{i}{:});
    else
        fprintf(log_fid, '  %15s:  %.2f\n', names{i}, vals{i});
    end
end

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

function handles = getTraces_Callback(hObject, eventdata, handles)

%----- Find peak locations from total intensity

% Do nothing if no file has been loaded. This function may be triggered by
% changing the settings fields before a file is open.
if ~isfield(handles, 'stkfile')
    return;
end

% Locate single molecules
stkData = getappdata(handles.figure1,'stkData');
[stkData,peaks] = gettraces( stkData, handles.params );

% The alignment may involve shifting (or distorting) the fields to get a
% registered donor+acceptor field. Show this distorted imaged so the user
% can see what the algorithm is doing.
% axes(handles.axTotal);
val = get(handles.scaleSlider,'value');
minimum = get(handles.scaleSlider,'min');
val = max(val,minimum+1);

imshow( stkData.total_t, [minimum*2 val*2], 'Parent',handles.axTotal );
colormap(handles.colortable);  zoom on;  title('Total');


% If no alignment data given (for example in single-channel recordings),
% don't display any status messages.
set( handles.txtAlignWarning, 'Visible','off' );

if ~isfield(stkData,'alignStatus') || isempty(stkData.alignStatus),
    set( handles.txtAlignStatus, 'String', '' );
    
% Display alignment status to inform user if realignment may be needed.
% Format: translation deviation (x, y), absolute deviation (x, y)
else
    a = stkData.alignStatus;
    
    if isfield(a,'quality'),
        text = 'Alignment applied';
    else
        text = 'Alignment deviation';
    end
    
    text = sprintf('%s:\n%0.1f (x), %0.1f (y), %0.1f° (rot),\n%0.1f (abs)', ...
                    text, [a.dx a.dy a.theta a.abs_dev] );
                
    if isfield(a,'residual_dev')
        text = [text sprintf(', %0.1f (res)',a.residual_dev)];
    end
    
    set( handles.txtAlignStatus, 'String', text );

    if a.abs_dev>0.25,
        set( handles.txtAlignStatus, 'ForegroundColor', [(3/2)*min(2/3,a.abs_dev) 0 0] );
    else
        set( handles.txtAlignStatus, 'ForegroundColor', [0 0 0] );
    end
    
    % Total misalignment (no corresponding peaks) can give a relatively low
    % alignment since the algorithm can only search a 1px neighborhood. So
    % when things get close to 1, we need a big warning. 
    if a.abs_dev>=0.7,
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
handles.num = numel(handles.x)/handles.params.geometry;
set(handles.nummoles,'String',num2str(handles.num));

set(handles.saveTraces,'Enable','on');

guidata(hObject,handles);

% end function



function highlightPeaks(handles)
% Draw circles around each selected fluorescence spot, defined in handles.x
% and handles.y
%

    
% Parameters for drawing circles around each detected molecule.
style1 = {'LineStyle','none','marker','o','color','w','EraseMode','background'};
style2 = {'LineStyle','none','marker','o','color','y','EraseMode','background'};

% Parameters for overlapped molecules that won't actually be recorded.
% These are darker to de-emphasize them, but yet still be visible.
style1b = {'LineStyle','none','marker','o','color',[0.4,0.4,0.4],'EraseMode','background'};
style2b = {'LineStyle','none','marker','o','color',[0.4,0.4,0.0],'EraseMode','background'};


[nrow,ncol] = size(handles.stk_top);

% Clear any existing selection markers from previous calls.
delete(findobj(gcf,'type','line'));


if handles.params.geometry==1, %single-channel
    % Draw markers on selection points (total intensity composite image)
    axes(handles.axTotal);
    line(handles.x,handles.y, style2{:});
    line(handles.rx,handles.ry, style2b{:});
    
elseif handles.params.geometry==2, %dual-channel
    indD = 1:2:numel(handles.x);
    indA = 2:2:numel(handles.x);
        
    % Draw markers on selection points (donor side)
    axes(handles.axDonor);
    line(handles.x(indD),handles.y(indD), style1{:});
    line(handles.rx(1:2:end),handles.ry(1:2:end), style1b{:});

    % Draw markers on selection points (acceptor side)
    axes(handles.axAcceptor);%default 3-color channel assignments.
    line(handles.x(indA)-(ncol/2),handles.y(indA), style1{:});
    line(handles.rx(2:2:end)-(ncol/2),handles.ry(2:2:end), style1b{:});

    % Draw markers on selection points (total intensity composite image)
    axes(handles.axTotal);
    line(handles.total_x,handles.total_y, style2{:});
    line(handles.rtotal_x,handles.rtotal_y, style2b{:});
    
elseif handles.params.geometry>2, %quad-channel
    % Assign each spot to a quadrant.
    indUL = find( handles.x<=(ncol/2) & handles.y<=(nrow/2) );
    indLL = find( handles.x<=(ncol/2) & handles.y> (nrow/2) );
    indLR = find( handles.x> (ncol/2) & handles.y> (nrow/2) );
    indUR = find( handles.x> (ncol/2) & handles.y<=(nrow/2) );
    
    axes(handles.axUL);
    line(handles.x(indUL),handles.y(indUL), style1{:});
    
    axes(handles.axLL);
    line(handles.x(indLL),handles.y(indLL)-(nrow/2), style1{:});
    
    axes(handles.axLR);
    line(handles.x(indLR)-(ncol/2),handles.y(indLR)-(nrow/2), style1{:});
    
    axes(handles.axUR);
    line(handles.x(indUR)-(ncol/2),handles.y(indUR), style1{:});
    
    % Draw markers on selection points (total intensity composite image).
    axes(handles.axTotal);
    line(handles.total_x,handles.total_y, style2{:});
end


% end function highlightPeaks




% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

% --- Executes on button press in saveTraces.
function saveTraces_Callback(hObject, eventdata, handles)

% Integrate fluorophore point-spread functions, generate fluorescence
% traces, and save to file.
stkData = getappdata(handles.figure1,'stkData');
gettraces( stkData, handles.params, handles.stkfile );
clear stkData;






% --------------------- MISC. GUI CALLBACK FUNCTIONS --------------------- %

% --- Executes on slider movement.
function scaleSlider_Callback(hObject, eventdata, handles)
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
function txtIntensityThreshold_Callback(hObject, eventdata, handles)
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
function txtOverlap_Callback(hObject, eventdata, handles)
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
function txtIntegrationWindow_Callback(hObject, eventdata, handles)
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
function chkSaveLocations_Callback(hObject, eventdata, handles)
% Update gettraces parameters using specified values
handles.params.saveLocations = get(handles.chkSaveLocations,'Value');
guidata(hObject,handles);




function txtMaxIntensity_Callback(hObject, eventdata, handles)
% Update axes color limits from new slider value
val = str2num( get(hObject,'String') );
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
function cboGeometry_Callback(hObject, eventdata, handles)
%

handles.params.geometry = get(hObject,'Value');

% If a movie has already been loaded, reload movie with new setup
if isfield(handles,'stkfile'),
    handles = OpenStk( handles.stkfile, handles, hObject );
end

if handles.params.geometry==1, %Single-channel recordings
    set( handles.txtDACrosstalk,    'Enable','off', 'String','' );
    set( handles.chkAlignTranslate, 'Enable','off', 'Value',0   );
    set( handles.chkAlignRotate,    'Enable','off', 'Value',0   );
    
elseif handles.params.geometry==2, %Dual-channel recordings
    set( handles.txtDACrosstalk,    'Enable','on', 'String',num2str(handles.params.crosstalk) );
    set( handles.chkAlignTranslate, 'Enable','on', 'Value',handles.params.alignTranslate      );
    set( handles.chkAlignRotate,    'Enable','on', 'Value',handles.params.alignRotate  );
    
elseif handles.params.geometry>2, %Three-color recordings
    set( handles.txtDACrosstalk,    'Enable','on', 'String',num2str(handles.params.crosstalk) );
    set( handles.chkAlignTranslate, 'Enable','off', 'Value',0  );
    set( handles.chkAlignRotate,    'Enable','off', 'Value',0  );
end


guidata(hObject,handles);




function txtDACrosstalk_Callback(hObject, eventdata, handles)
% 
handles.params.crosstalk = str2num( get(hObject,'String') );
guidata(hObject,handles);





function txtPhotonConversion_Callback(hObject, eventdata, handles)
%
handles.params.photonConversion = str2num( get(hObject,'String') );
guidata(hObject,handles);





% --- Executes on button press in chkAlignTranslate.
function chkAlignTranslate_Callback(hObject, eventdata, handles)
%
handles.params.alignTranslate = get(hObject,'Value');

% Re-pick molecules with new settings.
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);



% --- Executes on button press in chkAlignRotate.
function chkAlignRotate_Callback(hObject, eventdata, handles)
%
handles.params.alignRotate = get(hObject,'Value');

% Re-pick molecules with new settings.
handles = getTraces_Callback( hObject, [], handles);
guidata(hObject,handles);




% --- Executes on button press in btnMetadata.
function btnMetadata_Callback(hObject, eventdata, handles)
% Display a simple diaglog with the MetaMorph metadata for the first frame.
% 

stkData = getappdata(handles.figure1,'stkData');

if isempty(stkData) || ~isfield(stkData.movie.stkHeader,'MM')
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




