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

% Last Modified by GUIDE v2.5 05-Aug-2009 18:13:45

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

% Choose default command line output for gettraces
handles.output = hObject;

% Setup initial values for parameter values
% params.don_thresh = 0; %not specified = auto pick
params.overlap_thresh = 2.5;
params.nPixelsToSum   = 4;

set( handles.txtIntensityThreshold,'String','' );
set( handles.txtOverlap,'String',num2str(params.overlap_thresh) );
set( handles.txtIntegrationWindow,'String',num2str(params.nPixelsToSum) );

% Update handles structure
handles.params = params;
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
    '*.stk;*.stk.bz2;*.movie','Choose a stk file');
if datafile==0, return; end

handles.stkfile = strcat(datapath,datafile);

% Load the movie
handles = OpenStk( handles.stkfile, handles, hObject );

% Update GUI now that data has been loaded.
set(handles.getTraces,'Enable','on');
guidata(hObject,handles);



% --------------------- OPEN SINGLE MOVIE --------------------- %
function handles = OpenStk(filename, handles, hObject)


set(handles.txtFilename,'String',filename);

% Clear the original stack to save memory
setappdata(handles.figure1,'stk',[]);

% Load colormap for image viewer
fid=fopen('colortable.txt','r');
colortable=fscanf(fid,'%d',[3 256]);
colortable=colortable'/255;
fclose(fid);

% Load movie data
[stkData] = gettraces( filename );
handles.stk_top = stkData.stk_top;

% Since the image stack is very large, it is stored in ApplicationData
% instead of GUIData for memory efficiency
setappdata(handles.figure1,'stkData', stkData);


% Setup slider bar (adjusting maximum value in image, initially 2x max)
low = min(min(handles.stk_top));
high = max(max(handles.stk_top));
high = ceil(high*2);

set(handles.scaleSlider,'min',low);
set(handles.scaleSlider,'max',high);
set(handles.scaleSlider,'value', (low+high)/2);

%
image_t    = handles.stk_top-stkData.background+mean2(stkData.background);
[nrow,ncol] = size(image_t);

donor_t    = image_t(:,1:ncol/2);
acceptor_t = image_t(:,(ncol/2)+1:end);
total_t    = donor_t+acceptor_t;

% Show donor image
axes(handles.axDonor);
imshow( donor_t, [low (high+low)/2] );
colormap(colortable);  zoom on;

% Show acceptor image
axes(handles.axAcceptor);
imshow( acceptor_t, [low (high+low)/2] );
colormap(colortable);  zoom on;

% Show total intensity image
axes(handles.axTotal);
imshow( total_t, [low*2 (high+low)] );
colormap(colortable);  zoom on;

linkaxes( [handles.axDonor handles.axAcceptor handles.axTotal] );


% Finish up
guidata(hObject,handles);





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
    movieFilenames  = rdir([direct filesep '**' filesep '*.stk']);
    movieFilenames  = [movieFilenames; rdir([direct filesep '**' filesep '*.stk.bz2'])];
    movieFilenames  = [movieFilenames; rdir([direct filesep '**' filesep '*.movie'])];
else
    movieFilenames  = rdir([direct filesep '*.stk']);
    movieFilenames  = [movieFilenames; rdir([direct filesep '*.stk.bz2'])];
    movieFilenames  = [movieFilenames; rdir([direct filesep '*.movie'])];
end

nFiles = length(movieFilenames);



% ---- For each file in the user-selected directory

nTraces  = zeros(nFiles,1); % number of peaks found in file X
existing = zeros(nFiles,1); % true if file was skipped (traces file exists)

% Show progress information
h = waitbar(0,'Extracting traces from movies...');
set(handles.txtProgress,'String','Creating traces, please wait...');

% For each file...
for i=1:nFiles
    stk_fname = movieFilenames(i).name;
    handles.stkfile = stk_fname;
    
    % Skip if previously processed (.traces file exists)
    stk_fname = strrep(stk_fname,'.bz2','');
    [p,name] = fileparts(stk_fname);
    traceFname = [p filesep name '.traces'];
    
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
    waitbar(i/nFiles);
end
close(h);



% ----- Create log file with results
log_fid = fopen( [direct filesep 'gettraces.log'], 'w' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');

names = fieldnames(  handles.params );
vals  = struct2cell( handles.params );

for i=1:numel(names),
    fprintf(log_fid, '  %15s:  %.2f\n', names{i}, vals{i});
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

% Locate single molecules
stkData = getappdata(handles.figure1,'stkData');
[stkData,peaks] = gettraces( stkData, handles.params );

% Update guidata with peak selection coordinates
handles.x = peaks(:,1);
handles.y = peaks(:,2);
handles.num = numel(handles.x)/2;


%----- Graphically show peak centers

ncol = stkData.stkX;
clear stkData;

% Clear selection markers
delete(findobj(gcf,'type','line'));

% Draw markers on selection points (donor side)
l = 1:2:numel(handles.x);

axes(handles.axDonor);
line(handles.x(l),handles.y(l),'LineStyle','none','marker','o','color','y','EraseMode','background');

% Draw markers on selection points (acceptor side)
ll = 2:2:numel(handles.x);

axes(handles.axAcceptor);
line(handles.x(ll)-(ncol/2),handles.y(ll),'LineStyle','none','marker','o','color','w','EraseMode','background');

% Draw markers on selection points (total intensity composite image)
axes(handles.axTotal);
line(handles.x(l),handles.y(l),'LineStyle','none','marker','o','color','w','EraseMode','background');

% Update GUI controls
set(handles.nummoles,'String',num2str(handles.num));
set(handles.saveTraces,'Enable','on');

guidata(hObject,handles);




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
set( handles.axDonor,    'CLim',[get(hObject,'min') get(hObject,'value')] );
set( handles.axAcceptor, 'CLim',[get(hObject,'min') get(hObject,'value')] );
set( handles.axTotal,    'CLim',[get(hObject,'min')*2 get(hObject,'value')*2] );

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
guidata(hObject,handles);


% --- Overlap rejection threshold specification
function txtOverlap_Callback(hObject, eventdata, handles)
% Update gettraces parameters using specified values
text = get(hObject,'String');
if ~isempty( text )
    handles.params.overlap_thresh = str2double(text);
elseif isfield(handles.params,'overlap_thresh');
    handles.params = rmfield( handles.params,'overlap_thresh' );
end
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
guidata(hObject,handles);


