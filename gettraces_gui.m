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
%      threshold and overlap rejection. Then select batch mode and choose a
%      directory. All the .stk files in that directory will be converted to
%      separate .traces files using the established threshold and overlap
%      rejection. Batch mode is not recursive.
%
%               -JBM 12/06
%
% Changes:
% 
% 0. Much faster run time, less memory usage (AppData), comments
% 1. Overlap distances now calculated for centroid of fluor peak.
% 2. Full bg correction used for peak picking
% 3. Peak search uses plus-shaped window instead of 3x3
% 
%           -DT 10/22

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help gettraces

% Last Modified by GUIDE v2.5 20-Aug-2008 15:45:44

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

% Update handles structure
guidata(hObject, handles);

set(handles.saveTraces,'Enable','off');
set(handles.getTraces,'Enable','off');
% set(handles.getTracesCy5,'Enable','off');
% set(handles.batchmode,'Enable','off');


% --- Outputs from this function are returned to the command line.
function varargout = gettraces_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------------- OPEN STK MOVIE (CALLBACK) --------------------- %

% --- Executes on button press in openstk.
function openstk_Callback(hObject, eventdata, handles)
% hObject    handle to openstk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load stk file
[datafile,datapath]=uigetfile( ...
    '*.stk;*.stk.bz2;*.movie','Choose a stk file');
if datafile==0, return; end

handles.stkfile = strcat(datapath,datafile);

% Load the movie
handles = OpenStk( handles.stkfile, handles, hObject );


set(handles.getTraces,'Enable','on');
guidata(hObject,handles);



% --------------------- OPEN STK MOVIE --------------------- %
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
handles.background = stkData.background;
handles.time = stkData.time;
handles.endBG = stkData.endBackground;

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
image_t    = handles.stk_top-handles.background+mean2(handles.background);
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

skipExisting = get(handles.chkOverwrite,'Value');
don_thresh = str2double(get(handles.donthresh,'String'));
overlap_thresh = str2double(get(handles.overlap,'String'));

%iptsetpref('imshowInitialMagnification','fit');

direct=uigetdir('','Choose directory:');
if direct==0, return; end

% Create header for log file
log_fid = fopen( [direct filesep 'gettraces.log'], 'w' );
fprintf(log_fid, 'Donor Thresh = %.1f\n',don_thresh);
fprintf(log_fid, 'Overlap = %.1f\n',overlap_thresh);
% fprintf(log_fid, 'Acceptor Thresh = %.1f\n',acc_thresh);

fprintf(log_fid,'\n%s\n\n%s\n%s\n\n%s\n',date,'DIRECTORY',direct,'FILES');

disp(direct);

% Create progress bar at bottom of window
set(handles.txtProgress,'String','Creating traces, please wait...');

% Get list of files in current directory (option: and all subdirectories)
recursive = get(handles.chkRecursive,'Value');

if recursive
    stk_files  = rdir([direct filesep '**' filesep '*.stk*']);
else
    stk_files  = rdir([direct filesep '*.stk*']);
end

h = waitbar(0,'Extracting traces from movies...');
tic;
% For each file in the user-selected directory
i = 1;
nFiles = length(stk_files);
for file = stk_files'
    handles.stkfile = file.name;
    
    % Skip if previously processed (.traces file exists)
    stk_fname = strrep(file.name,'.bz2','');
    [p,name]=fileparts(stk_fname);
    traceFname = [p filesep name '.traces'];
    
    if skipExisting && exist(traceFname,'file'),
        disp( ['Skipping (already processed): ' stk_fname] );
        continue;
    end
    
    % Load STK file
    handles = OpenStk(handles.stkfile,handles, hObject);
    
    % Pick molecules using default parameter values
    handles = getTraces_Callback(hObject, [], handles);
    
    % Save the traces to file
    saveTraces_Callback(hObject, [], handles);
    
    % Save entry in log file
    fprintf(log_fid, '%.0f %s\n', handles.num, file.name);
    
    
%     text = sprintf('Creating traces: %.0f%%', 100*(i/size(stk_files,1)) );
%     set(handles.txtProgress,'String',text);
%     guidata(hObject,handles);
    waitbar(i/nFiles);
    i = i+1;
end

close(h);

set(handles.txtProgress,'String','Finished.');
fclose(log_fid);
toc
guidata(hObject,handles);







% --------------- PICK MOLECULES CALLBACKS --------------- %


%------------- Pick Cy3 spots (CALLBACK) ----------------- 
% --- Executes on button press in getTraces.
function handles = getTraces_Callback(hObject, eventdata, handles)
% hObject    handle to getTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% NOTE: ignores 3 pixels from all edges, only looks at left side (cy3)


%----- Load picking parameters
overlap_thresh = str2double(get(handles.overlap,'String'));

% Donor threshold is relative to background level.
% If not chosen by user, is 10x standard dev. of background
don_thresh = str2double(get(handles.donthresh,'String'));

image_t=handles.stk_top-handles.background;
[nrow ncol] = size(image_t);

params.overlap_thresh = overlap_thresh;
params.don_thresh = don_thresh;

%----- Find peak locations from total intensity

stkData = getappdata(handles.figure1,'stkData');

% [handles.x,handles.y] = getPeaks( image_t, don_thresh, overlap_thresh,handles );
[stkData,peaks] = gettraces( stkData, params );

handles.x = peaks(:,1);
handles.y = peaks(:,2);
handles.num = numel(handles.x)/2;

clear stkData;

%----- GUI stuff

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

% Draw markers on selection points (acceptor side)
axes(handles.axTotal);
line(handles.x(l),handles.y(l),'LineStyle','none','marker','o','color','w','EraseMode','background');

% Update GUI controls
set(handles.nummoles,'String',num2str(handles.num));
set(handles.saveTraces,'Enable','on');

guidata(hObject,handles);




% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

% --- Executes on button press in saveTraces.
function saveTraces_Callback(hObject, eventdata, handles)
% hObject    handle to saveTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stkData = getappdata(handles.figure1,'stkData');

params.overlap_thresh = str2double(get(handles.overlap,'String'));
params.don_thresh = str2double(get(handles.donthresh,'String'));

% integrateAndSave( stkData.stk, stk_top, peaks, handles.stkfile, handles.time );
gettraces( stkData, params, handles.stkfile );

clear stkData;

%guidata(hObject,handles);


% --------------------- MISC. GUI CALLBACK FUNCTIONS --------------------- %

% --- Executes on slider movement.
function scaleSlider_Callback(hObject, eventdata, handles)
% hObject    handle to scaleSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set( handles.axDonor,    'CLim',[get(hObject,'min') get(hObject,'value')] );
set( handles.axAcceptor, 'CLim',[get(hObject,'min') get(hObject,'value')] );
set( handles.axTotal,    'CLim',[get(hObject,'min')*2 get(hObject,'value')*2] );




