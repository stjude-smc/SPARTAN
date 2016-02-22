function varargout = simulate_gui(varargin)
% SIMULATE_GUI M-file for simulate_gui.fig
%      SIMULATE_GUI, by itself, creates a new SIMULATE_GUI or raises the existing
%      singleton*.
%
%      H = SIMULATE_GUI returns the handle to a new SIMULATE_GUI or the handle to
%      the existing singleton*.
%
%      SIMULATE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULATE_GUI.M with the given input arguments.
%
%      SIMULATE_GUI('Property','Value',...) creates a new SIMULATE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simulate_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to simulate_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% TODO:
%  - WARNING: if photoressurection rates are set to be equal, then
%    all states are equally likely after a blink!  This effectively brings
%    the system out of equilibrium!!!

% Last Modified by GUIDE v2.5 23-Nov-2015 13:23:26


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulate_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @simulate_gui_OutputFcn, ...
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


% --- Executes just before simulate_gui is made visible.
function simulate_gui_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simulate_gui (see VARARGIN)

updateSpartan; %check for updates

% Choose default command line output for simulate_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

constants = cascadeConstants;
set( handles.figure1, 'Name', ['simulate - ' constants.software] );

% UIWAIT makes simulate_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = simulate_gui_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




%% -------------------  EXECUTION FUNCTIONS  -------------------

% --- Executes on button press in btnSimulate.
function btnSimulate_Callback(~, ~, handles)   %#ok<DEFNU>
% hObject    handle to btnSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

m = handles.model;
options = handles.options;


%----- Prompt user for location to save file in...
[f,p] = uiputfile('*.traces','Save simulated data as...');
if f==0, return; end  %user pressed "cancel"
fname_output = [p f];


%----- If requested, prompt user for background movie location:
simMovies   = get(handles.chkSimulateMovies,'Value');

if simMovies
    % Get movie location from user.
    bgMovieDir = uigetdir(pwd,'Select directory where background movies are stored');
    if ~ischar(bgMovieDir), return; end

    % Verify there are background movies to use.
    d = regexpdir(bgMovieDir,'^.*\.(tiff?|stk)$',false);
    bgMovieFilenames = {d.name};
    
    if numel(bgMovieFilenames)<1,
        error('No movies found');
    end
end

set(handles.figure1, 'pointer', 'watch'); drawnow;


%----- Gather parameter values from GUI, prepare for calling simulate()
nTraces  = options.nTraces;
traceLen = options.traceLen;
sampling = options.sampling/1000;

options.stdBackground = options.totalIntensity/(sqrt(2)*options.snr);
options.kBleach = 1/options.totalTimeOn;


%----- Run simulate and save results

% If background movies are being used, shot noise and background noise are added
% at that stage instead of in simulate(). This alters the input to simulate().
options2 = options;
if simMovies
    options2.stdBackground = 0;
    options2.shotNoise = false;
    options2.emNoise = false;
end

% Ensure model dimensions for saving
m.mu = reshape( m.mu, [numel(m.mu),1] );
m.sigma = reshape( m.sigma, [numel(m.sigma),1] );

% Simulate FRET and fluorescence traces
dataSize = [nTraces traceLen];
fretModel = [m.mu m.sigma];

data = TracesFret(nTraces,traceLen);

[dwt,data.fret,data.donor,data.acceptor] = simulate( dataSize, sampling, m, options2 );
data.time = 1000*sampling*( 0:(traceLen-1) );
data.fileMetadata(1).wavelengths = [532 640];

% Save resulting raw traces files
saveTraces( fname_output, data );

% Save the underlying state trajectory.
% Only traces passing autotrace are saved so that the "true" idealization
% can be directly compared to estimations.
% dwt = dwt(sel);
offsets = (0:(numel(dwt)-1))*traceLen * (sampling*1000);

[p,f] = fileparts(fname_output);
fname_idl = fullfile(p, [f '.sim.dwt']);
saveDWT( fname_idl, dwt, offsets, fretModel, 1 );


%----- Simulate wide-field fluroescence movies using simulated trajectories.
if simMovies
     simulateMovie( data, bgMovieFilenames, fullfile(p,[f '.tif']), options );
end


%----- Save a log file detailing the input parameters...
constants = cascadeConstants;
fid = fopen( fullfile(p,[f '.sim.log']), 'w' );

t = clock;
fprintf(fid, 'Run time:  %d/%d/%d %d:%d  (v%s)', t(1:5), constants.version);

% Save values of all other constants used
fprintf(fid, '\n\nSIMULATION PARAMETERS');

names = fieldnames(  options );
vals  = struct2cell( options );

for i=1:numel(names),
    if numel( vals{i} )==1,
        fprintf(fid, '\n  %22s:  %.2f', names{i}, vals{i} );
    elseif numel( vals{i} )>1,
        fprintf(fid, '\n  %22s:  %s', names{i}, mat2str(vals{i}) );
    end
end

fprintf(fid, '\n\nMODEL PARAMETERS');

names = fieldnames(  m );

for i=1:numel(names),
    if isstruct( m.(names{i}) ),
        fprintf(fid, '\n  %22s:  %s', names{i}, '(struct)' );
        continue;
    end
    
    if numel( m.(names{i}) )>1
        fprintf(fid, '\n  %22s:  %s', names{i}, mat2str(m.(names{i})) );
    else
        fprintf(fid, '\n  %22s:  %.2f', names{i}, m.(names{i}) );
    end
end

fprintf(fid,'\n\n');
fclose(fid);

set(handles.figure1, 'pointer', 'arrow');



%% -------------------  GUI CALLBACK FUNCTIONS  -------------------

function fieldChanged_Callback(hObject,fieldName)   %#ok<DEFNU>
% This function is called whenever a GUI textbox's value is updated.
% The value extracted from the current control and saved into the
% options structure, using the given fieldName.
% If the value in the GUI is invalid, an error will be raised in str2num.

handles = guidata(hObject);
style = get(hObject,'Style');

if strcmp(style,'checkbox'),
    handles.options.(fieldName) = get(hObject,'Value');
elseif strcmp(style,'edit')
    value = get(hObject,'String');
    handles.options.(fieldName) = str2double(value);
else
    error('Invalid field');
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in chkSimFluorescence.
function chkSimFluorescence_Callback(hObject, ~, handles)   %#ok<DEFNU>
controls = {'edTotalItensity','edStdTotalItensity','edSNR','edStdPhoton','edGamma','chkShotNoise','chkEMNoise'};
handleCheckbox( get(hObject,'Value'), controls, handles );

% Update handles structure
handles.options.simFluorscence = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in chkSimulateMovies.
function chkSimulateMovies_Callback(hObject, ~, handles)   %#ok<DEFNU>
controls = {'edDensity','edSigmaPSF','chkGrid','edAlignX','edAlignY','edAlignTheta','edADUPhoton'};  %,'edAlignMag'};
handleCheckbox( get(hObject,'Value'), controls, handles );

% Update handles structure
handles.options.simMovies = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in chkGrid.
function chkGrid_Callback(hObject, ~, handles)   %#ok<DEFNU>
% Save selection in parameter list.
if get(hObject,'Value'),
    handles.options.grid = 1;
else
    handles.options = rmfield(handles.options,'grid');
end

% Update handles structure
guidata(hObject, handles);


% ---
function handleCheckbox( val, controlNames, handles )
% If val=1, enable all controls listed in controlNames,
% if val=0, disable all controls listed.
if val, state='on'; else state='off'; end

 for i=1:numel(controlNames),
    set( handles.(controlNames{i}),'Enable',state );
 end

 

% --- Executes on button press in btnLoadModel.
function btnLoadModel_Callback(hObject, ~, handles)   %#ok<DEFNU>
%

% Try to send the user to the model location which may be far away from
% data.
persistent modelLoc;

if isempty(modelLoc),
    modelLoc = pwd;
end

[p] = fileparts(modelLoc);

% Ask the user for a filename
[fname,p] = uigetfile( [p filesep '*.qmf'], 'Load QuB model file' );
if fname==0, return; end
fname = [p fname];
modelLoc = fname;

% Load the model
handles.model = QubModel(fname);
set( handles.txtModelFilename, 'String',['...' fname(max(1,end-40):end)] );

% Show the model properties in the GUI. The model's properties are
% automatically updated whenever the model is modified in the GUI.
handles.model.showModel( handles.axModel );


% Enable saving the model
set( handles.btnSaveModel, 'Enable','on' );

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in btnSaveModel.
function btnSaveModel_Callback(~, ~, handles)   %#ok<DEFNU>
%  Save the currently loaded model to file.

if ~isfield(handles,'model') || isempty(handles.model),
    return;
end

fname = handles.model.filename;
[f,p] = uiputfile(fname,'Save model to file');

if f~=0,
    fname = fullfile(p,f);
    handles.model.save( fname );
    set( handles.txtModelFilename, 'String',['...' fname(max(1,end-50):end)] );
end
