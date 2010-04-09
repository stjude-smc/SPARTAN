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

% TODO:
%  - WARNING: if photoressurection rates are set to be equal, then
%    all states are equally likely after a blink!  This effectively brings
%    the system out of equilibrium!!!

% Edit the above text to modify the response to help simulate_gui

% Last Modified by GUIDE v2.5 09-Apr-2010 12:42:07

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
function simulate_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simulate_gui (see VARARGIN)


% Choose default command line output for simulate_gui
handles.output = hObject;


% Setup initial values for parameters
handles.hasModel = 0; %no model selected.

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes simulate_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = simulate_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




%% -------------------  EXECUTION FUNCTIONS  -------------------

% --- Executes on button press in btnSimulate.
function btnSimulate_Callback(hObject, eventdata, handles)


% Verify input parameters
if ~handles.hasModel,
    warning('Must load a model first!');
end

options = handles.options;
model = handles.model;
disp( options );


%----- Prompt user for location to save file in...
[f,p] = uiputfile('*.txt','Save simulated data as...');
if f==0, return; end  %user pressed "cancel"
fname_txt = [p f];


%----- If requested, prompt user for background movie location:
simMovies   = get(handles.chkSimulateMovies,'Value');

if simMovies
    % Get movie location from user.
    bgMovieDir = uigetdir;
    if ~ischar(bgMovieDir), return; end

    % Verify there are background movies to use.
    d = dir( [bgMovieDir filesep '*.stk'] );
    bgMovieFilenames = strcat( [bgMovieDir filesep], {d.name} );
    
    if numel(bgMovieFilenames)<1,
        error('No movies found');
    end
end


%----- Gather parameter values from GUI, prepare for calling simulate()
nTraces  = options.nTraces;
traceLen = options.traceLen;
nStates  = model.nStates;
sampling = options.sampling/1000;

options.stdBackground = options.totalIntensity/(sqrt(2)*options.snr);
options.kBleach = 1/options.totalTimeOn;

% if numel(options.sigma)==1,
%     options.sigma = repmat(options.sigma,nStates);
% elseif isempty(options.sigma),
%     % Use a default value (not used for simulation)
%     options.sigma = repmat(0.1,1,nStates);
% end

% Generate Q matrix from rates specified in GUI
% Q = zeros(nStates,nStates);
% idx = find( ~logical(eye(nStates)) );
% r = model.rates'; %is this still correct?
% Q(idx) = r(:);
% Q = Q';
Q = model.rates;



%----- Check parameter values for sanity?
% if numel(options.mu)~=nStates,
%     error('simulate_gui:btnSimulate_Callback','Inconsistent number of states...');
% end


%----- Run simulate and save results

% If background movies will be simulated, we don't want to
% add background noise twice.
if simMovies
    options.stdBackground = 0;
end

% Simulate FRET and fluorescence traces
dataSize = [nTraces traceLen];
fretModel = [model.mu(model.class)' model.sigma(model.class)'];
options.startProb = model.p0;

[dwt,fret,donor,acceptor] = simulate( dataSize, sampling, fretModel, Q, options );
time = 1000*sampling*( 0:(traceLen-1) );

% Generate ids for the traces
[p,name] = fileparts(fname_txt);
ids = cell(nTraces,1);
for j=1:nTraces;
    ids{j} = sprintf('%s_%d', name, j);
end

% Save resulting raw traces files
fname_trc = strrep(fname_txt, '.txt', '.traces');
saveTraces( fname_trc, 'traces', donor,acceptor, ids, time );

% Produce a .txt data file, as if it had been processed by autotrace/sorttraces.
% Traces without a descernable bleaching event are removed.
[d,a,f,ids,time] = loadTraces(fname_trc);

stats = traceStat(d,a,f);
sel = [stats.snr]>=1;
d = d(sel,:); a = a(sel,:); f = f(sel,:); ids = ids(sel);

saveTraces( fname_txt, 'txt', d,a,f,ids,time );

% Save the underlying state trajectory.
% Only traces passing autotrace are saved so that the "true" idealization
% can be directly compared to estimations.
dwt = dwt(sel);
offsets = (0:(numel(dwt)-1))*traceLen * (sampling*1000);

fname_idl = strrep(fname_txt, '.txt', '.sim.dwt');
saveDWT( fname_idl, dwt, offsets, fretModel, 1 );


%----- Simulate wide-field fluroescence movies using simulated trajectories.
if simMovies
    % Run simulateMovie.m with user-provided parameter values.
    % TODO?: Peak locations are saved in simulateMovie for accuracy evaluation.
     simulateMovie( fname_trc, bgMovieFilenames, options );
end


%----- Save a log file detailing the input parameters...
constants = cascadeConstants;

logname = strrep( fname_txt, '.txt','.log' );
fid = fopen(logname,'w');

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

fprintf(fid,'\n\n');
fclose(fid);



%% -------------------  GUI CALLBACK FUNCTIONS  -------------------

function fieldChanged_Callback(hObject,fieldName)
% This function is called whenever a GUI textbox's value is updated.
% The value extracted from the current control and saved into the
% options structure, using the given fieldName.
% If the value in the GUI is invalid, an error will be raised in str2num.

handles = guidata(hObject);

value = get(hObject,'String');
handles.options.(fieldName) = str2num(value);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in chkSimFluorescence.
function chkSimFluorescence_Callback(hObject, eventdata, handles)
controls = {'edTotalItensity','edStdTotalItensity','edSNR','edStdPhoton','edGamma'};
handleCheckbox( get(hObject,'Value'), controls, handles );

% Update handles structure
handles.options.simFluorscence = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in chkSimulateMovies.
function chkSimulateMovies_Callback(hObject, eventdata, handles)
controls = {'edDensity','edSigmaPSF','chkGrid'};
handleCheckbox( get(hObject,'Value'), controls, handles );

% Update handles structure
handles.options.simMovies = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in chkGrid.
function chkGrid_Callback(hObject, eventdata, handles)
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
 


% --- Executes on button press in btnSaveParameters.
function btnSaveParameters_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnLoadParameters.
function btnLoadParameters_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edModelFilename_Callback(hObject, eventdata, handles)

% Make sure the model filename is valid.
filename = get(handles.edModelFilename,'String');
if exist(filename,'file')~=2,
    disp('Model file does not exist!');
    set(hObject,'ForegroundColor',[0.4 0.4 0.4]);
    return;
end

% Load model
handles.modelFilename = filename;
handles.model = qub_loadModel( filename );
handles.hasModel = 1;

% Update GUI
set(hObject,'ForegroundColor',[0 0 0]);
guidata(hObject, handles);


% --- Executes on button press in btnBrowseModel.
function btnBrowseModel_Callback(hObject, eventdata, handles)

% Set startup location for selecting a model
filename = get(handles.edModelFilename,'String');
if isempty(filename) || filename(1)=='(',
    constants = cascadeConstants;
    filename = [constants.modelLocation filesep '*.qmf'];
end

% Get filename of a QuB model from user.
[f,p] = uigetfile( {'*.qmf','QuB Model File (*.qmf)';'*.*','All Files'}, ...
                   'Select a QuB model file...',filename);
if f==0, return; end
filename = [p f]

set(handles.edModelFilename,'String',filename);
set(handles.edModelFilename,'ForegroundColor',[0 0 0]);

% Load model data.
handles.modelFilename = filename;
handles.model = qub_loadModel( filename );
handles.hasModel = 1;

guidata(hObject, handles);


% --- Executes on button press in btnEditModel.
function btnEditModel_Callback(hObject, eventdata, handles)
% 
if handles.hasModel
    [hFigure,hText,hRate] = drawModel( handles.modelFilename );
end
