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

% Last Modified by GUIDE v2.5 27-Aug-2009 11:22:53

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
% hObject    handle to btnSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp( handles.options );


%----- Prompt use for location to save file in...
[f,p] = uiputfile('*.txt','Save simulated data as...');
if f==0, return; end  %user pressed "cancel"
fname_txt = [p f];


%----- Gather parameter values from GUI, prepare for calling simulate()
options = handles.options;

nTraces  = options.nTraces;
traceLen = options.traceLen;
nStates  = options.nStates;
sampling = options.sampling/1000;

options.stdBackground = options.totalIntensity/(sqrt(2)*options.snr);
options.kBleach = 1/options.totalTimeOn;

if numel(options.sigma)==1,
    options.sigma = repmat(options.sigma,nStates);
elseif isempty(options.sigma),
    % Use a default value (not used for simulation)
    options.sigma = repmat(0.1,1,nStates);
end

% Generate Q matrix from rates specified in GUI
Q = zeros(nStates,nStates);
idx = find( ~logical(eye(nStates)) );
r = options.rates';
Q(idx) = r(:);
Q = Q';



%----- Check parameter values for sanity?
if numel(options.mu)~=nStates,
    error('simulate_gui:btnSimulate_Callback','Inconsistent number of states...');
end


%----- Run simulate and save results
dataSize = [nTraces traceLen];
fretModel = [options.mu' options.sigma'];

[dwt,fret,donor,acceptor] = simulate( dataSize, sampling, fretModel, Q, options );
time = 1000*sampling*( 0:(traceLen-1) );

% Save resulting raw traces files
fname_trc = strrep(fname_txt, '.txt', '.traces');
saveTraces( fname_trc, 'traces', donor,acceptor, [], time );

% Produce a .txt data file, as if it had been processed by autotrace/sorttraces.
% Traces without a descernable bleaching event are removed.
[d,a,f,ids,time] = loadTraces(fname_trc);

stats = traceStat(d,a,f);
sel = [stats.snr]>=1;
d = d(sel,:);
a = a(sel,:);
f = f(sel,:);
ids = ids(sel);

saveTraces( fname_txt, 'txt', d,a,f,ids,time );

% Save the underlying state trajectory.
% Only traces passing autotrace are saved so that the "true" idealization
% can be directly compared to estimations.
dwt = dwt(sel);
offsets = (0:(numel(dwt)-1))*traceLen * (sampling*1000);

fname_idl = strrep(fname_txt, '.txt', '.sim.dwt');
saveDWT( fname_idl, dwt, offsets, fretModel, 1 );


%----- Save a log file detailing the input parameters...
logname = strrep( fname_txt, '.txt','.log' );
fid = fopen(logname,'w');

t = clock;
fprintf(fid, 'Run time:  %d/%d/%d %d:%d', t(1:5));

fprintf(fid, '\n\nSimulating %d traces of length %d at %d sec/frame.', ...
             nTraces, traceLen, sampling);

fprintf(fid, '\n\nFRET mean: ');
fprintf(fid, '%0.2f ', fretModel(:,1) );
fprintf(fid, '\nFRET stdev: ');
fprintf(fid, '%0.2f ', fretModel(:,2) );

fprintf(fid, '\n\nTotal Intensity: %d ', options.totalIntensity );
fprintf(fid, '\nStdev Total Intensity: %d ', options.stdTotalIntensity );
fprintf(fid, '\nSNR: %d ', options.snr );
fprintf(fid, '\nStdev Fluorescence: %f ', options.stdPhoton );

fprintf(fid, '\n\nRate Matrix:\n');
fprintf(fid, '%s', get(handles.edRates,'String'));

fprintf(fid, '\n\nRandom seed: %d',options.randomSeed );

fclose(fid);



%% -------------------  CALLBACK FUNCTIONS  -------------------

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



% function edSNR_Callback(hObject, eventdata, handles)
% hObject    handle to edSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSNR as text
%        str2double(get(hObject,'String')) returns contents of edSNR as a double

% I  = str2double(get(handles.edTotalIntensity,'String'));
% IS = str2double(get(handles.edStdTotalIntensity,'String'));
% SNR  = str2double(get(handles.edSNR,'String'));
% 
% S = I/(sqrt(2)*SNR);
% 
% % text = sprintf('Noise = %0.1f ± %0.1f',I/S,IS/S );
% text = sprintf('Noise = %0.1f',S );
% set(handles.txtNoise, 'String', text);



% --- Executes on button press in chkSimFluorescence.
function chkSimFluorescence_Callback(hObject, eventdata, handles)
% hObject    handle to chkSimFluorescence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val  = get(hObject,'Value');

% If checking, enable fluorescence simulation parameter boxes
if val,
    set( handles.edTotalItensity,'Enable','on' );
    set( handles.edStdTotalItensity,'Enable','on' );
    set( handles.edSNR,'Enable','on' );
    set( handles.edStdPhoton,'Enable','on' );
    
    % Re-add these values from the options structure
    
% If unchecking, disable them.
else
    % Disable text boxes
    set( handles.edTotalItensity,'Enable','off' );
    set( handles.edStdTotalItensity,'Enable','off' );
    set( handles.edSNR,'Enable','off' );
    set( handles.edStdPhoton,'Enable','off' );
    
    % Remove these values from the options structure
end


% --- Executes on button press in chkSimulateMovies.
function chkSimulateMovies_Callback(hObject, eventdata, handles)
% hObject    handle to chkSimulateMovies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSimulateMovies



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


% --- Executes on button press in btnBrowseMovie.
function btnBrowseMovie_Callback(hObject, eventdata, handles)
% hObject    handle to btnBrowseMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


