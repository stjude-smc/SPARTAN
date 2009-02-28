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

% Last Modified by GUIDE v2.5 28-Feb-2009 15:52:14

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

% Prompt use for location to save file in...
[f,p] = uiputfile('*.txt','Save simulated data as...');
if f==0, return; end  %user pressed "cancel"
fname_txt = [p f];

% Gather parameter values from GUI
nStates   = str2double( get(handles.edNstates,  'String') );
traceLen  = str2double( get(handles.edNframes,  'String') );
nTraces   = str2double( get(handles.edNtraces,  'String') );
sampling = str2double( get(handles.edSampling,'String') );
sampling = sampling/1000;

mu        = eval( get(handles.edMu,   'String') );
sigma     = eval( get(handles.edSigma,'String') );
rates     = eval( get(handles.edRates,'String') );

totalIntensity    = str2double( get(handles.edTotalIntensity,   'String') );
stdTotalIntensity = str2double( get(handles.edStdTotalIntensity,'String') );
snr = str2double( get(handles.edSNR,  'String') );
stdFluorescence = totalIntensity/(sqrt(2)*snr);

% totalTon = str2double( get(handles.edTotalTon, 'String') );
% Ton      = str2double( get(handles.edTon,      'String') );
% Toff     = str2double( get(handles.edToff,     'String') );

randomSeed = str2double( get(handles.edRandomSeed,'String') );


%
if numel(sigma)==1,
    sigma = repmat(sigma,size(mu));
end


% Compile array of optional arguments to simulate...
options = { 'totalIntensity',totalIntensity, ...
            'stdTotalIntensity',stdTotalIntensity, ...
            'stdFluorescence',stdFluorescence, ...
            'randomSeed',randomSeed };


% Check parameter values for sanity?
if numel(mu)~=nStates || numel(sigma)~=nStates,
    error('simulate_gui:btnSimulate_Callback','Inconsistent number of states...');
end

% Generate Q matrix
Q = zeros(nStates,nStates);
idx = find( ~logical(eye(nStates)) );
r = rates';
Q(idx) = r(:);
Q = Q';

% Add blinking kinetic parameters to matrix...


% Run simulate
dataSize = [nTraces traceLen];
model = [mu' sigma'];

[dwt,fret,donor,acceptor] = simulate( dataSize, sampling, model, Q, options{:} );


% Save resulting data
fname_qub = strrep(fname_txt, '.txt', '.qub.txt');
fname_trc = strrep(fname_txt, '.txt', '.traces');

% saveTraces( fname_qub, 'qub', fret );
% saveTraces( fname_txt, 'txt', donor,acceptor,fret );
saveTraces( fname_trc, 'traces', donor,acceptor );

[d,a,f,ids] = loadTraces(fname_trc);
saveTraces( fname_txt, 'txt', d,a,f,ids );

% Save the underlying state trajectory
fname_idl = strrep(fname_txt, '.txt', '.sim.dwt');

tl = numel(dwt);
offsets = (0:(nTraces-1))*tl;
% dwt = cell(nTraces,1);
% 
% for i=1:nTraces,
%     dwt{i} = RLEncode( idl(i,:) );
% end

saveDWT( fname_idl, dwt, offsets, model, 1 );

% save a log file detailing the input parameters...
logname = strrep( fname_txt, '.txt','.log' );
fid = fopen(logname,'w');

t = clock;
fprintf(fid, 'Run time:  %d/%d/%d %d:%d', t(1:5));

fprintf(fid, '\n\nSimulating %d traces of length %d at %d sec/frame.', ...
             nTraces, traceLen, sampling);

fprintf(fid, '\n\nFRET mean: ');
fprintf(fid, '%0.2f ', model(:,1) );
fprintf(fid, '\nFRET stdev: ');
fprintf(fid, '%0.2f ', model(:,2) );

fprintf(fid, '\n\nTotal Intensity: %d ', totalIntensity );
fprintf(fid, '\nStdev Total Intensity: %d ', stdTotalIntensity );
fprintf(fid, '\nSNR: %d ', snr );
fprintf(fid, '\nStdev Fluorescence: %f ', stdFluorescence );

fprintf(fid, '\n\nRate Matrix:\n');
fprintf(fid, '%s', get(handles.edRates,'String'));

fprintf(fid, '\n\nRandom seed: %d',randomSeed );

fclose(fid);



%% -------------------  CALLBACK FUNCTIONS  -------------------

function edNstates_Callback(hObject, eventdata, handles)
% hObject    handle to edNstates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edNstates as text
%        str2double(get(hObject,'String')) returns contents of edNstates as a double



function edMu_Callback(hObject, eventdata, handles)
% hObject    handle to edMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edMu as text
%        str2double(get(hObject,'String')) returns contents of edMu as a double




function edSigma_Callback(hObject, eventdata, handles)
% hObject    handle to edSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSigma as text
%        str2double(get(hObject,'String')) returns contents of edSigma as a double




function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double




function edSampling_Callback(hObject, eventdata, handles)
% hObject    handle to edSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSampling as text
%        str2double(get(hObject,'String')) returns contents of edSampling as a double




function edRates_Callback(hObject, eventdata, handles)
% hObject    handle to edRates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRates as text
%        str2double(get(hObject,'String')) returns contents of edRates as a double



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double



function edNtraces_Callback(hObject, eventdata, handles)
% hObject    handle to edNtraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edNtraces as text
%        str2double(get(hObject,'String')) returns contents of edNtraces as a double





function edNframes_Callback(hObject, eventdata, handles)
% hObject    handle to edNframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edNframes as text
%        str2double(get(hObject,'String')) returns contents of edNframes as a double



function edTotalIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to edTotalIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edTotalIntensity as text
%        str2double(get(hObject,'String')) returns contents of edTotalIntensity as a double

I  = str2double(get(handles.edTotalIntensity,'String'));
IS = str2double(get(handles.edStdTotalIntensity,'String'));
SNR  = str2double(get(handles.edSNR,'String'));

S = I/(sqrt(2)*SNR);

% text = sprintf('Noise = %0.1f ± %0.1f',I/S,IS/S );
text = sprintf('Noise = %0.1f',S );
set(handles.txtNoise, 'String', text);


function edStdTotalIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to edStdTotalIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edStdTotalIntensity as text
%        str2double(get(hObject,'String')) returns contents of edStdTotalIntensity as a double
I  = str2double(get(handles.edTotalIntensity,'String'));
IS = str2double(get(handles.edStdTotalIntensity,'String'));
SNR  = str2double(get(handles.edSNR,'String'));

S = I/(sqrt(2)*SNR);

% text = sprintf('Noise = %0.1f ± %0.1f',I/S,IS/S );
text = sprintf('Noise = %0.1f',S );
set(handles.txtNoise, 'String', text);


function edSNR_Callback(hObject, eventdata, handles)
% hObject    handle to edSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSNR as text
%        str2double(get(hObject,'String')) returns contents of edSNR as a double

I  = str2double(get(handles.edTotalIntensity,'String'));
IS = str2double(get(handles.edStdTotalIntensity,'String'));
SNR  = str2double(get(handles.edSNR,'String'));

S = I/(sqrt(2)*SNR);

% text = sprintf('Noise = %0.1f ± %0.1f',I/S,IS/S );
text = sprintf('Noise = %0.1f',S );
set(handles.txtNoise, 'String', text);



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double



function edRandomSeed_Callback(hObject, eventdata, handles)
% hObject    handle to edRandomSeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRandomSeed as text
%        str2double(get(hObject,'String')) returns contents of edRandomSeed as a double



function edTotalTon_Callback(hObject, eventdata, handles)
% hObject    handle to edTotalTon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edTotalTon as text
%        str2double(get(hObject,'String')) returns contents of edTotalTon as a double



function edTon_Callback(hObject, eventdata, handles)
% hObject    handle to edTon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edTon as text
%        str2double(get(hObject,'String')) returns contents of edTon as a double



function edToff_Callback(hObject, eventdata, handles)
% hObject    handle to edToff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edToff as text
%        str2double(get(hObject,'String')) returns contents of edToff as a double



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edTotalIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edTotalIntensity as text
%        str2double(get(hObject,'String')) returns contents of edTotalIntensity as a double




function stdTotalIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to edStdTotalIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edStdTotalIntensity as text
%        str2double(get(hObject,'String')) returns contents of edStdTotalIntensity as a double




function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSNR as text
%        str2double(get(hObject,'String')) returns contents of edSNR as a double


