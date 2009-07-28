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

% Last Modified by GUIDE v2.5 24-Jul-2009 09:11:02

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
startProb = eval( get(handles.edStartProb,'String') );

totalTon = eval( get(handles.edTotalTon,'String') );
% Ton      = str2double( get(handles.edTon,      'String') );
% Toff     = str2double( get(handles.edToff,     'String') );

totalIntensity    = str2double( get(handles.edTotalIntensity,   'String') );
stdTotalIntensity = str2double( get(handles.edStdTotalIntensity,'String') );
snr = str2double( get(handles.edSNR,  'String') );
stdBackground = totalIntensity/(sqrt(2)*snr);
stdPhoton    = str2double( get(handles.edStdPhoton,   'String') );

randomSeed = str2double( get(handles.edRandomSeed,'String') );


%
if numel(sigma)==1,
    sigma = repmat(sigma,size(mu));
end


% Compile array of optional arguments to simulate...
options = { 'totalIntensity',totalIntensity, ...
            'stdTotalIntensity',stdTotalIntensity, ...
            'stdBackground',stdBackground, ...
            'stdPhoton',stdBackground, ...
            'randomSeed',randomSeed, ...
            'startProb',startProb, ...
            'kBleach', 1/totalTon };


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
time = 1000*sampling*[0:(traceLen-1)];

% Save resulting data
fname_qub = strrep(fname_txt, '.txt', '.qub.txt');
fname_trc = strrep(fname_txt, '.txt', '.traces');

% saveTraces( fname_qub, 'qub', fret );
% saveTraces( fname_txt, 'txt', donor,acceptor,fret );
saveTraces( fname_trc, 'traces', donor,acceptor, [], time );

[d,a,f,ids,time] = loadTraces(fname_trc);
saveTraces( fname_txt, 'txt', d,a,f,ids,time );

% Save the underlying state trajectory
fname_idl = strrep(fname_txt, '.txt', '.sim.dwt');

offsets = (0:(nTraces-1))*traceLen;
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



function edStdPhoton_Callback(hObject, eventdata, handles)
% hObject    handle to edStdPhoton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edStdPhoton as text
%        str2double(get(hObject,'String')) returns contents of edStdPhoton as a double


% --- Executes during object creation, after setting all properties.
function edStdPhoton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edStdPhoton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edGamma_Callback(hObject, eventdata, handles)
% hObject    handle to edGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edGamma as text
%        str2double(get(hObject,'String')) returns contents of edGamma as a double


% --- Executes during object creation, after setting all properties.
function edGamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edStdPhoton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edStdPhoton as text
%        str2double(get(hObject,'String')) returns contents of edStdPhoton as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edStdPhoton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnMoviesBrowse.
function btnMoviesBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to btnMoviesBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
    
% If unchecking, disable them.
else
    

end


% --- Executes on button press in chkSimulateMovies.
function chkSimulateMovies_Callback(hObject, eventdata, handles)
% hObject    handle to chkSimulateMovies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSimulateMovies


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


