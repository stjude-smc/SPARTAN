function varargout = bw_gui(varargin)
% BW_GUI M-file for bw_gui.fig
%      BW_GUI, by itself, creates a new BW_GUI or raises the existing
%      singleton*.
%
%      H = BW_GUI returns the handle to a new BW_GUI or the handle to
%      the existing singleton*.
%
%      BW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BW_GUI.M with the given input arguments.
%
%      BW_GUI('Property','Value',...) creates a new BW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bw_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bw_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% TODO:
%  - Add input for names of each state, so that rates.txt output files
%    have appropriate headings for each readability...
%  - WARNING: if photoressurection rates are set to be equal, then
%    all states are equally likely after a blink!  This effectively brings
%    the system out of equilibrium!!!

% Edit the above text to modify the response to help bw_gui

% Last Modified by GUIDE v2.5 28-Feb-2009 16:02:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bw_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @bw_gui_OutputFcn, ...
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


% --- Executes just before bw_gui is made visible.
function bw_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bw_gui (see VARARGIN)

% Choose default command line output for bw_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bw_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bw_gui_OutputFcn(hObject, eventdata, handles) 
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

tic;

clear handles.results;
clear handles.errorResults;

%------ Gather parameter values from GUI
nStates   = str2double( get(handles.edNstates,  'String') );
sampling = str2double( get(handles.edSampling,'String') );
sampling = sampling/1000; %to seconds

% useMu    = get(handles.chkMu,'Value');
% useSigma = get(handles.chkSigma,'Value');
useRates       = get(handles.chkRates,'Value');
% useInitialProb = get(handles.chkInitialProb,'Value');

maxItr = str2double( get(handles.edMaxItr,'String') );
LLConv = str2double( get(handles.edLLConv,'String') );
bootstrapN = str2double( get(handles.edBootstrapN,'String') );

reestMu    = get(handles.chkReestMu,   'Value');
reestSigma = get(handles.chkReestSigma,'Value');


% Generate Q matrix
if useRates,
    rates = eval( get(handles.edRates,'String') );
    
    Q = zeros(nStates,nStates);
    idx = find( ~logical(eye(nStates)) );
    r = rates';
    Q(idx) = r(:)+eps;
    Q = Q';

    A = Q*sampling;
end

% For now, these are required parameters...
% if useMu,
    mu = eval( get(handles.edMu, 'String') );
% end

% if useSigma,
    sigma = eval( get(handles.edSigma,'String') );
% end

if numel(sigma)==1,
    sigma = repmat(sigma,size(mu));
end


% Check parameter values for sanity?
% if numel(mu)~=nStates || numel(sigma)~=nStates,
%     error('simulate_gui:btnSimulate_Callback','Inconsistent number of states...');
% end



% Compile array of optional arguments to simulate...
fixMu    = repmat( 1-reestMu,    1,nStates );
fixSigma = repmat( 1-reestSigma, 1,nStates );

options = { 'maxItr',maxItr, 'LLConv',LLConv, ...
            'FixMu',fixMu, 'FixSigma',fixSigma, ...
            'bootstrapN',bootstrapN };



%------ Prompt use for location to save file in...
[f,p] = uigetfile('*.txt','Select datafile(s) to load','MultiSelect','on');
if p==0, return; end  %user pressed "cancel"
fname_txt = strcat(p,f);

% s = load('files.mat');
% fname_txt = s.files;


if ~iscell(fname_txt),
    fname_txt = {fname_txt};
end
nFiles = numel(fname_txt);

%----- Run BW optimizer on each file, combine results
BWparameters = { mu, sigma, options{:} };

h = waitbar(0,'Running Baum-Welch...');

for i=1:nFiles
    filename = fname_txt{i};
    
    sprintf('%d: %s', i,filename);
        
    % Estimate parameter values using Baum's method
    if bootstrapN>1
        [results(i),errorResults(i)] = runBW( ...
            filename, sampling, BWparameters);
    else
        results(i) = runBW( ...
            filename, sampling, BWparameters);
    end
    
    disp( results(i).mu );
    
    % Idealize data using Viterbi algorithm
    model = [results(i).mu' results(i).sigma'];
    p0 = results(i).p0;
    A = results(i).A;
    
    [d,a,traces] = loadTraces(filename);
    
    [dwt,offsets] = idealize( traces, model, p0, A );
    
    dwtFilename = strrep(filename,'.txt','.qub.dwt');
    saveDWT( dwtFilename, dwt, offsets, model, 1000*sampling );
    
    waitbar(i/nFiles,h);
end
close(h);

% Save results
handles.results = results;
save( [p 'results.mat'], 'results');

if bootstrapN>1
    handles.errorResults = errorResults;
    save( [p 'errorResults.mat'], 'errorResults');
end


% Enable GUI elements for saving results.
% Update handles structure
guidata(hObject, handles);

disp(toc);





function [results,errorResults] = runBW( ...
                filename, sampling, BWparameters )
%

[d,a,observations] = loadTraces(filename);
clear d; clear a;


%------ Run Baum-Welch optimization proceedure
[results,errorResults] = BWoptimize( observations, BWparameters{:} );

nStates = size(results.A,1);

framerate = 1/sampling;

% 1->[2,3,4], 2->[1,3,4], 3->[1,2,4], 4->[1,2,3]
Q = framerate.*results.A';
idx = find( ~logical(eye(nStates)) );
Q_parts = Q( idx );

Q_nonzero = Q(2:end,2:end);
idx_nonzero = find( ~logical(eye(nStates-1)) );
Q_nonzero = Q_nonzero(idx_nonzero);


% For now just print the whole thing to the console... see above
results.Q = Q;
results.rates = Q_parts;
results.nonzeroRates = Q_nonzero;

if ~isempty(errorResults),
    Q = framerate.*errorResults.A';
    Q_parts = Q( idx );
    errorResults.rates = Q_parts;
    errorResults.Q = Q;
    
    Q_nonzero = Q(2:end,2:end);
    Q_nonzero = Q_nonzero(idx_nonzero);
    errorResults.nonzeroRates = Q_nonzero;
end

[p,f,e] = fileparts(filename);
results.fname = strrep(f,'_auto','');





%% -------------------  CALLBACK FUNCTIONS  -------------------


% --- Executes on button press in btnSaveRates.
function btnSaveRates_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveRates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 

saveFilename = uiputfile('rates.txt','Save rate matrix as...');
if saveFilename==0, return; end


nonZeroOnly = get(handles.chkNonZeroOnly,'Value');
if nonZeroOnly
    start=2;
else
    start=1;
end


% Extract names for each datafile/condition
results = handles.results;
filenames = {results.fname};

for i=1:length(filenames)
    [p,f] = fileparts(filenames{i});
    filenames{i} = strrep(f,'_auto','');
end

nStates = numel( results(1).mu );

% Extract rates from results matrix...
if nonZeroOnly
    rates = [results.nonzeroRates]';
else
    rates = [results.rates]';
end
[nFiles,nRates] = size(rates);

% if isfield(handles,'errorResults'),
%     errorResults = handles.errorResults;
%     rateErrors = [errorResults.nonzeroRates]';
% else
    rateErrors = zeros( size(rates) );
% end

% Combine rates and errors into a single matrix
output = zeros(nFiles,nRates*2);
output(:,1:2:end) = rates;
output(:,2:2:end) = rateErrors;

% Construct names for each of the rates
rateNames = cell(0,1);

for i=start:nStates,
    for j=start:nStates
        if i==j, continue; end
        rateNames{end+1} = sprintf('k%d->%d', i,j );
        rateNames{end+1} = ['d' rateNames{end}];
    end
end


% Save results to file with headers
fid = fopen(saveFilename,'w');

header = {'Dataset',rateNames{:}};
fprintf(fid, '%s\t', header{:});

for i=1:nFiles,
    fprintf(fid, '\n%s',   filenames{i});
    fprintf(fid, '\t%.4f', output(i,:));
end


fclose(fid);







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



function edMaxItr_Callback(hObject, eventdata, handles)
% hObject    handle to edMaxItr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edMaxItr as text
%        str2double(get(hObject,'String')) returns contents of edMaxItr as a double


function edStdTotalIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to edStdTotalIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edStdTotalIntensity as text
%        str2double(get(hObject,'String')) returns contents of edStdTotalIntensity as a double



function edLLConv_Callback(hObject, eventdata, handles)
% hObject    handle to edLLConv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edLLConv as text
%        str2double(get(hObject,'String')) returns contents of edLLConv as a double



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
% hObject    handle to edMaxItr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edMaxItr as text
%        str2double(get(hObject,'String')) returns contents of edMaxItr as a double




function stdTotalIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to stdTotalIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stdTotalIntensity as text
%        str2double(get(hObject,'String')) returns contents of stdTotalIntensity as a double




function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edLLConv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edLLConv as text
%        str2double(get(hObject,'String')) returns contents of edLLConv as a double


% --- Executes on button press in chkReestMu.
function chkReestMu_Callback(hObject, eventdata, handles)
% hObject    handle to chkReestMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkReestMu


% --- Executes on button press in chkReestSigma.
function chkReestSigma_Callback(hObject, eventdata, handles)
% hObject    handle to chkReestSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkReestSigma




% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edSteadyState_Callback(hObject, eventdata, handles)
% hObject    handle to edSteadyState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSteadyState as text
%        str2double(get(hObject,'String')) returns contents of edSteadyState as a double



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edSampling as text
%        str2double(get(hObject,'String')) returns contents of edSampling as a double



function edFinalLL_Callback(hObject, eventdata, handles)
% hObject    handle to edFinalLL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edFinalLL as text
%        str2double(get(hObject,'String')) returns contents of edFinalLL as a double



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double



function edBootstrapN_Callback(hObject, eventdata, handles)
% hObject    handle to edBootstrapN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edBootstrapN as text
%        str2double(get(hObject,'String')) returns contents of edBootstrapN as a double


% --- Executes on button press in chkMu.
function chkMu_Callback(hObject, eventdata, handles)
% hObject    handle to chkMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkMu


% --- Executes on button press in chkSigma.
function chkSigma_Callback(hObject, eventdata, handles)
% hObject    handle to chkSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkSigma


% --- Executes on button press in chkRates.
function chkRates_Callback(hObject, eventdata, handles)
% hObject    handle to chkRates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkRates


% --- Executes on button press in chkInitialProb.
function chkInitialProb_Callback(hObject, eventdata, handles)
% hObject    handle to chkInitialProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkInitialProb


% --- Executes on button press in chkNonZeroOnly.
function chkNonZeroOnly_Callback(hObject, eventdata, handles)
% hObject    handle to chkNonZeroOnly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkNonZeroOnly


