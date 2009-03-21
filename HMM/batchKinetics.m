function varargout = batchKinetics(varargin)
% BATCHKINETICS M-file for batchKinetics.fig
%      BATCHKINETICS, by itself, creates a new BATCHKINETICS or raises the existing
%      singleton*.
%
%      H = BATCHKINETICS returns the handle to a new BATCHKINETICS or the handle to
%      the existing singleton*.
%
%      BATCHKINETICS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BATCHKINETICS.M with the given input arguments.
%
%      BATCHKINETICS('Property','Value',...) creates a new BATCHKINETICS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before batchKinetics_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to batchKinetics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help batchKinetics

% Last Modified by GUIDE v2.5 13-Mar-2009 01:47:25


%% GUI Callbacks

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @batchKinetics_OpeningFcn, ...
                   'gui_OutputFcn',  @batchKinetics_OutputFcn, ...
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


% --- Executes just before batchKinetics is made visible.
function batchKinetics_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to batchKinetics (see VARARGIN)

% Choose default command line output for batchKinetics
handles.output = hObject;

handles.hasModel = 0;
handles.hasData = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes batchKinetics wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = batchKinetics_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edSampling_Callback(hObject, eventdata, handles)
% hObject    handle to edSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.sampling = str2double( get(handles.edSampling,'String') );
guidata(hObject, handles);


function edBootstrapN_Callback(hObject, eventdata, handles)
% hObject    handle to edBootstrapN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% bootstrapN = str2double( get(hObject,'String') );
% if boostrapN<1,
%     warning('
%     set(hObject,'String',1'


%% CALLBACKS: Loading Data...

% --- Executes on button press in btnLoadData.
function btnLoadData_Callback(hObject, eventdata, handles)

%------ Prompt use for location to save file in...
% FIXME: use a custom dialog that allows selecting many files
%  from disparate locations...
[f,p] = uigetfile('*.txt','Select datafile(s) to analyze','MultiSelect','on');
if p==0, return; end  %user pressed "cancel"
fname_txt = strcat(p,f);

if ~iscell(fname_txt),
    fname_txt = {fname_txt};
end
fname_txt = sort(fname_txt)
handles.dataFilenames = fname_txt;
handles.dataPath = p;
handles.hasData = 1;

% If a model is already loaded, enable the Execute button
if handles.hasModel,
    set(handles.btnExecute,'Enable','on');
end

guidata(hObject, handles);


%% CALLBACKS: Loading Model...

% --- Executes on button press in btnBrowseModel.
function btnBrowseModel_Callback(hObject, eventdata, handles)
% hObject    handle to btnBrowseModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get list of data files from user
fname = get(handles.edModelFilename,'String');
if isempty(fname)
    fname = '*.qmf';
end

[f,p] = uigetfile(fname,'Select a QuB model file...');
if f==0, return; end
fname = [p f];

set(handles.edModelFilename,'String',fname);

handles = loadModel(handles,fname);
guidata(hObject, handles);



function handles = loadModel(handles,filename)

model = qub_loadModel( filename );
handles.model = model;
handles.hasModel = 1;

% Update model status text -- FIXME
info = sprintf('FRET values: %s',mat2str(model.mu));
set(handles.txtModelInfo,'String',info);

% If data are already loaded, enable the Execute button
if handles.hasData,
    set(handles.btnExecute,'Enable','on');
end



function edModelFilename_Callback(hObject, eventdata, handles)
% hObject    handle to edModelFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = get(handles.edModelFilename,'String');
if ~exist(filename,'file'),
    warning('Model file does not exist!');
else
    handles = loadModel(handles,filename);
end

guidata(hObject, handles);



%% CALLBACKS: Execution...


% --- Executes on button press in btnExecute.
function btnExecute_Callback(hObject, eventdata, handles)
% hObject    handle to btnExecute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.hasModel || ~handles.hasData,
    set(handles.btnExecute,'Enable','off');
    warning('Missing model or data');
    return;
end

% Get location to save results...
f = [handles.dataPath 'result.qrf'];
[f,p] = uiputfile(f,'Save results as...');
if f==0, return; end
resultFilename = [p f];

% Get sampling interval of data...
sampling = str2double( get(handles.edSampling,'String') );
sampling = sampling/1000; %to seconds

% Get optimizer options
bootstrapN = floor( str2double(get(handles.edBootstrapN,'String')) );
options.bootstrapN = max(1,bootstrapN);

% Run the optimizer...
set(handles.btnExecute,'Enable','off');
set(handles.btnStop,'Enable','on');

% t0 = clock;
delete('resultTree.mat');
delete('mil_result.qtr');
delete('result.qmf');
delete('result.qrf');
resultTree = runParamOptimizer(handles.model,handles.dataFilenames, ...
                                        sampling,options);

% Save results
% Users can see visually plot the results using another program. FIXME
save('resultTree.mat','resultTree');
qub_saveTree(resultTree,resultFilename);
qub_saveTree(resultTree.milResults(1).ModelFile,'result.qmf','ModelFile');

% Finish up...
% et = etime(t0,clock);
% disp(sprintf('Run time: %d:%d:%d (h:m:s)',et(4:6)));

set(handles.btnStop,'Enable','off');
set(handles.btnExecute,'Enable','on');

disp('Finished!');




% --- Executes on button press in btnStop.
function btnStop_Callback(hObject, eventdata, handles)
% hObject    handle to btnStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set a variable somewhere to tell algorithm its time to stop....
% FIXME!!!

% Allow the user to start again...
set(handles.btnExecute,'Enable','on');





%% Parameter optimization engine...

function resultTree = runParamOptimizer( model,dataFilenames,sampling,options)

nFiles = numel(dataFilenames);
tic;

%----- STEP 1: Run Baum-Welch to optimize FRET parameters and estimate kinetics
h = waitbar(0,'Optimizing FRET model using Baum-Welch...');

bwOptions.maxItr = 100;
bwOptions.convLL = 1e-2;
    
for i=1:nFiles
    filename = dataFilenames{i};
    sprintf('%d: %s', i,filename);
    
    % Estimate parameter values using Baum's method
    bwResults(i) = runBW( filename, sampling, model, bwOptions );
    
    waitbar(0.33*i/nFiles,h);
end

toc

%----- STEP 2: Idealize data using Baum-Welch parameters, save to disk
h = waitbar(0.33,h,'Idealizing data...');

for i=1:nFiles
    filename = dataFilenames{i};
    
    % Idealize data using Viterbi algorithm
    fretModel = [bwResults(i).mu' bwResults(i).sigma'];
    p0 = bwResults(i).p0;
    A = bwResults(i).A;
    
    [d,a,traces] = loadTraces(filename);
    
    [dwt,offsets] = idealize( traces, fretModel, p0, A );
    
    dwtFilename{i} = strrep(filename,'.txt','.qub.dwt');
    saveDWT( dwtFilename{i}, dwt, offsets, fretModel, 1000*sampling );
    
    waitbar(0.33+0.33*i/nFiles,h);
end

toc

%----- STEP 3: Refine kinetic parameter estimates using MIL
waitbar(0.66,h,'Refining kinetic model using QuB...');

% Eventually, we want to use the Baum-Welch optimized model
% as a starting point for QuB... FIXME...
% Save starting model to temporary location...
% mfname = [tempname '.qmf'];
mfname = 'bwmodel.qmf';
delete(mfname);
qub_saveTree( model.qubTree, mfname, 'ModelFile' );
% qub_saveModel(model,mfname);

for i=1:nFiles
    disp( sprintf('%d: %s', i,dwtFilename{i}) );
    
    % Save optimized BW results as a model for MIL
%     md = bwResults(i);
%     md.qubTree = model.qubTree;
%     md.rates = md.Q;  % HACK
%     md = qub_saveModel(md,mfname);
    
    % Run MIL
    milResults(i) = qub_milOptimize(dwtFilename{i},mfname);
    
    waitbar(0.66+0.33*i/nFiles,h);
end

%delete(mfname);

toc

%----- STEP 4: Compile results into a qubtree
waitbar(1,h,'Saving results...');

resultTree.bwResults = bwResults;
resultTree.milResults = milResults;


close(h)


















function [results,errorResults] = runBW( ...
                filename, sampling, model, BWparameters )
%

[d,a,observations] = loadTraces(filename);
clear d; clear a;


%------ Run Baum-Welch optimization proceedure
[results,errorResults] = BWoptimize( ...
                            observations, sampling, model, BWparameters );

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

[p,f] = fileparts(filename);
results.fname = strrep(f,'_auto','');









