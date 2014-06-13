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

% Last Modified by GUIDE v2.5 15-Dec-2009 15:43:18


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

% FIXME: these should either be SET in the GUI or obtained from the GUI so
% the .fig file always matches these settings when batchKinetics is first
% loaded.
options.sampling = 0.04;
options.bootstrapN = 1;
options.deadTime = 0.5;
options.idealizeMethod = 'Segmental k-means';
options.seperately = 1; %SKM: analyze each trace individually
options.kineticsMethod = 'MIL Together';
handles.options = options;

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




%% CALLBACKS: Loading Data...

% --- Executes on button press in btnLoadData.
function btnLoadData_Callback(hObject, eventdata, handles)

%------ Prompt use for location to save file in...
handles.dataFilenames = getFiles([],'Select traces files to analyze');
handles.dataPath = pwd;
handles.hasData = 1;

% If a model is loaded, enable the Execute button & update GUI
if handles.hasModel,
    set(handles.btnExecute,'Enable','on');
    
    text = sprintf('%d dataset files loaded.',numel(handles.dataFilenames));
    set(handles.txtFileInfo,'String',text);
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
    constants = cascadeConstants;
    fname = [constants.modelLocation filesep '*.qmf'];
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
info = sprintf('FRET values: %s\n',mat2str(model.mu(model.class)));
info = [ info sprintf('FRET stdev: %s\n',mat2str(model.sigma(model.class))) ];
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
% Run the data analysis pipeline with user-specified data & model.

% Verify data and model have been specified by user in GUI.
if ~handles.hasModel || ~handles.hasData,
    set(handles.btnExecute,'Enable','off');
    warning('Missing model or data');
    return;
end

% Process analysis parameters from GUI
options  = handles.options;
model = handles.model;

if isfield(options,'fixFret'),
    assert( all(options.fixFret<=model.nStates) );
    model.fixMu = zeros(model.nStates,1);
    model.fixMu( options.fixFret ) = 1;
end
if isfield(options,'fixStdev'),
    assert( all(options.fixStdev<=model.nStates) );
    model.fixSigma = zeros(model.nStates,1);
    model.fixSigma( options.fixStdev ) = 1;
end

% Update GUI for "Running" status.
set(handles.btnExecute,'Enable','off');
set(handles.btnStop,'Enable','on');

% Run the analysis algorithms...
resultTree = runParamOptimizer(model,handles.dataFilenames,options);

% Save results to file for later processing by the user.
save('resultTree.mat','resultTree');
% qub_saveTree(resultTree,resultFilename);
% qub_saveTree(resultTree.milResults(1).ModelFile,'result.qmf','ModelFile');

% Update GUI for finished status.
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
%  ========================================================================

function resultTree = runParamOptimizer( model,dataFilenames,options)

resultTree = struct([]);

% disp(options);
% disp(model);
% 
% return;


h = waitbar(0,'Initializing...');

% ...
sampling = options.sampling;

nFiles = numel(dataFilenames);


% Setup algorithm settings
skmOptions.maxItr = 40;
skmOptions.convLL = 1e-4;
skmOptions.seperately = options.seperately;
if skmOptions.seperately,
    skmOptions.quiet = 1;
end

bwOptions = skmOptions;
thresholdOptions = struct([]);


% Remove intermediate files from previous runs.
warning('off','MATLAB:DELETE:FileNotFound');
delete('resultTree.mat');
delete('mil_result.qtr');
delete('result.qmf');
delete('result.qrf');
delete('bwmodel.qmf');


if ~strcmp(options.idealizeMethod,'Do Nothing'),
    %----- STEP 1: Optimize params using SKM and idealize
    waitbar(0,h,'Idealizing data using Segmental k-means...');
    
    skmLL = zeros(nFiles,1);

    for i=1:nFiles
        filename = dataFilenames{i};
        sprintf('%d: %s', i,filename);

        % Load data
        d = loadTraces(dataFilenames{i});
        data = d.fret;

        % Idealize data using user-specified algorithm...
        if strcmp(options.idealizeMethod,'Segmental k-means'),
            
            [dwt,optModel,LL,offsets] = skm( data, sampling, model, skmOptions );
            skmLL(i) = LL(end);
            
        elseif strcmp(options.idealizeMethod,'Baum-Welch'),
            
            result = BWoptimize( data, sampling, model, bwOptions );
            fretModel = [model.mu' model.sigma'];
            optModel = model;
            % TODO: update optModel with optimized parameter values.
            
            [dwt,idl,offsets,LL] = idealize( ...
                    data, fretModel, result.p0, result.A );
            skmLL(i) = mean(LL);
            
        elseif  strcmp(options.idealizeMethod,'Thresholding'),
            
            [dwt,offsets,optModel] = tIdealize( data, model, thresholdOptions );
            skmLL(i) = 0;
            
        end
        
        % Remove final zero-state dwell, if it exists.
        % These dwells confuse MIL.
        keep = ones(numel(dwt),1);
        for j=1:numel(dwt),
            states = dwt{j}(:,1);
            times  = dwt{j}(:,2);
            
            % Remove last dwell if in dark state. These dwells result from the
            % photobleached state, which is not considered in kinetic analysis.
            if numel(states)>0 && states(end)==1,
                states = states(1:end-1);
                times  = times(1:end-1);
            end
            
            % Remove last dwell, which is cut short due to photobleaching.
            % This prevents bias in kinetic parameter estimation because
            % the last dwell is frequently an artificial "step" from
            % high FRET to the dark state.
            if numel(states)<1,
                keep(j) = 0;
            else
                if times(end)<=1
                    states = states(1:end-1);
                    times  = times(1:end-1);
                end
            end
            
            % Save changes, marking empty traces for removal (keep=0)
            if numel(states)==0,
                keep(j) = 0;
            else
                dwt{j} = [states times];
            end
        end
        dwt = dwt( logical(keep) );
        offsets = offsets( logical(keep) );
        
        % Save results.
        % NOTE: there is no one result when each trace is optimized seperately.
        % An "average" model is extracted instead.
        if skmOptions.seperately,
            skmModels(i).mu    = mean( [optModel.mu], 2 );
            skmModels(i).sigma = mean( [optModel.sigma], 2 );
            skmModels(i).LL    = mean( LL, 2 );
        else
            skmModels(i) = optModel(1);
        end

        % Save the idealization
        [p,n] = fileparts(filename);
        dwtFilename = fullfile( p, [n '.qub.dwt'] );
        fretModel = [skmModels(i).mu skmModels(i).sigma];
        saveDWT( dwtFilename, dwt, offsets, fretModel, 1000*sampling );

        waitbar(0.33*i/nFiles,h);
        drawnow;
    end

    % Save results
    resultTree(1).skmModels = skmModels;
    resultTree.skmLL = skmLL;
end


% If no further action is neccessary, exit.
if strcmp(options.kineticsMethod,'Do Nothing'),
    close(h);
    return;
end




%----- STEP 2: Refine kinetic parameter estimates using MIL
waitbar(0.33,h,'Refining kinetic model using QuB...');

% Eventually, we want to use the Baum-Welch optimized model
% as a starting point for QuB... FIXME...
% Save starting model to temporary location...
% mfname = [tempname '.qmf'];
mfname = 'bwmodel.qmf';
delete(mfname);
qub_saveTree( model.qubTree, mfname, 'ModelFile' );

avgRates = cell(nFiles,1);
avgStdRates = cell(nFiles,1);

for i=1:nFiles
    [p,n] = fileparts(dataFilenames{i});
    dwtFilename = fullfile( p, [n '.qub.dwt'] );
    disp( sprintf('%d: %s', i,dwtFilename) );
    
    % Verify the idealization has been performed
    if ~exist(dwtFilename,'file'),
        error('File not idealized. Select an idealization method.');
    end

    % Run MIL
    waitbarBounds = 0.33+0.66*[(i-1) i]/nFiles;
    [avgRates{i},avgStdRates{i},milResults(i)] = bootstrapMIL( ...
                dwtFilename, mfname, options.bootstrapN );

    waitbar(waitbarBounds(end),h);
    drawnow;
end

% Save results
resultTree(1).milResults = milResults;


%----- STEP 3: Compile results into a qubtree
waitbar(1,h,'Saving results...');


rates = [];
stdRates = [];
nStates = model.nStates;

for i=1:nFiles
%     modelTree = milResults(i).ModelFile;
%     model = qub_loadModel( modelTree );
%     Q = model.rates';
    Q    = avgRates{i}';
    Qstd = avgStdRates{i}';

    idx         = find( ~logical(eye(nStates)) );
    idx_nonzero = find( ~logical(eye(nStates-1)) );

    % 1->[2,3,4], 2->[1,3,4], 3->[1,2,4], 4->[1,2,3]
    Q_parts   = Q(idx)';
    Q_nonzero = Q(2:end,2:end);
    Q_nonzero = Q_nonzero(idx_nonzero)';

    Qstd_parts   = Qstd(idx)';
    Qstd_nonzero = Qstd(2:end,2:end);
    Qstd_nonzero = Qstd_nonzero(idx_nonzero)';

    rates    = [rates    ; Q_nonzero   ];
    stdRates = [stdRates ; Qstd_nonzero];
end
nRates = size(rates,2);

% Combine rates and errors into a single matrix
output = zeros(nFiles,nRates*2);
output(:,1:2:end) = rates;
output(:,2:2:end) = stdRates;

% Construct names for each of the rates
rateNames = cell(0,1);
for i=2:nStates,
    for j=2:nStates
        if i==j, continue; end
        rateNames{end+1} = sprintf('k%d->%d', i,j );
        rateNames{end+1} = ['d' rateNames{end}];
    end
end

% Save results to file with headers
fid = fopen('rates.txt','w');

header = {'Dataset',rateNames{:}};
fprintf(fid, '%s\t', header{:});

for i=1:nFiles,
    [p,fname] = fileparts( dataFilenames{i} );
    fprintf(fid, '\n%s',  fname );
    fprintf(fid, '\t%.4f', output(i,:));
end
fclose(fid);

resultTree.rates = rates;

close(h);


%-------------------- END FUNCTION runParamOptimizer ---------------------









%%
function [avgRates,stdRates,firstResult] = bootstrapMIL( ...
                          dwtFilename, modelFilename, nBootstrap )

if nargin < 3,
    nBootstrap = 1;
end


% Load the dwell-time data and the model.
[dwt,sampling,offsets,fret_model] = loadDWT( dwtFilename );
nTraces = numel(dwt);

nStates = numel(fret_model)/2;
X = eye(nStates);
X(1,:) = 1; X(:,1) = 1;
idx_nonzero = find( ~logical(X) );


% Construct a set of bootstrap samples of the dwell-time data.
% These samples are then processed by MIL in parallel using the
% jobQueue function.
% idxBootstrap = cell(nBootstrap,1);
tempDwtNames = cell(nBootstrap,1);

for i=1:nBootstrap,
    % The first bootstrap sample is all data (no sampling) - this insures
    % that the "mean" result is always the same, as expected.
    if i == 1,
        idxBootstrap = 1:nTraces;
    else
        idxBootstrap = floor(rand(nTraces,1)*nTraces)+1;
    end
    
    tempDwtNames{i} = [tempname '.dwt'];
    saveDWT( tempDwtNames{i}, dwt(idxBootstrap), ...
             offsets(idxBootstrap), fret_model, sampling );
end

% Run MIL
result = qub_milOptimize( tempDwtNames, modelFilename );

% Process the results.
bootstrapRates = [];

for i=1:nBootstrap,
    
    % Save the first result as "the" result; others are only informative of
    % the error in the analysis.
    if i == 1,
        firstResult = result(i);
    end
    
    % Process the results for averaging
    modelTree = result(i).ModelFile;
    model = qub_loadModel( modelTree );
    rates = model.rates(:);
    
    disp( rates(:)' );
    if any( rates(idx_nonzero) < 10e-9 )
        warning( 'Key rate estimated as zero. Ignoring result.' );
        continue;
    end

    if any( rates(~logical(eye(nStates))) > 2*1000/sampling )
        warning( 'Rate estimate way out of range. Ignoring.' );
        continue;
    end
    
    bootstrapRates(i,:) = rates(:);
end

% Calculate average rates
avgRates = bootstrapRates(1,:);
% avgRates = mean(bootstrapRates,1);
stdRates = std(bootstrapRates,0,1);

if nBootstrap==1,
    stdRates = zeros( size(avgRates) );
end

avgRates = reshape( avgRates, size(model.rates) );
stdRates = reshape( stdRates, size(model.rates) );

save( [dwtFilename '.ratedata.txt'], 'bootstrapRates', '-ASCII' );


% Delete temporary data to save disk space.
for i=1:numel(tempDwtNames),
    delete( tempDwtNames{i} );
end
delete('.milresult*');




%% Other GUI Callbacks
%  ========================================================================



% --- Executes on button press in chkFixFret.
function chkFixFret_Callback(hObject, eventdata, handles)
% Idealization options: Fix FRET values checkbox

state = get(hObject,'Value');

if ~state,
    % unchecked - remove setting
    set(handles.edFixFret,'Enable','off');
    if isfield(handles.options,'fixFret')
        handles.options = rmfield( handles.options, 'fixFret' );
    end
else
    % checked - restore previous settings
    set(handles.edFixFret,'Enable','on');
    text = get(handles.edFixFret,'String');
    handles.options.fixFret = str2num(text);
end
guidata(hObject, handles);


% --- Executes on button press in chkFixStdev.
function chkFixStdev_Callback(hObject, eventdata, handles)
% Idealization options: Fix FRET standard deviations checkbox

state = get(hObject,'Value');

if ~state,
    % unchecked - remove setting
    set(handles.edFixStdev,'Enable','off');
    if isfield(handles.options,'fixStdev')
        handles.options = rmfield( handles.options, 'fixStdev' );
    end
else
    % checked - restore previous settings
    set(handles.edFixStdev,'Enable','on');
    text = get(handles.edFixStdev,'String');
    handles.options.fixStdev = str2num(text);
end
guidata(hObject, handles);


% --- Executes on selection change in cboIdealizationMethod.
function cboIdealizationMethod_Callback(hObject, eventdata, handles)
% Idealization options: idealization method combo box

% Update method to use for idealization
text = get(hObject,'String');
handles.options.idealizeMethod = text{get(hObject,'Value')};
guidata(hObject, handles);

% If user selected "Do Nothing", disable idealization option controls.
if get(hObject,'Value')==1,
    set(handles.edFixFret,'Enable','off');
    set(handles.edFixStdev,'Enable','off');
    set(handles.chkFixFret,'Enable','off');
    set(handles.chkFixStdev,'Enable','off');
else
    set(handles.chkFixFret,'Enable','on');
    set(handles.chkFixStdev,'Enable','on');
    % Update text box states
    chkFixFret_Callback( handles.chkFixFret, [], handles );
    chkFixStdev_Callback( handles.chkFixStdev, [], handles );
end
    
    
% --- Executes on selection change in cboKineticsMethod.
function cboKineticsMethod_Callback(hObject, eventdata, handles)
% Idealization options: idealization method combo box

% Update method to use for kinetic parameter estimation.
text = get(hObject,'String');
handles.options.kineticsMethod = text{get(hObject,'Value')};
guidata(hObject, handles);

% If user selected "Do Nothing", disable idealization option controls.
if get(hObject,'Value')==1,
    set(handles.edBootstrapN,'Enable','off');
    set(handles.edDeadTime,'Enable','off');
else
    set(handles.edBootstrapN,'Enable','on');
    set(handles.edDeadTime,'Enable','on');
end




function edSampling_Callback(hObject, eventdata, handles)
% Input data options: experimental sampling interval.
sampling = str2double( get(handles.edSampling,'String') );
handles.options.sampling = sampling/1000; %convert to seconds.
guidata(hObject, handles);


function edFixFret_Callback(hObject, eventdata, handles)
% Idealization options: fix FRET values to specified values.
handles.options.fixFret = str2num( get(hObject,'String') );
guidata(hObject, handles);


function edFixStdev_Callback(hObject, eventdata, handles)
% Idealization options: fix FRET standard deviationsto specified values.
handles.options.fixStdev = str2num( get(hObject,'String') );
guidata(hObject, handles);


function edBootstrapN_Callback(hObject, eventdata, handles)
% Kinetic estimation options: number of bootstrap samples to test

N = floor(str2double( get(hObject,'String') ));
N = max(1,N); %must be at least =1.
set(hObject,'String',num2str(N));
handles.options.bootstrapN = N;
guidata(hObject, handles);


function edDeadTime_Callback(hObject, eventdata, handles)
% Idealization options: fix FRET standard deviationsto specified values.
handles.options.deadTime = str2double( get(hObject,'String') );
guidata(hObject, handles);


% --- Executes on button press in chkIdealizeSeperately.
function chkIdealizeSeperately_Callback(hObject, eventdata, handles)
% Idealization options: fix FRET standard deviationsto specified values.
handles.options.seperately = get(hObject,'Value');
guidata(hObject, handles);
