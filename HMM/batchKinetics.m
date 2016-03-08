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

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% Last Modified by GUIDE v2.5 14-Aug-2015 17:46:22


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
function batchKinetics_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to batchKinetics (see VARARGIN)

updateSpartan; %check for updates

% Choose default command line output for batchKinetics
handles.output = hObject;

% Set initial internal state of the program
handles.modelFilename = fullfile(pwd,'*.qmf');
handles.model = [];
handles.dataFilenames = {};
handles.dwtFilenames  = {};

% Set default analysis settings. FIXME: put these in cascadeConstants?
options.bootstrapN = 1;
options.deadTime = 0.5;
options.seperately = 1; %SKM: analyze each trace individually
options.maxItr = 100;

handles.options = options;
guidata(hObject, handles);

% Update GUI to reflect these default settings.
set( handles.cboIdealizationMethod, 'Value',2 );  %SKM
cboIdealizationMethod_Callback(handles.cboIdealizationMethod,[],handles);

handles = guidata(hObject);
set( handles.cboKineticsMethod, 'Value',1 );  %Do nothing
cboKineticsMethod_Callback(handles.cboKineticsMethod,[],handles);

set( handles.edBootstrapN,          'String', options.bootstrapN );
set( handles.edDeadTime,            'String', options.deadTime   );
set( handles.chkIdealizeSeperately, 'Value',  options.seperately );
set( handles.edMaxIterations,       'String', options.maxItr     );
% set( handles.tblFixFret, 'Data', num2cell(false(3,2)) );

constants = cascadeConstants;
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );


% UIWAIT makes batchKinetics wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = batchKinetics_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




%% CALLBACKS: Loading Data...

% --- Executes on button press in btnLoadData.
function btnLoadData_Callback(hObject, ~, handles) %#ok<DEFNU>

%------ Prompt use for location to save file in...
handles.dataFilenames = getFiles([],'Select traces files to analyze');
handles.dataPath = pwd;

% If a model is loaded, enable the Execute button & update GUI
if ~isempty(handles.model),
    set(handles.btnExecute,'Enable','on');
end

text = sprintf('%d files loaded.',numel(handles.dataFilenames));
set(handles.txtFileInfo,'String',text);
set(handles.btnMakeplots,'Enable','on');

guidata(hObject, handles);


%% CALLBACKS: Loading Model...

% --- Executes on button press in btnLoadModel.
function btnLoadModel_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to btnLoadModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Ask the user for a filename
[fname,p] = uigetfile( handles.modelFilename, 'Select a QuB model file...' );
if fname==0, return; end
fname = fullfile(p,fname);
set( handles.edModelFilename, 'String',['...' fname(max(1,end-45):end)] );

% Load the model and show the model properties in the GUI.
% The model's properties are automatically updated whenever the model is
% modified in the GUI.
handles.model = QubModel(fname);
handles.model.showModel( handles.axModel );

% Enable saving the model
set( handles.btnSaveModel, 'Enable','on' );

% If data are already loaded, enable the Execute button
if ~isempty(handles.dataFilenames),
    set(handles.btnExecute,'Enable','on');
end

% Automatically update the parameter table when the model is altered.
handles.modelUpdateListener = addlistener(handles.model, ...
                            {'mu','fixMu','sigma','fixSigma'}, 'PostSet', ...
                            @(s,e)modelUpdate_Callback(handles.tblFixFret,e) );
handles.model.mu = handles.model.mu;  %trigger table update

% Update handles structure
guidata(hObject, handles);




%% CALLBACKS: Execution...


% --- Executes on button press in btnExecute.
function btnExecute_Callback(hObject, ~, handles) %#ok<DEFNU>
% Run the data analysis pipeline with user-specified data & model.

% Verify data and model have been specified by user in GUI.
if isempty(handles.model) || isempty(handles.dataFilenames),
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
[resultTree,handles.dwtFilenames] = runParamOptimizer(model,handles.dataFilenames,options); %#ok<ASGLU>

% Save results to file for later processing by the user.
save('resultTree.mat','resultTree');
% qub_saveTree(resultTree,resultFilename);
% qub_saveTree(resultTree.milResults(1).ModelFile,'result.qmf','ModelFile');

% Update GUI for finished status.
set(handles.btnStop,'Enable','off');
set(handles.btnExecute,'Enable','on');
set(handles.btnLifetimeExp,'Enable','on');
disp('Finished!');

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in btnStop.
function btnStop_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to btnStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set a variable somewhere to tell algorithm its time to stop....
% FIXME!!!

% Allow the user to start again...
set(handles.btnExecute,'Enable','on');





%% Parameter optimization engine...
%  ========================================================================

function [resultTree,dwtFilenames] = runParamOptimizer( model,dataFilenames,options)

resultTree = struct([]);

% disp(options);
% disp(model);
% 
% return;


h = waitbar(0,'Initializing...');

nFiles = numel(dataFilenames);
dwtFilenames = cell( size(dataFilenames) );

% Setup algorithm settings
skmOptions.maxItr = options.maxItr;
skmOptions.convLL = 1e-4;
skmOptions.seperately = options.seperately;
if skmOptions.seperately,
    skmOptions.quiet = 1;
end

bwOptions = skmOptions;
thresholdOptions = struct([]);


% Remove intermediate files from previous runs.
warning('off','MATLAB:DELETE:FileNotFound');
delete('resultTree.mat','mil_result.qtr','result.qmf','result.qrf','bwmodel.qmf');


if ~strcmp(options.idealizeMethod,'Do Nothing'),
    %----- STEP 1: Optimize params using SKM and idealize
    waitbar( 0, h, ['Idealizing data using ' options.idealizeMethod '...'] );
    
    skmLL = zeros(nFiles,1);

    for i=1:nFiles
        filename = dataFilenames{i};
        sprintf('%d: %s', i,filename);

        % Load data
        d = loadTraces(dataFilenames{i});
        data = d.fret;
        sampling = d.sampling;

        % Idealize data using user-specified algorithm...
        if strcmp(options.idealizeMethod,'Segmental k-means'),
            
            [dwt,optModel,LL,offsets] = skm( data, sampling, model, skmOptions );
            skmLL(i) = LL(end);
            
        elseif strcmp(options.idealizeMethod,'Baum-Welch'),
            
            result = BWoptimize( data, sampling, model, bwOptions );
            fretModel = [to_col(model.mu) to_col(model.sigma)];
            optModel = model;
            % TODO: update optModel with optimized parameter values.
            
            [dwt,~,offsets,LL] = idealize( ...
                    data, fretModel, result.p0, result.A );
            skmLL(i) = mean(LL);
            
        elseif  strcmp(options.idealizeMethod,'Thresholding'),
            
            [dwt,offsets] = tIdealize( data, model, thresholdOptions );
            optModel = model;
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
        dwtFilenames{i} = fullfile( p, [n '.qub.dwt'] );
        fretModel = [to_col(skmModels(i).mu) to_col(skmModels(i).sigma)];
        saveDWT( dwtFilenames{i}, dwt, offsets, fretModel, sampling );

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
    [~,fname] = fileparts( dataFilenames{i} );
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
[dwt,data.sampling,offsets,fret_model] = loadDWT( dwtFilename );
nTraces = numel(dwt);

nStates = numel(fret_model)/2;
X = eye(nStates);
X(1,:) = 1; X(:,1) = 1;
idx_nonzero = ~logical(X);


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
             offsets(idxBootstrap), fret_model, data.sampling );
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

    if any( rates(~logical(eye(nStates))) > 2*1000/data.sampling )
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

% --- Executes on selection change in cboIdealizationMethod.
function cboIdealizationMethod_Callback(hObject, ~, handles)
% Idealization options: idealization method combo box

% Update method to use for idealization
text = get(hObject,'String');
handles.options.idealizeMethod = text{get(hObject,'Value')};
guidata(hObject, handles);

% If user selected "Do Nothing", disable idealization option controls.
names = {'tblFixFret','chkIdealizeSeperately','edMaxIterations'};
if get(hObject,'Value')==1,
    for i=1:numel(names),
        set(handles.(names{i}),'Enable','off');
    end
else
    for i=1:numel(names),
        set(handles.(names{i}),'Enable','on');
    end
end
    
    
% --- Executes on selection change in cboKineticsMethod.
function cboKineticsMethod_Callback(hObject, ~, handles)
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


function edBootstrapN_Callback(hObject, ~, handles) %#ok<DEFNU>
% Kinetic estimation options: number of bootstrap samples to test

N = floor(str2double( get(hObject,'String') ));
N = max(1,N); %must be at least =1.
set(hObject,'String',num2str(N));
handles.options.bootstrapN = N;
guidata(hObject, handles);


function edDeadTime_Callback(hObject, ~, handles) %#ok<DEFNU>
% Idealization options: fix FRET standard deviationsto specified values.
handles.options.deadTime = str2double( get(hObject,'String') );
guidata(hObject, handles);


% --- Executes on button press in chkIdealizeSeperately.
function chkIdealizeSeperately_Callback(hObject, ~, handles) %#ok<DEFNU>
% Idealization options: fix FRET standard deviationsto specified values.
handles.options.seperately = get(hObject,'Value');
guidata(hObject, handles);



function edMaxIterations_Callback(hObject, ~, handles) %#ok<DEFNU>
% Maximum number of iterations boxed changed.
handles.options.maxItr = str2double( get(hObject,'String') );
guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in tblFixFret.
function tblFixFret_CellEditCallback(hObject, ~, handles) %#ok<DEFNU>
% Update QubModel object with new settings from the table

handles.modelUpdateListener.Enabled = false;

data = get(hObject,'Data');
handles.model.mu       = [data{:,1}];
handles.model.fixMu    = [data{:,2}];
handles.model.sigma    = [data{:,3}];
handles.model.fixSigma = [data{:,4}];

handles.modelUpdateListener.Enabled = true;

% END FUNCTION tblFixFret_CellEditCallback


% --- Executes when the QubModel object is altered.
function modelUpdate_Callback(tblFixFret,event)
% Update tblFixFret to reflect current model parameters

model = event.AffectedObject;

celldata = num2cell(false(model.nClasses,4));
celldata(:,1) = num2cell(model.mu);
celldata(:,2) = num2cell(model.fixMu);
celldata(:,3) = num2cell(model.sigma);
celldata(:,4) = num2cell(model.fixSigma);
set( tblFixFret, 'Data', celldata );

% END FUNCTION modelUpdate_Callback


% --- Executes on button press in btnSaveModel.
function btnSaveModel_Callback(~, ~, handles) %#ok<DEFNU>
% Save current model to file

if ~isfield(handles,'model') || isempty(handles.model),
    return;
end

fname = handles.model.filename;
[f,p] = uiputfile(fname,'Save model to file');

if f~=0,
    fname = fullfile(p,f);
    handles.model.save( fname );
    set( handles.edModelFilename, 'String',['...' fname(max(1,end-45):end)] );
end


% --- Executes on button press in btnMakeplots.
function btnMakeplots_Callback(~, ~, handles) %#ok<DEFNU>
% Display ensemble plots with the currently-loaded data.
if ~isempty(handles.dataFilenames),
    makeplots( handles.dataFilenames );
end


% --- Executes on button press in btnLifetimeExp.
function btnLifetimeExp_Callback(~, ~, handles) %#ok<DEFNU>
% Compare state lifetimes for all loaded dwell-time files.
% FIXME: should look for both .dwt and .qub.dwt
if ~isempty(handles.dwtFilenames),
    lifetime_exp( handles.dwtFilenames );
end
