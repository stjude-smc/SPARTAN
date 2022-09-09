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
%      unrecognized property name or invalid value makesroperty application
%      stop.  All inputs are passed to batchKinetics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%   Copyright 2007-2017 Cornell University All Rights Reserved.

% Last Modified by GUIDE v2.5 07-Sep-2022 17:59:56


%% ----------------------  GUIDE INITIALIZATION  ---------------------- %%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
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

updateSpartan; %check for updates

% Choose default command line output for batchKinetics
handles.output = hObject;

% Set initial internal state of the program
[handles.model] = deal([]);
[handles.dataFilenames,handles.dwtFilenames] = deal({});

% Set default analysis settings. FIXME: put these in cascadeConstants?
% FIXME: these really should be setup per method...
options.updateModel = true;  %update model viewer during optimization
options.seperately = true; %SKM: analyze each trace individually
options.maxItr = 100;
options.minStates = 1;
options.maxStates = 5;
options.maxRestarts = 10;
options.threshold = 1e-5;
% HMJP parameters
options.alpha    = 2;   % Uniformization; determines further refinements of jump times within a frame period
options.beta     = 10;      % higher number = slow rates. 1/(beta*eta)=peak escape rate in prior
options.eta      = 2;      % gamma distribution shape parameter: 4=peaked prior, 2=exp prior.
options.HMC_eps  = 0.01;   % Hamiltonian Monte Carlo integration step size
options.HMC_L    = 50;     % Hamiltonian Monte Carlo number of Leap-frog integration steps.
handles.options = options;

% Update GUI to reflect these default settings. MIL not supported on Macs
methods = {'Segmental k-Means','Baum-Welch','ebFRET','MIL (Together)','MIL (Separately)','MPL','HMJP'};
set( handles.cboIdealizationMethod, 'String',methods, 'Value',1 );  %SKM
handles = cboIdealizationMethod_Callback(handles.cboIdealizationMethod,[],handles);

constants = cascadeConstants;
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );

handles.traceViewer = TraceListViewer(handles.panTraces);
% hold(handles.axTraces,'on');
% box(handles.axTraces,'on');
guidata(hObject,handles);

if nargin>=4
    btnLoadData_Callback( hObject, [], handles, varargin{:} );
end

% END FUNCTION batchKinetics_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = batchKinetics_OutputFcn(~, ~, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;
% END FUNCTION batchKinetics_OutputFcn






%% ----------------------  LOAD/SAVE DATA CALLBACKS  ---------------------- %%

function btnLoadData_Callback(hObject, ~, handles, files)
% Executes on button press in btnLoadData.

if nargin<4
    % Prompt use for location to save file in...
    filter = {'*.traces','Binary Traces Files (*.traces)'; ...
              '*.*','All Files (*.*)'};
    newFiles = getFiles(filter,'Select traces files to analyze',false);
    if isempty(newFiles), return; end  %user hit cancel.

    handles.dataFilenames = [newFiles handles.dataFilenames];
    handles.dwtFilenames = [findDwt(newFiles) handles.dwtFilenames];
else
    if ischar(files), files={files}; end
    handles.dataFilenames = files;
    handles.dwtFilenames = findDwt(files);
end

[~,names] = cellfun(@fileparts, handles.dataFilenames, 'UniformOutput',false);
set(handles.lbFiles, 'String',names, 'Value',1 );
guidata(hObject,handles);

% Update GUI, showing the first file
lbFiles_Callback(handles.lbFiles, [], handles);
enableControls(handles);

% END FUNCTION btnLoadData_Callback



function enableControls(handles)
% Enable or disable toolbar buttons and menus according to current state.

hasData = ~isempty(handles.dataFilenames);
set( [handles.btnMakeplots handles.mnuViewMakeplots handles.btnSorttraces ...
      handles.mnuSorttraces handles.mnuSimMovie handles.btnFrethist ...
      handles.mnuFrethist], 'Enable',onoff(hasData) );
set( allchild(handles.mnuFileList), 'Enable',onoff(hasData) );

hasModel = ~isempty(handles.model);
set( [handles.btnSaveModel handles.tblFixFret handles.btnSim handles.mnuSim ...
      handles.mnuSimPhoton handles.btnSaveModel], 'Enable',onoff(hasModel) );
set( [handles.btnExecute handles.mnuExecute], 'Enable',onoff(hasData&hasModel) );
  
isIdealized = any( ~cellfun(@isempty,handles.dwtFilenames) );
set( [handles.btnDwellhist handles.btnDwellhist2 handles.mnuDwellhist ...
      handles.mnuDwellhist2 handles.btnPT handles.mnuViewPercentTime ...
      handles.mnuViewTPS handles.btnViewTPS handles.btnOccTime ...
      handles.mnuViewOccTime], 'Enable',onoff(isIdealized));

isMIL = strcmpi(handles.options.idealizeMethod(1:3),'MIL');
set( [handles.mnuExecuteAll handles.btnExecuteAll], 'Enable',...
                                         onoff(hasData & hasModel & ~isMIL) );
set( [handles.chkUpdateModel handles.tblFixFret], 'Enable',...
                                         onoff(hasModel & ~isMIL) );
% END FUNCTION enableControls



function mnuLoadIdl_Callback(hObject, ~, handles) %#ok<DEFNU>
% Load an alternate idealization from file (with conversion from vbFRET, etc.)

idxfile = get(handles.lbFiles,'Value');

dwtfname = getFile( {'*.dwt','QuB Idealization (*.dwt)'; ...
                  '*.mat','vbFRET Idealization (*.mat)'}, 'Load Idealization' );

% Save changes and update GUI.
if ~isempty(dwtfname)
    handles.dwtFilenames{idxfile} = dwtfname;
    guidata(hObject, handles);
    lbFiles_Callback(handles.lbFiles, [], handles);
end

% END FUNCTION mnuLoadIdl_Callback



function mnuSaveIdl_Callback(~, ~, handles) %#ok<DEFNU>
% Save currently-loaded idealization to an alternate file.

idxfile = get(handles.lbFiles,'Value');

[f,p] = uiputfile( {'*.dwt','QuB Idealization (*.dwt)'}, 'Save idealization', ...
                   handles.dwtFilenames{idxfile} );
if isequal(f,0), return; end  %user hit cancel
dwtfname = fullfile(p,f);

% Copy the current idealization file to the new location.
% (idealizations are never stored only in memory in batchKinetics).
if ~strcmp( handles.dwtFilenames{idxfile}, dwtfname )
    copyfile( handles.dwtFilenames{idxfile}, dwtfname );
    fprintf('Saved idealization to %s\n',dwtfname);
end

% END FUNCTION mnuSaveIdl_Callback






%% ------------------  MODEL LOAD/SAVE/EDIT CALLBACKS  ------------------ %%

function btnLoadModel_Callback(hObject, ~, handles, filename)
% Executes on button press in btnLoadModel.
% FIXME: allow QubModelViewer to be updated, rather than recreated.

if nargin<4
    % Ask the user for a filename
    filter = {'*.model','SPARTAN model files (*.model)'; ...
              '*.qmf','QuB format model files (*.qmf)'; ...
              '*.*','All Files (*.*)'};
    [fname,p] = uigetfile(filter, 'Load Model');
    if isequal(fname,0), return; end
    filename = fullfile(p,fname);
end

% Load the model and show the model properties in the GUI.
% The model's properties are automatically updated whenever the model is
% modified in the GUI.
handles.model = QubModel(filename);
handles.modelViewer = QubModelViewer(handles.model, handles.axModel);
handles.traceViewer.model = handles.model;

% Save in most recent list
addRecent(handles.mnuRecentModels, filename);

% Automatically update the parameter table when the model is altered.
handles.modelUpdateListener = addlistener(handles.model,'UpdateModel', ...
                        @(s,e)modelUpdate_Callback(handles.tblFixFret,e) );
handles.model.mu = handles.model.mu;  %trigger listener updates

enableControls(handles);
guidata(hObject, handles);

% END FUNCTION btnLoadModel_Callback



function mnuRecentModels_Callback(hObject,~)
btnLoadModel_Callback( hObject, [], guidata(hObject), get(hObject,'UserData') );
% END FUNCTION mnuRecentModels_Callback



function addRecent(hMenu, filename)
% Add a newly loaded/saved model file to the "Recent" menu list.
recent = get( findobj('Parent',hMenu), 'UserData' );
if ~iscell(recent), recent={recent}; end

if ~any(  cellfun( @(x)strcmp(x,filename), recent)  )
    [~,f,e] = fileparts(filename);
    uimenu(hMenu, 'Label',[f e], 'UserData',filename, ...
               'Callback',@mnuRecentModels_Callback);
end
set(hMenu, 'Enable','on');
%end function addRecent



function btnNewModel_Callback(hObject, ~, handles) %#ok<DEFNU>
% Create a new model with two states and display it. See btnLoadModel_Callback.
% FIXME: allow QubModelViewer to be updated, rather than recreated.

handles.model = QubModel(2);
handles.modelViewer = QubModelViewer(handles.model, handles.axModel);
handles.traceViewer.model = handles.model;

% Automatically update the parameter table when the model is altered.
handles.modelUpdateListener = addlistener(handles.model,'UpdateModel', ...
                        @(s,e)modelUpdate_Callback(handles.tblFixFret,e) );
handles.model.mu = handles.model.mu;  %trigger listener updates

enableControls(handles);
guidata(hObject, handles);

% END FUNCTION btnNewModel_Callback



function btnSaveModel_Callback(~, ~, handles) %#ok<DEFNU>
% Save current model to file
if isfield(handles,'model') && ~isempty(handles.model),
    handles.modelViewer.save_callback();
    addRecent(handles.mnuRecentModels, handles.model.filename);
end
% END FUNCTION btnSaveModel_Callback



function tblFixFret_CellEditCallback(hObject, ~, handles) %#ok<DEFNU>
% tblFixFret was altered. Update current QubModel to match.
enableListener(handles.modelUpdateListener, false);

data = get(hObject,'Data');
handles.model.mu       = [data{:,1}];
handles.model.fixMu    = [data{:,2}];
handles.model.sigma    = [data{:,3}];
handles.model.fixSigma = [data{:,4}];

enableListener(handles.modelUpdateListener, true);
handles.traceViewer.showModelLines();
guidata(hObject,handles);

% END FUNCTION tblFixFret_CellEditCallback



function modelUpdate_Callback(tblFixFret,event)
% Called whenever current QubModel is altered. Updates tblFixFret to match.
model = event.Source;

celldata = num2cell(false(model.nClasses,4));
celldata(:,1) = num2cell(model.mu);
celldata(:,2) = num2cell(model.fixMu);
celldata(:,3) = num2cell(model.sigma);
celldata(:,4) = num2cell(model.fixSigma);
set( tblFixFret, 'Data', celldata );

handles = guidata(tblFixFret);
handles.traceViewer.showModelLines();

% END FUNCTION modelUpdate_Callback






%% -------------------------  EXECUTE ANALYSIS  ------------------------- %%

function handles = btnExecute_Callback(hObject, ~, handles)
% Run the data analysis pipeline with user-specified data & model.

% Verify data and model have been specified by user in GUI.
idxfile  = get(handles.lbFiles,'Value');

% Get options from traceViewer
options = handles.options;
options.dataField = handles.traceViewer.dataField;
options.updateModel = get(handles.chkUpdateModel,'Value');
dwtfname = handles.dwtFilenames{idxfile};

if ~handles.traceViewer.data.isChannel( options.dataField )
    options = settingdlg( options, {'dataField'}, {'Data field to analyze'}, ...
                          {handles.traceViewer.data.channelNames} );
    if isempty(options), return; end  %user hit cancel
end

set(handles.figure1,'pointer','watch');
set(handles.txtStatus,'String','Analyzing...'); drawnow;

try
    % Run MIL rate optimizer using only dwell-times as input
    if strcmpi(options.idealizeMethod(1:3),'MIL')

        % Load dwell-time information
        assert( ~isempty(handles.traceViewer.idl), 'Traces must be idealized before running MIL');
        dwt = idlToDwt( handles.traceViewer.idl(~handles.traceViewer.exclude,:) );
        dt = handles.traceViewer.data.sampling/1000;

        % Run MIL, only updating model rates.
        % NOTE: optModel will have the .qubTree model values, which only reflect 
        % the model as originally loaded from file. FIXME.
        if strcmpi( options.idealizeMethod, 'MIL (Together)' )
            options.updateModel = true;
            optModel = milOptimize(dwt, dt, handles.model, options);
            handles.model.rates = optModel.rates;
        else
            rates = milOptimizeSeparately(dwt, dt, handles.model);
            ratehist(rates);
        end

    % Run any other method that uses FRET data as input.
    else
        handles.traceViewer.loadIdealization();  %clear idealization
        
        % If no .dwt file name was previously specified, construct one now.
        if isempty( dwtfname )
            [p,f] = fileparts( handles.dataFilenames{idxfile} );
            dwtfname = fullfile( p, [f '.qub.dwt'] );
            handles.dwtFilenames{idxfile} = dwtfname;
        end

        % Run optimization algorithm.
        options.exclude = handles.traceViewer.exclude;
        [idl,optModel] = runParamOptimizer( handles.traceViewer.data, ...
                                            dwtfname, handles.model, options );
        
        % Update model and idealization display
        handles.traceViewer.loadIdealization( idl, optModel.mu );
        if get(handles.chkUpdateModel,'Value'),
            handles.model.copyValuesFrom( optModel );
        end
    end
    
    set(handles.txtStatus,'String','Finished');
    guidata(hObject,handles);
    enableControls(handles);

catch e
    if strcmpi(e.identifier,'spartan:op_cancelled')
        set(handles.txtStatus,'String','Operation cancelled');
    else
        errordlg(['Error: ' e.message]);
        set(handles.txtStatus,'String','Operation failed');
    end
    return;
end

set(handles.figure1,'pointer','arrow');
% END FUNCTION btnExecute_Callback



function btnExecuteAll_Callback(~, ~, handles) %#ok<DEFNU>
% Analyize each loaded file in sequence (batch mode).
for i=1:numel(handles.dataFilenames),
    set(handles.lbFiles,'Value',i);
    lbFiles_Callback(handles.lbFiles, [], handles);
    handles = btnExecute_Callback(handles.btnExecute, [], handles);
end
% END FUNCTION btnExecuteAll_Callback



function btnStop_Callback(~, ~, handles) %#ok<DEFNU>
% Executes on button press in btnStop.
% FIXME: this should stop any task mid-execution.
set(handles.btnExecute,'Enable','on');
% END FUNCTION btnStop_Callback




%% -------------------------  SIMULATE DATA  ------------------------- %%

function mnuSim_Callback(hObject, ~, handles) %#ok<DEFNU>
% Simulate traces using current model.

if isempty(handles.model), return; end  %model required.

% Get simulation settings.
persistent opt;
if isempty(opt)
    opt = struct('nTraces',1000, 'nFrames',2000, 'sampling',40, ...
                 'snr',30, 'shotNoise',true, 'gamma',1, ...
                 'totalIntensity',500, 'stdTotalIntensity',0, ...
                 'accLife',20, 'donLife',30, 'simExpBleach',true );
end
prompt = {'Traces',   'Frames',     'Integration Time (ms)', ... 
          'Signal:background noise ratio', 'Shot noise', 'Apparent gamma', ...
          'Intensity (photons)', 'Intensity stdev', ...
          'Acceptor Lifetime (s)', 'Donor Lifetime (s)', ...
          'Exponential photobleaching'};
newOpt = settingdlg(opt, fieldnames(opt), prompt); %, 'Simulation parameters');
if isempty(newOpt), return; end
opt = newOpt;

% Get output filename from user.
[f,p] = uiputfile('sim.traces','Save simulated data as...');
if f==0, return; end  %user pressed "cancel"


% Simulate new data.
set(handles.figure1,'pointer','watch');
set(handles.txtStatus,'String','Simulating...'); drawnow;

try
    newOpt = rmfield(newOpt, {'nTraces','nFrames','sampling'});
    data = simulate( opt.nTraces,opt.nFrames, opt.sampling/1000, handles.model, newOpt );
    saveTraces( fullfile(p,f), data );
catch e
    if ~strcmpi(e.identifier,'spartan:op_cancelled')
        errordlg( ['Error: ' e.message], mfilename );
    end
    set(handles.txtStatus,'String',['Error: ' e.message]);
    set(handles.figure1,'pointer','arrow');
    return;
end

% Load the new simulated file and clear any others loaded.
handles.dataFilenames = { fullfile(p,f) };
handles.dwtFilenames = cell( size(handles.dataFilenames) );  %findDwt(handles.dataFilenames);  %FIXME?

[~,names] = cellfun(@fileparts, handles.dataFilenames, 'UniformOutput',false);
set(handles.lbFiles, 'Value',1, 'String',names);

% Update GUI
guidata(hObject,handles);
enableControls(handles);
lbFiles_Callback(handles.lbFiles, [], handles);

set(handles.figure1,'pointer','arrow');
set(handles.txtStatus,'String','Finished.'); drawnow;

% END FUNCTION mnuSim_Callback





% ~
function mnuSimPhoton_Callback(hObject, ~, handles) %#ok<DEFNU>
% Simulate fluorescence traces one photon at a time using a full photophysical
% model (Jablonski diagram)


if isempty(handles.model), return; end  %model required.

% Get simulation settings.
persistent opt;
if isempty(opt)
    opt = struct('nTraces',1000, 'nFrames',2000, 'sampling',40, 'snr',30, ...
                 'detection',22);
end
prompt = {'Traces', 'Frames', 'Sampling (ms)', 'Signal:background noise ratio', ...
          'Photon detection efficiency (%)'};
newOpt = settingdlg(opt, fieldnames(opt), prompt); %, 'Simulation parameters');
if isempty(newOpt), return; end  %user hit cancel

% Get output filename from user.
[f,p] = uiputfile('sim.traces','Save simulated data as...');
if f==0, return; end  %user pressed "cancel"


% Simulate new data.
% FIXME: simulate.m should return a valid traces object.
set(handles.figure1,'pointer','watch');
set(handles.txtStatus,'String','Simulating...'); drawnow;

opt = newOpt;
% newOpt.stdBackground = opt.totalIntensity/(sqrt(2)*opt.snr);

try
    data = simphotons( [opt.nTraces,opt.nFrames], opt.sampling/1000, handles.model, newOpt );
    saveTraces( fullfile(p,f), data );
catch e
    if ~strcmpi(e.identifier,'spartan:op_cancelled')
        errordlg( ['Error: ' e.message], mfilename );
    end
    set(handles.txtStatus,'String',['Error: ' e.message]);
    set(handles.figure1,'pointer','arrow');
    return;
end

% Load the new simulated file and clear any others loaded.
handles.dataFilenames = { fullfile(p,f) };
handles.dwtFilenames = cell( size(handles.dataFilenames) );  %findDwt(handles.dataFilenames);  %FIXME?

[~,names] = cellfun(@fileparts, handles.dataFilenames, 'UniformOutput',false);
set(handles.lbFiles, 'Value',1, 'String',names);

% Update GUI
guidata(hObject,handles);
enableControls(handles);
lbFiles_Callback(handles.lbFiles, [], handles);

set(handles.figure1,'pointer','arrow');
set(handles.txtStatus,'String','Finished.'); drawnow;

% END FUNCTION mnuSimPhoton_Callback




function mnuSimMovie_Callback(~, ~, handles)  %#ok<DEFNU>
% Executed when the user selects the "Actions->Simulate Movie" menu.
% Simulates a wide-field fluorescence movie by distributing the fluorescence
% from each trace in the currently loaded file across pixels in simulated
% point-spread functions.
% Here we assume the user wants to use all of the particles!

% Get simulation settings.
% FIXME: get these from cascadeConstants.
persistent opt;
if isempty(opt)
    opt = struct('sigmaPSF',0.8, 'density',handles.traceViewer.data.nTraces, ...
                 'aduPhoton',2, 'grid',false, 'alignX',0, 'alignY',0, ...
                 'alignTheta',0, 'alignScale',1);
end
prompt = {'PSF Size (px stdev):', 'Particles to simulate:', ...
          'ADU to photon conversion', 'Use a regular grid', 'X deviation (px):', ...
          'Y deviation (px):', 'Rotation (degrees):', 'Scaling factor'};
newOpt = settingdlg(opt, fieldnames(opt), prompt); %, 'Simulation parameters');
if isempty(newOpt), return; end
opt = newOpt;

newOpt.density = min(handles.traceViewer.data.nTraces, newOpt.density);

% Simulate the movie (simulateMovie.m will ask for the movie filenames)
set(handles.figure1,'pointer','watch');
set(handles.txtStatus,'String','Simulating...'); drawnow;

try
    simulateMovie(handles.traceViewer.data, [],[], newOpt);
catch e
    if ~strcmpi(e.identifier,'spartan:op_cancelled')
        errordlg(['Error: ' e.message]);
    end
    set(handles.txtStatus,'String',['Error: ' e.message]);
end
set(handles.txtStatus,'String','Finished.'); drawnow;
set(handles.figure1,'pointer','arrow');

% END FUNCTION mnuSimMovie_Callback






%% ---------------------   FILE LIST CALLBACKS   --------------------- %%

function handles = lbFiles_Callback(hObject, ~, handles)
% Callback function for changes to lbFiles data file list.
% Draws traces from current file in the trace viewer panel.

if isempty(handles.dataFilenames),
    % No files loaded. Clear trace viewer.
    handles.traceViewer.loadTraceData();
else
    idxFile = get(hObject,'Value');
    handles.traceViewer.loadTraceData( handles.dataFilenames{idxFile} );
    handles.traceViewer.loadIdealization( handles.dwtFilenames{idxFile} );
end

% END FUNCTION lbFiles_Callback



function mnuFileRemove_Callback(hObject, ~, handles) %#ok<DEFNU>
% Close currently-selected file in GUI list.

idxfile = get(handles.lbFiles,'Value');
handles.dataFilenames(idxfile) = [];
handles.dwtFilenames(idxfile) = [];
guidata(hObject,handles);

names = get(handles.lbFiles,'String');
names(idxfile) = [];
set(handles.lbFiles,'String',names, 'Value',max(1,idxfile-1));

enableControls(handles);
lbFiles_Callback(handles.lbFiles, [], handles);

% END FUNCTION mnuFileRemove_Callback



function mnuFileRemoveAll_Callback(hObject, ~, handles) %#ok<DEFNU>
% Close all open files.

handles.dataFilenames = {};
handles.dwtFilenames = {};
set(handles.lbFiles, 'String',{}, 'Value',1);
guidata(hObject,handles);

enableControls(handles);
lbFiles_Callback(handles.lbFiles, [], handles);

% END FUNCTION mnuFileRemoveAll_Callback



function mnuFileUp_Callback(hObject, ~, handles, inc) %#ok<DEFNU>
% Move currently-selected file up in the GUI list.
% The last parameter specifies the direction to move (+1 up, -1 down).

names   = get(handles.lbFiles,'String');
idxfile = get(handles.lbFiles,'Value');  %currently selected file.
idxnew  = max(1,  min( numel(names), idxfile+inc)  );

set( handles.lbFiles, 'String',shiftvec(names, idxfile, idxnew), 'Value',idxnew );
handles.dataFilenames = shiftvec( handles.dataFilenames, idxfile, idxnew );
handles.dwtFilenames  = shiftvec( handles.dwtFilenames,  idxfile, idxnew );

guidata(hObject,handles);

% END FUNCTION mnuFileUp_Callback


function vector = shiftvec( vector, idx, idxfinal )
% Move the element in the index IDX of VECTOR to the final index IDXFINAL.
vector = to_col(vector);
temp = vector(idx);
vector(idx) = [];
vector = [ vector(1:idxfinal-1); temp; vector(idxfinal:end) ];
% END shiftvec





%% ------------------------  PLOTTING FUNCTIONS  ------------------------ %%
% Executed when plotting menu or toolbar buttons are clicked.

function mnuSorttraces_Callback(~, ~, handles) %#ok<DEFNU>
% FIXME: create class method instead.
idxFile  = get(handles.lbFiles,   'Value');
idxTrace = handles.traceViewer.idxShown();
sorttraces( 0, handles.dataFilenames{idxFile}, idxTrace(1) );
% END FUNCTION



function mnuDwellhist_Callback(~, ~, handles, showFits) %#ok<DEFNU>
% Draw dwell-time distributions, with model fits.
if nargin<4, showFits=false; end

if ~isempty(handles.model) && showFits
    params.model = handles.model;
    idxFile = get(handles.lbFiles, 'Value');
    dwellhist(handles.dwtFilenames(idxFile), params);
else
    params.model = [];
    dwellhist(handles.dwtFilenames, params);
end
% END FUNCTION







%% ------------------------  SETTINGS DIALOGS  ------------------------ %%

function handles = cboIdealizationMethod_Callback(hObject, ~, handles)
% Update method to use for idealization
% FIXME: consider getting this value only when needed (execution).

text = get(hObject,'String');
handles.options.idealizeMethod = text{get(hObject,'Value')};
guidata(hObject, handles);

enableControls(handles);

% END FUNCTION cboIdealizationMethod_Callback



function mnuIdlSettings_Callback(hObject, ~, handles) %#ok<DEFNU>
% Change idealization settings

switch upper(handles.options.idealizeMethod(1:3))  %#ok<*MULCC>
    case {'SEG','BAU'}  %SKM, Baum-Welch
        prompt = {'Max iterations','Analyze traces individually:'}; %'LL Convergence:', 'Grad. Convergence:'
        fields = {'maxItr','seperately'};  %'gradLL', 'gradConv'
        
    case {'VBF','EBF'}  %vb/ebFRET
        prompt = {'Min states','Max states','Max restarts:','Convergence'};
        fields = {'minStates', 'maxStates', 'maxRestarts',  'threshold'};
    
    case 'HMJ'
        prompt = {'Alpha (uniformization)','Beta (rate prior distribution scale parameter)',...
                  'Eta (rate prior distribution shape parameter)','HMC integration step size',...
                  'HMC leap-frog steps','Max iterations'};
        fields = {'alpha','beta','eta','HMC_eps','HMC_L','maxItr'};
        
%     case {'MI','MP'}  %MIL or MPL
%         prompt = {'LL threshold','Gradient threshold'};
%         fields = {'convLL',      'convGrad'};

    otherwise
        return;
end

options = settingdlg(handles.options, fields, prompt);
if ~isempty(options),
    handles.options = options;
    guidata(hObject,handles);
end

% END FUNCTION mnuIdlSettings_Callback



% --------------------------------------------------------------------
function mnuSelRates_Callback(~, ~, handles) %#ok<DEFNU>
% Select traces using ranges of rate constants.
% FIXME: assumes number of fitted rates == number of traces.
handles.traceViewer.mnuSelRates_Callback();
% END FUNCTION mnuSelRates_Callback


% --------------------------------------------------------------------
function mnuSelDwells_Callback(~, ~, handles) %#ok<DEFNU>
% Select traces by total number number of dwells in any state
handles.traceViewer.mnuSelByDwells_Callback();
% END FUNCTION mnuSelDwells_Callback


% --------------------------------------------------------------------
function mnuSelOccupancy_Callback(~, ~, handles) %#ok<DEFNU>
% Select traces by state occupancy
handles.traceViewer.mnuSelOccupancy_Callback();
% END FUNCTION mnuSelOccupancy_Callback


% --------------------------------------------------------------------
function btnFrethist_ClickedCallback(~, ~, handles) %#ok<DEFNU>
if ~isempty( handles.traceViewer.data )
    frethistComparison( handles.dataFilenames );
end
% END FUNCTION
