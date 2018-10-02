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

% Last Modified by GUIDE v2.5 09-Mar-2018 16:58:35


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
options.updateModel = true;  %update model viewer during optimization
options.seperately = true; %SKM: analyze each trace individually
options.maxItr = 100;
options.minStates = 1;
options.maxStates = 5;
options.maxRestarts = 10;
options.threshold = 1e-5;
options.dataField = 'fret';  %which 
handles.options = options;

% Update GUI to reflect these default settings. MIL not supported on Macs
methods = {'Segmental k-Means','Baum-Welch','ebFRET','MIL','MPL'};
set( handles.cboIdealizationMethod, 'String',methods, 'Value',1 );  %SKM
handles = cboIdealizationMethod_Callback(handles.cboIdealizationMethod,[],handles);

constants = cascadeConstants;
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );

handles.traceViewer = TraceListViewer(handles.axTraces, handles.sldTraces, handles.sldTracesX);
hold(handles.axTraces,'on');
box(handles.axTraces,'on');
guidata(hObject,handles);

% END FUNCTION batchKinetics_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = batchKinetics_OutputFcn(~, ~, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;
% END FUNCTION batchKinetics_OutputFcn






%% ----------------------  LOAD/SAVE DATA CALLBACKS  ---------------------- %%

function btnLoadData_Callback(hObject, ~, handles) %#ok<DEFNU>
% Executes on button press in btnLoadData.

% Prompt use for location to save file in...
handles.dataFilenames = getFiles([],'Select traces files to analyze');
if isempty(handles.dataFilenames), return; end  %user hit cancel.

[~,names] = cellfun(@fileparts, handles.dataFilenames, 'UniformOutput',false);
set(handles.lbFiles, 'Value',1, 'String',names);

% Look for .dwt files if data were already analyzed.
handles.dwtFilenames = findDwt(handles.dataFilenames);
guidata(hObject,handles);

% Update GUI, showing the first file
lbFiles_Callback(handles.lbFiles, [], handles);
enableControls(handles);

% END FUNCTION btnLoadData_Callback



function enableControls(handles)
% Enable or disable toolbar buttons and menus according to current state.

hasData = ~isempty(handles.dataFilenames);
set( [handles.btnMakeplots handles.mnuViewMakeplots handles.btnSorttraces ...
      handles.mnuSorttraces handles.mnuSimMovie], 'Enable',onoff(hasData) );
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
    filter = {'*.model;*.qmf','All model files (*.model;*.qmf)'; ...
              '*.model','SPARTAN model files (*.model)'; ...
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
trcfile  = handles.dataFilenames{idxfile};

% Get options from traceViewer
options = handles.options;
options.dataField = handles.traceViewer.dataField;
options.exclude = handles.traceViewer.exclude;
options.updateModel = get(handles.chkUpdateModel,'Value');

% Verify external modules installed
if strcmpi(options.idealizeMethod,'ebFRET') && isempty(which('ebfret.analysis.hmm.vbayes'))
    errordlg('ebFRET not found. Check your path.',mfilename);
    disp('Go to https://ebfret.github.io/ to download ebFRET, then add to the MATLAB path.');
    return;
end

% Run the analysis algorithms...
set(handles.figure1,'pointer','watch');
set(handles.txtStatus,'String','Analyzing...'); drawnow;

if strcmpi(options.idealizeMethod(1:3),'MIL')
    options.updateModel = true;
    dwtfname = handles.dwtFilenames{idxfile};
    
    if isempty(dwtfname) || ~exist(dwtfname,'file'),
        errordlg('Traces must be idealized before running MIL');
        set(handles.figure1,'pointer','arrow');
        return;
    end
    [dwt,sampling] = loadDWT(dwtfname);
    
    % Run MIL, only updating model rates.
    % NOTE: optModel will have the .qubTree model values, which only reflect 
    % the model as originally loaded from file. FIXME.
    try
        optModel = milOptimize(dwt, sampling/1000, handles.model, options);
    catch e
        if ~strcmpi(e.identifier,'spartan:op_cancelled')
            errordlg(['Error: ' e.message]);
        end
        set(handles.txtStatus,'String',['Error: ' e.message]);
        set(handles.figure1,'pointer','arrow');
        return;
    end
    handles.model.rates = optModel.rates;
    handles.modelViewer.redraw();
else
    % Clear current idealization (FIXME: also delete .dwt?)
    handles.traceViewer.idl = [];
    handles.traceViewer.showTraces();
    
    try
        [dwtfname,optModel] = runParamOptimizer(handles.model, trcfile, options);
    catch e
        if ~strcmpi(e.identifier,'spartan:op_cancelled')
            errordlg(['Error: ' e.message]);
        end
        set(handles.txtStatus,'String',['Error: ' e.message]);
        set(handles.figure1,'pointer','arrow');
        return;
    end
    
    handles.dwtFilenames{idxfile} = dwtfname;
    handles.traceViewer.idl = loadIdl(dwtfname, handles.traceViewer.data);
    handles.traceViewer.showTraces();

    if get(handles.chkUpdateModel,'Value'),
        handles.model.muteListeners = true;
        handles.model.rates = optModel.rates;
        handles.model.mu    = optModel.mu;
        handles.model.sigma = optModel.sigma;
        handles.model.p0    = optModel.p0;
        handles.model.muteListeners = false;
        handles.modelViewer.redraw();
        
        handles.model.mu = optModel.mu; %trigger update
        
        handles.traceViewer.showModelLines();
    end
end

% Display summary plots of the results to the GUI
statehist( handles.axResult1, dwtfname, trcfile );
tplot( handles.axResult2, tdplot(dwtfname,trcfile) );  ylabel( handles.axResult2, '');
dhparam.model = optModel;
dwellhist( handles.axResult3, dwtfname, dhparam );

guidata(hObject,handles);
enableControls(handles);
set(handles.figure1,'pointer','arrow');
set(handles.txtStatus,'String','Finished');

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
    handles.traceViewer.data = [];
    handles.traceViewer.idl  = [];
    handles.traceViewer.dataFilename = '';
else
    set(handles.figure1,'pointer','watch');
    
    idxFile = get(hObject,'Value');
    datafname = handles.dataFilenames{idxFile};
    dwtfname  = handles.dwtFilenames{idxFile};
    data = loadTraces( datafname );
    
    handles.traceViewer.dataFilename = datafname;
    handles.traceViewer.data = data;
    handles.traceViewer.idl = loadIdl(dwtfname, data);
    
    set(handles.figure1,'pointer','arrow');
end

handles.traceViewer.redraw();

% END FUNCTION lbFiles_Callback



function idl = loadIdl(dwtfname, data)
% Returns the idealization for the currently selected file

idl = [];
try
    if ~isempty(dwtfname) && exist(dwtfname,'file'),
        [dwt,~,offsets,model] = loadDWT(dwtfname);
        idl = dwtToIdl(dwt, offsets, data.nFrames, data.nTraces);

        assert( size(model,2)==2 );
        fretValues = [NaN; model(:,1)];
        idl = fretValues( idl+1 );
    end
catch e
    disp(['Idealization failed to load: ' e.message])
end

% END FUNCTION loadIdl



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
idxFile  = get(handles.lbFiles,   'Value');
idxTrace = get(handles.sldTraces,'Max')-floor(get(handles.sldTraces,'Value'));
sorttraces( 0, handles.dataFilenames{idxFile}, idxTrace );
% END FUNCTION



function mnuDwellhist_Callback(~, ~, handles, showFits) %#ok<DEFNU>
% Draw dwell-time distributions, with model fits.
if nargin<4, showFits=false; end

if ~isempty(handles.model) && showFits
    params.model = handles.model;
    dwellhist(handles.dwtFilenames(1), params);
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

% Parameters common to all methods
prompt = {'Max iterations'};
fields = {'maxItr'};

switch upper(handles.options.idealizeMethod(1:2))  %#ok<*MULCC>
    case {'SE','BA'}  %SKM, Baum-Welch
        prompt = [prompt, {'Analyze traces individually:'}]; %'LL Convergence:', 'Grad. Convergence:'
        fields = [fields, {'seperately'}];  %'gradLL', 'gradConv'
        
    case {'VB','EB'}  %vb/ebFRET
        prompt = {'Min states','Max states','Max restarts:','Convergence'};
        fields = {'minStates', 'maxStates', 'maxRestarts',  'threshold'};
    
%     case {'MI','MP'}  %MIL or MPL
%         prompt = {'LL threshold','Gradient threshold'};
%         fields = {'convLL',      'convGrad'};
end

options = settingdlg(handles.options, fields, prompt);
if ~isempty(options),
    handles.options = options;
    guidata(hObject,handles);
end

% END FUNCTION mnuIdlSettings_Callback
