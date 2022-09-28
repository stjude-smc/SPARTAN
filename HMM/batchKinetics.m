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

% Last Modified by GUIDE v2.5 28-Sep-2022 15:47:10


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
[handles.model,handles.rates] = deal([]);
[handles.dataFilenames,handles.dwtFilenames] = deal({});

% Update GUI to reflect these default settings. MIL not supported on Macs
methods = {'Segmental k-Means','Baum-Welch','ebFRET','MIL (Together)','MIL (Separately)','MPL','HMJP'};
set( handles.cboIdealizationMethod, 'String',methods, 'Value',1 );  %SKM
handles = cboIdealizationMethod_Callback(handles.cboIdealizationMethod,[],handles);

constants = cascadeConstants;
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );

handles.traceViewer = TraceListViewer(handles.panTraces);
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
assert( ~isempty(names) );
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

idealizeMethod = getIdealizeMethod(handles);
isMIL = strcmpi(idealizeMethod, 'MIL');
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
tv = handles.traceViewer;
saveDWT( dwtfname, idlToDwt(tv.idl), tv.idlValues(2:end), handles.data.sampling );

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

data = handles.traceViewer.data;
idl = handles.traceViewer.idl;

% Verify data and model have been specified by user in GUI.
idxfile  = get(handles.lbFiles,'Value');

% Get options from traceViewer
idealizeMethod = getIdealizeMethod(handles);
options = hmmopt(idealizeMethod);
options.dataField = handles.traceViewer.dataField;
options.updateModel = get(handles.chkUpdateModel,'Value');

if ~data.isChannel( options.dataField )
    options = settingdlg( options, {'dataField'}, {'Data field to analyze'}, ...
                                                   {data.channelNames} );
    if isempty(options), return; end  %user hit cancel
end

set(handles.figure1,'pointer','watch');
set(handles.txtStatus,'String','Analyzing...'); drawnow;

try
    %handles.traceViewer.loadIdealization();  %clear idealization
    %drawnow;

    % Run optimization algorithm.
    options.idealizeMethod = idealizeMethod;
    options.exclude = handles.traceViewer.exclude;
    options.truncate = handles.traceViewer.truncate;
    [idl,optModel,handles.rates] = runParamOptimizer( data, idl, handles.model, options );
    
    % Save the idealization.
    if ~strcmpi(options.idealizeMethod(1:3),'MIL')
        dwtfname = handles.dwtFilenames{idxfile};
        if isempty( dwtfname )
            [p,f] = fileparts( handles.dataFilenames{idxfile} );
            dwtfname = fullfile( p, [f '.qub.dwt'] );
            handles.dwtFilenames{idxfile} = dwtfname;
        end
        saveDWT( dwtfname, idl, optModel, data.sampling );
        handles.traceViewer.loadIdealization( idl, optModel.mu );
    end
    
    % Update model and idealization display
    if get(handles.chkUpdateModel,'Value'),
        handles.model.copyValuesFrom( optModel );
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
    assert( ~isempty(idxFile) );
    handles.traceViewer.loadTraceData( handles.dataFilenames{idxFile} );
    handles.traceViewer.loadIdealization( handles.dwtFilenames{idxFile} );
end

handles.rates = [];

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


function btnFrethist_ClickedCallback(~, ~, handles) %#ok<DEFNU>
if ~isempty( handles.traceViewer.data )
    frethistComparison( handles.dataFilenames );
end
% END FUNCTION



% --------------------------------------------------------------------
function mnuDwellCorr_Callback(~, ~, handles) %#ok<DEFNU>
% 
idl = handles.traceViewer.idl;
if ~isempty(idl)
    idl = idl( ~handles.traceViewer.exclude, : );
    memtrace_JBM(idl);
end





%% ------------------------  SETTINGS DIALOGS  ------------------------ %%

function result = getIdealizeMethod(handles)
% Get the current idealization method from cboIdealizationMethod
text = get(handles.cboIdealizationMethod,'String');
result = text{get(handles.cboIdealizationMethod,'Value')};


function handles = cboIdealizationMethod_Callback(~, ~, handles)
% Update method to use for idealization
enableControls(handles);
% END FUNCTION cboIdealizationMethod_Callback


function mnuIdlSettings_Callback(~, ~, handles) %#ok<DEFNU>
% Change idealization settings
hmmopt( getIdealizeMethod(handles), true );
% END FUNCTION mnuIdlSettings_Callback



% --------------------------------------------------------------------
function mnuSelRates_Callback(~, ~, handles) %#ok<DEFNU>
% Select traces using ranges of rate constants.
% FIXME: assumes number of fitted rates == number of traces.

    persistent defaults;
    ex = handles.traceViewer.exclude;
    rates = handles.rates;
    
    if isempty(rates)
        errordlg('No rate info available. Please run optimizer.');
        return;
    end
    if numel(ex)~=size(rates,3)
        errordlg('Rate matrix size mismatch. Please run optimizer again.');
        return;
    end

    % Set order of state pairs that describe each rate constant
    [src,dst] = find( nanmin(rates,[],3)>0 );  %& ~model.fixRates;
    [src,idx] = sort(src);
    dst = dst(idx);

    prompt = cell( 2*numel(src), 1 );
    for i=1:numel(src)
        j = (i-1)*2 +1;
        prompt{j}   = sprintf('k%d,%d >', src(i), dst(i) );
        prompt{j+1} = sprintf('k%d,%d <', src(i), dst(i) );
    end

    if numel(defaults) ~= numel(prompt)
        defaults = repmat( {''}, [2*numel(src) 1] );
    end
    answer = inputdlg( prompt, 'Select traces by fitted rate constants', ...
                       1, defaults );
    if isempty(answer), return; end

    % Update exclusion list and update display.
    for i=1:numel(src)
        j = (i-1)*2 +1;
        values = squeeze(  rates( src(i), dst(i), : )  );
        lb = str2double( answer{j} );
        ub = str2double( answer{j+1} );

        if ~isnan(lb)
            ex( values <= lb ) = true;
        end
        if ~isnan(ub)
            ex( values >= ub ) = true;
        end
        ex( isnan(values) ) = true;
    end
    
    defaults = answer;
    handles.traceViewer.exclude = ex;
    handles.traceViewer.showTraces();
    fprintf('%d of %d selected (%.0f%%)\n', sum(~ex),numel(ex),100*sum(~ex)/numel(ex) );

% END FUNCTION mnuSelRates_Callback


% --------------------------------------------------------------------
function mnuSelDwells_Callback(~, ~, handles) %#ok<DEFNU>
% Select traces by total number number of dwells in any state

    persistent defaults;
    if isempty(handles.traceViewer.idl), return; end
    
    % Prompt user for trace selection criteria
    if isempty(defaults), defaults={'',''};  end
    answer = inputdlg( {'Minimum:','Maximum:'}, 'Select traces by number of dwells', ...
                       1, defaults );
    if isempty(answer), return; end
    bounds = cellfun( @str2double, answer );
    if any( isnan(bounds) & ~cellfun(@isempty,answer) )
        errordlg('Invalid input value');
        return;
    end
    if isnan(bounds(1)), bounds(1)=-Inf; end
    if isnan(bounds(2)), bounds(2)=Inf; end
    
    % Calculate number of dwells in each trace
    nDwells = cellfun( @numel, idlToDwt(handles.traceViewer.idl,true) );
    
    % TESTING
    %figure;
    %hist( nDwells, 0:10:100 );
    keep = nDwells>=bounds(1) & nDwells<=bounds(2);
    fprintf('Selected %d of %d traces (%.0f%%) \n', ...
        sum(keep), numel(keep), 100*sum(keep)/numel(keep) );
    
    % Update exclusion list and update display.
    defaults = answer;
    handles.traceViewer.exclude( nDwells<bounds(1) | nDwells>bounds(2) ) = true;
    handles.traceViewer.showTraces();

% END FUNCTION mnuSelDwells_Callback


% --------------------------------------------------------------------
function mnuSelOccupancy_Callback(~, ~, handles) %#ok<DEFNU>
% Select traces by state occupancy

    persistent defaults;
    if isempty(handles.traceViewer.idl), return; end
    ex = handles.traceViewer.exclude;

    nClass = numel(handles.traceViewer.idlValues)-1;
    prompt = cell( nClass, 1 );
    for i=1:nClass
        prompt{i} = sprintf('Minimum frames in class %d', i );
    end

    % Prompt user for minumum number of frames in a state
    if numel(defaults) ~= numel(prompt)
        defaults = repmat( {''}, nClass );
    end
    answer = inputdlg( prompt, 'Select traces by state occupancy', 1, defaults );
    if isempty(answer), return; end

    bounds = cellfun( @str2double, answer );
    if any( isnan(bounds) & ~cellfun(@isempty,answer) )
        errordlg('Invalid input value');
        return;
    end
    bounds( isnan(bounds) ) = 0;

    % Update exclusion list and update display.
    for i=1:nClass
        if ~isnan(bounds(i))
            occupancy = sum( handles.traceViewer.idl==i, 2 );
            ex = ex | occupancy<bounds(i);
        end
    end
    defaults = answer;
    handles.traceViewer.exclude = ex;
    handles.traceViewer.showTraces();

% END FUNCTION mnuSelOccupancy_Callback



% --------------------------------------------------------------------
function mnuFretSep_Callback(~, ~, handles) %#ok<DEFNU>
% Select traces with state FRET values close to ensemble average.
    tv = handles.traceViewer;
    if isempty(tv.idl), return; end
    tv.exclude = tv.exclude | ~fret_model_filter( tv.data, tv.idl );
    tv.showTraces();
% END FUNCTION mnuFretSep_Callback


% --------------------------------------------------------------------
function mnuTrunc1_Callback(~, ~, handles) %#ok<DEFNU>
% Truncate traces to donor photobleaching
data = handles.traceViewer.data;
handles.traceViewer.truncate = zeros( 1, data.nTraces );
for i=1:data.nTraces
    x = find( data.fret(i,:)~=0, 1, 'last' );
    if ~isempty(x)
        handles.traceViewer.truncate(i) = x;
    end
end
handles.traceViewer.showTraces();


% --------------------------------------------------------------------
function mnuTrunc2_Callback(~, ~, handles) %#ok<DEFNU>
% Truncate traces to first donor blink
data = handles.traceViewer.data;
handles.traceViewer.truncate = repmat( data.nFrames, [1 data.nTraces] );
for i=1:data.nTraces
    x = find( data.fret(i,:)==0, 1, 'first' );
    if ~isempty(x)
        handles.traceViewer.truncate(i) = x-1;
    end
end
handles.traceViewer.showTraces();


% --------------------------------------------------------------------
function mnuTrunc3_Callback(~, ~, handles) %#ok<DEFNU>
% Truncate traces by FRET value
% TODO: prompt user for a specific FRET value.
data = handles.traceViewer.data;
handles.traceViewer.truncate = repmat( data.nFrames, [1 data.nTraces] );
for i=1:data.nTraces
    x = find( rleFilter(data.fret(i,:)>=0.1,5), 1, 'last' );
    if ~isempty(x)
        handles.traceViewer.truncate(i) = x;
    end
end
handles.traceViewer.showTraces();


% --------------------------------------------------------------------
function mnuTrunc4_Callback(~, ~, handles) %#ok<DEFNU>
% Truncate traces by idealized state

data = handles.traceViewer.data;
idl  = handles.traceViewer.idl;
if isempty(idl)
    errordlg('Traces must be idealized');
    return;
end

% Prompt user for a state number and minimum number of frames.
% persistent options
% if isempty(options)
%     options = struct('state',1, 'minFrames',1);
% end
% nStates = max(idl(:));
% states = strsplit( num2str(1:nStates) );
% settingdlg( options, fieldnames(options), {'State #:',''} , {states,@isfinite} );
% if ~ismember(options.state,unique(idl(:))) || isnan(options.minFrames)
%     errordlg('Invalid parameter values');
%     return;
% end

% Get truncation point for each trace
handles.traceViewer.truncate = repmat( data.nFrames, [1 data.nTraces] );
for i=1:data.nTraces
    x = find( idl(i,:)<=1, 1, 'first' );
    if ~isempty(x)
        handles.traceViewer.truncate(i) = x-1;
    end
end
handles.traceViewer.showTraces();


% --------------------------------------------------------------------
function mnuRemoveDarkState_Callback(~, ~, handles) %#ok<DEFNU>
% Remove the lowest-FRET class from model and idealization.

tv = handles.traceViewer;
if isempty(tv.idl) || isempty(handles.model)
    errordlg('Idealization and model must be loaded first.');
    return;
end

% Truncate idealizations to first dark-state dwell (if any)
for i=1:size(tv.idl,1)
    x = find( tv.idl(i,:)<=1, 1, 'first' );
    if ~isempty(x)
        handles.traceViewer.truncate(i) = x;
        tv.idl(i,x:end) = 0;
    end
end
tv.idl = max(0,tv.idl-1);
tv.idlValues = [NaN; tv.idlValues(3:end)];
tv.showTraces();

handles.model.removeClass(1);
handles.modelViewer.redraw();  %FIXME: need a new notifier for this?





