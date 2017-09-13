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

% Last Modified by GUIDE v2.5 19-Apr-2017 20:16:23


%% GUI Callbacks

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
[handles.modelFilename,handles.model,handles.data,handles.idl] = deal([]);
[handles.dataFilenames,handles.dwtFilenames] = deal({});
handles.nTracesToShow = 6;  %number displayed in trace display panel
handles.showStateMarkers = true;  %show dotted lines for model FRET values

% Set default analysis settings. FIXME: put these in cascadeConstants?
options.bootstrapN = 1;
options.deadTime = 0.5;
options.seperately = 1; %SKM: analyze each trace individually
options.maxItr = 100;
options.minStates = 1;
options.maxStates = 5;
options.maxRestarts = 10;
options.threshold = 1e-5;
options.dataField = 'fret';  %which 
handles.options = options;

% Update GUI to reflect these default settings. MIL not supported on Macs
methods = {'Segmental k-Means','Baum-Welch','ebFRET','MIL (Rate Optimizer)'};
if isempty(which('ebfret.analysis.hmm.vbayes')), methods(3)=[]; end
if ismac, methods(end)=[]; end

set( handles.cboIdealizationMethod, 'String',methods, 'Value',1 );  %SKM
handles = cboIdealizationMethod_Callback(handles.cboIdealizationMethod,[],handles);

constants = cascadeConstants;
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );

% Trace viewer pane callbacks
handles.sldTracesListener(1) = addlistener( handles.sldTraces, 'Value', ...
          'PostSet',@(h,e)showTraces(guidata(e.AffectedObject))  );

handles.sldTracesListener(2) = addlistener( handles.sldTracesX, 'Value', ...
          'PostSet',@(h,e)sldTracesX_Callback(h,e,guidata(e.AffectedObject))  );

hold(handles.axTraces,'on');
box(handles.axTraces,'on');
guidata(hObject,handles);

% END FUNCTION batchKinetics_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = batchKinetics_OutputFcn(~, ~, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;
% END FUNCTION batchKinetics_OutputFcn




%% =======================  CALLBACK FUNCTIONS  ======================= %%

function btnLoadData_Callback(~, ~, handles) %#ok<DEFNU>
% Executes on button press in btnLoadData.

% Prompt use for location to save file in...
handles.dataFilenames = getFiles([],'Select traces files to analyze');
if isempty(handles.dataFilenames), return; end  %user hit cancel.

[~,names] = cellfun(@fileparts, handles.dataFilenames, 'UniformOutput',false);
set(handles.lbFiles, 'Value',1, 'String',names);

% Look for .dwt files if data were already analyzed.
handles.dwtFilenames = findDwt(handles.dataFilenames);

% Update GUI, showing the first file
lbFiles_Callback(handles.lbFiles, [], handles);
enableControls(handles);

% END FUNCTION btnLoadData_Callback


function enableControls(handles)
% Enable or disable toolbar buttons and menus according to current state.

hasData = ~isempty(handles.dataFilenames);
set( [handles.btnMakeplots handles.mnuViewMakeplots handles.btnSorttraces ...
      handles.mnuSorttraces handles.mnuSimMovie], 'Enable',onoff(hasData) );

hasModel = ~isempty(handles.model);
set( [handles.btnSaveModel handles.tblFixFret handles.btnSim handles.mnuSim ...
      handles.btnSaveModel], 'Enable',onoff(hasModel) );
set( [handles.btnExecute handles.mnuExecute], 'Enable',onoff(hasData&hasModel) );
  
isIdealized = any( ~cellfun(@isempty,handles.dwtFilenames) );
set( [handles.btnDwellhist handles.mnuDwellhist handles.btnPT ...
      handles.mnuViewPercentTime handles.mnuViewTPS handles.btnViewTPS...
      handles.btnOccTime handles.mnuViewOccTime], 'Enable',onoff(isIdealized));

isMIL = strcmpi(handles.options.idealizeMethod(1:3),'MIL');
set( [handles.mnuExecuteAll handles.btnExecuteAll], 'Enable',...
                                         onoff(hasData & hasModel & ~isMIL) );
set( [handles.chkUpdateModel handles.tblFixFret], 'Enable',...
                                         onoff(hasModel & ~isMIL) );

% END FUNCTION enableControls


function btnLoadModel_Callback(hObject, ~, handles, filename)
% Executes on button press in btnLoadModel.

if nargin<4,
    if isempty(handles.modelFilename),
        handles.modelFilename = fullfile(pwd,'*.qmf');
    end

    % Ask the user for a filename
    [fname,p] = uigetfile( handles.modelFilename, 'Select a QuB model file...' );
    if fname==0, return; end
    filename = fullfile(p,fname);
end

% Load the model and show the model properties in the GUI.
% The model's properties are automatically updated whenever the model is
% modified in the GUI.
handles.model = QubModel(filename);
handles.modelViewer = QubModelViewer(handles.model, handles.axModel);

% Enable relevant GUI controls
enableControls(handles);

% Save in most recent list
recent = get( findobj('Parent',handles.mnuRecentModels), 'UserData' );
if ~iscell(recent), recent={recent}; end

if ~any(  cellfun( @(x)strcmp(x,filename), recent)  )
    uimenu(handles.mnuRecentModels, 'Label',fname, 'UserData',filename, ...
               'Callback',@mnuRecentModels_Callback);
end
set(handles.mnuRecentModels, 'Enable','on');

% Automatically update the parameter table when the model is altered.
handles.modelUpdateListener = addlistener(handles.model,'UpdateModel', ...
                        @(s,e)modelUpdate_Callback(handles.tblFixFret,e) );
handles.model.mu = handles.model.mu;  %trigger table update

guidata(hObject, handles);

lbFiles_Callback(handles.lbFiles, [], handles);

% END FUNCTION btnLoadModel_Callback



function mnuRecentModels_Callback(hObject,~)
btnLoadModel_Callback( hObject, [], guidata(hObject), get(hObject,'UserData') );
% END FUNCTION mnuRecentModels_Callback




% --------------------------------------------------------------------
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


% --------------------------------------------------------------------
function mnuSaveIdl_Callback(~, ~, handles) %#ok<DEFNU>
% Save currently-loaded idealization to an alternate file.

idxfile = get(handles.lbFiles,'Value');

[f,p] = uiputfile( {'*.dwt','QuB Idealization (*.dwt)'}, 'Save idealization', ...
                   handles.dwtFilenames{idxfile} );
if isequal(f,0), return; end  %user hit cancel
dwtfname = fullfile(p,f);

% Copy the current idealization file to the new location.
% (idealizations are never stored only in memory in batchKinetics).
copyfile( handles.dwtFilenames{idxfile}, dwtfname );
fprintf('Saved idealization to %s\n',dwtfname);

% END FUNCTION mnuSaveIdl_Callback



function handles = btnExecute_Callback(hObject, ~, handles)
% Run the data analysis pipeline with user-specified data & model.

% Verify data and model have been specified by user in GUI.
idxfile = get(handles.lbFiles,'Value');
trcfile  = handles.dataFilenames{idxfile};
dwtfname = handles.dwtFilenames{idxfile};

% Verify external modules installed
if strcmpi(handles.options.idealizeMethod,'ebFRET') && isempty(which('ebfret.analysis.hmm.vbayes'))
    errordlg('ebFRET not found. Check your path.',mfilename);
    disp('Go to https://ebfret.github.io/ to download ebFRET, then add to the MATLAB path.');
    return;
end

% Run the analysis algorithms...
% FIXME: ideally we want idl (or dwt) returned directly for speed.
set(handles.figure1,'pointer','watch');
set(handles.txtStatus,'String','Analyzing...'); drawnow;

if strcmpi(handles.options.idealizeMethod(1:3),'MIL')
    if isempty(dwtfname) || ~exist(dwtfname,'file'),
        errordlg('Traces must be idealized before running MIL');
        set(handles.figure1,'pointer','arrow');
        return;
    end
    
    optModel = milOptimize(dwtfname, handles.model, handles.options);
    handles.model.rates = optModel.rates;
    handles.modelViewer.redraw();
else
    % Clear current idealization (FIXME also delete .dwt?)
    handles.idl = [];
    set(handles.hIdlLine,'Visible','off');
    
    try
        [handles.dwtFilenames{idxfile},optModel] = runParamOptimizer(...
                                 handles.model, trcfile, handles.options);
    catch e
        set(handles.figure1,'pointer','arrow');
        set(handles.txtStatus,'String',['Error: ' e.message]);
        return;
    end
end

if get(handles.chkUpdateModel,'Value'),
    handles.model.rates = optModel.rates;
    handles.model.mu    = optModel.mu;
    handles.model.sigma = optModel.sigma;
    handles.model.p0    = optModel.p0;
    handles.modelViewer.redraw();
end

% Save results to file for later processing by the user.
% save('resultTree.mat','resultTree');
% qub_saveTree(resultTree,resultFilename);
% qub_saveTree(resultTree.milResults(1).ModelFile,'result.qmf','ModelFile');

% Load and draw idealization, show traces, and update toolbar/menu state.
handles.idl = loadIdl(handles);
guidata(hObject,handles);
showTraces(handles);
set(handles.hIdlLine,'Visible','on');

enableControls(handles);
set(handles.figure1,'pointer','arrow');
set(handles.txtStatus,'String','Finished'); %drawnow;

% END FUNCTION btnExecute_Callback


function btnExecuteAll_Callback(~, ~, handles) %#ok<DEFNU>
% Analyize each loaded file in sequence (batch mode).
for i=1:numel(handles.dataFilenames),
    set(handles.lbFiles,'Value',i);
    handles = lbFiles_Callback(handles.lbFiles, [], handles);
    handles = btnExecute_Callback(handles.btnExecute, [], handles);
end
% END FUNCTION btnExecuteAll_Callback


function idl = loadIdl(handles)
% Returns the idealization for the currently selected file

dwtfname = handles.dwtFilenames{ get(handles.lbFiles,'Value') };

if ~isempty(dwtfname) && exist(dwtfname,'file'),
    [dwt,~,offsets,model] = loadDWT(dwtfname);
    idl = dwtToIdl(dwt, offsets, handles.data.nFrames, handles.data.nTraces);
    
    assert( size(model,2)==2 );
    fretValues = [NaN; model(:,1)];
    idl = fretValues( idl+1 );
else
    idl = [];
end

% END FUNCTION loadIdl


function btnStop_Callback(~, ~, handles) %#ok<DEFNU>
% Executes on button press in btnStop.
% FIXME: this should stop any task mid-execution.
set(handles.btnExecute,'Enable','on');
% END FUNCTION btnStop_Callback





%% Other GUI Callbacks
%  ========================================================================

% --- Executes on selection change in cboIdealizationMethod.
function handles = cboIdealizationMethod_Callback(hObject, ~, handles)
% Update method to use for idealization
% FIXME: consider getting this value only when needed (execution).

text = get(hObject,'String');
handles.options.idealizeMethod = text{get(hObject,'Value')};
guidata(hObject, handles);

enableControls(handles);

% END FUNCTION cboIdealizationMethod_Callback


% --- Executes when entered data in editable cell(s) in tblFixFret.
function tblFixFret_CellEditCallback(hObject, ~, handles) %#ok<DEFNU>
% Update QubModel object with new settings from the table.
% FIXME: avoid triggering many calls to modelUpdate_Callback.

enableListener(handles.sldTracesListener, false);

data = get(hObject,'Data');
handles.model.mu       = [data{:,1}];
handles.model.fixMu    = [data{:,2}];
handles.model.sigma    = [data{:,3}];
handles.model.fixSigma = [data{:,4}];

enableListener(handles.sldTracesListener, true);
guidata(hObject,handles);

% END FUNCTION tblFixFret_CellEditCallback


% --- Executes when the QubModel object is altered.
function modelUpdate_Callback(tblFixFret,event)
% Update tblFixFret to reflect current model parameters.
% FIXME: This gets called four times above -- one for each model change.

model = event.Source;

celldata = num2cell(false(model.nClasses,4));
celldata(:,1) = num2cell(model.mu);
celldata(:,2) = num2cell(model.fixMu);
celldata(:,3) = num2cell(model.sigma);
celldata(:,4) = num2cell(model.fixSigma);
set( tblFixFret, 'Data', celldata );
showModelLines(guidata(gcbf));

% END FUNCTION modelUpdate_Callback


% --- Executes on button press in btnSaveModel.
function btnSaveModel_Callback(~, ~, handles) %#ok<DEFNU>
% Save current model to file
if isfield(handles,'model') && ~isempty(handles.model),
    handles.modelViewer.save_callback();
end
% END FUNCTION btnSaveModel_Callback


% --- Executes on button press in btnSaveModel.
function btnNewModel_Callback(hObject, ~, handles) %#ok<DEFNU>
% Create a new model object.
% FIXME: this shares some code with btnLoadModel.

% Create a new model with two states and display it.
handles.model = QubModel(2);
handles.modelViewer = QubModelViewer(handles.model, handles.axModel);

% Enable relevant GUI controls
enableControls(handles);

% Automatically update the parameter table when the model is altered.
handles.modelUpdateListener = addlistener(handles.model,'UpdateModel', ...
                        @(s,e)modelUpdate_Callback(handles.tblFixFret,e) );
handles.model.mu = handles.model.mu;  %trigger table update

handles.modelFilename = [];
guidata(hObject, handles);
showModelLines(handles);

% END FUNCTION btnSaveModel_Callback


% ========================  PLOTTING FUNCTIONS  ======================== %
% Executed when plotting menu or toolbar buttons are clicked.

function mnuSorttraces_Callback(~, ~, handles) %#ok<DEFNU>
idxFile  = get(handles.lbFiles,   'Value');
idxTrace = get(handles.sldTraces,'Max')-floor(get(handles.sldTraces,'Value'));
sorttraces( 0, handles.dataFilenames{idxFile}, idxTrace );



% =========================  SETTINGS DIALOGS  ========================= %

function mnuIdlSettings_Callback(hObject, ~, handles) %#ok<DEFNU>
% Change idealization settings

switch upper(handles.options.idealizeMethod(1:2))  %#ok<*MULCC>
    case {'SE','BA'}  %SKM, Baum-Welch
        prompt = {'Analyze traces individually:', 'Max iterations:'}; %'LL Convergence:', 'Grad. Convergence:'
        fields = {'seperately', 'maxItr'};  %'gradLL', 'gradConv'
    
    case {'VB','EB'}  %vb/ebFRET
        prompt = {'Min states','Max states','Max restarts:','Max iterations:','Convergence'};
        fields = {'minStates', 'maxStates', 'maxRestarts',  'maxItr',         'threshold'};
    
    case 'MI'  %MIL
        prompt = {'Max iterations:'};  %'Dead time (frames):'
        fields = {'maxIter'};  %'deadTime'
        
    otherwise
        return;
end

options = settingdlg(handles.options, fields, prompt);
if ~isempty(options),
    handles.options = options;
    guidata(hObject,handles);
end

% END FUNCTION mnuIdlSettings_Callback



% ========================  TRACE VIEWER PANEL  ======================== %

function mnuDisplaySettings_Callback(~, ~, handles) %#ok<DEFNU>
% Change display settings (e.g., number of traces displayed).

opt = struct('dataField',handles.options.dataField,'nTracesToShow',handles.nTracesToShow, 'showStateMarkers',handles.showStateMarkers);
prompt = {'Data field', 'Number of traces to show', 'Show model FRET values over traces'};
types = { handles.data.channelNames(handles.data.idxFret), @(x)(x==round(x)&x>0) };  %'total' field inaccessable...
opt = settingdlg(opt, fieldnames(opt), prompt, types);

if ~isempty(opt),
    handles.options.dataField = opt.dataField;
    handles.nTracesToShow = opt.nTracesToShow;
    handles.showStateMarkers = opt.showStateMarkers;
    lbFiles_Callback(handles.lbFiles, [], handles);
end

% END FUNCTION mnuDisplaySettings_Callback


function sldTracesX_Callback(~, ~, handles)
% User adjusted the trace view slider -- show a different subset of traces.

xlimit = floor(get(handles.sldTracesX,'Value'));
set( handles.axTraces, 'XLim',[0 handles.data.time(xlimit)/1000] );

for i=1:handles.nTracesToShow,
    p = get(handles.hTraceLabel(i), 'Position');
    p(1) = 0.98*handles.data.time(xlimit)/1000;
    set( handles.hTraceLabel(i), 'Position',p );
end

% END FUNCTION sldTraces_Callback


% --------------------------------------------------------------------
function handles = lbFiles_Callback(hObject, ~, handles)
% Callback function for changes to lbFiles data file list.
% Draws traces from current file in the trace viewer panel.
% Creates line objects for data plotting, which are updated when the user
% scrolls or makes other changes by calling showTraces() below.
% 
% sldTraces: scrolls vertical list of traces. Value is trace shown at top.
%            Max means the Value showing the first traces (starting with 1).
% sldTracesX: Value means the truncation length.

if isempty(handles.dataFilenames), return; end  %no data loaded.

set(handles.figure1,'pointer','watch');

idxFile = get(hObject,'Value');
data = loadTraces( handles.dataFilenames{idxFile} );
handles.data = data;
handles.idl = loadIdl(handles);

% Disable the trace scroll bar when it won't be valid (prevents errors).
set( handles.sldTraces, 'Enable',onoff(data.nTraces>handles.nTracesToShow) );

enableListener(handles.sldTracesListener, false);
set(handles.sldTraces,  'Min',min(data.nTraces,handles.nTracesToShow),  ...
                        'Max',data.nTraces, 'Value',data.nTraces);
set(handles.sldTracesX, 'Min',10, 'Max',data.nFrames, 'Value',data.nFrames);
enableListener(handles.sldTracesListener, true);

% Setup axes for plotting traces.
cla(handles.axTraces);
[handles.hFretLine, handles.hIdlLine, handles.hTraceLabel] = deal([]);

time = handles.data.time/1000;
xlim( handles.axTraces, [0 time(end)] );
dt = handles.data.sampling/2/1000; %put idl transitions between FRET datapoints.

for i=1:handles.nTracesToShow,
    y_offset = 1.18*(handles.nTracesToShow-i) +0.2;
           
    plot( handles.axTraces, time([1,end]), y_offset+[0 0], 'k:' );  %baseline marker
        
    % Determine colors for data display
    dataField = handles.options.dataField;
    
    if strcmpi(dataField,'fret'),
        traceColor='b';
    elseif strcmpi(dataField,'fret2')
        traceColor='m';
    else
        idxch = find(strcmpi(dataField,data.channelNames) );
        if ~isempty(idxch) && any(idxch==handles.data.idxFluor) && ...
                              isfield(data.fileMetadata,'wavelengths')
            traceColor = Wavelength_to_RGB( data.fileMetadata.wavelengths(idxch) );
        else
            traceColor='k';
        end
    end
    
    % Draw traces, idealizations, and trace number labeles.
    handles.hFretLine(i) = plot( handles.axTraces, time, ...
                          y_offset+zeros(1,data.nFrames), 'Color',traceColor );

    handles.hIdlLine(i)  = stairs( handles.axTraces, time-dt, ...
                          y_offset+zeros(1,data.nFrames), 'r-' );
    
    handles.hTraceLabel(i) = text( 0.98*time(end),y_offset+0.1, '', ...
               'Parent',handles.axTraces, 'BackgroundColor','w', ...
               'HorizontalAlignment','right', 'VerticalAlignment','bottom' );
end

ylim(handles.axTraces,[0 1.2*handles.nTracesToShow]);
xlabel(handles.axTraces,'Time (s)');
set(handles.figure1,'pointer','arrow');

guidata(hObject,handles);
showTraces(handles);
showModelLines(handles);  %Draw state markers underneath data

% END FUNCTION lbFiles_Callback



function showTraces(handles)
% Update axTraces to show the current subset -- called by sldTraces.
% i is the index into the lines in the viewer.
% idx is the index into the traces in the whole file.

idxStart = get(handles.sldTraces,'Max')-floor(get(handles.sldTraces,'Value'));
nToShow = min(handles.nTracesToShow, handles.data.nTraces);

for i=1:nToShow
    idx = i+idxStart;
    y_offset = 1.18*(handles.nTracesToShow-i) +0.2;
    
    % Redraw FRET and idealization traces
    ydata = handles.data.(handles.options.dataField)(idx,:);
    ydata = y_offset + min(1.15, max(-0.15,ydata) );  %clip outliers, position w/i viewer
    set( handles.hFretLine(i), 'YData',ydata );
    
    if ~isempty(handles.idl)
        idldata = y_offset + handles.idl(idx,:);
        set( handles.hIdlLine(i), 'YData',idldata );
    end
    
    set( handles.hTraceLabel(i), 'String',sprintf('%d',idx) );
end

% Clear final lines if there isn't enough data to fill them.
set( handles.hFretLine(nToShow+1:end), 'YData',nan(1,handles.data.nFrames) );
set( handles.hIdlLine(nToShow+1:end),  'YData',nan(1,handles.data.nFrames) );
set( handles.hTraceLabel(nToShow+1:end), 'String','' );

% END FUNCTION showTraces



function handles = showModelLines(handles)
% Draw dotted lines to indicate model FRET values in trace viewer panel.

if ~handles.showStateMarkers || isempty(handles.dataFilenames) || isempty(handles.model)
    return;
end

% Remove any existing markers without completely clearing the graph.
delete( findall(handles.axTraces,'Tag','ModelMarker') );

nToShow = min(handles.nTracesToShow, handles.data.nTraces);
colors = handles.modelViewer.colors;
time = handles.data.time([1,end])/1000;
mu = handles.model.mu;
h = [];

% Redraw model FRET value markers.
for i=1:nToShow
    y_offset = 1.18*(handles.nTracesToShow-i) +0.2;
    
    for k=2:handles.model.nClasses
        h(end+1) = plot( handles.axTraces, time, y_offset+mu([k k]), ':', ...
                         'Color',colors{k}, 'Tag','ModelMarker' );
    end
end

uistack(h,'bottom');  %draw markers below data.

% END FUNCTION showModelLines



% --------------------------------------------------------------------
function wheelScroll_callback(~, eventData, handles) %#ok<DEFNU>
% Mouse wheel scrolling moves the trace viewer pane up and down.
% The event is triggered at the figure level.

loc = get(handles.sldTraces, 'Value')-3*eventData.VerticalScrollCount;
loc = min( loc, get(handles.sldTraces,'Max') );
loc = max( loc, get(handles.sldTraces,'Min') );
set(handles.sldTraces, 'Value', loc);  %triggers listener, updating viewer.

% END FUNCTION wheelScroll_callback



% --------------------------------------------------------------------
function mnuSim_Callback(hObject, ~, handles) %#ok<DEFNU>
% Simulate traces using current model.

if isempty(handles.model), return; end  %model required.

% Get simulation settings.
persistent opt;
if isempty(opt)
    opt = struct('nTraces',1000, 'nFrames',2000, 'sampling',40, ...
                 'snr',30, 'shotNoise',true, 'gamma',1, ...
                 'totalIntensity',500, 'stdTotalIntensity',100, ...
                 'stdPhoton',0, 'totalTimeOn',20 );
end
prompt = {'Traces',   'Frames',     'Sampling (ms)', ...
          'Signal:background noise ratio', 'Shot noise', 'Apparent gamma', ...
          'Intensity (photons)', 'Intensity stdev', ...
          'Excess noise stdev',  'FRET Lifetime (s)'};
newOpt = settingdlg(opt, fieldnames(opt), prompt); %, 'Simulation parameters');
if isempty(newOpt), return; end

% Get output filename from user.
[f,p] = uiputfile('sim.traces','Save simulated data as...');
if f==0, return; end  %user pressed "cancel"


% Simulate new data.
% FIXME: simulate.m should return a valid traces object.
set(handles.figure1,'pointer','watch');
set(handles.txtStatus,'String','Simulating...'); drawnow;

opt = newOpt;
newOpt.stdBackground = opt.totalIntensity/(sqrt(2)*opt.snr);
newOpt.kBleach = 1/opt.totalTimeOn;

try
    [~,data] = simulate( [opt.nTraces,opt.nFrames], opt.sampling/1000, handles.model, newOpt );
    saveTraces( fullfile(p,f), data );
catch e
    set(handles.figure1,'pointer','arrow');
    set(handles.txtStatus,'String',['Error: ' e.message]);
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


% --------------------------------------------------------------------
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
    opt = struct('sigmaPSF',0.8, 'density',handles.data.nTraces, ...
                 'aduPhoton',2, 'grid',false, 'alignX',0, 'alignY',0, ...
                 'alignTheta',0, 'alignScale',1);
end
prompt = {'PSF Size (px stdev):', 'Particles to simulate:', ...
          'ADU to photon conversion', 'Use a regular grid', 'X deviation (px):', ...
          'Y deviation (px):', 'Rotation (degrees):', 'Scaling factor'};
newOpt = settingdlg(opt, fieldnames(opt), prompt); %, 'Simulation parameters');
if isempty(newOpt), return; end
opt = newOpt;

newOpt.density = min(handles.data.nTraces, newOpt.density);

% Simulate the movie.
% (simulateMovie.m will ask for the movie filenames)
simulateMovie(handles.data, [],[], newOpt);

% END FUNCTION mnuSimMovie_Callback



% --------------------------------------------------------------------
function mnuFileRemove_Callback(hObject, ~, handles) %#ok<DEFNU>
% Close currently-selected file in GUI list.

idxfile = get(handles.lbFiles,'Value');
handles.dataFilenames(idxfile) = [];
handles.dwtFilenames(idxfile) = [];
guidata(hObject,handles);

names = get(handles.lbFiles,'String');
names(idxfile) = [];
set(handles.lbFiles,'String',names, 'Value',max(1,idxfile-1));

% Update GUI.
if numel(names)>0
    lbFiles_Callback(handles.lbFiles, [], handles);
else
    cla(handles.axTraces);
    enableControls(handles);
end

% END FUNCTION mnuFileRemove_Callback


% --------------------------------------------------------------------
function mnuFileRemoveAll_Callback(hObject, ~, handles) %#ok<DEFNU>
% Close all open files.

handles.dataFilenames = {};
handles.dwtFilenames = {};
set(handles.lbFiles, 'String',{}, 'Value',1);
guidata(hObject,handles);

cla(handles.axTraces);
enableControls(handles);

% END FUNCTION mnuFileRemoveAll_Callback


% --------------------------------------------------------------------
function mnuFileUp_Callback(hObject, ~, handles, inc) %#ok<DEFNU>
% Move currently-selected file up in the GUI list.
% The last parameter specifies the direction to move (+1 up, -1 down).

% Reorder in handles objects
names   = get(handles.lbFiles,'String');
idxfile = get(handles.lbFiles,'Value');  %currently selected file.
idxnew  = idxfile+inc;  %new position after move.
idxnew  = max(1, min(numel(names),idxnew) );  %can't move past the end

handles.dataFilenames = swap( handles.dataFilenames, idxfile, idxnew );
handles.dwtFilenames  = swap( handles.dwtFilenames,  idxfile, idxnew );
guidata(hObject,handles);

% Reorder within listbox control
names = swap( names, idxfile, idxnew );
set(handles.lbFiles, 'String',names, 'Value',idxnew);

% END FUNCTION mnuFileUp_Callback



function out = swap(in, idx1, idx2)
% Swap elements in a matrix by index.
out = in;
out(idx1) = in(idx2);
out(idx2) = in(idx1);
%end
