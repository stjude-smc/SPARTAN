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

% Last Modified by GUIDE v2.5 15-Nov-2016 20:20:43


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
[handles.modelFilename,handles.model,handles.idl] = deal([]);
[handles.dataFilenames,handles.dwtFilenames] = deal({});
handles.nTracesToShow = 6;  %number displayed in trace display panel

% Set default analysis settings. FIXME: put these in cascadeConstants?
options.bootstrapN = 1;
options.deadTime = 0.5;
options.seperately = 1; %SKM: analyze each trace individually
options.maxItr = 100;
options.minStates = 1;
options.maxStates = 5;
options.maxRestarts = 10;
options.threshold = 1e-5;
handles.options = options;

% Update GUI to reflect these default settings. MIL not supported on Macs
methods = {'Segmental k-Means','Baum-Welch','ebFRET','MIL (Rate Optimizer)'};
if isempty(which('ebfret.analysis.hmm.vbayes')), methods(3)=[]; end
if ismac, methods(end)=[]; end

set( handles.cboIdealizationMethod, 'String',methods);
set( handles.cboIdealizationMethod, 'Value',1 );  %SKM
handles = cboIdealizationMethod_Callback(handles.cboIdealizationMethod,[],handles);

constants = cascadeConstants;
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );

% Trace viewer pane callbacks
handles.sldTracesListener(1) = addlistener( handles.sldTraces, 'Value', ...
          'PostSet',@(h,e)showTraces(guidata(e.AffectedObject))  );

handles.sldTracesListener(2) = addlistener( handles.sldTracesX, 'Value', ...
          'PostSet',@(h,e)sldTracesX_Callback(h,e,guidata(e.AffectedObject))  );

hold(handles.axTraces,'on');
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
handles.dataPath = pwd;

% If a model is loaded, enable the Execute button & update GUI
if ~isempty(handles.model),
    set(handles.btnExecute,'Enable','on');
end
set([handles.btnMakeplots handles.mnuViewMakeplots], 'Enable','on');

[~,names] = cellfun(@fileparts, handles.dataFilenames, 'UniformOutput',false);
set(handles.lbFiles, 'Value',1, 'String',names);

% Look for .dwt files if data were already analyzed.
handles.dwtFilenames = findDwt(handles.dataFilenames);

if ~any( cellfun(@isempty,handles.dwtFilenames) )
    set( [handles.btnDwellhist handles.mnuDwellhist handles.btnPT ...
          handles.mnuViewPercentTime handles.mnuViewTPS ...
          handles.btnOccTime handles.mnuViewOccTime], 'Enable','on');
end

% Show the first file.
lbFiles_Callback(handles.lbFiles, [], handles);

% END FUNCTION btnLoadData_Callback



function btnLoadModel_Callback(hObject, ~, handles) %#ok<DEFNU>
% Executes on button press in btnLoadModel.

if isempty(handles.modelFilename),
    handles.modelFilename = fullfile(pwd,'*.qmf');
end

% Ask the user for a filename
[fname,p] = uigetfile( handles.modelFilename, 'Select a QuB model file...' );
if fname==0, return; end

% Load the model and show the model properties in the GUI.
% The model's properties are automatically updated whenever the model is
% modified in the GUI.
handles.model = QubModel( fullfile(p,fname) );
handles.modelViewer = QubModelViewer(handles.model, handles.axModel);
title(handles.axModel, fname, 'interpreter','none');

% Enable relevant GUI controls
set([handles.btnSaveModel handles.tblFixFret], 'Enable','on');
set(handles.btnExecute,'Enable',onoff(~isempty(handles.dataFilenames)));

% Automatically update the parameter table when the model is altered.
handles.modelUpdateListener = addlistener(handles.model,'UpdateModel', ...
                        @(s,e)modelUpdate_Callback(handles.tblFixFret,e) );
handles.model.mu = handles.model.mu;  %trigger table update

guidata(hObject, handles);

% END FUNCTION btnLoadModel_Callback


function btnExecute_Callback(hObject, ~, handles) %#ok<DEFNU>
% Run the data analysis pipeline with user-specified data & model.

% Verify data and model have been specified by user in GUI.
idxfile = get(handles.lbFiles,'Value');

if isempty(handles.model) || isempty(idxfile),
    set(handles.btnExecute,'Enable','off');
    warning('Missing model or data');
    return;
end

trcfile  = handles.dataFilenames{idxfile};
dwtfname = handles.dwtFilenames{idxfile};

% Verify external modules installed
if strcmpi(handles.options.idealizeMethod,'ebFRET') && isempty(which('ebfret.analysis.hmm.vbayes'))
    errordlg('ebFRET not found. Check your path.',mfilename);
    disp('Go to https://ebfret.github.io/ to download ebFRET, then add to the MATLAB path.');
    return;
end

% Update GUI for "Running" status.
% set(handles.btnExecute,'Enable','off');
% set(handles.btnStop,'Enable','on');

% Run the analysis algorithms...
% FIXME: ideally we want idl (or dwt) returned directly for speed.
if strcmpi(handles.options.idealizeMethod(1:3),'MIL')
    if isempty(dwtfname) || ~exist(dwtfname,'file'),
        errordlg('Traces must be idealized before running MIL');
        return;
    end
    
    optModel = milOptimize(dwtfname, handles.model, handles.options);
    handles.model.rates = optModel.rates;
else
    % Clear current idealization (FIXME also delete .dwt?)
    handles.idl = [];
    set(handles.hIdlLine,'Visible','off');
    
    [handles.dwtFilenames{idxfile},optModel] = runParamOptimizer(...
                                 handles.model, trcfile, handles.options);
end

if get(handles.chkUpdateModel,'Value'),
    handles.model.rates = optModel.rates;
    handles.model.mu    = optModel.mu;
    handles.model.sigma = optModel.sigma;
    handles.model.p0    = optModel.p0;
end

handles.modelViewer.redraw();

% Save results to file for later processing by the user.
% save('resultTree.mat','resultTree');
% qub_saveTree(resultTree,resultFilename);
% qub_saveTree(resultTree.milResults(1).ModelFile,'result.qmf','ModelFile');

% Update GUI for finished status.
% set(handles.btnStop,'Enable','off');
set( [handles.btnExecute handles.btnDwellhist handles.btnMakeplots ...
      handles.mnuDwellhist handles.mnuViewPercentTime handles.btnPT ...
      handles.mnuViewTPS handles.mnuViewOccTime handles.btnOccTime], 'Enable','on');
disp('Finished!');

% Reload and draw idealization.
handles.idl = loadIdl(handles);
guidata(hObject,handles);
showTraces(handles);

% END FUNCTION btnExecute_Callback


function idl = loadIdl(handles)
% Returns the idealization for the currently selected file

dwtfname = handles.dwtFilenames{ get(handles.lbFiles,'Value') };

if ~isempty(dwtfname)
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
% Idealization options: idealization method combo box

% Update method to use for idealization
text = get(hObject,'String');
handles.options.idealizeMethod = text{get(hObject,'Value')};
guidata(hObject, handles);

if strcmpi(handles.options.idealizeMethod(1:3),'MIL')
    set(handles.chkUpdateModel,'Enable','off');
else
    set(handles.chkUpdateModel,'Enable','on');
end

% END FUNCTION cboIdealizationMethod_Callback


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

model = event.Source;

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
    handles.model.save( fullfile(p,f) );
    title(handles.axModel, f, 'interpreter', 'none');
end

% END FUNCTION btnSaveModel_Callback


% ========================  PLOTTING FUNCTIONS  ======================== %
% Executed when plotting menu or toolbar buttons are clicked.

function btnMakeplots_Callback(~, ~, handles) %#ok<DEFNU>
makeplots(handles.dataFilenames);

function btnDwellhist_ClickedCallback(~, ~, handles) %#ok<DEFNU>
if ~isempty(handles.dwtFilenames),
    dwellhist(handles.dwtFilenames);
end

function btnPT_ClickedCallback(~, ~, handles) %#ok<DEFNU>
if ~isempty(handles.dwtFilenames),
    percentTime(handles.dwtFilenames);
end

function btnOccTime_ClickedCallback(~, ~, handles) %#ok<DEFNU>
if ~isempty(handles.dwtFilenames),
    occtime(handles.dwtFilenames);
end

function mnuViewTPS_Callback(~, ~, handles) %#ok<DEFNU>
if ~isempty(handles.dwtFilenames),
    transitionsPerSecond(handles.dwtFilenames);
end



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

handles.options = settingsDialog(handles.options, fields, prompt);
guidata(hObject,handles);
% END FUNCTION mnuIdlSettings_Callback



% ========================  TRACE VIEWER PANEL  ======================== %

% --------------------------------------------------------------------
function sldTraces_Callback(~, ~, handles) %#ok<DEFNU>
% User adjusted the trace view slider -- show a different subset of traces.
showTraces(handles);
% END FUNCTION sldTraces_Callback

function sldTracesX_Callback(~, ~, handles)
% User adjusted the trace view slider -- show a different subset of traces.

xlimit = floor(get(handles.sldTracesX,'Value'));
set( handles.axTraces, 'XLim',[0 handles.data.time(xlimit)] );

for i=1:handles.nTracesToShow,
    p = get(handles.hTraceLabel(i), 'Position');
    p(1) = 0.98*handles.data.time(xlimit);
    set( handles.hTraceLabel(i), 'Position',p );
end

% END FUNCTION sldTraces_Callback


% --------------------------------------------------------------------
function lbFiles_Callback(hObject, ~, handles)
% User selected a file. Show traces in the trace viewer panel.
% FIXME: could be somewhat faster if plotting one long trace rather than
% many line objects...

idxFile = get(hObject,'Value');
data = loadTraces( handles.dataFilenames{idxFile} );
handles.data = data;
handles.idl = loadIdl(handles);

[handles.sldTracesListener.Enabled] = deal(false);
set(handles.sldTraces,  'Min',0,  'Max',data.nTraces-handles.nTracesToShow, 'Value',data.nTraces-handles.nTracesToShow);
set(handles.sldTracesX, 'Min',10, 'Max',data.nFrames,    'Value',data.nFrames);
[handles.sldTracesListener.Enabled] = deal(true);

% Setup axes for plotting traces.
% Some code duplication with showTraces().
cla(handles.axTraces);

xlimit = floor(get(handles.sldTracesX,'Value'));
time = handles.data.time(1:xlimit);
xlim( handles.axTraces, [0 time(end)] );

for i=1:handles.nTracesToShow,
    y_offset = 1.18*(handles.nTracesToShow-i) +0.2;
           
    plot( handles.axTraces, time([1,end]), y_offset+[0 0], 'k:' );  %baseline marker
    
    handles.hFretLine(i) = plot( handles.axTraces, time, ...
                              y_offset+zeros(1,xlimit), 'b-' );

    handles.hIdlLine(i)  = plot( handles.axTraces, time, ...
                              y_offset+zeros(1,xlimit), 'r-' );
    
    handles.hTraceLabel(i) = text( 0.98*time(end),y_offset+0.1, '', ...
               'Parent',handles.axTraces, 'BackgroundColor','w', ...
               'HorizontalAlignment','right', 'VerticalAlignment','bottom' );
end

ylim(handles.axTraces,[0 1.2*handles.nTracesToShow]);

guidata(hObject,handles);
showTraces(handles);

% END FUNCTION lbFiles_Callback



function showTraces(handles)
% Update axTraces to show the current subset -- called by sldTraces.
% FIXME: this will crash if there are less than N traces in the file!

idxStart = get(handles.sldTraces,'Max')-floor(get(handles.sldTraces,'Value'));

for i=1:handles.nTracesToShow,
    idx = i+idxStart;
    y_offset = 1.18*(handles.nTracesToShow-i) +0.2;
    set( handles.hFretLine(i), 'YData',y_offset+handles.data.fret(idx,:) );
    
    if ~isempty(handles.idl)
        set( handles.hIdlLine(i), 'YData', y_offset+handles.idl(idx,:) );
    end
    
    set( handles.hTraceLabel(i), 'String',sprintf('%d',idx) );
end

set(handles.hIdlLine, 'Visible',onoff(~isempty(handles.idl)) );

% END FUNCTION function



function wheelScroll_callback(~, eventData, handles) %#ok<DEFNU>
% Mouse wheel scrolling moves the trace viewer pane up and down.
% The event is triggered at the figure level.

loc = get(handles.sldTraces, 'Value')-3*eventData.VerticalScrollCount;
loc = min( loc, get(handles.sldTraces,'Max') );
loc = max( loc, get(handles.sldTraces,'Min') );
set(handles.sldTraces, 'Value', loc);

% The trace viewer is automatically updated by the Value property listener.

% END FUNCTION wheelScroll_callback
