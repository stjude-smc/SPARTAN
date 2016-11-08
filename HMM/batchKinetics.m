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

% Last Modified by GUIDE v2.5 08-Nov-2016 15:21:07


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
handles.modelFilename = [];
handles.model = [];
handles.dataFilenames = {};
handles.dwtFilenames  = {};
handles.idl = [];

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

constants = cascadeConstants;
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );

% Trace viewer pane callbacks
addlistener( [handles.sldTraces handles.sldTracesX], 'Value', ...
             'PostSet',@(h,e)showTraces(guidata(e.AffectedObject))  );

ylim(handles.axTraces,[0 12]);
% xlim(handles.axTraces,[0 xx]);
hold(handles.axTraces,'on');

% END FUNCTION batchKinetics_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = batchKinetics_OutputFcn(~, ~, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;
% END FUNCTION batchKinetics_OutputFcn




%% =======================  CALLBACK FUNCTIONS  ======================= %%

function btnLoadData_Callback(hObject, ~, handles) %#ok<DEFNU>
% Executes on button press in btnLoadData.

% Prompt use for location to save file in...
handles.dataFilenames = getFiles([],'Select traces files to analyze');
handles.dataPath = pwd;

% If a model is loaded, enable the Execute button & update GUI
if ~isempty(handles.model),
    set(handles.btnExecute,'Enable','on');
end
set([handles.btnMakeplots handles.mnuViewMakeplots], 'Enable','on');

[~,names] = cellfun(@fileparts, handles.dataFilenames, 'UniformOutput',false);
set(handles.lbFiles, 'Value',[], 'String',names);

% Look for .dwt files if data were already analyzed.
try
    handles.dwtFilenames = findDwt(handles.dataFilenames);
  
    set( [handles.btnDwellhist handles.mnuDwellhist handles.btnPT ...
          handles.mnuViewPercentTime handles.mnuViewTPS ...
          handles.btnOccTime handles.mnuViewOccTime], 'Enable','on');
catch
end
guidata(hObject, handles);
cla(handles.axTraces);

% END FUNCTION btnLoadData_Callback



function btnLoadModel_Callback(hObject, ~, handles) %#ok<DEFNU>
% Executes on button press in btnLoadModel.

if isempty(handles.modelFilename),
    handles.modelFilename = fullfile(pwd,'*.qmf');
end

% Ask the user for a filename
[fname,p] = uigetfile( handles.modelFilename, 'Select a QuB model file...' );
if fname==0, return; end
fname = fullfile(p,fname);

% Load the model and show the model properties in the GUI.
% The model's properties are automatically updated whenever the model is
% modified in the GUI.
handles.model = QubModel(fname);
handles.modelViewer = QubModelViewer(handles.model, handles.axModel);
title(handles.axModel, ['...' fname(max(1,end-40):end)], 'interpreter','none');

% Enable relevant GUI controls
set([handles.btnSaveModel handles.tblFixFret], 'Enable','on');
set(handles.btnExecute,'Enable',onoff(~isempty(handles.dataFilenames)));

% Automatically update the parameter table when the model is altered.
handles.modelUpdateListener = addlistener(handles.model, ...
                        {'mu','fixMu','sigma','fixSigma'}, 'PostSet', ...
                        @(s,e)modelUpdate_Callback(handles.tblFixFret,e) );
handles.model.mu = handles.model.mu;  %trigger table update

% Update handles structure
guidata(hObject, handles);

% END FUNCTION btnLoadModel_Callback


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
% set(handles.btnExecute,'Enable','off');
% set(handles.btnStop,'Enable','on');

% Run the analysis algorithms...
[resultTree,handles.dwtFilenames] = runParamOptimizer(model,handles.dataFilenames,options); %#ok<ASGLU>

% Save results to file for later processing by the user.
save('resultTree.mat','resultTree');
% qub_saveTree(resultTree,resultFilename);
% qub_saveTree(resultTree.milResults(1).ModelFile,'result.qmf','ModelFile');

% Update GUI for finished status.
% set(handles.btnStop,'Enable','off');
set( [handles.btnExecute handles.btnDwellhist handles.btnMakeplots ...
      handles.mnuDwellhist handles.mnuViewPercentTime handles.btnPT ...
      handles.mnuViewTPS handles.mnuViewOccTime handles.btnOccTime], 'Enable','on');
disp('Finished!');

% Update handles structure
guidata(hObject, handles);

lbFiles_Callback(handles.lblFiles, [], handles);

% END FUNCTION btnExecute_Callback



function btnStop_Callback(~, ~, handles) %#ok<DEFNU>
% Executes on button press in btnStop.
% FIXME: this should stop any task mid-execution.
set(handles.btnExecute,'Enable','on');






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
enable = onoff( get(hObject,'Value')>1 );
set(handles.mnuIdlSettings, 'Enable',enable);

% END FUNCTION cboIdealizationMethod_Callback
    
    
% --- Executes on selection change in cboKineticsMethod.
function cboKineticsMethod_Callback(hObject, ~, handles)
% Idealization options: idealization method combo box

% Update method to use for kinetic parameter estimation.
text = get(hObject,'String');
handles.options.kineticsMethod = text{get(hObject,'Value')};
guidata(hObject, handles);

% If user selected "Do Nothing", disable idealization option controls.
enable = onoff( get(hObject,'Value')>1 );
set(handles.mnuKineticsSettings, 'Enable',enable);

% END FUNCTION cboKineticsMethod_Callback



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
    title(handles.axModel, ['...' fname(max(1,end-40):end)], 'interpreter', 'none');
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
% Change idealization (SKM) settings

prompt = {'Idealize traces individually:', 'Max iterations:'}; %'LL Convergence:', 'Grad. Convergence:'
fields = {'seperately', 'maxItr'};  %'gradLL', 'gradConv'
handles.options = settingsDialog(handles.options, fields, prompt);
guidata(hObject,handles);

% END FUNCTION mnuIdlSettings_Callback


function mnuKineticsSettings_Callback(hObject, ~, handles) %#ok<DEFNU>
% Change kinetic parameter estimation (MIL) settings

prompt = {'Bootstrap samples:', 'Dead time (frames):'};
fields = {'bootstrapN', 'deadTime'};
handles.options = settingsDialog(handles.options, fields, prompt);
guidata(hObject,handles);

% END FUNCTION mnuIdlSettings_Callback



% --------------------------------------------------------------------
function sldTraces_Callback(~, ~, handles) %#ok<DEFNU>
% User adjusted the trace view slider -- show a different subset of traces.
showTraces(handles);
% END FUNCTION sldTraces_Callback


% --------------------------------------------------------------------
function lbFiles_Callback(hObject, ~, handles)
% User selected a file. Show traces in the trace viewer panel

idxFile   = get(hObject,'Value');
data = loadTraces( handles.dataFilenames{idxFile} );
handles.data = data;

if ~isempty(handles.dwtFilenames{idxFile})
    [dwt,~,offsets,model] = loadDWT( handles.dwtFilenames{idxFile} );
    handles.idl = dwtToIdl(dwt, offsets, data.nFrames, data.nTraces);
    
%     for i=1:data.nTraces,
%         fretValues = [NaN; model(:,1)];
%         handles.idl(i,:) = fretValues( handles.idl(i,:)+1 );
%     end
    assert( size(model,2)==2 );
    fretValues = [NaN; model(:,1)];
    handles.idl = fretValues( handles.idl+1 );
else
    handles.idl = [];
end

set(handles.sldTraces,'Min',0,'Max',handles.data.nTraces-10,'Value',0);
set(handles.sldTracesX, 'Min',10, 'Max',handles.data.nFrames, 'Value',handles.data.nFrames);
set([handles.sldTraces handles.sldTracesX],'SliderStep',[0.01 0.1]); %move to GUIDE?

guidata(hObject,handles);
showTraces(handles);

% END FUNCTION lbFiles_Callback


function showTraces(handles)
% Update trace viewer panel

cla(handles.axTraces);

idxTraceStart = floor(get(handles.sldTraces,'Value'));
xlimit = floor(get(handles.sldTracesX,'Value'));
time = handles.data.time(1:xlimit);

for i=1:10,
    idx = i+idxTraceStart;
    y_offset = 1.18*(i-1) +0.2;
    plot( handles.axTraces, time([1,end]), y_offset+[0 0], 'k:' ); %baseline marker
    plot( handles.axTraces, time, y_offset+handles.data.fret(idx,1:xlimit), 'b-' );
    
    if ~isempty(handles.idl)
        plot( handles.axTraces, time, y_offset+handles.idl(idx,1:xlimit), 'r-' );
    end
    
    % TODO: draw trace number.
end

xlim(handles.axTraces, [0 time(xlimit)]); %FIXME: slow and only required when X slider is changed.

% END FUNCTION function



function wheelScroll_callback(~, eventData, handles) %#ok<DEFNU>
% Mouse wheel scrolling moves the trace viewer pane up and down.
% The event is triggered at the figure level.

% disp(get(gcbo,'type'));
% if ~strcmp(get(gcbo,'type'),'axes'), return; end %&& gco==handles.axTraces

% Update scroll bar position when scrolling with the mouse wheel.
loc = get(handles.sldTraces, 'Value')-2*eventData.VerticalScrollCount;
loc = min( loc, get(handles.sldTraces,'Max') );
loc = max( loc, get(handles.sldTraces,'Min') );
set(handles.sldTraces, 'Value', loc);

% The trace viewer is automatically updated by the Value property listener.

% END FUNCTION wheelScroll_callback

