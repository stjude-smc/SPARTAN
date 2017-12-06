function varargout = autotrace(varargin)
% AUTOTRACE  Select traces according to user-defined criteria
%
%   Autotrace is a GUI that displays histograms of statistics calculated
%   for each trace in loaded file(s). These are calculated in traceStat.m.
%   The user can then select traces according to a defined set of criteria
%   and save that subset to a new file, typically ending in "_auto.traces".
%   A .log file is also saved that includes the criteria and other details.
%
%   NOTE: selection can introduce bias in the data analysis process.
%   Ensure that the selected subset reflects the whole ensemble of traces
%   and use consistent criteria when comparing datasets.
%
%   See also traceStat, pickTraces, loadPickSaveTraces.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


% Last Modified by GUIDE v2.5 14-Dec-2016 10:33:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @autotrace_OpeningFcn, ...
    'gui_OutputFcn',  @autotrace_OutputFcn, ...
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






%#########################################################################
%------------------------- INITIALIZATION (GUI) -------------------------%
%#########################################################################


%----------INITIALIZATION OF THE GUI----------%
function autotrace_OpeningFcn(hObject, ~, handles, varargin)
% Executes just before autotrace is made visible.

updateSpartan; %check for updates
constants = cascadeConstants();
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );

% Standard criteria control names and corresponding criteria.
% See PickTraces_Callback() and setCriteria() functions.
handles.chkNames = {'chk_acclife','chk_corr','chk_corr','chk_snr','chk_bg','chk_ncross','chk_maxTotalSigma'};
handles.txtNames = {'ed_acclife','ed_min_corr','ed_max_corr','ed_snr','ed_bg','ed_ncross','ed_maxTotalSigma'};
handles.criteriaNames = {'min_acclife','min_corr','max_corr','min_snr','max_bg','max_ncross','maxTotalSigma'};

%---- Setup drop-down boxes listing trace statistics.
% Get names of trace statistics
ln = traceStat;  %get long statistic names
handles.statLongNames = ln;
longNames  = struct2cell(ln);
shortNames = fieldnames(ln);

% Populate statistic dropdowns above histograms and set defaults.
handles.cboNames = strcat('cboStat',{'1','2','3','4'});
statNames = {'t','donorlife','corr','snr'};
for id=1:length(handles.cboNames),
    set( handles.(handles.cboNames{id}), 'String',longNames, ...
             'Value',find(strcmp(statNames{id},shortNames)) );
end

% Setup special criteria selection drop-down boxes.
handles.nCriteriaBoxes = 7;
criteriaNames = [{''}; longNames];

for id=1:handles.nCriteriaBoxes
    set( handles.(['cboCriteria' num2str(id)]), 'String', criteriaNames );
end

%---- Add context menus to the plots to launch curve fitting or copy data.
zoom off;
hZoom = zoom(handles.figure1);
set(hZoom, 'UIContextMenu', handles.mnuStat);
zoom on;

% Choose default command line output for autotrace
handles.output=hObject;
guidata(hObject,handles);

% Set controls to default selection criteria
mnuCriteriaReset_Callback(hObject, [], handles);

% END FUNCTION autotrace_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = autotrace_OutputFcn(~, ~, handles)
% Get default command line output from handles structure
varargout{1}=handles.output;
% END FUNCTION autotrace_OutputFcn






%#########################################################################
%----------------------- LOAD, FILTER, SAVE TRACES ----------------------%
%#########################################################################


%----------OPENS INDIVIDUAL TRACES FILES----------%

% --- Executes on button press in OpenTracesFile.
function OpenTracesFile_Callback(hObject, ~, handles) %#ok<*DEFNU>
% This method is called when the user clicks the "Open Traces File.."
% button in the GUI.  The trace is parsed,

% Prompt user to select traces files to load.
filter = {'*.rawtraces;*.traces','All traces files (*.rawtraces,*.traces)'; ...
          '*.rawtraces','Raw traces files (*.rawtraces)'; ...
          '*.*','All files (*.*)'  };
      
[datafile,datapath] = uigetfile( filter,'Choose a traces file:', ...
                                                'MultiSelect','on');
if datapath==0, return; end %user hit cancel

% Convert filename list into a cell array
if ~iscell(datafile), datafile = {datafile}; end
filename = strcat(datapath,datafile);

handles.inputdir = datapath;
handles.inputfiles = filename;

%--- Update GUI listing of number of files and file types
for i=1:numel(handles.inputfiles),
    disp( handles.inputfiles{i} );
end

if numel(filename) == 1,
    fileDisplayText = handles.inputfiles{1};
else
    % Find common substring within files.
    tmp = char(filename{:});
    first = find( sum(abs(bsxfun(@minus,tmp(1,:),tmp)))~=0, 1,'first' ); %find first differing character
    fileDisplayText = sprintf('%s... (%d files)', tmp(1,1:first-1), numel(filename));
end

% Shorten display name by removing initial part of the path if too long.
if length(fileDisplayText)>100,
    fileDisplayText = ['...' fileDisplayText(end-100:end)];
end
set(handles.editFilename,'String', fileDisplayText);

% Load the traces files.
OpenTracesBatch( hObject, handles );

% END FUNCTION OpenTracesFile_Callback





%----------OPEN ALL TRACES FILES WITHIN THE SAME DIRECTORY----------%
% --- Executes on button press in btnOpenDirectory.
function btnOpenDirectory_Callback(hObject, ~, handles)
% Load all .rawtraces files in the selected directory

% Select directory by user interface.
datapath=uigetdir;
if datapath==0, return; end %user hit cancel.

% Create list of .rawtraces files in the directory.
target = [datapath filesep '*.rawtraces'];
traces_files = dir(target);

if numel(traces_files) == 0
    disp('No files in this directory!');
    return;
end

handles.inputdir = datapath;
handles.inputfiles = strcat( [datapath filesep], {traces_files.name} );

% Update the GUI with the new data location name.
fileDisplayText = sprintf('%s (%d files)', target, numel(traces_files));
set(handles.editFilename, 'String',fileDisplayText);

% Load the traces files.
OpenTracesBatch( hObject, handles );

% END FUNCTION btnOpenDirectory_Callback





%----------BATCH ANALYSIS----------%
function btnBatchMode_Callback(hObject, ~, handles, datapath)
% For each .rawtraces file in the current directory, and optionally all
% sub-directories, select traces according to the current criteria and save
% a corresponding _auto.traces file.

skipExisting = strcmpi( get(handles.mnuBatchOverwrite,'Checked'), 'on');
recursive = strcmpi( get(handles.mnuBatchRecursive,'Checked'), 'on');

% Select directory by user interface.
% A path is given if called from the timer for automatic processing.
if nargin>=4 && exist(datapath,'dir'),
    auto = true;
else
    % User clicked the "Batch Mode" button directly. Ask for a location.
    datapath=uigetdir;
    if datapath==0, return; end
    auto = false;
end

% Get a list of all raw traces files under the current directory.
trace_files = regexpdir(datapath,'^.*\.rawtraces$', recursive);

if numel(trace_files)==0,
    return;
end
trace_files = {trace_files.name};

% For each file in the user-selected directory
for i=1:numel(trace_files)
    
    % If running in the automatic mode, skip already-processed files.
    [p,f] = fileparts(trace_files{i});
    if skipExisting || auto,
        auto_name = fullfile(p,[f '_auto.traces']);
        if exist(auto_name,'file'), continue; end
    end
    
    handles.inputdir = p;
    handles.inputfiles = trace_files(i);
    set(handles.editFilename,'String',trace_files{i});

    % Load traces and calculate statistics.
    % (this calls PickTraces_Callback, which saves GUI data).
    handles = OpenTracesBatch( hObject, handles );

    % Save picked traces to a new _auto.txt file.
    if handles.picked_mols > 0
        disp( handles.outfile );
        SaveTraces( handles.outfile, handles );
    else
        disp('No files to save, skipping');
    end
end

% END FUNCTION btnGo_Callback





function handles = OpenTracesBatch( hObject, handles )
% Calculate and display trace statistics for the current list of files.

set(handles.figure1, 'pointer','watch'); drawnow;

% Clear out old data to save memory.
if isappdata(handles.figure1,'infoStruct') %if data previously loaded.
    rmappdata(handles.figure1,'infoStruct');
end

% Determine default filename to use when saving.
[p,f] = fileparts( handles.inputfiles{1} );
handles.outfile = fullfile(p, [f '_auto.traces']);

% Calculate trace stats
try
    [infoStruct,nTracesPerFile] = traceStat(handles.inputfiles);
catch e
    if ~strcmpi(e.identifier,'parfor_progressbar:cancelled')
        errordlg( ['Error: ' e.message], mfilename );
    end
    set(handles.figure1, 'pointer','arrow'); drawnow;
    return;
end

handles.nTraces = numel( infoStruct );
handles.nTracesPerFile = nTracesPerFile;
                    
% Save the trace properties values to application data
setappdata(handles.figure1,'infoStruct', infoStruct);
clear infoStruct;

% Select traces according to the current criteria.
handles = PickTraces_Callback(hObject,handles);

set(handles.figure1, 'pointer', 'arrow'); drawnow;

% END FUNCTION OpenTracesBatch



%---------------  SAVE PICKED TRACES TO FILE (CALLBACK) ---------------%
% --- Executes on button press in SaveTraces.
function outfile = SaveTraces_Callback(hObject, ~, handles)
% Save the currently selected traces to a new _auto.traces file.

% Create a name for the output file
[f,p] = uiputfile('.traces','Save picked traces as:',handles.outfile);
if f==0,
    outfile = [];
    return;
else
    % Save picked traces to a new _auto.txt file.
    outfile = fullfile(p,f);
    SaveTraces( outfile, handles );
end

% Disabling the button as a means of confirming operation success.
handles.outfile = outfile;
guidata(hObject,handles);




%--------------------  SAVE PICKED TRACES TO FILE --------------------%
function SaveTraces( filename, handles )

set(handles.figure1, 'pointer', 'watch'); drawnow;

% Build list of trace indexes in each file
picks = handles.inds_picked;
nTracesPerFile = handles.nTracesPerFile;

idxStart = cumsum([0; nTracesPerFile(1:end-1)]);
nFiles = numel(nTracesPerFile);
picksByFile = cell( nFiles,1 );

for i=1:nFiles,
    picksByFile{i} = picks(  picks>idxStart(i) & picks<=idxStart(i)+nTracesPerFile(i)  ) ...
                   - idxStart(i);
end

% Save selected traces to a new file
options.indexes = picksByFile;
options.stats = getappdata(handles.figure1,'infoStruct');
options.outFilename = filename;

loadPickSaveTraces( handles.inputfiles, handles.criteria, options );

set(handles.figure1, 'pointer', 'arrow');
set([handles.mnuFileSave handles.tbFileSave],'Enable','off');
drawnow;


% END FUNCTION SaveTraces_Callback



%----------------  SAVE MOLECULE PROPERTIES TO FILE ----------------%

% --- Executes on button press in btnSaveProperties.
function btnSaveProperties_Callback(hObject, ~, handles)

% Get the output filename fromt he user
[p,f] = fileparts( handles.outfile );
filename = fullfile(p, [f '_prop.txt']);

[f,p] = uiputfile('.txt','Save statistics as:',filename);
if f==0,  return;  end
filename = fullfile(p,f);

% Retrieve trace statistics, extract to matrix
stats = getappdata(handles.figure1,'infoStruct');
stats = stats(handles.inds_picked);

data = cellfun( @double, struct2cell(stats) );
data = squeeze(data);
names = fieldnames(stats);

% Write column headers
fid = fopen( filename, 'w' );
fprintf( fid, '%s\t', names{:} );
fprintf( fid, '\n' );
fclose(fid);

% Write trace statistics
dlmwrite( filename, data', '-append', 'delimiter','\t' );

% Disable button as a means of confirming operation success.
set(handles.mnuFileSaveProp,'Enable','off');
guidata(hObject,handles);

% END FUNCTION SaveProperties_Callback
      




% --------------- VIEW TRACES FUNCTIONALITY --------------- %

% --- Executes on button press in ViewPickedTraces.
function ViewPickedTraces_Callback(hObject, ~, handles)

% If not already, save the selected traces to file.
if strcmpi( get(handles.mnuFileSave,'Enable'), 'on' );
    outfile = SaveTraces_Callback(hObject, [], handles);
else
    outfile = handles.outfile;
end

% Run sorttraces interface so traces can be viewed.
if ~isempty(outfile),
    sorttraces(0, outfile);
end





%----------APPLIES PICKING CRITERIA TO TRACES----------
% --- Executes on button press in PickTraces.
function handles = PickTraces_Callback(hObject, handles)
%

criteria = struct();

% Get criteria associated with the "free-form" drop-down boxes.
shortNames = fieldnames(handles.statLongNames);
equalityText = {'min_','max_','eq_'};

for id=1:handles.nCriteriaBoxes
    selection = get( handles.(['cboCriteria' num2str(id)]), 'Value' );
    equality  = get( handles.(['cboEquality' num2str(id)]), 'Value' );
    edText    = get( handles.(['edCriteria'  num2str(id)]), 'String');
    
    if selection>1 && equality>1 && ~isempty(edText),
        name = [equalityText{equality-1} shortNames{selection-1}];
        criteria.(name) = str2double(edText);
    end
end

% Get criteria values for all the fixed GUI elements.
% chk_corr is listed twice because it applies to two distinct criteria.
% Overlap must be handled seperately since it doesn't have an associated textbox.
for i=1:numel(handles.chkNames),
    hchk = handles.chkNames{i};
    htxt = handles.txtNames{i};
    name = handles.criteriaNames{i};
    
    if get(handles.(hchk),'Value'),
        criteria.(name) = str2double( get(handles.(htxt),'String') );
    end
end

if get(handles.chk_overlap,'Value'),
    criteria.eq_overlap = 0;  %zero means remove overlapped traces.
end

handles.criteria = criteria;
guidata(hObject,handles);  %for saving criteria before traces loaded.

% Find which molecules pass the selection criteria
stats = getappdata(handles.figure1,'infoStruct');
if isempty(stats), return; end %no data loaded, nothing to do.

handles.inds_picked = pickTraces( stats, criteria );
handles.picked_mols = numel(handles.inds_picked);
guidata(hObject,handles);

set(handles.MoleculesPicked,'String', ...
            sprintf('%d of %d',[handles.picked_mols,handles.nTraces]));

% If at least one trace is picked, turn some buttons on.
set( [handles.mnuFileSave handles.tbFileSave handles.mnuViewPlots ...
      handles.tbViewPlots handles.mnuViewTraces handles.tbViewTraces ...
      handles.mnuFileSaveProp], 'Enable',onoff(handles.picked_mols>0) );

% Update trace statistic histograms for traces passing selection criteria.
for i=1:length(handles.cboNames),
    cboStat_Callback(handles.(handles.cboNames{i}), [], handles);
end

% END FUNCTION PickTraces_Callback




%#########################################################################
%---------------------- HISTOGRAMS & CONTOUR PLOTS ----------------------%
%#########################################################################

%----------MAKE CONTOUR PLOT----------
% --- Executes on button press in MakeContourPlot.
function MakeContourPlot_Callback(hObject, ~, handles)
% Display FRET contour plot of currently selected traces.

% If not already, save the currently selected traces to file.
if strcmpi( get(handles.mnuFileSave,'Enable'), 'on' );
    outfile = SaveTraces_Callback(hObject, [], handles);
else
    outfile = handles.outfile;
end

% Run makeplots to display FRET contour plot
if ~isempty(outfile),
    makeplots(outfile);
end

% END FUNCTION MakeContourPlot_Callback




%#########################################################################
%------------------------- INPUT FIELD HANDLING -------------------------%
%#########################################################################

%----------MENU OF CONTOUR PLOT SETTINGS----------

% --- Executes on selection change in any of the drop down boxes above each
%     of the trace statistic histogram plots.
function cboStat_Callback(hObject, ~, handles)
% 

% Get ID of this combo control
id = get(hObject,'UserData');

% Get trace statistics
stats = getappdata(handles.figure1,'infoStruct');
stats = stats(handles.inds_picked);

% Get user selection
selected   = get(hObject,'Value');
statNames = fieldnames(handles.statLongNames);
statToPlot = statNames{selected};

% Make sure it is a recognized stat
if ~ismember( statToPlot, fieldnames(stats) ),
    error('Selected trace statistic is unknown');
end

% Define bin positions for certain parameters.
if strcmp(statToPlot,'corr'),
    bins = -1:0.1:1;
    statDecimals = '%.2f';
else
    bins = 40; %let hist choose N bin centers.
    statDecimals = '%.1f';
end

% Plot the distribution of the statistic
statData = [stats.(statToPlot)];
if any(isnan(statData))
    disp( 'warning: NaN values found' );
%     statData( isnan(statData) ) = 0;
end

% SNRs data sometimes goes to infinity and the histograms cannot be
% plotted. Avoid this by setting an absolute maximum
if strcmp(statToPlot,'snr_s') || strcmp(statToPlot,'snr'),
    statData(statData>1000) = 1000;
end

[data,binCenters] = hist( statData,bins );
data = 100*data/sum(data);  %normalize the histograms

ax = handles.(['axStat' num2str(id)]);
bar( ax, binCenters, data, 1 );
grid(ax,'on');

if id==1,
    ylabel(handles.axStat1, 'Count (%)');
end

% Save histogram data in plot for launching cftool
if length(stats)>=1,
    set( ax, 'UserData', [binCenters;data] );
end

% Display a mean value for easier interpretation.
set(  handles.(['txtStat' num2str(id)]), 'String', ...
             sprintf(['Mean: ' statDecimals],nanmean([stats.(statToPlot)]))  );

% END FUNCTION cboStat_Callback



function launchFitTool_Callback(ax)
% Callback for context menu for trace statistics plots.
% Launches Curve Fitting Tool using the data in the selected plot.

% Get ID of this combo control
histData = get(ax,'UserData');
if isempty(histData) || numel(histData)<2, return; end

binCenters = histData(1,:);
data = histData(2,:);

cftool(binCenters,data);


function copyPlotData_Callback(ax)
% Callback for context menu for trace statistics plots.
% Launches Curve Fitting Tool using the data in the selected plot.

% Get ID of this combo control
histData = get(ax,'UserData');
if isempty(histData) || numel(histData)<2, return; end

% Copy to clipboard
clipboard('copy', sprintf('%f %f\n',histData) );




% --- Executes on button press in chkAutoBatch.
function chkAutoBatch_Callback(hObject, ~, handles)
% Creates a timer to look for and automatically process .rawtraces files.
 
% If another timer is running, stop it.
fileTimer = timerfind('Name','autotrace_fileTimer');
isRunning = strcmpi(get(handles.mnuAutoBatch,'Checked'), 'off');

if ~isempty(fileTimer),
    stop(fileTimer);
    delete(fileTimer);
end

% Start a new timer if requested
if isRunning,
    % Ask the user for a directory location
    targetDir = uigetdir('','Choose directory:');
    if targetDir==0,  %user hit cancel
        return;
    end
    disp(targetDir);
    
    % Start a thread that will periodically check for new .rawtraces files.
    fileTimer = timer('ExecutionMode','fixedSpacing','StartDelay',1, ...
                      'Name','autotrace_fileTimer',...
                      'TimerFcn', {@updateFileTimer,hObject,targetDir}, ...
                      'StopFcn',{@stopFileTimer,hObject}, ...
                      'Period',2.0, 'BusyMode','drop');
    start(fileTimer);
    set(handles.mnuAutoBatch,'Checked','on');
else
    set(handles.mnuAutoBatch,'Checked','off');
end %if

% END FUNCTION chkAutoBatch_Callback



function stopFileTimer(~,~,hObject)
% This function is called when there is an error during the timer callback
% or when the timer is stopped.
if ishandle(hObject)
    handles = guidata(hObject);
    set(handles.mnuAutoBatch,'Checked','off');
end

% END FUNCTION stopFileTimer


function updateFileTimer(~,~,hObject,targetDir)
% This function runs each time the timer is fired, looking for any new
% movies that may have appeared on the path.

if ~ishandle(hObject),
    fileTimer = timerfind('Name','autotrace_fileTimer');
    stop(fileTimer);
    delete(fileTimer);
    return;
end
handles = guidata(hObject);

% Kill the timer if the directory is inaccessible
if ~exist(targetDir,'dir'),
    disp('Autotrace: stopping batch mode: directory was moved or is inaccessible');
    set(handles.mnuAutoBatch,'Checked','off');
    
    fileTimer = timerfind('Name','autotrace_fileTimer');
    stop(fileTimer);
    delete(fileTimer);
end

btnBatchMode_Callback( handles.figure1, [], handles, targetDir );

% END FUNCTION updateFileTimer



%#########################################################################
%------------------------ LOAD AND SAVE CRITERIA ------------------------%
%#########################################################################

% --------------------------------------------------------------------
function mnuCriteriaLoad_Callback(hObject, ~, handles)
% Load criteria from .mat file saved by autotrace.
fname = getFile('*.mat','Load Criteria');
if isempty(fname), return; end  %user hit cancel
input = load(fname);

if ~isfield(input,'criteria'),
    errordlg('Invalid criteria file');
else
    handles.criteria = input.criteria;
    guidata(hObject,handles);
    setCriteria(handles);
end


% --------------------------------------------------------------------
function mnuCriteriaSave_Callback(~, ~, handles)
% Save criteria struct to a .mat file.
[f,p] = uiputfile('criteria.mat','Save criteria:');
if isequal(f,0), return; end  %user hit cancel
criteria = handles.criteria; %#ok<NASGU>
save(fullfile(p,f),'criteria');


% --------------------------------------------------------------------
function mnuCriteriaReset_Callback(hObject, ~, handles)
% Reset criteria to factory defaults in cascadeConstants.
constants = cascadeConstants;
criteria = constants.defaultAutotraceCriteria;
handles.criteria = criteria;
guidata(hObject,handles);
setCriteria(handles);


% --------------------------------------------------------------------
function setCriteria(handles)
% Update GUI controls to reflect current criteria values.
% Called on load and by mnuCriteriaReset_Callback().

criteria = handles.criteria;

% Set standard criteria
for i=1:numel(handles.chkNames),
    hchk = handles.chkNames{i};
    htxt = handles.txtNames{i};
    name = handles.criteriaNames{i};
    
    set( handles.(hchk), 'Value',isfield(criteria,name) );
    if isfield(criteria,name),
        set( handles.(htxt), 'String', num2str(criteria.(name)) );
    end
end
overlap = isfield(criteria, 'eq_overlap') && criteria.eq_overlap==0;
set( handles.chk_overlap, 'Value',overlap);

% Clear non-standard fields
for i=1:handles.nCriteriaBoxes,
    set( handles.(['cboCriteria' num2str(i)]), 'Value', 1  );
    set( handles.(['cboEquality' num2str(i)]), 'Value', 1  );
    set( handles.(['edCriteria'  num2str(i)]), 'String','' );
end

% Set non-standard criteria.
fields = setdiff( fieldnames(criteria), [handles.criteriaNames 'eq_overlap'] );
assert( numel(fields)<=handles.nCriteriaBoxes, 'Too many criteria' );

shortNames = fieldnames(handles.statLongNames);
equalityText = {'min','max','eq'};

for i=1:numel(fields),
    temp = strsplit(fields{i},'_');
    [equality,name] = temp{:};
    
    set( handles.(['cboCriteria' num2str(i)]), 'Value',  find(strcmp(name,shortNames))+1       );
    set( handles.(['cboEquality' num2str(i)]), 'Value',  find(strcmp(equality,equalityText))+1 );
    set( handles.(['edCriteria'  num2str(i)]), 'String', num2str(criteria.(fields{i}))         );
end

% Update trace selections
PickTraces_Callback(handles.figure1, handles);

% END FUNCTION setCriteria


