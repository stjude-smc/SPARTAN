function varargout = autotrace(varargin)
% AUTOTRACE  Select traces according to user-defined criteria
%
%   Autotrace is a GUI that displays histograms of statistics calculated
%   for each trace in loaded file(s). These are calculated in traceStat.m.
%   The user can then select traces according to a defined set of criteria
%   and save that subse to a new file, typically ending in "_auto.traces".
%   A .log file is also saved that includes the criteria and other details.
%
%   "Batch Mode": For each .rawtraces file in the current directory and
%   sub-directories, load the traces, calculate statistics, select traces
%   according to the current criteria, and saved as an "_auto.traces" file.
%
%   NOTE that selection can bias the data. Always use the same selection
%   criteria when comparing datasets.
%
%   See also traceStat, pickTraces, loadPickSaveTraces.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Last Modified by GUIDE v2.5 23-Oct-2015 18:27:13

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
% --- Executes just before autotrace is made visible.
function autotrace_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to autotrace (see VARARGIN)


%---- PROGRAM CONSTANTS
constants = cascadeConstants();
handles.constants = constants;

set( handles.figure1, 'Name', ['autotrace (version ' constants.version ')'] );

%---- INITIAL VALUES FOR PICKING CRITERIA
criteria = constants.defaultAutotraceCriteria;
handles.criteria = criteria;

%---- OTHER USER-TUNABLE PARAMETERS
handles.nHistBins=40; % histogram bins size


%---- Initialize input fields with values defined above.
if isfield(criteria,'eq_overlap') && criteria.eq_overlap==0,
    set(handles.chk_overlap,'Value',1);
else
    set(handles.chk_overlap,'Value',0);
end

set( handles.ed_min_corr, 'String',num2str(criteria.min_corr)    );
set( handles.ed_max_corr, 'String',num2str(criteria.max_corr)    );
set( handles.ed_snr,      'String',num2str(criteria.min_snr)     );
set( handles.ed_bg,       'String',num2str(criteria.max_bg)      );
set( handles.ed_ncross,   'String',num2str(criteria.max_ncross)  );
set( handles.ed_acclife,  'String',num2str(criteria.min_acclife) );


%---- Setup drop-down boxes listing trace statistics.

% Get names of trace statistics
ln = traceStat;  %get long statistic names
handles.statLongNames = ln;
longNames  = struct2cell(ln);
shortNames = fieldnames(ln);

% Add trace statistic names to dropdown boxes above histogram axes.
handles.cboNames = strcat('cboStat',{'1','2','3','4','5'});

for id=1:length(handles.cboNames),
    % Also set options in combobox
    set( handles.(['cboStat' num2str(id)]), 'String', longNames );
end

% Set default selections for the drop-down boxes.
set( handles.cboStat1, 'Value', find(strcmp('t',shortNames))  );
set( handles.cboStat2, 'Value', find(strcmp('donorlife',shortNames))  );
set( handles.cboStat3, 'Value', find(strcmp('corr',shortNames))  );
set( handles.cboStat4, 'Value', find(strcmp('snr',shortNames))   );
set( handles.cboStat5, 'Value', find(strcmp('bg',shortNames))    );

% Setup special criteria selection drop-down boxes.
handles.nCriteriaBoxes = 7;
criteriaNames = [{''}; longNames];

for id=1:handles.nCriteriaBoxes
    set( handles.(['cboCriteria' num2str(id)]), 'String', criteriaNames );
end

% FIXME: in future, default criteria values may include some that will show
% up in the "special" boxes. Code is needed here to set these up.


%---- Add context menus to the plots to launch curve fitting or copy data.
% Ideally, you should be able to drop-in any variant of the interface,
% with variable numbers of histogram boxes, etc and have it still work.
menu = uicontextmenu;

% Context menu for launching Matlab's curve fitting tool
uimenu( menu, 'Label','Curve Fitting...', 'Callback',...
       'autotrace(''launchFitTool_Callback'',gca)' );

% Context menu for copying the raw statistic data
% for fitting in other programs (like Origin)
uimenu( menu, 'Label','Copy data', 'Callback',...
       'autotrace(''copyPlotData_Callback'',gca)' );

% Add context menu to all axes (histograms of trace statistics).
zoom off;
hZoom = zoom(handles.figure1);
set(hZoom, 'UIContextMenu', menu );
zoom on;

% Choose default command line output for autotrace
handles.output=hObject;

% Update handles structure
guidata(hObject,handles);

% END FUNCTION autotrace_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = autotrace_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

%--- Open file with user-interface.
filter = {'*.rawtraces;*.traces','All traces files (*.rawtraces,*.traces)'; ...
          '*.rawtraces','Raw traces files (*.rawtraces)'; ...
          '*.txt','Old format traces files (*.txt)'; ...
          '*.*','All files (*.*)'  };
      
[datafile,datapath] = uigetfile( filter,'Choose a traces file:', ...
                                                'MultiSelect','on');

if datapath==0, return; end %user hit cancel


%--- Convert filename list into a cell array
if ~iscell(datafile), datafile = {datafile}; end
filename = strcat(datapath,datafile);

%--- Save file list for later use.
handles.inputdir = datapath;
handles.inputfiles = filename;

%--- Update GUI listing of number of files and file types
disp('Loading files...');
for i=1:numel(handles.inputfiles),
    disp( handles.inputfiles{i} );
end

if numel(filename) == 1,
    fileDisplayText = handles.inputfiles{1};
else
    fileDisplayText = handles.inputdir;
end

% Shorten display name by removing initial part of the path if too long.
if length(fileDisplayText)>70,
    fileDisplayText = ['...' fileDisplayText(end-70:end)];
end

set(handles.editFilename,'String', fileDisplayText);

%--- Load the traces files.
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
traces_files = dir( [datapath filesep '*.rawtraces'] );

if numel(traces_files) == 0
    disp('No files in this directory!');
    return;
end

handles.inputdir = datapath;
handles.inputfiles = strcat( [datapath filesep], {traces_files.name} );

% Update the GUI with the new data location name.
disp(handles.inputdir);
set(handles.editFilename,'String',handles.inputdir);

% Load the traces files.
OpenTracesBatch( hObject, handles );

% END FUNCTION btnOpenDirectory_Callback





%----------BATCH ANALYSIS----------%
function btnBatchMode_Callback(hObject, ~, handles, datapath)
% For each .rawtraces file in the current directory, and optionally all
% sub-directories, select traces according to the current criteria and save
% a corresponding _auto.traces file.

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
% FIXME: for now keeping this not recursive by default.
trace_files = regexpdir(datapath,'^.*\.rawtraces$');

if numel(trace_files)==0,
    return;
end

trace_files = {trace_files.name};


% For each file in the user-selected directory
% FIXME: consider avoiding passing handles this way.
% wbh = waitbar(0,'Batch mode progress');

for i=1:numel(trace_files)
    
    % If running in the automatic mode, skip already-processed files.
    [p,f] = fileparts(trace_files{i});
    if auto,
        auto_name = fullfile(p,[f '_auto.traces']);
        if exist(auto_name,'file'), continue; end
    end
    
    handles.inputdir = p;
    handles.inputfiles = trace_files(i);
    set(handles.editFilename,'String',trace_files{i});

    % Load traces and calculate statistics.
    % (this calls PickTraces_Callback, which saves GUI data).
    handles = OpenTracesBatch( hObject, handles );
    
    %waitbar((i-0.5)/numel(trace_files),wbh);

    % Save picked traces to a new _auto.txt file.
    disp( handles.outfile );
    SaveTraces( handles.outfile, handles );

    %waitbar(i/numel(trace_files),wbh);
end

% close(wbh);

% END FUNCTION btnGo_Callback





function handles = OpenTracesBatch( hObject, handles )
% Calculate and display trace statistics for the current list of files.

set(handles.figure1, 'pointer', 'watch')
drawnow;

% Clear out old data to save memory.
if isappdata(handles.figure1,'infoStruct') %if data previously loaded.
    rmappdata(handles.figure1,'infoStruct');
end

% Determine default filename to use when saving.
[p,f] = fileparts( handles.inputfiles{1} );
handles.outfile = fullfile(p, [f '_auto.traces']);

% Calculate trace stats
[infoStruct,nTracesPerFile] = traceStat( handles.inputfiles, handles.constants );
handles.nTraces = numel( infoStruct );
handles.nTracesPerFile = nTracesPerFile;
                    
% Save the trace properties values to application data
setappdata(handles.figure1,'infoStruct', infoStruct);
clear infoStruct;

% Select traces according to the current criteria.
handles = PickTraces_Callback(hObject,handles);

set(handles.figure1, 'pointer', 'arrow');

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

set(handles.figure1, 'pointer', 'watch')
drawnow;

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
set(handles.SaveTraces,'Enable','off');


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
set(handles.btnSaveProperties,'Enable','off');
guidata(hObject,handles);

% END FUNCTION SaveProperties_Callback
      




% --------------- VIEW TRACES FUNCTIONALITY --------------- %

% --- Executes on button press in ViewPickedTraces.
function ViewPickedTraces_Callback(hObject, ~, handles)

% If not already, save the selected traces to file.
if strcmpi( get(handles.SaveTraces,'Enable'), 'on' );
    outfile = SaveTraces_Callback(hObject, [], handles);
else
    outfile = handles.outfile;
end

% Run sorttraces interface so traces can be viewed
% Could pass description/title as third param (currently empty!)
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
    if selection==1, continue; end %no selection
    
    equality = get(handles.(['cboEquality' num2str(id)]),'Value');
    if equality==1, continue; end %no inequality selected
    
    criteriaName = [equalityText{equality-1} shortNames{selection-1}];
    edText = get(handles.(['edCriteria' num2str(id)]),'String');
    if isempty(edText), continue; end %no criteria value.
    
    criteria.(criteriaName) = ...
        str2double( get(handles.(['edCriteria' num2str(id)]),'String') );
end


% Get criteria values for all the fixed GUI elements.
% chk_corr is listed twice because it applies to two distinct criteria.
% Overlap must be handled seperately since it doesn't have an associated textbox.
chkNames = {'chk_acclife','chk_corr','chk_corr','chk_snr','chk_bg','chk_ncross','chk_maxTotalSigma'};
txtNames = {'ed_acclife','ed_min_corr','ed_max_corr','ed_snr','ed_bg','ed_ncross','ed_maxTotalSigma'};
criteriaNames = {'min_acclife','min_corr','max_corr','min_snr','max_bg','max_ncross','maxTotalSigma'};

for i=1:numel(chkNames),
    if get(handles.(chkNames{i}),'Value'),
        criteria.(criteriaNames{i}) = str2double( get(handles.(txtNames{i}),'String') );
    end
end

% Handle overlap checkbox special case.
if get(handles.chk_overlap,'Value'),
    criteria.eq_overlap = 0;  %zero means remove overlapped traces.
end

handles.criteria = criteria;


% Find which molecules pass the selection criteria
stats = getappdata(handles.figure1,'infoStruct');
if isempty(stats), return; end %no data loaded, nothing to do.

picks = pickTraces( stats, criteria );

% The number of traces picked.
handles.inds_picked = picks;
handles.picked_mols = numel(handles.inds_picked);

% If at least one trace is picked, turn some buttons on.
if handles.picked_mols > 0
    set(handles.SaveTraces,'Enable','on');
    set(handles.MakeContourPlot,'Enable','on');
    set(handles.ViewPickedTraces,'Enable','on');
    set(handles.btnSaveProperties,'Enable','on');
end

% Turn some other buttons on/off.
set(handles.MoleculesPicked,'String', ...
            sprintf('%d of %d',[handles.picked_mols,handles.nTraces]));

        
% Update trace statistic histograms for traces passing selection criteria.
nAxes = length( handles.cboNames );

for i=1:nAxes,
    cboStat_Callback(handles.(handles.cboNames{i}), [], handles);
end


guidata(hObject,handles);

% END FUNCTION PickTraces_Callback




%#########################################################################
%---------------------- HISTOGRAMS & CONTOUR PLOTS ----------------------%
%#########################################################################

%----------MAKE CONTOUR PLOT----------
% --- Executes on button press in MakeContourPlot.
function MakeContourPlot_Callback(hObject, ~, handles)
% Display FRET contour plot of currently selected traces.

% If not already, save the currently selected traces to file.
if strcmpi( get(handles.SaveTraces,'Enable'), 'on' );
    outfile = SaveTraces_Callback(hObject, [], handles);
else
    outfile = handles.outfile;
end

assert( exist(outfile,'file')~=0 );

% Run makeplots to display FRET contour plot
if ~isempty(outfile),
    [~,title] = fileparts(outfile);
    title = strrep( title,'_',' ' );
    makeplots( outfile, title );
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
    bins = handles.nHistBins; %let hist choose N bin centers.
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
    ylabel( handles.axStat1, 'Number of Traces (%)' );
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
function chkAutoBatch_Callback(hObject, ~, ~)
% Creates a timer to look for and automatically process .rawtraces files.
 
% If another timer is running, stop it.
fileTimer = timerfind('Name','autotrace_fileTimer');

if ~isempty(fileTimer),
    stop(fileTimer);
    delete(fileTimer);
end


% Start a new timer if requested
if get(hObject,'Value') == get(hObject,'Max'),
    % Ask the user for a directory location
    targetDir = uigetdir('','Choose directory:');
    if targetDir==0,  %user hit cancel
        set( hObject, 'Value', get(hObject,'Min') );
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
end %if

% END FUNCTION chkAutoBatch_Callback



function stopFileTimer(~,~,hObject)
% This function is called when there is an error during the timer callback
% or when the timer is stopped.
handles = guidata(hObject);
set(handles.chkAutoBatch,'Value',0);

% END FUNCTION stopFileTimer


function updateFileTimer(~,~,hObject,targetDir)
% This function runs each time the timer is fired, looking for any new
% movies that may have appeared on the path.

handles = guidata(hObject);

% Kill the timer if the directory is inaccessible
if ~exist(targetDir,'dir'),
    disp('Autotrace: stopping batch mode: directory was moved or is inaccessible');
    set(handles.chkAutoBatch,'Value',0);
    
    fileTimer = timerfind('Name','autotrace_fileTimer');
    stop(fileTimer);
    delete(fileTimer);
end

btnBatchMode_Callback( handles.figure1, [], handles, targetDir );


% END FUNCTION updateFileTimer
