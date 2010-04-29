function varargout = realtimeAnalysis(varargin)
% REALTIMEANALYSIS  Trace processing and filtering
%


% Last Modified by GUIDE v2.5 08-Apr-2010 17:25:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @realtimeAnalysis_OpeningFcn, ...
    'gui_OutputFcn',  @realtimeAnalysis_OutputFcn, ...
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





%%
%#########################################################################
%------------------------- INITIALIZATION (GUI) -------------------------%
%#########################################################################


%----------INITIALIZATION OF THE GUI----------%
% --- Executes just before realtimeAnalysis is made visible.
function realtimeAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to realtimeAnalysis (see VARARGIN)

% Leave everything alone if the program is already running.
% This initialization proceedure will confuse the program state.
if isfield(handles,'criteria'),
    disp('realtimeAnalysis GUI is already running!');
    return;
end


%---- Generate settings dialog box, hide from view
handles.hSettings = realtimeAnalysis_settings( ...
                            @realtimeAnalysis_Notify, hObject );
% set(handles.hSettings,'Visible','off');
                        
% Load settings from the dialog...
settingsHandle = guidata( handles.hSettings );
handles.criteria = settingsHandle.criteria;
handles.gettracesOptions = settingsHandle.options;


%---- PROGRAM CONSTANTS
handles.constants = cascadeConstants();

handles.nHistBins=40; % histogram bins size

makeplotsOptions.contour_bin_size = 0.035;
makeplotsOptions.no_tdp = 1;
makeplotsOptions.pophist_offset = 0;
% makeplotsOptions.hideText = true;
handles.makeplotsOptions = makeplotsOptions;


%---- SETUP GUI TRACE STATISTIC CONTROLS

% Names of trace statistics -- these should be defined somewhere else!
ln = traceStat;  %get long statistic names
handles.statLongNames = ln;
longNames  = struct2cell(ln);
shortNames = fieldnames(ln);

% Add context menus to the plots to launch curve fitting
% Ideally, you should be able to drop-in any variant of the interface,
% with variable numbers of histogram boxes, etc and have it still work.
handles.cboNames = strcat('cboStat',{'1','2'});
handles.nPlots = length(handles.cboNames);

for id=1:handles.nPlots,
    set( handles.(['cboStat' num2str(id)]), 'String', longNames );
end

% Set default selections
set( handles.cboStat1, 'Value', find(strcmp('t',shortNames))  );
set( handles.cboStat2, 'Value', find(strcmp('acclife',shortNames))  );


%---- INITIALIZE DATA STORAGE VARIABLES
handles = resetGUI(handles);


%---- Other setup:
handles.inputdir = pwd;
set(handles.txtDirectory,'String',pwd);

% Choose default command line output for realtimeAnalysis
handles.output=hObject;

% Update handles structure
guidata(hObject,handles);

% END FUNCTION realtimeAnalysis_OpeningFcn



% Callback for when settings dialog is updated...
function realtimeAnalysis_Notify(result, hObject)

% Load data in both figures
handles = guidata(hObject);

% Update selection criteria
handles.criteria = result.criteria;
handles.gettracesOptions = result.gettracesOptions;

disp(handles.criteria);
disp(handles.gettracesOptions);

% Update handles structure
guidata(hObject,handles);

% Update analysis.
updateGUI( hObject, handles, 1 );



% --- Outputs from this function are returned to the command line.
function varargout = realtimeAnalysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=handles.output;

% END FUNCTION realtimeAnalysis_OutputFcn



%%
%#########################################################################
%----------------------- LOAD, FILTER, SAVE TRACES ----------------------%
%#########################################################################


%----------BATCH ANALYSIS----------%
function updateGUI( hObject, handles, force )
% This allows for realtimeAnalysis to run in batch mode. A directory is selected,
% and all the traces files are combined, and a single output file is
% generated. NOTE: be careful not to combine data from different
% experiments. Store different data sets in separate folders.

% Process input parameters
if nargin<3, force = 0; end

% Only run one analysis instance at once.
if isfield(handles,'isExecuting') && handles.isExecuting,
    %disp('Analysis is already executing!');
    return;
end

% Get user-selected directory name from GUI.
datapath = get(handles.txtDirectory,'String');
handles.inputdir = datapath;
if isempty(handles.inputdir), return; end

% Prevent GUI interactions from triggering this callback while running.
handles.isExecuting = 1;
guidata(hObject,handles);
needUpdate = 0; %true if GUI needs to be updated.


%---- Extract traces from all movies in current directory.
set( handles.txtStatus, 'String', 'Extracting traces from movies...' ); drawnow;
gettracesOptions = handles.gettracesOptions;
gettracesOptions.skipExisting = 1;
gettracesOptions.quiet = 1;
gettraces( datapath, gettracesOptions );


%---- Load traces files that have not already been loaded.
tracesLoaded = getappdata(handles.figure1,'tracesLoaded');
stats = getappdata(handles.figure1,'infoStruct');

% Create list of .traces files in the selected directory.
tracesToLoad = dir( [datapath filesep '*.traces'] );
inputfiles = strcat( [datapath filesep], {tracesToLoad.name} );


% If a file (and the corresponding .stk) were removed, trace data is removed
% from the cache.
[filesToKeep,idxToKeep] = intersect( handles.filesLoaded, inputfiles );

if numel(idxToKeep)<numel(handles.filesLoaded),
    idxTraces = handles.idxTraces(idxToKeep,:);
    tracesToKeep = mergeIndexes(idxTraces);
    
    tracesLoaded.d   = tracesLoaded.d(tracesToKeep,:);
    tracesLoaded.a   = tracesLoaded.a(tracesToKeep,:);
    tracesLoaded.f   = tracesLoaded.f(tracesToKeep,:);
    tracesLoaded.ids = tracesLoaded.ids(tracesToKeep);
    setappdata(handles.figure1,'tracesLoaded',tracesLoaded);
    
    stats = stats(tracesToKeep);
    setappdata(handles.figure1,'infoStruct',stats);
    
    handles.nTraces = size( tracesLoaded.d,1 );
    handles.idxTraces = squashIndexes( idxTraces );
    handles.filesLoaded = filesToKeep;
    
    needUpdate = 1;
end


% Load new files into the cache.
filesToLoad = setdiff( inputfiles, handles.filesLoaded );

if ~isempty(filesToLoad)
    % Load trace data from new files
    set( handles.txtStatus, 'String', 'Loading traces files...' ); drawnow;
    [traceData,indexes] = loadTracesBatch( filesToLoad );
    if isempty(handles.timeAxis),
        handles.timeAxis = traceData.time;
    end
    handles.filesLoaded = [handles.filesLoaded filesToLoad];

    % Merge new data into existing cache
    handles.idxTraces = [handles.idxTraces ; indexes+handles.nTraces];
    
    tracesLoaded(1).d   = [tracesLoaded.d ;  traceData.d  ];
    tracesLoaded.a   = [tracesLoaded.a ;  traceData.a  ];
    tracesLoaded.f   = [tracesLoaded.f ;  traceData.f  ];
    tracesLoaded.ids = [tracesLoaded.ids traceData.ids];
    setappdata(handles.figure1,'tracesLoaded',tracesLoaded);
    clear traceData;

    handles.nTraces = size( tracesLoaded.d,1 );
    handles.nFiles  = numel( inputfiles );
    
    % Calculate trace statistics for selection.
    % TODO: only run this for the subset of data being loaded, then
    %       merge with existing (cached) results!!!
    stats = traceStat(tracesLoaded.d,tracesLoaded.a,tracesLoaded.f);
    setappdata(handles.figure1,'infoStruct',stats);
    
    needUpdate = 1;
end

% If there are no new/lost files, no need to update (unless settings changed).
if ~needUpdate && ~force,
    set( handles.txtStatus, 'String','Finished.' );
    handles.isExecuting = 0;
    guidata(hObject,handles);
    return;
end


%---- Select traces based on user-defined criteria.
set( handles.txtStatus, 'String', 'Selecting traces with user-defined criteria...' );

picks = pickTraces( stats, handles.criteria );
handles.inds_picked = picks; %just picked traces.
handles.picked_mols = numel(picks);


%---- Update contour plots
set( handles.txtStatus, 'String','Updating plots...' ); drawnow;

options = handles.makeplotsOptions;
options.targetAxes = { handles.axFretContourAll };
frethist = makecplot( tracesLoaded.f(picks,:), options );
makeplots( frethist, 'All movies', options );

% Also show only traces from the last movie.
selectedPicks = picks( picks>=handles.idxTraces(end,1) & picks<=handles.idxTraces(end,2) );
options.targetAxes = { handles.axFretContourSelected };
frethist = makecplot( tracesLoaded.f(selectedPicks,:), options );
makeplots( frethist, 'Last movie', options );




%---- Show population statistics in GUI
% All traces
handles = updateStats( stats, handles );

% Just traces from the last (most recent) movie.
lastMovie = 1:handles.nTraces;
lastMovie = lastMovie>=handles.idxTraces(end,1) & lastMovie<=handles.idxTraces(end,2);
handles = updateStats( stats(lastMovie), handles, 'Last' );
set( handles.edAcceptanceLast,'String','' ); %not working yet.


% Select traces with quantifiable stats.
% selectedIdx = idx( idx>=handles.idxTraces(end,1) & idx<=handles.idxTraces(end,2) );
% nSelected = diff( handles.idxTraces(end,:) );



%---- Update GUI
set( handles.txtStatus, 'String','Finished.' ); drawnow;

handles.isExecuting = 0;
guidata(hObject,handles);



% Update trace statistic histograms for traces passing selection criteria.
nAxes = length( handles.cboNames );

for i=1:nAxes,
    cboStat_Callback(handles.(handles.cboNames{i}), handles);
end


% END FUNCTION OpenTracesBatch




function handles = updateStats( stats, handles, postfix )

if nargin<3, postfix=''; end

dt = diff(handles.timeAxis(1:2));

%---- Show population statistics in GUI

clear basicCriteria;
basicCriteria.min_snr = 2;
basicCriteria.min_lifetime = 10;
basicCriteria.overlap = 1;
idx = pickTraces( stats, basicCriteria );

t = [stats.t];
snr   = [stats.snr];
snr_s = [stats.snr_s];
acceptorLifetime = [stats.acclife];
donorLifetime = [stats.lifetime];

% Calculate breakdown of traces removed during selection.
nTraces = numel(stats); %all data.
nSingle = sum( snr>0 & ~[stats.overlap] );
nFRET   = sum( acceptorLifetime>1 );
nPicked = handles.picked_mols;

% Calculate important trace statistics.
avgI    = median( t(idx) );
avgSNR  = median( snr(idx)   );
avgSNRs = median( snr_s(idx) );

% Calculate donor bleaching rate...
[donorDist,donorAxes] = hist( donorLifetime(idx), 40 );
donorDist = 1 - [0 cumsum(donorDist)]/sum(donorDist);
donorAxes = [0 donorAxes];
f = fit( donorAxes',donorDist','exp1' );
% figure; plot( f ); hold on; bar( donorAxes',donorDist' );

if handles.timeAxis(1)==0
    donorT = sprintf('%.1f sec',-1/f.b/dt);
else
    donorT = sprintf('%.1f frames',-1/f.b);
end

% Calculate acceptor bleaching rate
acclifeCriteria = basicCriteria;
acclifeCriteria.min_acclife = 1;
idxAcc = pickTraces( stats, acclifeCriteria );

[accDist,accAxes] = hist( acceptorLifetime(idxAcc), 40 );
accDist = 1 - [0 cumsum(accDist)]/sum(accDist);
accAxes = [0 accAxes];
f = fit( accAxes',accDist','exp1' );

if handles.timeAxis(1)==0
    acceptorT = sprintf('%.1f sec',-1/f.b/dt);
else
    acceptorT = sprintf('%.1f frames',-1/f.b);
end


% Create data table for displaying states in GUI
tableData = { sprintf('%.0f%% (%d)',100*nSingle/nTraces,nSingle); ...
              sprintf('%.0f%% (%d)',100*nFRET/nTraces,  nFRET  ); ...
              sprintf('%.0f%% (%d)',100*nPicked/nTraces,nPicked); ...
              %' '; ... %spacer row
              sprintf('%.1f',        avgI           ); ...
              sprintf('%.1f (%.1f)', avgSNR,avgSNRs ); ...
              donorT ; acceptorT  };
          
% set( handles.tblStats, 'Data',tableData );

textboxNames = {'edSingleDonor','edHasFret','edAcceptance',...
                'edIntensity', 'edSNR', 'edLTDonor','edLTAcceptor'};
setString( handles, strcat(textboxNames,postfix), tableData );

    
% END FUNCTION updateStats




function setString( handles, fields, values )

for i=1:numel(fields),
    set( handles.(fields{i}),'String',values{i} );
end

% END FUNCTION setString








%%
%#########################################################################
%------------------------- INPUT FIELD HANDLING -------------------------%
%#########################################################################


% --- Executes on button press in btnBrowse.
function btnBrowse_Callback(hObject, eventdata, handles)
% CALLED: when the user clicked "Browse..."
% ACTION: Get location of data to process
 
datadir = uigetdir(handles.inputdir, 'Select a directory with all data to process');

if datadir~=0,
    % If this is a new location, reset the GUI.
    if ~strcmp(datadir,handles.inputdir),
        handles = resetGUI( handles );
    end
    
    % Update GUI with new data location
    handles.inputdir = datadir;
    set(handles.txtDirectory,'String',datadir);
end

% Update handles structure
guidata(hObject,handles);


function txtDirectory_Callback(hObject, eventdata, handles)
% CALLED: when the user modifies the data location textbox
% ACTION: reset GUI and prepare to load data from the new location.

datadir = get(hObject,'String');

if ~exist(datadir,'dir'),
    warning('realtimeAnalysis: directory doesn''t exist!');
    set(hObject,'String',handles.inputdir);
else
    % If this is a new location, reset the GUI.
    if ~strcmp(datadir,handles.inputdir),
        handles.inputdir = datadir;
        handles = resetGUI( handles );
    end
end

% Update handles structure
guidata(hObject,handles);


function handles = resetGUI( handles )
% Resets all GUI elements and clears background variables.

% Clear trace stat histograms
for id=1:handles.nPlots,
    cla( handles.(['axStat' num2str(id)]) );
end

% Clear contour plots
cla( handles.axFretContourSelected  );
cla( handles.axFretContourAll  );

% Clear trace stat table
% Data = cell(8,2);
% [Data{:}] = deal(' ');
% set( handles.tblStats, 'Data',Data );

% Delete saved trace data and reset
handles.filesLoaded = {};
handles.idxTraces = zeros(0,2);
handles.nFiles  = 0;
handles.nTraces = 0;
handles.inds_picked = [];
handles.picked_mols = 0;

tracesLoaded = struct( 'd',[], 'a',[], 'f',[], 'ids',{} );
handles.timeAxis = [];
setappdata(handles.figure1,'tracesLoaded',tracesLoaded);
setappdata(handles.figure1,'infoStruct',[]);

% END FUNCTION resetGUI



% --- Executes on button press in btnSettings.
function btnSettings_Callback(hObject, eventdata, handles)
% If the realtime analysis window has not been launched, do so now.
% If it has been launched, simply make it visible.
set( handles.hSettings, 'Visible','on' );




function close_Callback(hObject)
% Set figure as hidden (not Visible) so its properties can be
% accessed without the window being in the way.
handles = guidata(hObject);

% Also kill the settings dialog and timer. We don't want to hanging around
% after the main app has closed!
if isfield(handles,'fileTimer')
    disp('Timer deleted.');
    delete(handles.fileTimer);
end
delete(handles.hSettings);
delete( hObject );




% --- Executes on button press in chkAutoUpdate.
function btnGo_Callback(hObject, eventdata, handles)

handles.autoUpdate = get(hObject,'Value');

% If turned on, create a timer object to check the current directory
if handles.autoUpdate    
    disp('Timer started.');
    handles.fileTimer = timer('ExecutionMode','fixedSpacing','StartDelay',1,...
                              'TimerFcn',{@timer_Callback,handles.figure1},...
                              'Period',3.0,'BusyMode','drop');
    start(handles.fileTimer);

% If turned off, disable the current timer
else
    
    if isfield(handles,'fileTimer')
        disp('Timer deleted.');
        stop(handles.fileTimer);
        delete(handles.fileTimer);
    end
    
    % Reset the isExecuting flag. Need to have some way of reseting this
    % variable in case of a crash...
    set( handles.txtStatus, 'String','IDLE' );
    handles.isExecuting = 0;

end

guidata(hObject,handles);



function timer_Callback(timerObject,eventdata,hObject)

% disp('Timer tick');
updateGUI( hObject, guidata(hObject) )

% END FUNCTION timer_Callback






%%
%#########################################################################
%------------------------ TRACE STATISTIC DISPLAY -----------------------%
%#########################################################################


% --- Executes on selection change in any of the drop down boxes above each
%     of the trace statistic histogram plots.
function cboStat_Callback(hObject, handles)

% Get ID of this combo control and the handle of the corresponding axes.
id = get(hObject,'UserData');
hAx = handles.(['axStat' num2str(id)]);

% Get trace statistics
stats = getappdata(handles.figure1,'infoStruct');
% picks = find([stats.snr]>0);
picks = handles.inds_picked;

if isempty(picks),  cla(hAx); return;  end


% Get user selection
selected   = get(hObject,'Value');
statNames = fieldnames(handles.statLongNames);
statToPlot = statNames{selected};

assert( ismember( statToPlot, fieldnames(stats) ), ...
        'Selected trace statistic is unknown' );

% Also show only traces from the last movie.
selectedPicks = picks( picks>=handles.idxTraces(end,1) & picks<=handles.idxTraces(end,2) );
otherPicks = setdiff( picks, selectedPicks );

% Plot the distribution of the statistic
statData = [stats.(statToPlot)];

[~,binCenters] = hist( statData(picks), handles.nHistBins);

[dataOther,binCenters] = hist( statData(otherPicks), binCenters);
dataOther = 100*dataOther/numel(statData);  %normalize the histograms

[dataSelected,binCenters] = hist( statData(selectedPicks), binCenters);
dataSelected = 100*dataSelected/numel(statData); %normalize the histograms

dataOther = reshape( dataOther, size(dataSelected) );


% Plot histograms
hBar = bar( hAx, binCenters, [dataSelected' dataOther'], 1, 'stacked' );
% zoom on;
% grid on;
set(hBar(1),'facecolor',[1 0 0]);
set(hBar(2),'facecolor',[251 251 196]/256);

% if id==1,
%     ylabel( handles.axStat1, 'Number of Traces (%)' );
%     legend( {'Last movie','All movies'} );
% end





% --- Executes on button press in btnSaveTraces.
function btnSaveTraces_Callback(hObject, eventdata, handles)

% Prevent callbacks from running.
handles.isExecuting = 1;
guidata(hObject,handles);

% Create a name for the output file
handles.outfile = strrep(handles.filesLoaded{1}, '.traces', '_auto.txt');
handles.outfile = strrep(handles.outfile, '_01_auto.txt', '_auto.txt');

[inputfile inputpath]=...
    uiputfile('.txt','Save picked traces as:',handles.outfile);

if inputfile~=0,
    handles.outfile=[inputpath inputfile];
    [p,n,ext] = fileparts( handles.outfile );
    assert( strcmp(ext,'.txt'), 'must be .txt' );


    %---- Save selected traces files.
    tracesLoaded = getappdata(handles.figure1,'tracesLoaded');
    picks = handles.inds_picked; %just picked traces.

    saveTraces( handles.outfile, 'txt', tracesLoaded.d(picks,:), ...
                tracesLoaded.a(picks,:), tracesLoaded.f(picks,:), ...
                tracesLoaded.ids(picks), handles.timeAxis );

    qubFile = strrep(handles.outfile,'.txt','.qub.txt');
    saveTraces( qubFile, 'qub', tracesLoaded.f(picks,:) );
end
        
% Clean up
handles.isExecuting = 0;
guidata(hObject,handles);
        
        
