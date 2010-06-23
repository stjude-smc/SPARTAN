function varargout = autotrace(varargin)
% AUTOTRACE  Trace processing and filtering
%
%   Loads fluorescence traces files (produced using gettraces), makes
%   corrections for crosstalk, background intensity, and sets FRET=0
%   where the donor dye is dark.  Descriptive properties of each traces
%   is then calculated (SNR, D/A correlation, etc).  Using defined
%   criteria, the user can then select a portion of the dataset.
%
%   Use "Batch Mode" to load all traces within a directory.
%   Be careful not to mix multiple experiments in the same directory!
%
%   Use "Save Traces" to save the resulting corrected+filtered traces.
%   This will also produce a log file with useful information.
%
%   The typical criteria used are:
%     FRET-lifetime    > 15 frames       (trace must show FRET)
%     D/A correlation  < 0.5             (remove aggregates)
%     Signal-to-Noise  > 8               (sufficient resolution)
%     Background noise < 1500            (background drift)
%     N. Donor Blinks  < 3               (remove aggregates)
%     Remove traces with multiple dyes = YES
%
%   NOTE that many of these criteria can bias the data, especially
%   correlation.  When comparing datasets, use the same criteria.

% Created by: James Munro, Daniel Terry (Scott Blanchard Lab)
% Cascade smFRET Analysis Pipeline 1.3, Copyright (C) 2008 Scott Blanchard
% Date Created: Oct 11, 2007

%   A gamma (sensitivity/quantum yield ratio) correction is used in
%   calculating total intensity and SNR.  The value comes from
%   cascadeConstants and was calculated for our equipment with ribosome
%   samples.  The value will vary based on equipment and sample studied.
%   
%   Fluorescence data are no longer stored in handles because loading
%   too many traces at once cause out-of-memory errors.  They are now
%   oaded on-demand through GetPickedTraces.  Post-synchronization of traces
%   is no longer implemented!
%   8/2007  -DT
%
%   FRET lifetime, N. Donor Blinks, signal overlap detection criteria all
%   added from original version by JBM.
%   4/2008  -DT


% Last Modified by GUIDE v2.5 26-Mar-2010 16:16:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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
function autotrace_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to autotrace (see VARARGIN)

% Leave everything alone if the program is already running.
% This initialization proceedure will confuse the program state.
if isfield(handles,'criteria'),
    disp('Autotrace is already running!');
    return;
end

%---- PROGRAM CONSTANTS
constants = cascadeConstants();
handles.constants = constants;


%---- INITIAL VALUES FOR PICKING CRITERIA
criteria.overlap = 1; % Remove overlapping molecules

% These can be changed while running the program.
criteria.min_corr=-1.1;    % 
criteria.max_corr=0.5;     % D/A correlation < 0.5
criteria.min_snr=8;       % SNR over background
criteria.max_bg=1500;      % Background noise
criteria.max_ncross = 4;   % donor blinking events
criteria.min_acclife = 15; % FRET lifetime

handles.criteria = criteria;


%---- OTHER USER-TUNABLE PARAMETERS
handles.nHistBins=40; % histogram bins size
handles.contour_bin_size=0.035;
handles.sync='n';



%---- Initialize input fields with values defined above.
set(handles.lowThresh,'String',num2str(criteria.min_corr));
set(handles.highThresh,'String',num2str(criteria.max_corr));
set(handles.SignalNoiseThresh,'String',num2str(criteria.min_snr));
set(handles.BackgroundNoiseThresh,'String',num2str(criteria.max_bg));
set(handles.editNCross,'String',num2str(criteria.max_ncross));
set(handles.editAccLife,'String',num2str(criteria.min_acclife));

set(handles.FRETBinSize,'String',num2str(handles.contour_bin_size));


%---- Setup drop-down boxes listing trace statistics.

% Get names of trace statistics
ln = traceStat;  %get long statistic names
handles.statLongNames = ln;
longNames  = struct2cell(ln);
shortNames = fieldnames(ln);

% Add trace statistic names to dropdown boxes above histogram axes.
handles.cboNames = strcat('cboStat',{'1','2','3','4','5'});
handles.nPlots = length(handles.cboNames);

for id=1:handles.nPlots,
    % Also set options in combobox
    set( handles.(['cboStat' num2str(id)]), 'String', longNames );
end

% Set default selections for the drop-down boxes.
set( handles.cboStat1, 'Value', find(strcmp('t',shortNames))  );
set( handles.cboStat2, 'Value', find(strcmp('lifetime',shortNames))  );
set( handles.cboStat3, 'Value', find(strcmp('corr',shortNames))  );
set( handles.cboStat4, 'Value', find(strcmp('snr',shortNames))   );
set( handles.cboStat5, 'Value', find(strcmp('bg',shortNames))    );

% Setup special criteria selection drop-down boxes.
handles.nCriteriaBoxes = 7;
criteriaNames = [{''}; longNames];

for id=1:handles.nCriteriaBoxes
    set( handles.(['cboCriteria' num2str(id)]), 'String', criteriaNames );
end


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



handles.isBatchMode = 0;


%warning off MATLAB:divideByZero

% Choose default command line output for autotrace
handles.output=hObject;

% Update handles structure
guidata(hObject,handles);

% END FUNCTION autotrace_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = autotrace_OutputFcn(hObject, eventdata, handles)
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
% This method is called when the user clicks the "Open Traces File.."
% button in the GUI.  The trace is parsed,
% --- Executes on button press in OpenTracesFile.
function OpenTracesFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenTracesFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Open file with user-interface.
filter = {'*.traces*;*.traces.old;*.txt','All traces files (*.txt, *.traces)'; ...
          '*.txt','Autotrace Files (*.txt)'; ...
          '*.traces*','Raw traces files (*.traces)'};
      
[datafile,datapath] = uigetfile( filter,'Choose a traces file:', ...
                                                'MultiSelect','on');

if datapath==0, return; end %user hit cancel


% Convert filename list into a cell array
if ~iscell(datafile), datafile = {datafile}; end
filename = strcat(datapath,datafile);

% Save file list for later use.
handles.inputdir = datapath;
handles.inputfiles = filename;

% Update GUI listing of number of files and file types
handles.inputfiles{:}
if numel(filename) == 1,
    set(handles.editFilename,'String',handles.inputfiles{1});
else
    set(handles.editFilename,'String',handles.inputdir);
end

% Load the traces files.
OpenTracesBatch( hObject, handles );


% END FUNCTION OpenTracesFile_Callback





%----------OPEN ALL TRACES FILES WITHIN THE SAME DIRECTORY----------%
% --- Executes on button press in btnOpenDirectory.
function btnOpenDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpenDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This allows for autotrace to run in batch mode. A directory is selected,
% and all the traces files are combined, and a single output file is
% generated. NOTE: be careful not to combine data from different
% experiments. Store different data sets in separate folders.

% Select directory by user interface.
datapath=uigetdir;
if datapath==0, return; end %user hit cancel.

% Create list of .traces files in the directory.
traces_files = dir( [datapath filesep '*.traces'] );

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
OpenTracesBatch( hObject, handles )

% END FUNCTION btnOpenDirectory_Callback





%----------BATCH ANALYSIS----------%
% --- Executes on button press in btnGo.
function btnBatchMode_Callback(hObject, eventdata, handles)
% hObject    handle to btnGo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This allows for autotrace2 to run in batch mode. A directory is selected,
% and all the traces files are combined, and a single output file is
% generated. NOTE: be careful not to combine data from different
% experiments. Store different data sets in separate folders.

% Select directory by user interface.
datapath=uigetdir;
if datapath==0, return; end

handles.isBatchMode = 1;


% Get a list of all traces files under the current directory
trace_files  = rdir([datapath filesep '**' filesep '*.traces']);

% Pool these files into
data_dirs = {};

for i=1:numel(trace_files),
    % Extract path of the file
    f = fileparts(trace_files(i).name);
    
    % If not already in the list, insert it
    nMatches = sum( cellfun( @(arg)strcmp(f,arg), data_dirs ) );
    if nMatches==0,
        data_dirs{end+1} = f;
    end
end


% For each file in the user-selected directory
wb = waitbar(0,'Processing data....');

for i=1:numel(data_dirs)
    
    datapath = data_dirs{i};
    
    % Create list of .traces files in the directory.
    traces_files = dir( [datapath filesep '*.traces'] );
    handles.nFiles = numel(traces_files);

    if handles.nFiles == 0
        disp('No files in this directory!');
        return;
    end

    handles.inputdir = datapath;
    handles.inputfiles = strcat( [datapath filesep], {traces_files.name} );

    % Update GUI listing of number of files and file types
    disp(handles.inputdir);
    set(handles.editFilename,'String',handles.inputdir);

    % Load all traces in the current directory
    OpenTracesBatch( hObject, handles )
    handles = guidata(hObject);
    
    waitbar((i-0.5)/numel(data_dirs),wb);

    % Save picked data to handles.outfile
    SaveTraces( handles.outfile, handles );
    handles = guidata(hObject);

    waitbar(i/numel(data_dirs),wb);
end
close(wb);

handles.isBatchMode = 0;
guidata(hObject,handles);

% END FUNCTION btnGo_Callback










function OpenTracesBatch( hObject, handles )

% Clear out old data to save memory.
if isappdata(handles.figure1,'infoStruct') %if data previously loaded.
    rmappdata(handles.figure1,'infoStruct');
    rmappdata(handles.figure1,'tracedata');
    rmappdata(handles.figure1,'ids');
end

% Determine default filename to use when saving.
handles.outfile = strrep(handles.inputfiles{1}, '.traces', '_auto.txt');
handles.outfile = strrep(handles.outfile, '_01_auto.txt', '_auto.txt');

% Load all data into one large dataset.
[tracedata,ids] = loadTracesBatch( handles.inputfiles );
handles.timeAxis = tracedata.time;
[handles.nTraces,handles.len] = size( tracedata.d );

% Calculate trace stats
infoStruct = traceStat( tracedata.d, tracedata.a, tracedata.f, ...
                        handles.constants );

                    
% Save the trace properties values to application data
setappdata(handles.figure1,'infoStruct', infoStruct);
setappdata(handles.figure1,'tracedata', tracedata);
setappdata(handles.figure1,'ids', ids);
clear infoStruct; clear tracedata;


% Initialize a variable for storing the number of molecules picked.
handles.picked_mols=0;

guidata(hObject,handles);

% Automatically run Pick Traces
PickTraces_Callback(hObject,handles);


% END FUNCTION OpenTracesBatch



%---------------  SAVE PICKED TRACES TO FILE (CALLBACK) ---------------%
% --- Executes on button press in SaveTraces.
function SaveTraces_Callback(hObject, eventdata, handles)
% hObject    handle to SaveTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Create a name for the output file
[inputfile inputpath]=...
    uiputfile('.txt','Save picked traces as:',handles.outfile);
if inputfile==0, return; end

handles.outfile=[inputpath inputfile];
[p,n,ext] = fileparts( handles.outfile );
assert( strcmp(ext,'.txt'), 'must be .txt' );

% Save picked data to handles.outfile
SaveTraces( handles.outfile, handles );

% Update GUI controls:
% Disabling the button as a means of confirming operation success.
set(handles.SaveTraces,'Enable','off');
guidata(hObject,handles);




%--------------------  SAVE PICKED TRACES TO FILE --------------------%
function SaveTraces( filename, handles, txtOnly )

% Get trace data
data = getappdata(handles.figure1,'tracedata');
ids = getappdata(handles.figure1,'ids');

%---- Save only selected traces to disk.
inds = handles.inds_picked;
fret = data.f(inds,:);

saveTraces( filename, 'txt', data.d(inds,:), data.a(inds,:), fret, ...
                             ids(inds), handles.timeAxis );

if nargin<3 || ~txtOnly,
    qub_fname = strrep( filename, '.txt', '.qub.txt' );
    saveTraces( qub_fname, 'qub', fret );
end
clear data;


%---- Generate log file containing informtion about how the traces were picked.
logfile=strrep(filename,'.txt','.log');
fid=fopen(logfile,'w');

fprintf(fid,'%s\n\n%s\n',date,'DIRECTORY');
fprintf(fid,'  %s\n\n%s\n',handles.inputdir,'FILES');

if ~iscell(handles.inputfiles)
    fprintf(fid,'%s\n',handles.inputfiles);
else
    for k=1:numel(handles.inputfiles)
        [p,fname,ext] = fileparts( handles.inputfiles{k} );
        fprintf(fid,'  %s%s\n',  fname,ext);
    end
end


fprintf(fid,'\nMolecules Picked:\t%d of %d (%.1f%%)\n\n\n', ...
            handles.picked_mols, handles.nTraces, ...
            100*handles.picked_mols/handles.nTraces );  

        
% Descriptive statistics about dataset
stats = getappdata(handles.figure1,'infoStruct');

total = handles.nTraces;
isMolecule      = sum( [stats.snr]>0 );
singleMolecule  = sum( [stats.snr]>0 & [stats.overlap]==0 );
hasFRET         = sum( [stats.snr]>0 & [stats.overlap]==0 & [stats.acclife]>=5 );
other           = handles.picked_mols;

fprintf(fid,'PICKING RESULTS\n');
fprintf(fid, '  %20s:  %-5d (%.1f%%)\n', 'Donor photobleaches', isMolecule,     100*isMolecule/total);
fprintf(fid, '  %20s:  %-5d (%.1f%%)\n', 'Single donor',        singleMolecule, 100*singleMolecule/isMolecule);
fprintf(fid, '  %20s:  %-5d (%.1f%%)\n', 'Have FRET',           hasFRET, 100*hasFRET/singleMolecule);
fprintf(fid, '  %20s:  %-5d (%.1f%%)\n', 'Pass other criteria', other,   100*other/hasFRET);
fprintf(fid, '\n\n');


% Save picking criteria used
% NOTE: this does not include the specialized criteria!!!!
fprintf(fid,'PICKING CRITERIA\n');

criteria = getSpecialCriteria( handles );
names = fieldnames(  criteria );
vals  = struct2cell( criteria );

for i=1:numel(names),
    if isempty( vals{i} ), continue; end  %skip unchecked criteria
    fprintf(fid, '  %22s:  %.2f\n', names{i}, vals{i});
end

% Save values of all other constants used
fprintf(fid, '\n\nCONSTANTS\n');

names = fieldnames(  handles.constants );
vals  = struct2cell( handles.constants );

for i=1:numel(names),
    if isstruct( vals{i} ) || numel( vals{i} )>1, continue; end
    
    fprintf(fid, '  %22s:  %.2f\n', names{i}, vals{i});
end


fprintf(fid,'\n\n');
fclose(fid);


% END FUNCTION SaveTraces_Callback



%----------------  SAVE MOLECULE PROPERTIES TO FILE ----------------%

% --- Executes on button press in btnSaveProperties.
function btnSaveProperties_Callback(hObject, eventdata, handles)

% Retrieve trace statistics, extract to matrix
stats = getappdata(handles.figure1,'infoStruct');
stats = stats(handles.inds_picked);

data = cellfun( @double, struct2cell(stats) );
data = squeeze(data);
names = fieldnames(stats);

% Write column headers
filename = strrep( handles.outfile, '.txt', '_prop.txt' );
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
function ViewPickedTraces_Callback(hObject, eventdata, handles)

[inputfile inputpath]=...
    uiputfile('.txt','Save picked traces as:',handles.outfile);
if inputfile==0, return; end

handles.outfile=[inputpath inputfile];
guidata(hObject,handles);


% Save picked data to handles.outfile
SaveTraces( handles.outfile, handles, 1 );

% Update GUI controls:
% Disabling the button as a means of confirming operation success.
set(handles.SaveTraces,'Enable','off');
set(handles.ViewPickedTraces,'Enable','off');
guidata(hObject,handles);

% Run sorttraces interface so traces can be viewed
% Could pass description/title as third param (currently empty!)
sorttraces(0, handles.outfile);





function criteria = getSpecialCriteria( handles )
% Generate a selection criteria structure (see pickTraces.m) using the
% "Specialized Selection Criteria" drop-down boxes.

criteria = handles.criteria;

% Update criteria for combo-box selections
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



%----------APPLIES PICKING CRITERIA TO TRACES----------
% --- Executes on button press in PickTraces.
function PickTraces_Callback(hObject, handles)

criteria = getSpecialCriteria( handles );

% Find which molecules pass the selection criteria
stats = getappdata(handles.figure1,'infoStruct');
if isempty(stats), return; end %no data loaded, nothing to do.

picks = pickTraces( stats, criteria );
clear stats;

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
set(handles.SaveContourPlot,'Enable','off');

% Save data in handles object.
guidata(hObject,handles);


% Update trace statistic histograms for traces passing selection criteria.
nAxes = length( handles.cboNames );

for i=1:nAxes,
    cboStat_Callback(handles.(handles.cboNames{i}), [], handles);
end

% END FUNCTION PickTraces_Callback




%#########################################################################
%---------------------- HISTOGRAMS & CONTOUR PLOTS ----------------------%
%#########################################################################

%----------MAKE CONTOUR PLOT----------
% --- Executes on button press in MakeContourPlot.
function MakeContourPlot_Callback(hObject, eventdata, handles)
% Builds and displays contour plot.

% Get FRET data for selected traces.
inds = handles.inds_picked;
data = getappdata(handles.figure1,'tracedata');

% Make contour plots using makeplots.m
options.contour_bin_size = handles.contour_bin_size;
options.pophist_offset = 0;
frethist = makecplot( data.f(inds,:), options );
clear data;

[p,title] = fileparts(handles.outfile);
title = strrep( title,'_',' ' );
makeplots( frethist, title );


% Clean up and save FRET histogram data for saving.
set(handles.SaveContourPlot,'Enable','on');
guidata(hObject,handles);


% END FUNCTION MakeContourPlot_Callback

 
%----------SAVE CONTOUR PLOT----------
% --- Executes on button press in SaveContourPlot.
function SaveContourPlot_Callback(hObject, eventdata, handles)

% Get FRET data for selected traces.
inds = handles.inds_picked;
data = getappdata(handles.figure1,'tracedata');

% Make contour plots using makeplots.m
options.contour_bin_size = handles.contour_bin_size;
options.pophist_offset = 0;
frethist = makecplot( data.f(inds,:), options );
clear data;

% Write original file
histfile=strrep(handles.outfile,'.txt','_hist.txt');
dlmwrite(histfile,frethist,' ');

% GUI stuff
set(hObject,'Enable','off');
guidata(hObject,handles);



%#########################################################################
%------------------------- INPUT FIELD HANDLING -------------------------%
%#########################################################################



%----------MENU OF CONTOUR PLOT SETTINGS----------
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
opt=get(hObject,'Value');
choices = 'nsyi';
handles.sync = choices(opt);
guidata(hObject,handles);


function FRETBinSize_Callback(hObject, eventdata, handles)
handles.contour_bin_size=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes on button press in CorrelationCoefficientBox.
function CorrelationCoefficientBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.min_corr=str2double(get(handles.lowThresh,'String'));
    handles.criteria.max_corr=str2double(get(handles.highThresh,'String'));
    set(handles.lowThresh,'Enable','on');
    set(handles.highThresh,'Enable','on');
else
    handles.criteria.min_corr=[];
    handles.criteria.max_corr=[];
    set(handles.lowThresh,'Enable','off');
    set(handles.highThresh,'Enable','off');
end
guidata(hObject,handles);
PickTraces_Callback(hObject,handles);


% --- CALLED when "Remove overlapping traces" checkbox setting is changed.
% Updates the selection criteria automatically, w/o running PickTraces().
function chkOverlap_Callback(hObject, eventdata, handles)
handles.criteria.overlap = get(hObject,'Value');
guidata(hObject,handles);
PickTraces_Callback(hObject,handles);


% --- CALLED when any of the primary criteria textboxes is changed.
% Updates the selection criteria automatically, w/o running PickTraces().
function updateCriteria_Callback( hObject, handles, criteriaName )
handles.criteria.(criteriaName) = str2double(get(hObject,'String'));
guidata(hObject,handles);
PickTraces_Callback(hObject,handles);



% --- CALLED when any of the primary selection criteria checkbox setting is
% changed. Updates the selection criteria automatically,
% w/o running PickTraces().
function criteriaCheckbox_Callback(hObject, handles, criteriaName, textboxName)
textbox = handles.(textboxName);

if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.(criteriaName) = str2double(get(textbox,'String'));
    set(textbox,'Enable','on');
else
    handles.criteria.(criteriaName)=[];
    set(textbox,'Enable','off');
end

guidata(hObject,handles);
PickTraces_Callback(hObject,handles);




% --- Executes on selection change in any of the drop down boxes above each
%     of the trace statistic histogram plots.
function cboStat_Callback(hObject, eventdata, handles)
% hObject    handle to cboStat5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get ID of this combo control
id = get(hObject,'UserData');

% Get trace statistics
stats = getappdata(handles.figure1,'infoStruct');
stats = stats(handles.inds_picked);

% Get user selection
% statNames  = get(hObject,'String');
selected   = get(hObject,'Value');
% statToPlot = statNames{selected};
statNames = fieldnames(handles.statLongNames);
statToPlot = statNames{selected};

% Make sure it is a recognized stat
if ~ismember( statToPlot, fieldnames(stats) ),
    error('Selected trace statistic is unknown');
end

% Plot the distribution of the statistic
statData = [stats.(statToPlot)];
if any(isnan(statData))
    warning( 'NaN values found' );
%     statData( isnan(statData) ) = 0;
end

[data,binCenters] = hist( statData, handles.nHistBins);
data = 100*data/sum(data);  %normalize the histograms

axes( handles.(['axStat' num2str(id)]) );
bar( binCenters, data, 1 );
% xlabel(statToPlot);
zoom on;
grid on;

if id==1,
    ylabel( handles.axStat1, 'Number of Traces (%)' );
end

% Save histogram data in plot for launching cftool
if length(stats)>=1,
    set( handles.(['axStat' num2str(id)]), 'UserData', [binCenters;data] );
end


guidata(hObject,handles);






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
y = num2str(histData);
clipboard('copy', sprintf([y(1,:) '\n' y(2,:)]) );



