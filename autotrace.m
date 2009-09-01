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


% Last Modified by GUIDE v2.5 26-Jun-2009 17:17:51

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
criteria.minCorrelation=-1.1;
criteria.maxCorrelation=0.5;
criteria.minSNR=8;
criteria.maxBackground=1500;
criteria.maxDonorBlinks = 4;
criteria.minFretLifetime = 15;

handles.criteria = criteria;


%---- OTHER USER-TUNABLE PARAMETERS
handles.nHistBins=40; % histogram bins size
handles.contour_bin_size=0.035;
handles.sync='n';



% Initialize input fields with values defined above.
set(handles.lowThresh,'String',num2str(criteria.minCorrelation));
set(handles.highThresh,'String',num2str(criteria.maxCorrelation));
set(handles.SignalNoiseThresh,'String',num2str(criteria.minSNR));
set(handles.BackgroundNoiseThresh,'String',num2str(criteria.maxBackground));
set(handles.editNCross,'String',num2str(criteria.maxDonorBlinks));
set(handles.editAccLife,'String',num2str(criteria.minFretLifetime));

set(handles.FRETBinSize,'String',num2str(handles.contour_bin_size));



% Names of trace statistics -- these should be defined somewhere else!
ln = traceStat;  %get long statistic names
handles.statLongNames = ln;
longNames  = struct2cell(ln);
shortNames = fieldnames(ln);


% Add context menus to the plots to launch curve fitting
% Ideally, you should be able to drop-in any variant of the interface,
% with variable numbers of histogram boxes, etc and have it still work.
handles.cboNames = strcat('cboStat',{'1','2','3','4','5'});
handles.nPlots = length(handles.cboNames);

for id=1:handles.nPlots,
    menu = uicontextmenu;
    
    % Context menu for launching Matlab's curve fitting tool
    uimenu( menu, 'Label','Curve Fitting...', 'Callback',...
           ['autotrace(''launchFitTool_Callback'',gcbo,' num2str(id) ',guidata(gcbo))'] );
    
    % Context menu for copying the raw statistic data
    % for fitting in other programs (like Origin)
    uimenu( menu, 'Label','Copy data', 'Callback',...
           ['autotrace(''copyPlotData_Callback'',gcbo,' num2str(id) ',guidata(gcbo))'] );
   
    % Add the context menu to the axes
    set( handles.(['axStat' num2str(id)]), 'UIContextMenu', menu );
    
    % Also set options in combobox
    set( handles.(['cboStat' num2str(id)]), 'String', longNames );
end

% Set default selections
set( handles.cboStat1, 'Value', find(strcmp('t',shortNames))  );
set( handles.cboStat2, 'Value', find(strcmp('lifetime',shortNames))  );
set( handles.cboStat3, 'Value', find(strcmp('corr',shortNames))  );
set( handles.cboStat4, 'Value', find(strcmp('snr',shortNames))   );
set( handles.cboStat5, 'Value', find(strcmp('bg',shortNames))    );


%
handles.nCriteriaBoxes = 7;
criteriaNames = [{''}; longNames];

for id=1:handles.nCriteriaBoxes
    
    % Set options in selection crtieria comboboxes
    set( handles.(['cboCriteria' num2str(id)]), 'String', criteriaNames );

end


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
[datafile,datapath] = uigetfile( {'*.traces';'*.txt'},'Choose a traces file:', ...
                                 'MultiSelect','on');
if datapath==0, return; end
if ~iscell(datafile), datafile = {datafile}; end
filename = strcat(datapath,datafile);

if isempty(filename),  return;  end
if ~iscell(filename),  filename = {filename};  end

handles.inputdir = datapath;
handles.inputfiles = filename;
handles.nFiles = numel( filename );

disp(handles.inputfiles);
if numel( handles.inputfiles ) == 1,
    set(handles.editFilename,'String',handles.inputfiles{1});
else
    set(handles.editFilename,'String',handles.inputdir);
end

handles.outfile = strrep(filename{1}, '.traces', '_auto.txt');
handles.outfile = strrep(handles.outfile, '_01_auto.txt', '_auto.txt');

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
if datapath==0, return; end

% Create list of .traces files in the directory.
traces_files = dir( [datapath filesep '*.traces'] );
handles.nFiles = numel(traces_files);

if handles.nFiles == 0
    disp('No files in this directory!');
    return;
end

handles.inputdir = datapath;
handles.inputfiles = strcat( [datapath filesep], {traces_files.name} );


disp(handles.inputdir);
set(handles.editFilename,'String',handles.inputdir);

handles.outfile = strrep(handles.inputfiles{1}, '.traces', '_auto.txt');
handles.outfile = strrep(handles.outfile, '_01_auto.txt', '_auto.txt');



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

% Automatically run gettraces if only STKs are available.
% Uses automatic threshold picking, which may not work as well as a
% manual value!
% traces_files = dir( [datapath filesep '*.traces'] );
% stk_files    = dir( [datapath filesep '*.stk'] );
% stk_files    = [  stk_files  ;  dir([datapath filesep '*.stk.bz2'])  ];
% 
% Extract traces from unanalyzed movies, if requested
% if numel(traces_files)>0 && numel(stk_files)>0, 
%     answer = questdlg( ...
%         'Traces files already created. Reprocess and overwrite?', ...
%         'autotrace2: Reprocess movies?', 'Yes','No','No');
%     if strcmp(answer,'Yes')
%         gettraces_backend( datapath );
%     end
% else
%     gettraces_backend( datapath );
% end



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


    disp(handles.inputdir);
    set(handles.editFilename,'String',handles.inputdir);

    handles.outfile = strrep(handles.inputfiles{1}, '.traces', '_auto.txt');
    handles.outfile = strrep(handles.outfile, '_01_auto.txt', '_auto.txt');

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

handles.ids = cell(0);  % trace names (name_file#_trace#)
handles.nTracesPerFile = zeros(handles.nFiles,1);


% Open each file of traces and build the raw data array. Works the same as
% above, but loops through each file in the directory.
if ~handles.isBatchMode, wb=waitbar(0,'Loading traces...'); end

for k=1:handles.nFiles  % for each file...    
    
    % Load the traces file.
    % If raw data, corrections for background and crosstalk are made
    [donor,acceptor,fret,ids,timeAxis] = loadTraces( ...
                handles.inputfiles{k}, handles.constants);
    
    % Make sure all movies have the same number of frames
    if exist('len','var') && len~=size(donor,2),
        error('Trace lengths do not match! Use resizeTraces');
    end
    
    % Calculate lifetimes, and average amplitudes, etc for current file
    ss = traceStat(donor,acceptor,fret, handles.constants);
    
    % Add to values for all files
    if ~exist('infoStruct','var')
        infoStruct = ss;
    else
        infoStruct = cat(2, infoStruct, ss  );
    end
    
    % Add new data to handles
    [Ntraces,len] = size(donor);
    handles.nTracesPerFile(k) = Ntraces;
    handles.ids = [handles.ids ids];
    
    if ~handles.isBatchMode, waitbar(k/handles.nFiles,wb); end

end 
if ~handles.isBatchMode, close(wb); end

handles.len = len;
handles.Ntraces = sum(handles.nTracesPerFile);
handles.timeAxis = timeAxis;


% Save the trace properties values to application data
setappdata(handles.figure1,'infoStruct', infoStruct);
clear infoStruct;


% Initialize a variable for storing the number of molecules picked.
handles.picked_mols=0;

% Turn on/off buttons, and initialize the stat outputs.
set(handles.PickTraces,'Enable','on');

guidata(hObject,handles);

% Automatically run Pick Traces
PickTraces_Callback(handles.PickTraces,[],guidata(handles.PickTraces));


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

guidata(hObject,handles);



%--------------------  SAVE PICKED TRACES TO FILE --------------------%
function SaveTraces( filename, handles )

fid=fopen(filename,'w');
disp( ['Saving to ' filename] );

% Save forQUB
qub_fname = strrep( filename, '.txt', '.qub.txt' );
qubfid=fopen(qub_fname,'w');
disp( ['Saving to ' qub_fname] );

% Write time markers (first row)
fprintf(fid,'%d ', handles.timeAxis);
fprintf(fid,'\n');

if ~handles.isBatchMode, wb=waitbar(0,'Saving traces...'); end

pick_offset = [0; cumsum(handles.nTracesPerFile)];

for index = 1:handles.nFiles  %for each file in batch...
    
    %---- Load trace data from file, make corrections
    % inds_picked is indexes as if all the traces data were in one huge array.
    % This is translating into an offset at the start of this particular file
    Ntraces = handles.nTracesPerFile(index);
    picks = handles.inds_picked - pick_offset(index);
    picks = picks( picks>0 & picks<=Ntraces );
    
    if numel(picks)==0, continue; end
    
    [donor,acceptor,fret] = LoadTraces( handles.inputfiles{index}, ...
                                        handles.constants, picks );
    
    
    %--- Write fluorescence data: {name} {datapoints...}
    % 3 lines per molecule: donor, acceptor, fret
    indexes = picks + pick_offset(index);  %indexes into whole dataset
    
    for j=1:size(donor,1)  %for each molecule in file...
        
        % output name
        name = '';
        if ~isempty(handles.ids)
            name = sprintf('%s ',handles.ids{indexes(j)});
        end

        % output fluorescence data
        fprintf(fid,'%s', name);
        fprintf(fid,'%g ', donor(j,:));
        assert( ~isnan(donor(j,1)) );
        fprintf(fid,'\n');
        
        fprintf(fid,'%s', name);
        fprintf(fid,'%g ', acceptor(j,:));
        fprintf(fid,'\n');

        fprintf(fid,'%s', name);
        fprintf(fid,'%g ', fret(j,:));
        fprintf(fid,'\n');
        
        % Write FRET data
        fprintf(qubfid, '%f\n', fret(j,:));
        
    end % for each molecule
    
    
    if ~handles.isBatchMode, waitbar(index/handles.nFiles,wb); end
    
end % for each file
fclose(fid);
fclose(qubfid);
if ~handles.isBatchMode, close(wb); end


% Generate log file containing informtion about how the traces were picked.
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
            handles.picked_mols, handles.Ntraces, ...
            100*handles.picked_mols/handles.Ntraces );  

        
% Descriptive statistics about dataset
stats = getappdata(handles.figure1,'infoStruct');
% [picks,values] = pickTraces( stats, handles.criteria, handles.constants );

total = handles.Ntraces;
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

% --- Executes on button press in SaveProperties.
function SaveProperties_Callback(hObject, eventdata, handles)

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
fclose(fid);

% Write trace statistics
dlmwrite( filename, data', '-append', 'delimiter','\t' );


% END FUNCTION SaveProperties_Callback
      




% --------------- VIEW TRACES FUNCTIONALITY --------------- %

% --- Executes on button press in ViewPickedTraces.
function ViewPickedTraces_Callback(hObject, eventdata, handles)

[inputfile inputpath]=...
    uiputfile('.txt','Save picked traces as:',handles.outfile);
if inputfile==0, return; end

handles.outfile=[inputpath inputfile];


% Save picked data to handles.outfile
SaveTraces( handles.outfile, handles );

guidata(hObject,handles);


% Run sorttraces interface so traces can be viewed
% Could pass description/title as third param (currently empty!)
sorttraces(0, handles.outfile);





function criteria = getSpecialCriteria( handles )

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
    value    = str2double( get(handles.(['edCriteria' num2str(id)]),'String') );
    
    criteria.(criteriaName) = value;
end



%----------APPLIES PICKING CRITERIA TO TRACES----------
% --- Executes on button press in PickTraces.
function PickTraces_Callback(hObject, eventdata, handles)

criteria = getSpecialCriteria( handles );


% Find which molecules pass the selection criteria
stats = getappdata(handles.figure1,'infoStruct');
[picks,values] = pickTraces( stats, criteria );
clear stats;

% The number of traces picked.
handles.inds_picked = picks;
handles.picked_mols = numel(handles.inds_picked);

% If at least one trace is picked, turn some buttons on.
if handles.picked_mols > 0
    set(handles.SaveTraces,'Enable','on');
    set(handles.MakeContourPlot,'Enable','on');
    set(handles.ViewPickedTraces,'Enable','on');
end

% Turn some other buttons on/off.
set(handles.MoleculesPicked,'String', ...
            sprintf('%d of %d',[handles.picked_mols,handles.Ntraces]));
set(handles.SaveContourPlot,'Enable','off');




guidata(hObject,handles);


% Replot
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

if 0,  % simpler way to do it...but slower.
    % Save current selections as a temporary file
    filename = [tempname '_auto.txt'];
    SaveTraces( filename, handles );

    % Make contour plots using makeplots.m
    [p,title] = fileparts(handles.outfile);
    title = strrep( title,'_',' ' );
    makeplots( {filename}, {title} );

    return;
end

constants = handles.constants;


h2=waitbar(0,'Making contour plot...');

pick_offset = [0; cumsum(handles.nTracesPerFile)];
    
% Axes for histogram:
time_axis=1:handles.len;  %sd(2)
fret_axis=-0.1:handles.contour_bin_size:1.0;

% Initialize histogram array, setting the time step in the first row,
% and the FRET bins in the first column. This is done for import into
% Origin.
handles.frethist = zeros( length(fret_axis)+1, length(time_axis)+1 );
handles.frethist(1,2:end) = time_axis;
handles.frethist(2:end,1) = fret_axis';

for index = 1:handles.nFiles  %for each file in batch...

    % inds_picked is indexes as if all the traces data were in one huge array.
    % This is translating into an offset at the start of this particular file
    Ntraces = handles.nTracesPerFile(index);
    picks = handles.inds_picked - pick_offset(index);
    picks = picks( picks>0 & picks<=Ntraces );  %indexes just into this file
    
    if numel(picks)==0, continue; end

    filename = handles.inputfiles{index};

    if handles.sync == 'n',
        [donor,acceptor,fret] = loadTraces(filename,constants,picks);

    % If a synchronization method is selected, apply it here...
    elseif handles.sync == 's',
        fret = autosort( filename, picks );
        nTraces = size(fret,1);
        
        for i=1:nTraces
            s = find(fret(i,:)>=0.125,1);
            s = max(s-8,1);
            fret(i,:) = [fret(i,s:end) zeros(1,s-1)];
        end
        
        if size(fret,1) == 0, continue; end
    end

    assert( size(fret,2) == handles.len, ...
            'Trace lengths do not match! Use resizeTraces.' );

    % hist has different behavior for a 1xN vector than MxN, so we
    % pad zeros to make it do the same thing in both cases.
    if size(fret,1) == 1
        fret = [fret ; repmat(NaN,1,handles.len)];
    end

    handles.frethist(2:end,2:end) = handles.frethist(2:end,2:end) + ...
                                    hist( fret, fret_axis  );

    waitbar(index/handles.nFiles,h2);
end %for each file of traces

% Plot the histograms
[p,title] = fileparts(handles.outfile);
title = strrep( title,'_',' ' );
makeplots( handles.frethist, title );

set(handles.SaveContourPlot,'Enable','on');



close(h2);  %close waitbar
guidata(hObject,handles);



 
%----------SAVE CONTOUR PLOT----------
% --- Executes on button press in SaveContourPlot.
function SaveContourPlot_Callback(hObject, eventdata, handles)

% Write original file
histfile=strrep(handles.outfile,'.txt','_hist.txt');
dlmwrite(histfile,handles.frethist,' ');

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
switch opt
    case 1 % No post-synchronization. This is the default.
        handles.sync='n';
    case 2 % Seperate events and synchronize
        handles.sync='s';
    case 3 % Post-synchronize each trace to the first point above a threshold.
        handles.sync='y';
    case 4 % Post-synchronize each FRET event which crosses a threshold.
        handles.sync='i';
end
guidata(hObject,handles);


function FRETBinSize_Callback(hObject, eventdata, handles)
handles.criteria.contour_bin_size=str2double(get(hObject,'String'));
guidata(hObject,handles);



%----------CHECKBOX FOR MEAN TOTAL INTENSITY CRITERIA----------
% --- Executes on button press in MeanTotalIntensityBox.
function MeanTotalIntensityBox_Callback(hObject, eventdata, handles)
% If box is checked, get the values input by the user.
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minTotalIntensity=str2double(get(handles.MeanTotalIntensityLow,'String'));
    handles.criteria.maxTotalIntensity=str2double(get(handles.MeanTotalIntensityHigh,'String'));
    set(handles.MeanTotalIntensityLow,'Enable','on');
    set(handles.MeanTotalIntensityHigh,'Enable','on');
else
    % If the box is unchecked, use values which will nullify the criterium.
    handles.criteria.minTotalIntensity=[];
    handles.criteria.maxTotalIntensity=[];
    set(handles.MeanTotalIntensityLow,'Enable','off');
    set(handles.MeanTotalIntensityHigh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR FRET THRESHOLD----------
% --- Executes on button press in fretThreshold.
function fretThreshold_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minFret=str2double(get(handles.fretSlopeThresh,'String'));
    set(handles.fretSlopeThresh,'Enable','on');
else
    handles.criteria.minFret=[];
    set(handles.fretSlopeThresh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR FLUORESCENCE LIFETIME CRITERIA----------
% --- Executes on button press in FluorescenceLifetimeBox.
function FluorescenceLifetimeBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minTotalLifetime=str2double(get(handles.FluorescenceLifetime,'String'));
    handles.criteria.maxTotalLifetime=str2double(get(handles.FluorescenceLifetimeHigh,'String'));
    set(handles.FluorescenceLifetime,'Enable','on');
    set(handles.FluorescenceLifetimeHigh,'Enable','on');
else
    handles.criteria.minTotalLifetime=[];
    handles.criteria.maxTotalLifetime=[];
    set(handles.FluorescenceLifetime,'Enable','off');
    set(handles.FluorescenceLifetimeHigh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR CORRELATION COEFFICIENT CRITERIA----------
% --- Executes on button press in CorrelationCoefficientBox.
function CorrelationCoefficientBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minCorrelation=str2double(get(handles.lowThresh,'String'));
    handles.criteria.maxCorrelation=str2double(get(handles.highThresh,'String'));
    set(handles.lowThresh,'Enable','on');
    set(handles.highThresh,'Enable','on');
else
    handles.criteria.minCorrelation=[];
    handles.criteria.maxCorrelation=[];
    set(handles.lowThresh,'Enable','off');
    set(handles.highThresh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR SIGNAL-TO-NOISE CRITERIA----------
% --- Executes on button press in SignalNoiseBox.
function SignalNoiseBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minSNR=str2double(get(handles.SignalNoiseThresh,'String'));
    set(handles.SignalNoiseThresh,'Enable','on');
else
    handles.criteria.minSNR=[];
    set(handles.SignalNoiseThresh,'Enable','off');
end
guidata(hObject,handles);


%----------CHECKBOX FOR BACKGROUND NOISE CRITERIA----------
% --- Executes on button press in BackgroundNoiseBox.
function BackgroundNoiseBox_Callback(hObject, eventdata, handles)
% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.maxBackground=str2double(get(handles.BackgroundNoiseThresh,'String'));
    set(handles.BackgroundNoiseThresh,'Enable','on');
else
    handles.criteria.maxBackground=[];
    set(handles.BackgroundNoiseThresh,'Enable','off');
end
guidata(hObject,handles);


% --- Executes on button press in checkboxNCross.
function checkboxNCross_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkboxNCross

% Same as above...
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.maxDonorBlinks=str2double(get(handles.editNCross,'String'));
    set(handles.editNCross,'Enable','on');
else
    handles.criteria.maxDonorBlinks=[];
    set(handles.editNCross,'Enable','off');
end
guidata(hObject,handles);



% --- Executes on button press in checkboxAccLife.
function checkboxAccLife_Callback(hObject, eventdata, handles)
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minFretLifetime = str2double(get(handles.editAccLife,'String'));
    set(handles.editAccLife,'Enable','on');
else
    handles.criteria.minFretLifetime=[];
    set(handles.editAccLife,'Enable','off');
end
guidata(hObject,handles);


function chkOverlap_Callback(hObject, eventdata, handles)
handles.criteria.overlap = get(hObject,'Value');
guidata(hObject,handles);


function updateCriteria_Callback( hObject, handles, criteriaName )
handles.criteria.(criteriaName) = str2double(get(hObject,'String'));
guidata(hObject,handles);




% --- Executes on selection change in cboStat5.
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

if length(stats)<1, return; end

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
set( handles.(['axStat' num2str(id)]), 'UserData', [binCenters;data] );


guidata(hObject,handles);






function launchFitTool_Callback(hObject, id, handles)
% Callback for context menu for trace statistics plots.
% Launches Curve Fitting Tool using the data in the selected plot.

ax = handles.axStat5(id);

% Get ID of this combo control
histData = get(ax,'UserData');

binCenters = histData(1,:);
data = histData(2,:);

cftool(binCenters,data);


function copyPlotData_Callback(hObject, id, handles)
% Callback for context menu for trace statistics plots.
% Launches Curve Fitting Tool using the data in the selected plot.

ax = handles.axStat5(id);

% Get ID of this combo control
histData = get(ax,'UserData');

% Copy to clipboard
y = num2str(histData);
clipboard('copy', sprintf([y(1,:) '\n' y(2,:)]) );













% --- Executes on button press in chkIntSigma.
function chkIntSigma_Callback(hObject, eventdata, handles)
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.maxTotalSigma = str2double(get(handles.edIntSigma,'String'));
    set(handles.edIntSigma,'Enable','on');
else
    handles.criteria.maxTotalSigma=[];
    set(handles.edIntSigma,'Enable','off');
end
guidata(hObject,handles);



% --- Executes on button press in chkFretEvents.
function chkFretEvents_Callback(hObject, eventdata, handles)
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minFretEvents = str2double(get(handles.edFretEvents,'String'));
    set(handles.edFretEvents,'Enable','on');
else
    handles.criteria.minFretEvents=[];
    set(handles.edFretEvents,'Enable','off');
end
guidata(hObject,handles);



