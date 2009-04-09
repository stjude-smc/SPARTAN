function varargout = realtimeAnalysis(varargin)
% REALTIMEANALYSIS  Trace processing and filtering
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


% Last Modified by GUIDE v2.5 08-Apr-2009 17:23:53

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
                            {@realtimeAnalysis_Notify, hObject} );
set(handles.hSettings,'Visible','off');
                        
% Load settings from the dialog...
settingsHandles = guidata( handles.hSettings );
handles.criteria = settingsHandles.criteria;


%---- PROGRAM CONSTANTS
handles.constants = cascadeConstants();

% Choose default command line output for realtimeAnalysis
handles.output=hObject;

% Update handles structure
guidata(hObject,handles);

% END FUNCTION realtimeAnalysis_OpeningFcn



% Callback for when settings dialog is updated...
function realtimeAnalysis_Notify(hObject)

% Load data in both figures
handles = guidata(hObject);
settingsHandles = guidata(handles.hSettings);

% Update selection criteria
handles.criteria = settingsHandles.criteria;

disp(handles.criteria);

% Update handles structure
guidata(hObject,handles);

btnGo_Callback(hObject, [], handles);




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
% --- Executes on button press in btnGo.
function btnGo_Callback(hObject, eventdata, handles)
% hObject    handle to btnGo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This allows for realtimeAnalysis to run in batch mode. A directory is selected,
% and all the traces files are combined, and a single output file is
% generated. NOTE: be careful not to combine data from different
% experiments. Store different data sets in separate folders.


tic;
handles.isExecuting = 1;
guidata(hObject,handles);

set( handles.txtStatus, 'String', 'Loading traces...' );

disp( get(0,'CurrentFigure') )
settingsHandles = guidata(handles.hSettings);

% Update selection criteria
handles.criteria = settingsHandles.criteria;


datapath = get(handles.txtDirectory,'String');


% Run gettraces ...


% Create list of .traces files in the directory.
traces_files = dir( [datapath filesep '*.traces'] );
handles.nFiles = numel(traces_files);

if handles.nFiles == 0
    disp('No files in this directory!');
    return;
end

handles.inputdir = datapath;
handles.inputfiles = strcat( [datapath filesep], {traces_files.name} );
handles.inputstks = {};

disp(handles.inputdir);

handles.outfile = strrep(handles.inputfiles{1}, '.traces', '_auto.txt');
handles.outfile = strrep(handles.outfile, '_01_auto.txt', '_auto.txt');


OpenTracesBatch( hObject, handles )


set( handles.txtStatus, 'String','IDLE.' ); drawnow;
handles.isExecuting = 0;

guidata(hObject,handles);
toc

% END FUNCTION btnGo_Callback




function OpenTracesBatch( hObject, handles )

handles.ids = cell(0);  % trace names (name_file#_trace#)
handles.nTracesPerFile = zeros(handles.nFiles,1);


% Open each file of traces and build the raw data array. Works the same as
% above, but loops through each file in the directory.
for k=1:handles.nFiles  % for each file...    
    
    set( handles.txtStatus, 'String', ...
            sprintf('Loading traces (%d of %d)',k,handles.nFiles) );
    drawnow;
    
    % Load the traces file.
    % If raw data, corrections for background and crosstalk are made
    [donor,acceptor,fret,ids,time] = loadTraces( ...
                handles.inputfiles{k}, handles.constants);
            
    if size(donor,1)==0, continue; end  %skip empty files
    
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
    
end 
drawnow;

handles.len = len;
handles.Ntraces = sum(handles.nTracesPerFile);
if time(1)~=1,
    dt = time(2)-time(1);
end


% Save the trace properties values to application data
setappdata(handles.figure1,'infoStruct', infoStruct);
% clear infoStruct;


% Initialize a variable for storing the number of molecules picked.
handles.picked_mols=0;

guidata(hObject,handles);

% Automatically run Pick Traces
% stats = getappdata(handles.figure1,'infoStruct');
[picks,values] = pickTraces( infoStruct, handles.criteria );
% clear stats;

% The number of traces picked.
handles.inds_picked = picks;
handles.picked_mols = numel(handles.inds_picked);

guidata(hObject,handles);

% Save picked data to handles.outfile
SaveTraces( handles.outfile, handles );


% Run Hidden Markov Modeling analysis, if a model has been given

if isfield(handles,'model')
    skmOptions.maxItr = 100;
    skmOptions.convLL = 1e-4;
    handles.model.fixMu = ones( size(handles.model.mu) );
    handles.model.fixSigma = ones( size(handles.model.sigma) );


    set( handles.txtStatus, 'String','Idealizing data...' ); drawnow;

    filename = handles.outfile;
    [d,a,data,ids,time] = loadTraces( filename );
    if time(1)==1,
        f = inputdlg('What is the sampling interval (in ms) for this data?');
        sampling = str2double(f);
    else
        sampling = time(2)-time(1);
    end
    
    [dwt,newModel,LL,offsets] = skm( data, sampling, handles.model, skmOptions );

    mu = newModel.mu;
    sigma = newModel.sigma;
    fretModel = [mu sigma];
    
    % Save the idealization
    dwtFilename = strrep(filename,'.txt','.qub.dwt');
    saveDWT( dwtFilename, dwt, offsets, fretModel, sampling );
end


% Generate ensemble plots
targetAxes = { handles.axFretContour, handles.axFretHistogram };
set( handles.txtStatus, 'String','Making data plots...' ); drawnow;
makeplots( {handles.outfile}, {'auto'}, 'targetAxes',targetAxes );

% Show population statistics in GUI
set( handles.txtAcceptance,'String', ...
     sprintf('%.0f%% (%d of %d)',100*handles.picked_mols/handles.Ntraces, ...
                handles.picked_mols,handles.Ntraces) );

avgI = mean( [values.t] );
set( handles.txtIntensity,'String', ...
     sprintf('%.0f',avgI) );
 
avgSNR = mean( [values.snr] );
set( handles.txtSNR,'String', ...
     sprintf('%.1f',avgSNR) );

% Calculate donor bleaching rate...
[donorDist,donorAxes] = hist( [values.lifetime], 40 );
donorDist = 1 - [0 cumsum(donorDist)]/sum(donorDist);
donorAxes = [0 donorAxes];

f = fit( donorAxes',donorDist','exp1' );

if exist('dt','var')
    set( handles.txtLTDonor,'String', ...
         sprintf('%.1f sec',-1/f.b/dt) );
else
    set( handles.txtLTDonor,'String', ...
         sprintf('%.1f frames',-1/f.b) );
end

% Calculate acceptor bleaching rate
[accDist,accAxes] = hist( [values.acclife], 40 );
accDist = 1 - [0 cumsum(accDist)]/sum(accDist);
accAxes = [0 accAxes];

f = fit( accAxes',accDist','exp1' );

if exist('dt','var')
    set( handles.txtLTAcceptor,'String', ...
         sprintf('%.1f sec',-1/f.b/dt) );
else
    set( handles.txtLTAcceptor,'String', ...
         sprintf('%.1f frames',-1/f.b) );
end

% Calculate state occupancies
if isfield(handles,'model')
    info = sprintf('%.1f%% ', percentTime(dwtFilename));
    set(handles.edOccupancy,'String',info);
end

guidata(hObject,handles);



% END FUNCTION OpenTracesBatch


%--------------------  SAVE PICKED TRACES TO FILE --------------------%
function SaveTraces( filename, handles )

fid=fopen(filename,'w');
disp( ['Saving to ' filename] );

% Save forQUB
qub_fname = strrep( filename, '.txt', '.qub.txt' );
qubfid=fopen(qub_fname,'w');
disp( ['Saving to ' qub_fname] );


pick_offset = [0; cumsum(handles.nTracesPerFile)];

for index = 1:handles.nFiles  %for each file in batch...
    
    set( handles.txtStatus, 'String', ...
        sprintf('Saving traces (%.0f%%)',100*index/handles.nFiles) );
    drawnow;
    
    %---- Load trace data from file, make corrections
    % inds_picked is indexes as if all the traces data were in one huge array.
    % This is translating into an offset at the start of this particular file
    Ntraces = handles.nTracesPerFile(index);
    picks = handles.inds_picked - pick_offset(index);
    picks = picks( picks>0 & picks<=Ntraces );
    
    if numel(picks)==0, continue; end
    
    [donor,acceptor,fret,ids,time] = loadTraces( handles.inputfiles{index}, ...
                                        handles.constants, picks );
      
    % Write time markers (first row)
    if index==1,
        fprintf(fid,'%d ', time);
        fprintf(fid,'\n');
    end
    
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
    
    
    
end % for each file
fclose(fid);
fclose(qubfid);
drawnow;


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
fprintf(fid,'PICKING CRITERIA\n');

names = fieldnames(  handles.criteria );
vals  = struct2cell( handles.criteria );

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




%%
%#########################################################################
%------------------------- INPUT FIELD HANDLING -------------------------%
%#########################################################################


% --- Executes on button press in btnBrowse.
function btnBrowse_Callback(hObject, eventdata, handles)
% CALLED: when the user clicked "Browse..."
% ACTION: Get location of data to process
 
datadir = uigetdir(pwd, 'Select a directory with all data to process');

if datadir~=0,
    handles.inputdir = datadir;
    set(handles.txtDirectory,'String',datadir);

    handles.inputfiles = {};
    handles.nFiles = 0;
end

% Update handles structure
guidata(hObject,handles);


function txtDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to txtDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtDirectory as text
%        str2double(get(hObject,'String')) returns contents of txtDirectory as a double

loc = get(hObject,'String');
if ~exist(loc,'dir'),
    warning('realtimeAnalysis: directory doesn''t exist!');
    loc = handles.inputdir;
    set(hObject,'String',loc);
else
    handles.inputdir = loc;
end

handles.inputfiles = {};
handles.nFiles = 0;

guidata(hObject,handles);




% --- Executes on button press in btnSettings.
function btnSettings_Callback(hObject, eventdata, handles)
% If the realtime analysis window has not been launched, do so now.
% If it has been launched, simply make it visible.
if isfield(handles,'hSettings'),
    set( handles.hSettings, 'Visible','on' );
else
    handles.hSettings = realtimeAnalysis_settings( ...
                            {@realtimeAnalysis_Notify, hObject} );
    disp('loading settings gui');
end
guidata(hObject,handles);


% --- Executes on button press in chkAutoUpdate.
function chkAutoUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to chkAutoUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkAutoUpdate

handles.autoUpdate = get(hObject,'Value');

% If turned on, create a timer object to check the current directory
if handles.autoUpdate    
    disp('starting');
    handles.fileTimer = timer('ExecutionMode','fixedDelay','StartDelay',1,...
                              'TimerFcn',{@checkForFiles,handles.figure1},...
                              'Period',5.0,'BusyMode','drop');
    start(handles.fileTimer);

% If turned off, disable the current timer
else
    
    if isfield(handles,'fileTimer')
        disp('deleting');
        stop(handles.fileTimer);
        delete(handles.fileTimer);
    end

end

guidata(hObject,handles);



function checkForFiles(timerObject,event,figObj)

handles = guidata(figObj);

% If there is no current file list, exit
if ~isfield(handles,'inputfiles')
    return;
end

if isfield(handles,'isExecuting') && handles.isExecuting,
    disp('Already running, skipping test');
    return;
end

datapath = handles.inputdir;


% If there are new STKs, run gettraces
set( handles.txtStatus, 'String','Processing new movies...' );

params.overlap_thresh = 2.1;
params.don_thresh = 3000;
params.skipExisting = 1;
params.recursive = 0;
gettraces( datapath, params );

set( handles.txtStatus, 'String','IDLE.' ); drawnow;




% Check current list against one from previous polling
traces_files = dir( [datapath filesep '*.traces'] );
if numel(traces_files) == 0, return; end
inputfiles = strcat( [datapath filesep], {traces_files.name} );

% If new files exist, re-run analysis routine.
if numel(inputfiles)~=handles.nFiles || ...
   ~all( strcmp(handles.inputfiles,inputfiles) ),
%     handles.nFiles = numel(traces_files);
%     handles.inputfiles = inputfiles;
    disp('new file');
    figure(figObj);
    disp( get(0,'CurrentFigure') );
    realtimeAnalysis('btnGo_Callback',figObj,[],guidata(figObj))
%     btnGo_Callback( handles.btnGo,[],handles );
end













% --- Executes on button press in btnBrowseModel.
function btnBrowseModel_Callback(hObject, eventdata, handles)
% hObject    handle to btnBrowseModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get list of data files from user
fname = get(handles.edModelFilename,'String');
if isempty(fname)
    fname = '*.qmf';
end

[f,p] = uigetfile(fname,'Select a QuB model file...');
if f==0, return; end
fname = [p f];

set(handles.edModelFilename,'String',fname);

handles = loadModel(handles,fname);
guidata(hObject, handles);

% Update results -- rerun analysis pipeline
realtimeAnalysis('btnGo_Callback',hObject,[],handles);



function edModelFilename_Callback(hObject, eventdata, handles)
% hObject    handle to edModelFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = get(hObject,'String');
if ~exist(filename,'file'),
    warning('Model file does not exist!');
else
    handles = loadModel(handles,filename);
end

guidata(hObject, handles);

% Update results -- rerun analysis pipeline
realtimeAnalysis('btnGo_Callback',hObject,[],handles);



function handles = loadModel(handles,filename)

model = qub_loadModel( filename );
handles.model = model;
handles.hasModel = 1;

% Update model status text -- FIXME
info = sprintf('%s',mat2str(model.mu(model.class)));
set(handles.edFretValues,'String',info);

% If data are already loaded, enable the Execute button
% if handles.hasData,
%     set(handles.btnExecute,'Enable','on');
% end


