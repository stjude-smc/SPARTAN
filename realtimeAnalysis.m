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


% Last Modified by GUIDE v2.5 03-Mar-2009 15:52:55

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

handles.hSettings = realtimeAnalysis_settings( ...
                            {@realtimeAnalysis_Notify, hObject} );

% Load settings from the dialog...
settingsHandles = guidata( handles.hSettings );
handles.criteria = settingsHandles.criteria;

% Give a callback to the settings dialog so it can update the code here
% in the event of a change in criteria values


% Leave everything alone if the program is already running.
% This initialization proceedure will confuse the program state.
% if isfield(handles,'criteria'),
%     disp('Autotrace is already running!');
%     return;
% end

%---- PROGRAM CONSTANTS
constants = cascadeConstants();
handles.constants = constants;


% Choose default command line output for realtimeAnalysis
handles.output=hObject;

% Update handles structure
guidata(hObject,handles);

% END FUNCTION realtimeAnalysis_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = realtimeAnalysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=handles.output;

% END FUNCTION realtimeAnalysis_OutputFcn



function realtimeAnalysis_Notify(hObject)

% Load data in both figures
handles = guidata(hObject);
settingsHandles = guidata(handles.hSettings);

% Update selection criteria
handles.criteria = settingsHandles.criteria;

disp(handles.criteria);

% Update handles structure
guidata(hObject,handles);



%#########################################################################
%----------------------- LOAD, FILTER, SAVE TRACES ----------------------%
%#########################################################################


% --- Executes on button press in btnBrowse.
function btnBrowse_Callback(hObject, eventdata, handles)
% CALLED: when the user clicked "Browse..."
% ACTION: Get location of data to process
 
datadir = uigetdir(pwd, 'Select a directory with all data to process');

if datadir~=0,
    handles.directory = datadir;
    set(handles.txtDirectory,'String',datadir);
end

% Update handles structure
guidata(hObject,handles);



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
%         'realtimeAnalysis: Reprocess movies?', 'Yes','No','No');
%     if strcmp(answer,'Yes')
%         gettraces_backend( datapath );
%     end
% else
%     gettraces_backend( datapath );
% end


settingsHandles = guidata(handles.hSettings);

% Update selection criteria
handles.criteria = settingsHandles.criteria;


datapath = get(handles.txtDirectory,'String');


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
%     set(handles.editFilename,'String',handles.inputdir);

    handles.outfile = strrep(handles.inputfiles{1}, '.traces', '_auto.txt');
    handles.outfile = strrep(handles.outfile, '_01_auto.txt', '_auto.txt');


    OpenTracesBatch( hObject, handles )
end


% END FUNCTION btnGo_Callback




function OpenTracesBatch( hObject, handles )

handles.ids = cell(0);  % trace names (name_file#_trace#)
handles.nTracesPerFile = zeros(handles.nFiles,1);


% Open each file of traces and build the raw data array. Works the same as
% above, but loops through each file in the directory.
% wb=waitbar(0,'Loading traces...');
for k=1:handles.nFiles  % for each file...    
    
    % Load the traces file.
    % If raw data, corrections for background and crosstalk are made
    [donor,acceptor,fret,ids] = loadTraces( ...
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
    
%     waitbar(k/handles.nFiles,wb);
    
end 
% close(wb);

handles.len = len;
handles.Ntraces = sum(handles.nTracesPerFile);


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

% Write time markers (first row)
fprintf(fid,'%d ', 1:handles.len);
fprintf(fid,'\n');

% wb=waitbar(0,'Saving traces...');

pick_offset = [0; cumsum(handles.nTracesPerFile)];

for index = 1:handles.nFiles  %for each file in batch...
    
    %---- Load trace data from file, make corrections
    % inds_picked is indexes as if all the traces data were in one huge array.
    % This is translating into an offset at the start of this particular file
    Ntraces = handles.nTracesPerFile(index);
    picks = handles.inds_picked - pick_offset(index);
    picks = picks( picks>0 & picks<=Ntraces );
    
    if numel(picks)==0, continue; end
    
    [donor,acceptor,fret] = loadTraces( handles.inputfiles{index}, ...
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
    
    
%     waitbar(index/handles.nFiles,wb);
    
end % for each file
fclose(fid);
fclose(qubfid);
% close(wb);


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





%#########################################################################
%------------------------- INPUT FIELD HANDLING -------------------------%
%#########################################################################



%----------CHECKBOX FOR MEAN TOTAL INTENSITY CRITERIA----------
% --- Executes on button press in MeanTotalIntensityBox.
function MeanTotalIntensityBox_Callback(hObject, eventdata, handles)
% If box is checked, get the values input by the user.
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.criteria.minTotalIntesity=str2double(get(handles.MeanTotalIntensityLow,'String'));
    handles.criteria.maxTotalIntesity=str2double(get(handles.MeanTotalIntensityHigh,'String'));
    set(handles.MeanTotalIntensityLow,'Enable','on');
    set(handles.MeanTotalIntensityHigh,'Enable','on');
else
    % If the box is unchecked, use values which will nullify the criterium.
    handles.criteria.minTotalIntesity=[];
    handles.criteria.maxTotalIntesity=[];
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


function MeanTotalIntensityLow_Callback(hObject, eventdata, handles)
handles.criteria.minTotalIntesity=str2double(get(hObject,'String'));
guidata(hObject,handles);


function MeanTotalIntensityHigh_Callback(hObject, eventdata, handles)
handles.criteria.maxTotalIntesity=str2double(get(hObject,'String'));
guidata(hObject,handles);


function fretSlopeThresh_Callback(hObject, eventdata, handles)
handles.criteria.minFret=str2double(get(hObject,'String'));
guidata(hObject,handles);


function FluorescenceLifetime_Callback(hObject, eventdata, handles)
handles.criteria.minTotalLifetime=str2double(get(hObject,'String'));
guidata(hObject,handles);


function FluorescenceLifetimeHigh_Callback(hObject, eventdata, handles)
handles.criteria.maxTotalLifetime=str2double(get(hObject,'String'));
guidata(hObject,handles);


function lowThresh_Callback(hObject, eventdata, handles)
handles.criteria.minCorrelation=str2double(get(hObject,'String'));
guidata(hObject,handles);


function highThresh_Callback(hObject, eventdata, handles)
handles.criteria.maxCorrelation=str2double(get(hObject,'String'));
guidata(hObject,handles);


function SignalNoiseThresh_Callback(hObject, eventdata, handles)
handles.criteria.minSNR=str2double(get(hObject,'String'));
guidata(hObject,handles);



function BackgroundNoiseThresh_Callback(hObject, eventdata, handles)
handles.criteria.maxBackground=str2double(get(hObject,'String'));
guidata(hObject,handles);


function FRETBinSize_Callback(hObject, eventdata, handles)
handles.criteria.contour_bin_size=str2double(get(hObject,'String'));
guidata(hObject,handles);


function editNCross_Callback(hObject, eventdata, handles)
handles.criteria.maxDonorBlinks = str2double(get(hObject,'String'));
guidata(hObject,handles);


function editAccLife_Callback(hObject, eventdata, handles)
handles.criteria.minFretLifetime = str2double(get(hObject,'String'));
guidata(hObject,handles);


function chkOverlap_Callback(hObject, eventdata, handles)
handles.criteria.overlap = get(hObject,'Value');
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
    loc = handles.directory;
    set(hObject,'String',loc);
else
    handles.directory = loc;
end

guidata(hObject,handles);




% --- Executes on button press in btnSettings.
function btnSettings_Callback(hObject, eventdata, handles)
% If the realtime analysis window has not been launched, do so now.
% If it has been launched, simply make it visible.
if isfield('handles','hSettings'),
    set( handles.hSettings, 'Visible','on' );
else
    handles.hSettings = realtimeAnalysis;
end

