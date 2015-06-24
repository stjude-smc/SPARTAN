function varargout = sorttraces(varargin)
% SORTTRACES   GUI to manually select traces and make corrections.
%
%    sorttraces can be used to load a traces file; select individual
%    traces; make corrections to crosstalk, background subtraction, and
%    FRET calculation; and save these "binned" traces into new traces
%    files. By default, these are No, All, and Best FRET. "No" could be
%    unusable molecules with serious artifacts or no acceptor. "All" could
%    be all usable molecules. "Best" could be a few very high-quality,
%    repressentative traces used for publication.
%

% Depends on: sorttraces.fig, LoadTraces.m, CorrectTraces, cascadeConstants,
%    trace_stat (which requires: RLE_filter, CalcLifetime)

% Last Modified by GUIDE v2.5 15-Jun-2015 12:46:03


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @sorttraces_OpeningFcn, ...
    'gui_OutputFcn',  @sorttraces_OutputFcn, ...
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






% --- Executes just before sorttraces is made visible.
function sorttraces_OpeningFcn(hObject, ~, handles, varargin)
% Setup GUI controls in their default state (no file loaded).


% Initialize GUI if sorrtraces is being launched for the first time.
if ~isfield(handles,'constants')
    disp('Keyboard shortcuts:');
    disp('   Left arrow key  = Previous molecule');
    disp('   Right arrow key = Next molecule');
    disp('   a, s, d = Put molecule in No, All, or Best FRET, respectively');
    disp('   z       = Zoom in to trace. Press again to zoom out');    
    
    % Choose default command line output for sorttraces
    handles.output = hObject;
    handles.vals = [];

    % Initialize some variables, utilizing the handles structure
    handles.constants = cascadeConstants();
    set(handles.sldCrosstalk1, 'Value',  0);
    set(handles.sldThreshold, 'Value',  100);
    set(handles.edCrosstalk1,  'String', '0');

    % Link x-axes - zooming on one plot will automatically zoom on the other
    linkaxes([handles.axFluor handles.axTotal handles.axFret],'x');
    
    % SETUP AXES labels and settings.
    % Hold is needed so we don't lose the grid/zoom/etc settings, which are
    % lost when a new plot is generated on the axes...
    ylabel( handles.axFluor, 'Fluorescence' );
    grid( handles.axFluor, 'on' );
    zoom( handles.axFluor, 'on' );
    hold( handles.axFluor, 'on' );
    
    ylabel( handles.axTotal, 'Total Fluorescence' );
    grid( handles.axTotal, 'on' );
    zoom( handles.axTotal, 'on' );
    hold( handles.axTotal, 'on' );

    xlabel( handles.axFret, 'Frame Number' );
    ylabel( handles.axFret, 'FRET Efficiency' );
    ylim( handles.axFret, [-0.1 1] );
    grid( handles.axFret, 'on' );
    zoom( handles.axFret, 'on' );
    hold( handles.axFret, 'on' );
    
    set( zoom(handles.axFluor),'ActionPostCallback',@zoom_callback);
    set( zoom(handles.axTotal),'ActionPostCallback',@zoom_callback);
    set( zoom(handles.axFret), 'ActionPostCallback',@zoom_callback);
    
    handles.axFOV = [];
end

% Update handles structure
guidata(hObject, handles);


% If called by autotrace, 2nd argument is filename of traces output file.
if numel(varargin) > 1
    handles.filename = varargin{2};
    handles = OpenTracesFile( handles.filename, handles );
    guidata(hObject, handles);
        
    % If a trace number is also specified, jump to that trace.
    % TODO: Ideally, the file would be reloaded only if necessary.
    if numel(varargin) > 2,
        set( handles.editGoTo, 'String',num2str(varargin{3}) );
        editGoTo_Callback(handles.editGoTo, [], handles);
    end
end



%--- END FUNCTION sorttraces_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = sorttraces_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%--- END FUNCTION sorttraces_OutputFcn





%=========================================================================%
%=========================   LOAD TRACES FILES   =========================%


%----------"OPEN TRACES FILE" Button----------%
function btnOpen_Callback(~, ~, handles)
% Open a user-selected traces file.

% Get traces filename by menu driven input
filter = {'*.traces;*.rawtraces','Binary Traces Files (*.traces,*.rawtraces)'; ...
          '*.txt','Text Files (*.txt)'; ...
          '*.*','All Files (*.*)'};

[handles.datafile,handles.datapath] = uigetfile(filter,'Choose a traces file');

if handles.datafile==0, return; end
handles.filename = [handles.datapath handles.datafile];

% Load the file and fully initialize sorttraces.
OpenTracesFile( handles.filename, handles );


% guidata(hObject, handles);
% END FUNCTION btnOpen_Callback




%----------OPEN TRACES FROM FILE----------%
function handles = OpenTracesFile( filename, handles )
% Load a traces file, initialize GUI controls, and plot first trace.
% Called from btnOpen_Callback and sorttraces_OpeningFcn.
% Setting guidata is done there.


% Load the file
data = loadTraces( filename );

if isempty(data),
    warndlg('File is empty, so it cannot be loaded');
    return;
end

% Make sure time axis is in seconds (not frames)
if data.time(1)==1,
    f = inputdlg('What is the sampling interval (in ms) for this data?');
    if ~isempty(f),
        sampling = str2double(f);
        data.time = sampling.*(0:data.nFrames-1);
    end
end

% If this is multi-color FRET data and the metadata doesn't specify how FRET
% should be calculated, ask the user for clarification.
% NOTE: this is a temporary field and will be replaced by something more general
% in a future version.
if isChannel(data,'acceptor2') && ~isfield(data.fileMetadata,'isTandem3'),
    result = questdlg('Can you assume there is no donor->acceptor2 FRET?', ...
                    '3-color FRET calculation','Yes','No','Cancel','No');
    if strcmp(result,'Cancel'),  return;  end
    data.fileMetadata.isTandem3 = double( strcmp(result,'Yes') );
end


% Can't cancel after this point.
handles.filename = filename;
handles.data = data;


% Initialize array for tracking FRET donor-blinking threshold value.
% The default value of NaN is a marker that the value hasn't been
% calculated yet (but should be using traceStat).
% We don't calculate it here because then there would a long delay loading
% the file; small delays for each trace are not perceptible.
handles.fretThreshold = NaN( handles.data.nTraces, 1  );

% Set data correction starting values.
% The data is never modified. These parameters are used to adjust the data each
% time it is displayed or when it is ultimately saved. They are listed in the
% order they are applied to the trace. For now, background subtraction is the
% exception! The trace data is directly modified in that case.
% FIXME: intelligently choose size of these guys.
handles.adjusted   = false( handles.data.nTraces, 1 ); %bool if trace was changed
handles.background = zeros( handles.data.nTraces, 3 ); %not used (yet)
handles.gamma      = ones(  handles.data.nTraces, 3 );
handles.crosstalk  = zeros( handles.data.nTraces, 2 );

handles.binNames = {'No FRET', 'All FRET', 'Best FRET'};


% Trace indexes of binned molecules
handles.bins = cell(3,1);

% Check saved state file, including corrections and molecule selections.
% If no .mat file is found, look for .txt file (old version with only selections)
[p,fname] = fileparts( filename );
inds_fname = fullfile(p, [fname '_savedState.mat']);
if ~exist(inds_fname,'file'),
    inds_fname = fullfile(p, [fname '_picked_inds.txt']);
end

if exist(inds_fname,'file'),
    % FIXME: ask user whether to load corrections, selections, or both.
    text = 'These traces have been binned before.  ';
    text = [text 'Do you want to reload your selections?'];
    answer = questdlg(text, 'Load picking selections?','Yes','No','Yes');
    
    if strcmp(answer,'Yes'),
        handles = loadSavedState(handles, inds_fname);
    end
end

% Initialize picking boxes
for i=1:numel(handles.binNames),
    set( handles.(sprintf('editBin%d',i)), 'String',num2str(numel(handles.bins{i})) );
    set( handles.(sprintf('btnSelAll%d',i)), 'Enable','on' );
    set( handles.(sprintf('btnSelClear%d',i)), 'Enable','on' );
    set( handles.(sprintf('chkBin%d',i)), 'String',handles.binNames{i}, ...
         'Enable','on', 'Value', 0);
end

% Turn on other controls that can now be used now that a file is loaded.
if ismember('fret',handles.data.channelNames),
    isFret = 'on';
else
    isFret = 'off';
end
set(handles.edThreshold,  'Enable',isFret );
set(handles.sldThreshold, 'Enable',isFret );
set(handles.edCrosstalk1, 'Enable',isFret );
set(handles.sldCrosstalk1,'Enable',isFret );
set(handles.edGamma1,     'Enable',isFret,'String','1' );
set(handles.sldGamma1,    'Enable',isFret,'Value',1 );
set(handles.sldThreshold, 'min', 0, 'max', 200, 'sliderstep', [0.01 0.1] );

set(handles.btnPrint,    'Enable','on' );
set(handles.btnLoadDWT,  'Enable','on' );
set(handles.btnGettraces,'Enable','on' );

if isChannel(handles.data,'acceptor2') && ~isChannel(handles.data,'donor2')
    isThreeColor = 'on';
else
    isThreeColor = 'off';
end
set(handles.edCrosstalk2, 'Enable',isThreeColor,'String','0' );
set(handles.sldCrosstalk2,'Enable',isThreeColor,'Value',  0  );
set(handles.edGamma2,     'Enable',isThreeColor,'String','1' );
set(handles.sldGamma2,    'Enable',isThreeColor,'Value',  1  );

% Reset x-axis label to reflect time or frame-based.
time = handles.data.time;
if time(1)==1,
    xlabel( handles.axFret, 'Frame Number' );
else
    time = time/1000; %convert from ms to seconds
    xlabel( handles.axFret, 'Time (s)' );
end
xlim( handles.axFret, [time(1) time(end)] );

% Look for an corresponding idealization file and load it if found.
dwt_fname = fullfile(p, [fname '.qub.dwt']);
dwt_fname2 = fullfile(p, [fname '.dwt']);

if exist( dwt_fname, 'file' ),
    handles = loadDWT_ex( handles, dwt_fname );
elseif exist( dwt_fname2, 'file' ),
    handles = loadDWT_ex( handles, dwt_fname2 );
else
    handles.idl = [];
end


% Got to and plot first molecule.
% The "GoTo" callback ends up calling guidata() to save handles.
handles.molecule_no = 1;
set(handles.editGoTo,'Enable','on','String','1');
editGoTo_Callback(handles.editGoTo, [], handles);


% Add legends to the plotted traces.
% We want to do this once here because legend() is slow.
idxFluor = handles.data.idxFluor;

if sum(idxFluor)>1,
    legend( handles.axFluor, handles.data.channelNames{idxFluor} );
else
    legend( handles.axFluor, 'off' );
end

if isChannel( handles.data, 'fret2' ),
    legend( handles.axFret, {'fret1','fret2'} );
else
    legend( handles.axFret, 'off' );
end


% Adjust bottom axis, depending on the type of data.
if ismember('fret',handles.data.channelNames),
    if isa(data,'TracesFret4') && ~data.fileMetadata.isTandem3,
        ylabel(handles.axFret, 'Acceptor/Total');
    else
        ylabel(handles.axFret, 'FRET');
    end
    
    ylim(handles.axFret, [-0.1 1]);
else
    ylabel(handles.axFret, '');
    ylim(handles.axFret, 'auto');
end

% END FUNCTION OpenTracesFile




%----------LOAD SAVED STATE FROM FILE----------%
function handles = loadSavedState(handles, inds_fname)
% Load trace selections and adjustments from file.

[~,~,ext] = fileparts(inds_fname);

if strcmp('.mat',ext),
    savedState = load( inds_fname );
    requiredFields = {'bins','adjusted','crosstalk','gamma','background','fretThreshold'};

    if ~all( isfield(savedState,requiredFields) ),
        warning('Invalid sorttraces saved state file. Ignoring.');
        return;
    end

    % Load file settings into handles structure.
    for i=1:numel(requiredFields),
        handles.(requiredFields{i}) = savedState.(requiredFields{i});
    end

elseif strcmp('.txt',ext),
    % Load legacy saved state that only includes molecule selections.
    % These were made by versions 2.9 and earlier.
    fid = fopen(inds_fname,'r');
    for i=1:3,
        handles.bins{i} = sscanf( fgetl(fid), '%f' )';
    end
else
    warning('Invalid sorttraces saved state file. Ignoring.');
    return;
end

set(handles.btnSave,'Enable','on');
%TODO?: consider resetting scroll bars to fit range of loaded values.

%END FUNCTION loadSavedState




%=========================================================================%
%============================   NAVIGATION   =============================%


%----------GO TO MOLECULE----------%
function handles = editGoTo_Callback(hObject, ~, handles)
% Called when user changes the molecule number textbox. Jump to and plot
% the indicated trace.

% Get trace ID from GUI
mol=str2double( get(handles.editGoTo,'String') );

% If trace ID is invalid, reset it to what it was before.
if isnan(mol) || mol>handles.data.nTraces || mol<1,
    %disp('WARNING in sorttraces: Invalid trace number. Resetting.');
    set( hObject,'String',num2str(handles.molecule_no) );
    return;
else
    handles.molecule_no = mol;
end

% Make sure that the molecule selected actually exists.
if mol+1>handles.data.nTraces
    set(handles.btnNextBottom,'Enable','off');
else
    set(handles.btnNextBottom,'Enable','on');
end
if mol-1>=1
    set(handles.btnPrevBottom,'Enable','on');
else
    set(handles.btnPrevBottom,'Enable','off');
end

% Reset these values for the new trace.
fluorNames = handles.data.channelNames( handles.data.idxFluor );
handles.backgrounds = zeros( 1,numel(fluorNames) );

% If no value has been calculated for FRET threshold, do it now.
% FIXME: this calculation is duplicated from TracesFret*.
trace = handles.data.getSubset(mol);
handles.stats = traceStat(trace);
handles.trace = trace;
    
if trace.isChannel('fret') && isnan(handles.fretThreshold(mol)),
    total = zeros( size(trace) );
    
    % FIXME: this should only consider channels contributing to FRET (donor,
    % acceptor, acceptor2).
    for i=1:numel(fluorNames),
        total = total + trace.(fluorNames{i});
    end
    
    s = handles.stats.lifetime + 5;
    range = s:min(s+handles.constants.NBK,trace.nFrames);
    
    if numel(range)<10,
        handles.fretThreshold(mol) = 100; %arbitrary
    else
        handles.fretThreshold(mol) = handles.constants.blink_nstd * std(total(range));
    end
end

% Adjust scroll bar range if the new value falls outside of it.
sldMax = get( handles.sldThreshold, 'max' );
sldMax = max( sldMax, 2*handles.fretThreshold(mol) );
set( handles.sldThreshold, 'max', sldMax );

% Set bin checkboxes
for i=1:numel(handles.bins),
    chkName = sprintf('chkBin%d',i);
    set(handles.(chkName),'Value', any(handles.bins{i}==mol) );
end

% Re-initialize figure objects.
set( handles.edCrosstalk1,  'String', sprintf('%.3f',handles.crosstalk(mol,1)) );
set( handles.sldCrosstalk1, 'Value',  handles.crosstalk(mol,1) );

set( handles.edCrosstalk2,  'String', sprintf('%.3f',handles.crosstalk(mol,2)) );
set( handles.sldCrosstalk2, 'Value',  handles.crosstalk(mol,2) );

set( handles.edThreshold,  'String', sprintf('%.2f',handles.fretThreshold(mol)) );
set( handles.sldThreshold, 'Value',  handles.fretThreshold(mol) );

set( handles.edGamma1, 'String', sprintf('%.2f',handles.gamma(mol,2)) );
set( handles.sldGamma1, 'Value', handles.gamma(mol,2) );

set( handles.edGamma2, 'String', sprintf('%.2f',handles.gamma(mol,3)) );
set( handles.sldGamma2,'Value',  handles.gamma(mol,3) );

set(handles.btnSubDonor,    'Enable','on' );
set(handles.btnSubBoth,     'Enable','on' );
set(handles.btnSubAcceptor, 'Enable','on' );
set(handles.btnSubUndo,     'Enable','off');

guidata(hObject,handles);
plotter(handles);



%----------GO BACK TO PREVIOUS MOLECULE----------%
% --- Executes on button press in btnPrevTop.
function btnPrevTop_Callback(~, ~, handles)
% User clicked "previous molecule" button.
set( handles.editGoTo,'String',num2str(handles.molecule_no-1) );
editGoTo_Callback( handles.editGoTo, [], handles );



%----------GO TO NEXT MOLECULE----------%
% --- Executes on button press in btnNextTop - 'Next Molecule'.
function btnNextTop_Callback(~, ~, handles)
% User clicked "next molecule" button.
set( handles.editGoTo,'String',num2str(handles.molecule_no+1) );
editGoTo_Callback( handles.editGoTo, [], handles );




%=========================================================================%
%=========================   MOLECULE BINNING   ==========================%

% --- Executes on button press in chkBin1.
function addToBin_Callback(hObject, ~, handles, index)
% User clicked on one of the check boxes associated with each bin.
% The last parameter determines which bin was indicated.

mol = handles.molecule_no;
val = get(hObject,'Value');

if val==1,  %checking
    handles.bins{index} = [handles.bins{index} mol];
    
else  %unchecking
    handles.bins{index} = handles.bins{index}( handles.bins{index}~=mol );
end

% Update molecule counts for each bin in GUI.
for i=1:numel(handles.bins),
    edName = sprintf('editBin%d',i);
    set( handles.(edName), 'String',num2str(numel(handles.bins{i})) );
end

set(handles.btnSave,'Enable','on');
guidata(hObject,handles);


% --- Executes on button press in chkBin1.
function toggleBin(handles, index)
% User clicked on one of the check boxes associated with each bin.
% The last parameter determines which bin was indicated.

chkName = sprintf('chkBin%d',index);
val = ~get( handles.(chkName), 'Value' );
set( handles.(chkName), 'Value', val );

% Handle what normally happens after the box is checked.
addToBin_Callback(handles.(chkName),[],handles,index);




%----------SAVE TRACES----------%
% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, ~, handles)
% User clicked "Save Traces".
% Traces files are saved for each bin in which there are picked molecules.

[p,f]=fileparts(handles.filename);
baseFilename = fullfile(p,f);

%--- Save trace adjustments and indexes of picked molecules to file
% Saves a .mat file with a structure containing all the information necessary to
% recover the internal state of sorttraces as it is now.
savedState.version = handles.constants.version;

fields = {'binNames','bins','adjusted','crosstalk','gamma','background','fretThreshold'};
for i=1:numel(fields),
    savedState.(fields{i}) = handles.(fields{i});
end
save( [baseFilename '_savedState.mat'], '-struct', 'savedState' );


%--- Save files
binNames = lower( strrep(handles.binNames,' ','_') );  %FIXME: may need to remove other special characters.

for i=1:numel(handles.bins),
    filename = [baseFilename '_' binNames{i} '.traces'];
    savePickedTraces( handles, filename, handles.bins{i} );
end


% Finish up
set(hObject,'Enable','off');



function savePickedTraces( handles, filename, indexes )
% Save picked traces and idealizations to file.

% Sort indexes so they are in the same order as in the file, rather than in
% the order selected.
indexes = sort(indexes);

% If no traces remain, there is nothing to save. Rather than save an empty file,
% save nothing and delete any previously saved files.
if isempty(indexes),
    if exist(filename,'file'),
        delete(filename);
    end

    [p,f] = fileparts(filename);
    dwtfname = fullfile(p, [f '.qub.dwt']);
    if exist(dwtfname,'file'),
        delete(dwtfname);
    end

    return;
end

set(handles.figure1,'pointer','watch'); drawnow;

% Put together the subset of selected traces for saving, applying any
% adjustments made to the trace data.
data = adjustTraces(handles,indexes);  %creates a copy

% Save the trace data to file
[~,~,e] = fileparts(filename);
if strcmp(e,'.traces') || strcmp(e,'.rawtraces'),
    saveTraces( filename, data );
elseif strcmp(e,'.txt'),
    saveTraces( filename, 'txt', data );
else
    error('Unknown file format extension');
end

% Save idealizations of selected traces, if available.
% FIXME: dwt indexes may not match idl indexes if not all traces are
% idealized!!
if isfield(handles,'idl') && ~isempty(handles.idl),
    traceLen = numel(data.time);
    
    [p,f] = fileparts(filename);
    dwtFilename = fullfile(p, [f '.qub.dwt']);
    
    % Map DWT ID numbers to traces -- they may not be the same.
    offsets = (0:numel(indexes)-1)*traceLen;
    dwt_ids = handles.dwtIDs(indexes);
    
    % Remove traces from the DWT file without an idealization.
    offsets = offsets( dwt_ids>0 );
    dwt_ids = dwt_ids( dwt_ids>0 );
    
    saveDWT( dwtFilename, handles.dwt(dwt_ids), ...
             offsets, handles.dwtModel, handles.dwtSampling );
end

delete(data);  %clean up. not necessary, but fun.
set(handles.figure1,'pointer','arrow');

% end function savePickedTraces



% --- Executes on button press in btnSaveInPlace.
function btnSaveInPlace_Callback(~, ~, handles)
% Button to overwrite the current file with modifications, rather than
% saving results to selected traces.


% Ask the user for a target filename to save as, with the current file as the
% default.
[p,f,e] = fileparts(handles.filename);
filename = fullfile(p,[f '_adjusted' e]);

[inputfile,inputpath] = uiputfile('.traces','Save picked traces as:',filename);
if inputfile==0, return; end

filename=[inputpath inputfile];


% Put together the subset of selected traces for saving, applying any
% adjustments made to the trace data.
data = adjustTraces(handles);  %creates a copy

% Save the current data state to the file.
[~,~,e] = fileparts(filename);
if strcmp(e,'.traces') || strcmp(e,'.rawtraces'),
    saveTraces( filename, data );
elseif strcmp(e,'.txt'),
    saveTraces( filename, 'txt', data );
else
    error('Unknown file format extension');
end


% The idealization, if any, is NOT saved by this... FIXME?

set(handles.btnSaveInPlace,'Enable','off');


% END FUNCTION savePickedTraces




%=========================================================================%
%========================   TRACE CORRECTIONS   ==========================%


%----------HANDLE BACKGROUND SUBSTRACTION BUTTONS----------%
function btnSubBoth_Callback(hObject, ~, handles, mode)
% Subtract fluorescence background from the current x-axis region
% (presumably zoomed to a region after photobleaching). All of the
% subtraction buttons are handled with this one function. The mode
% parameter specifies which button was pressued and which function
% is to be performed.


m = handles.molecule_no; %selected molecule being viewed now.


% Get the current x-axis range to determine which area of the trace to use
% for background subtraction. Convert time (seconds) to frame number.
xlim = get(handles.axFluor,'XLim');
time = handles.data.time;

if time(1)~=1,  %time axis is in seconds if time(1)=0.
    dt = time(2)-time(1);
    xlim = floor(xlim./(dt/1000));
end

% Fix the x-axis range to be within the data range.
if xlim(1)<1, xlim(1)=1; end
if xlim(2)>handles.data.nFrames, xlim(2)=handles.data.nFrames; end
xrange = xlim(1):xlim(2);


% Get names/indexes for all of the fluorescence channels.
chNames = handles.data.channelNames;
idxFluor = cellfun( @isempty, strfind(chNames,'fret') );
fluorNames = chNames(idxFluor);

% Get the indexes of fluorescence channels the user wishes to subtract
idxSub = [];

if mode==1,      % Subtract DONOR background
    idxSub = find(  ~cellfun( @isempty, strfind(fluorNames,'donor') )  );
    
elseif mode==2,  % Subtract ACCEPTOR background
    idxSub = find(  ~cellfun( @isempty, strfind(fluorNames,'acceptor') )  );
    
elseif mode==3,  % Substrate background in ALL fluorescence channels
    idxSub = 1:numel(fluorNames);
end

% Subtract background from all selected channels
for i=1:numel(idxSub), %for every non-fret channel index
    ch = fluorNames{ idxSub(i) };
    bg = mean(  handles.data.(ch)(m,xrange)  );
    handles.backgrounds( idxSub(i) ) = bg;
    handles.data.(ch)(m,:) = handles.data.(ch)(m,:) - bg;
end
    
    
% UNDO background subtraction for all fluorescence channels
if mode==4,
    for i=1:numel(fluorNames), %for every non-fret channel
        ch = fluorNames{i};
        handles.data.(ch)(m,:) = handles.data.(ch)(m,:) + handles.backgrounds(i);
    end    
end


% Update GUI controls to allow undo if something has been changed.
if mode<4
    set(handles.btnSubUndo,'Enable','on');    %allow undo
else
    set(handles.btnSubUndo,'Enable','off');    %disable undo
end

handles = updateTraceData( handles );
guidata(hObject,handles);




%----------ADJUST CROSSTALK WITH SLIDER----------%
% --- Executes on slider movement.
function sldCrosstalk_Callback(hObject, ~, handles, ch )
% Called when user changes the scroll bar for specifying FRET threshold.
%

assert( ch==1 || ch==2 );
mol = handles.molecule_no;

% oldCrosstalk = handles.crosstalk(mol,ch);
handles.crosstalk(mol,ch) = get(hObject,'Value');
                         
% Save and display the result
name = sprintf('edCrosstalk%d',ch);
set( handles.(name), 'String',sprintf('%.3f',handles.crosstalk(mol,ch)) );

handles = updateTraceData( handles );
guidata(hObject,handles);



function edCrosstalk_Callback(hObject, eventdata, handles, ch )
% Called when user changes the text box for specifying FRET threshold.
%

assert( ch==1 || ch==2 );
crosstalk = str2double( get(hObject,'String') );

name = sprintf('sldCrosstalk%d',ch);  %slider control name

% Restrict value to the range of the slider to prevent errors.
% FIXME: changing the box should change the scale of the slider, within reason.
sldMax = get( handles.(name), 'max' );
sldMin = get( handles.(name), 'min' );
crosstalk = max(crosstalk,sldMin);
crosstalk = min(crosstalk,sldMax);

set( handles.(name), 'Value',crosstalk );

% Make the corrections and update plots.
sldCrosstalk_Callback( handles.(name), eventdata, handles, ch );



%---------------------------------------------------%

% --- Executes on slider movement.
function sldThreshold_Callback(hObject, ~, handles)
% Called when user changes the scroll bar for specifying FRET threshold.
%

mol = handles.molecule_no;

% Update edit box with new value
handles.fretThreshold(mol) = get(hObject,'Value');
set(handles.edThreshold,'String', sprintf('%.2f',handles.fretThreshold(mol)) );

% Plot data with new threshold value
handles = updateTraceData( handles );
guidata(hObject,handles);



function edThreshold_Callback(hObject, ~, handles)
% Called when user changes the text box for specifying FRET threshold.
%

mol = handles.molecule_no;
handles.fretThreshold(mol) = str2double( get(hObject,'String') );

% Verify the value is within the range of the slider bar. If not, it should
% be expanded to accomodate the new value or this will give an error.
sldMax = get( handles.sldThreshold,'max' );
sldMax = max( sldMax, 1.5*handles.fretThreshold(mol) );
set( handles.sldThreshold, 'max', sldMax );

% Update slider with new value
set( handles.sldThreshold, 'Value', handles.fretThreshold(mol) );

% Plot data with new threshold value
handles = updateTraceData( handles );
guidata(hObject,handles);




function edGamma_Callback(hObject, ~, handles)
% Text box for adjusting apparent gamma was changed. Scale acceptor1.
% A value of 1 means no adjustment. A value of 5 will multiply the acceptor
% by a factor of 5.

% Determine the acceptor channel number from control name
tag = get(hObject,'Tag');
ch = tag(end)-'0';
assert( ch==1 | ch==2, 'Invalid acceptor channel number' );
name = sprintf('sldGamma%d',ch);  %slider control name

gamma = str2double( get(hObject,'String') );

% Restrict value to the range of the slider to prevent errors.
% FIXME: changing the box should change the scale of the slider, within reason.
sldMax = get( handles.(name), 'max' );
sldMin = get( handles.(name), 'min' );
gamma = max(gamma,sldMin);
gamma = min(gamma,sldMax);

set( handles.(name), 'Value',gamma );

% Make the corrections and update plots.
sldGamma_Callback( handles.(name), [], handles );

% END FUNCTION edGamma1_Callback



% --- Executes on slider movement.
function sldGamma_Callback(hObject, ~, handles)
% 

% Determine the acceptor channel number from control name
tag = get(hObject,'Tag');
ch = tag(end)-'0';
assert( ch==1 | ch==2, 'Invalid acceptor channel number' );
name = sprintf('edGamma%d',ch);  %slider control name

mol = handles.molecule_no;
newGamma = get(hObject,'Value');
handles.gamma(mol,ch+1) = newGamma;

% Save and display the result
set( handles.(name), 'String',sprintf('%.2f',newGamma) );

handles = updateTraceData( handles );
guidata(hObject,handles);

% END FUNCTION sldGamma_Callback



% --- Executes on button press in btnResetAllCorrections.
function btnResetAllCorrections_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Clear trace corrections for ALL TRACES.
a = questdlg( 'This will reset corrections on ALL TRACES. Are you sure?', ...
              'Reset all corrections', 'OK','Cancel', 'OK' );

% Reset corrections variables to their defaults.
if strcmp(a,'OK'),
    handles.fretThreshold = NaN( handles.data.nTraces, 1  );
    handles.adjusted   = false( handles.data.nTraces, 1 );
    handles.background = zeros( handles.data.nTraces, 3 );
    handles.gamma      = ones(  handles.data.nTraces, 3 );
    handles.crosstalk  = zeros( handles.data.nTraces, 2 );
end

% Update GUI controls and redraw the trace.
guidata(hObject,handles);
editGoTo_Callback( hObject, [], handles );

% END FUNCTION btnResetAllCorrections_Callback



function handles = updateTraceData( handles )
% Recalculate FRET and stats. This is called following any changes made to
% the fluorescence traces (bg subtraction, crosstalk, etc) or the FRET
% threshold.

handles.adjusted(handles.molecule_no) = true;

% Since some kind of change was made, allow the user to save the file or
% selections in their new state.
if any( ~cellfun(@isempty,handles.bins) ),
    set(handles.btnSave,'Enable','on');
end
set(handles.btnSaveInPlace,'Enable','on');
plotter(handles);

% END FUNCTION updateTraceData




function displayData = adjustTraces( handles, indexes )
% Recalculate FRET and stats. This is called following any changes made to
% the fluorescence traces (bg subtraction, crosstalk, etc) or the FRET
% threshold.

if nargin<2,
    indexes = 1:handles.data.nTraces;
end

% Extract traces of interest, without modifying the input argument.
displayData = handles.data.getSubset(indexes);

% Make corrections for display. Crosstalk is only considered between
% neighboring channels.
chNames = displayData.channelNames(displayData.idxFluor);  %we assume these are in order of wavelength!!

for i=1:numel(indexes), 
    m = indexes(i); %i is index into data subset, m is index into all traces.
    
    if ~handles.adjusted(m), continue; end
    
    for j=2:numel(chNames),
        displayData.(chNames{j})(i,:) = displayData.(chNames{j})(i,:) - ...
                       handles.crosstalk(m,j-1)*displayData.(chNames{j-1})(i,:);
    end
    
    for j=1:numel(chNames),
        displayData.(chNames{j})(i,:) = displayData.(chNames{j})(i,:)*handles.gamma(m,j);
    end

end %for each trace
    
% Calculate total intensity and donor lifetime.
% FIXME: should thresholds be specified in metadata?
displayData.recalculateFret( handles.fretThreshold(indexes), handles.adjusted(indexes) );

% END FUNCTION adjustTraces





%=========================================================================%
%===========================   PLOT & PRINT   ============================%


%----------PRINT TRACE----------%
% --- Executes on button press in btnPrint.
function btnPrint_Callback(~, ~, handles)
% Opens the print system dialog to print the current trace.
printdlg(handles.figure1);



%----------PLOT TRACES----------%
function plotter(handles)
% Draw traces for current molecule and update displayed stats.

% Only "adjust" the data if some setting was changed. This way, the data appears
% exactly as it is in the file unless something is changed. This also makes just
% looking at traces somewhat faster.
m    = handles.molecule_no;

if handles.adjusted(m),
    data = adjustTraces(handles,m);
else
    data = handles.trace;
end
    

chNames = data.channelNames;


% If open, show the molecule location over the field image.
% Get the coordinates for all of the fluorescence channels.
if ishandle(handles.axFOV),
    traceMetadata = data.traceMetadata;

    % Verify that the currently-loaded movie matches the one current molecule.
    % They are not necessarily the same except for rawtraces files.
    % Assumes that files have properly formatted IDs. Not in old formats?
    output = split('#',traceMetadata.ids);
    [movieFilename,~] = deal( output{:} );
    
    if ~strcmp(handles.movieFilename,movieFilename),
        %btnGettraces_Callback(hObject, [], handles);
        close( get(handles.axFOV,'Parent') );
    
    else
        % Get x/y location of current molecule in each fluorescence channel.
        fields = fieldnames(traceMetadata);
        xs = find(  ~cellfun( @isempty, strfind(fields,'_x') )  );
        ys = find(  ~cellfun( @isempty, strfind(fields,'_y') )  );

        x = [];  y = [];
        for i=1:numel(xs),
            x = [ x ; traceMetadata.(fields{xs(i)}) ];
            y = [ y ; traceMetadata.(fields{ys(i)}) ];
        end
        disp([x y]);

        % Draw markers on selection points (total intensity composite image).
        % FIXME: try to draw a shape with scale dimensions so it gets bigger as
        % we zoom in. Otherwise it gets lost.
        axes(handles.axFOV);
        delete(findobj(handles.axFOV,'type','line'));
        line( x,y, 'LineStyle','none','marker','o','color','w' );

        % If available, draw a circle shape that scales with the image and is
        % easier to see. Requires 2014.
        if exist('viscircles','file'),
            viscircles( gca, [x y], repmat(3,numel(x),1), 'EdgeColor','w' );
        end

        figure(handles.figure1);  %return focus to main window.
    end
end


% Determine colors to user for plotting fluorescence.
fluorCh = chNames(handles.data.idxFluor);
nCh = numel(fluorCh);

if isfield(data.fileMetadata,'wavelengths'),
    chColors = Wavelength_to_RGB( data.fileMetadata.wavelengths );
elseif ismember('fret',data.channelNames),
    % For old FRET data (missing metadata), use the old standard colors.
    wavelengths = zeros(1,nCh);
    wavelengths( strcmp(chNames,'factor')    ) = 473;
    wavelengths( strcmp(chNames,'donor')     ) = 532;
    wavelengths( strcmp(chNames,'acceptor')  ) = 640;
    wavelengths( strcmp(chNames,'acceptor2') ) = 730;
    chColors = Wavelength_to_RGB(wavelengths);
else
    chColors = jet(nCh);  %color in order from blue to red as an approximation
end

% Get trace properties and reset GUI values with these results.
stats = handles.stats;
lt = stats.lifetime;
FRETlifetime = stats.acclife;

set(handles.editCorrelation,'String', sprintf('%.2f',stats.corr) );
set(handles.editSNR,'String', sprintf('%.1f, %.1f',stats.snr,stats.snr_s) );

if ismember('acceptor2',data.channelNames),  %isThreeColor,
    fret2Lifetime = stats.fret2Lifetime;
    set(handles.editLifetime,'String',  sprintf('%d, %d, %d', [FRETlifetime fret2Lifetime lt]));
else
    set(handles.editLifetime,'String',  sprintf('%d, %d', [FRETlifetime lt]));
end
set( handles.edZoomCorr, 'String','' );


[~,name,ext] = fileparts( handles.filename );
data_fname = [name ext];

% Plot fluorophore traces
time = data.time;
if time(1)~=1, %first time is 1 if in frame number (not ms)
    time = time/1000; %in seconds
end

cla( handles.axFluor );

for c=1:numel(fluorCh),
    trace = data.(fluorCh{c});
    plot( handles.axFluor, time,trace, 'Color',chColors(c,:) );
    
    if c==1
        total=trace;
    else
        total = total + trace;
    end
end

set( handles.txtTitle, 'String', [ 'Molecule ' num2str(m) ' of ' ...
                        num2str(handles.data.nTraces) ' of "' data_fname '"'] );
axis(handles.axFluor,'auto');

% Plot total fluorescence
cla( handles.axTotal );
plot( handles.axTotal, time,total,'k' );
axis(handles.axTotal,'auto');

% Draw lines representing donor (green) and acceptor (red) alive times
if ismember('fret',chNames),
    plot( handles.axTotal, time, repmat(handles.fretThreshold(m),1,data.nFrames), 'b-');
    
    mean_on_signal  = mean( total(1:lt) );
    mean_off_signal = mean( total(lt+5:end) );

    simplified_cy3=[mean_on_signal*ones(1,lt)...
            mean_off_signal*ones(1,data.nFrames-lt)];
    simplified_cy5=[mean_on_signal*ones(1,FRETlifetime)...
            mean_off_signal*ones(1,data.nFrames-FRETlifetime)];

    plot( handles.axTotal, time,simplified_cy5,'r' );
    plot( handles.axTotal, time,simplified_cy3,'g' );
end



% Plot FRET efficiency
cla( handles.axFret );
if ismember('fret',chNames),
    plot( handles.axFret, time,data.fret, 'b-');
end

if ismember('fret2',chNames),
    plot( handles.axFret, time,data.fret2, 'm-');
end

if isfield(handles,'idl') && ~isempty(handles.idl),
    stairs( handles.axFret, time, handles.idl(m,:), 'r-', 'LineWidth',1 );
end

xlim(handles.axFret, [time(1) time(end)]);
ylim(handles.axFret, [-0.1 1]);

drawnow;

% end function plotter.



function zoom_callback(hObject, ~)
% Called as ActionPostCallback for any of the axes objects upon zooming

handles = guidata(hObject);
m = handles.molecule_no;

if handles.adjusted(m),
    data = adjustTraces(handles,m);
else
    data = handles.trace;
end

% Get correlation over zoomed region.
if isChannel(data,'donor') && isChannel(data,'acceptor'),
    % Find the nearest datapoints in the zoomed region.
    lim = xlim(handles.axFluor);
    [~,idxLow]  = min( abs(data.time/1000-lim(1)) );
    [~,idxHigh] = min( abs(data.time/1000-lim(end)) );

    if idxLow==1 && idxHigh==data.nFrames,
        % Clear the box when fully zoomed out. This hints at the meaning of
        % the number -- it is only defined when zoomed in.
        set( handles.edZoomCorr, 'String','' );
    else        
        zcorr = corrcoef( data.donor(idxLow:idxHigh), data.acceptor(idxLow:idxHigh) );
        set( handles.edZoomCorr, 'String',sprintf('%.2f',zcorr(1,2)) );
    end
    
end;

%END FUNCTION




% --- Executes on button press in btnLoadDWT.
function btnLoadDWT_Callback(hObject, ~, handles)
% Loads an idealization (.dwt file) for later plotting in plotter().

handles = loadDWT_ex( handles );

% Save data to GUI
guidata(hObject,handles);

% Update GUI
plotter(handles);



function handles = loadDWT_ex( handles, filename)
% This function actually loads the dwell-time information.

time = handles.data.time;
sampling = time(2)-time(1);
nTraces  = handles.data.nTraces;
traceLen = handles.data.nFrames;

% Get filename for .dwt file from user and load it.
% FIXME: suggest a file name (or at least location) automatically).
if nargin>1,
    [dwt,dwtSampling,offsets,model] = loadDWT(filename);
else
    [dwt,dwtSampling,offsets,model] = loadDWT;
end

handles.dwt = dwt;
handles.dwtSampling = dwtSampling;
handles.dwtModel = model;

dwtSampling = double(dwtSampling);

handles.idl = [];

% Verify that the DWT matches the trace data.
if isempty(dwt), return; end

if sampling~=dwtSampling, 
    msgbox('Data and idealization sampling intervals do not match!','Error loading idealization','Error')

elseif offsets(end)>(nTraces*traceLen)
    msgbox('Idealization size does not match trace data.','Error loading idealization','Error');

else
    % Convert DWT to idealization (state sequence).
    [idl,handles.dwtIDs] = dwtToIdl( dwt, traceLen, offsets, nTraces );
    
    % Convert state sequence to idealized FRET trace.
    for i=1:size(idl,1),
        if iscell(model),
            m = model{i};
        else
            m = model;
        end
        
        fretValues = [NaN; m(:,1)];
        idl(i,:) = fretValues( idl(i,:)+1 );
    end
    
    % If the last few traces are not idealized, idl will be short.
    % Extend with NaN (no idealization marker) to avoid errors.
    handles.idl = [idl ; NaN(nTraces-size(idl,1), traceLen) ];
    
end %if errors

set(handles.btnClearIdl,'Enable','on');

% END FUNCTION loadDWT_ex



% --- Executes on button press in btnClearIdl.
function btnClearIdl_Callback(hObject, ~, handles)
% Clear the currently loaded idealization if any.

handles.dwt = [];
handles.idl = [];

% Save data and redraw the FRET plot without the idealization.
guidata(hObject,handles);
set(handles.btnClearIdl,'Enable','off');

plotter(handles);

% END FUNCTION btnClearIdl_Callback




% --- Executes on button press in btnSelAll.
function btnSelAll_Callback(hObject, ~, handles, index)
% User clicked the "select all" button above one of the bins.
% This is dangerous because all existing selections in that bin could be
% lost if this was accidental, so a warning dialog was added.

result = questdlg('Are you sure? All existing selections in this bin will be lost', ...
                            'Select all traces','OK','Cancel','Cancel');
if ~strcmp(result,'OK'),
    return;
end

handles.bins{index} = 1:handles.data.nTraces;
chkName = sprintf('chkBin%d',index);
set(handles.(chkName),'Value',1);

edName = sprintf('editBin%d',index);
set( handles.(edName), 'String',num2str(numel(handles.bins{index})) );

set(handles.btnSave,'Enable','on');
guidata(hObject,handles);

%end function btnSelAll_Callback



% --- Executes on button press in btnSelClear1.
function btnSelClear_Callback(hObject, ~, handles, index)   %#ok<DEFNU>
% User clicked the "clear selections" button above one of the bins.
% This is dangerous because all existing selections in that bin could be
% lost if this was accidental, so a warning dialog was added.

result = questdlg('This will clear ALL SELECTIONS in this bin. Are you sure?', ...
                            'Clear selections','OK','Cancel','Cancel');
if strcmp(result,'OK'),
    handles.bins{index} = [];
    chkName = sprintf('chkBin%d',index);
    set(handles.(chkName),'Value',0);

    edName = sprintf('editBin%d',index);
    set( handles.(edName), 'String',num2str(0) );

    set(handles.btnSave,'Enable','on');
    guidata(hObject,handles);
end

%end function btnSelClear_Callback





% --- Executes on button press in btnSelAll3.
function navKeyPress_Callback(hObject, eventdata, handles)
% Handles keyboard shortcut commands for moving through traces and putting
% them into bins. Called when keys are pressed when one of the navigation
% buttons has active focus.

ch = get(gcf,'CurrentCharacter');

switch ch
    case 28, %left arrow key
    btnPrevTop_Callback( hObject, eventdata, handles ); 
    
    case 29, %right arrow key
    btnNextTop_Callback( hObject, eventdata, handles );
    
    %case 30, %up arrow key
    %
    
    %case 31, %down arrow key
    %
    
    case 'a',
        toggleBin( handles, 1 );
    
    case 's',
        toggleBin( handles, 2 );
    
    case 'd',
        toggleBin( handles, 3 );
    
    case 'z',  %zoom in on trace
        dt  = diff(handles.data.time(1:2)) / 1000;  %time step in seconds
        lt  = handles.stats.lifetime*dt;
        lta = handles.stats.acclife*dt;
        x = xlim();  x = x(end);

        if x==lt+15*dt,
            % User already zoomed once; zoom in further.
            xlim( handles.axFluor, [0,lta+10*dt] );
        elseif x==lta+10*dt,
            % User already zoomed twice; zoom out.
            xlim( handles.axFluor, [0 handles.data.time(end)/1000] );
        else
            % Zoom in to show full trace.
            xlim( handles.axFluor, [0,lt+15*dt] );
        end
        
        zoom_callback(hObject,[]);
        
%     otherwise
%         disp( double(ch) );
end

%end function navKeyPress_Callback




% --- Executes on button press in btnGettraces.
function btnGettraces_Callback(hObject, ~, handles)
% Display an image of the field-of-view from the movie that the current trace
% came from and its physical location in each fluorescence channel. Iterating
% over traces will then update the molecule location.
%
% FIXME: this works correctly for showing the current trace, but we assume all
% traces are from the same movie, which is only the case for rawtraces files.
% Fixing this requres some work....


% Get the filename and movie coordinates of the selected trace.
m = handles.molecule_no;

if isfield(handles.data.traceMetadata,'ids'),
    id = handles.data.traceMetadata(m).ids;
else
    % TODO: disable the metadata button for visual feedback.
    disp('No gettraces metadata available for finding the original movie file.');
    return;
end

if any( id=='#' ),
    output = split('#',id);
    [movieFilename,~] = deal( output{:} );
    handles.movieFilename = movieFilename; %base name from IDs.
else
    disp('No gettraces metadata available. re-run gettraces!');
    return;
end

% Sometimes the traces filename is in the ID instead of the movie.
% Remove the file extension and add a guess.
[p,f,e] = fileparts(movieFilename);
if ~strcmp(e,'.stk') && ~strcmp(e,',tiff') && ~strcmp(e,'.tif'),
    movieFilename = fullfile(p,[f '.stk']);
end


% Find the movie data given in metadata, if it exists, plot an image from
% the first few frames, and indicate the location of the molecule.
% FIXME: need to save this image in the metadata rather than having to find
% the original movie file, which may not be around long.
% FIXME: automatically split the image up into each fluorescence channel.
if isempty(handles.axFOV) || ~ishandle(handles.axFOV),    
    
    % If the movie file doesn't exist, allow the user to look for it.
    if ~exist( movieFilename, 'file' ),
        idx = find( movieFilename=='\' | movieFilename=='/',1,'last' )+1;
        movieFilename = fullfile(pwd, movieFilename(idx:end));
        
        if ~exist( movieFilename, 'file' ),
            [f,p] = uigetfile( '*.stk', 'Manually find associated movie file', ...
                                    movieFilename );
            movieFilename = fullfile(p,f);

            % Verify the selected file exists.
            % FIXME: draw a blank background instead of cancelling.
            if ~ischar(f) || ~exist(movieFilename,'file'),
                return;
            end
        end
    end
    
    figure;
    handles.axFOV = gca;

    % Load colormap for image viewer
    fid=fopen('colortable.txt','r');
    colortable = fscanf(fid,'%d',[3 256]);
    colortable = colortable'/255;
    fclose(fid);

    stkData = gettraces( movieFilename );
    image_t = stkData.stk_top-stkData.background;

    % Display the field of view
    sort_px = sort(stkData.stk_top(:));
    val = sort_px( floor(0.98*numel(sort_px)) );
    imshow( image_t, [0 val], 'Parent',handles.axFOV );
    colormap(colortable);
    
    zoom on;
    guidata(hObject,handles);
end


% Show molecule location
plotter(handles);

% end function btnGettraces_Callback





% --- Executes when user attempts to close figure1.
function sorttraces_CloseRequestFcn(hObject, ~, handles)
% 

% Ask the user if the traces should be saved before closing.
if strcmpi(get(handles.btnSave,'Enable'),'on') || strcmpi(get(handles.btnSaveInPlace,'Enable'),'on')
    a = questdlg( 'Are you sure you want to exit? (You might lose changes)',...
                  'Exit without saving', 'OK','Cancel', 'OK' ); 
    if ~strcmp(a,'OK'),
        return;
    end
end

% Close the molecule location window if open
if ~isempty(handles.axFOV) && ishandle(handles.axFOV),
    close( get(handles.axFOV,'Parent') );
end

% Close the sorttraces window
delete(hObject);
