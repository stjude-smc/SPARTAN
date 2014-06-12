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

% Last Modified by GUIDE v2.5 26-May-2014 16:57:35


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
function sorttraces_OpeningFcn(hObject, eventdata, handles, varargin)
% Setup GUI controls in their default state (no file loaded).


% Initialize GUI if sorrtraces is being launched for the first time.
if ~isfield(handles,'constants')
    
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
    %ylim( handles.axFret, [-0.1 1] );
    grid( handles.axFret, 'on' );
    zoom( handles.axFret, 'on' );
    hold( handles.axFret, 'on' );
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



%--- END FUNCTIN sorttraces_OpeningFcn



% --- Outputs from this function are returned to the command line.
function varargout = sorttraces_OutputFcn(hObject, eventdata, handles)
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
function btnOpen_Callback(hObject, eventdata, handles)
% Open a user-selected traces file.

% Get traces filename by menu driven input
filter = {'*.traces;*.rawtraces','Binary Traces Files (*.traces,*.rawtraces)'; ...
          '*.txt','Text Files (*.txt)'; ...
          '*.*','All Files (*.*)'};

[handles.datafile,handles.datapath] = uigetfile(filter,'Choose a traces file');

if handles.datafile==0, return; end
handles.filename = [handles.datapath handles.datafile];

% Load the file and fully initialize sorttraces.
handles = OpenTracesFile( handles.filename, handles );


guidata(hObject, handles);
% END FUNCTION btnOpen_Callback




%----------OPEN TRACES FROM FILE----------%
function handles = OpenTracesFile( filename, handles )
% Load a traces file, initialize GUI controls, and plot first trace.
% Called from btnOpen_Callback and sorttraces_OpeningFcn.
% Setting guidata is done there.

handles.filename = filename;

% Load the file
handles.data = loadTraces( filename );

if isempty(handles.data),
    error('File is empty');
end

% Make sure time axis is in seconds (not frames)
if handles.data.time(1)==1,
    f = inputdlg('What is the sampling interval (in ms) for this data?');
    sampling = str2double(f)
    handles.data.time = sampling.*(0:handles.data.nFrames-1);
end


% Trace indexes of binned molecules
handles.NoFRETs_indexes=[];
handles.FRETs_indexes=[];
handles.Best_indexes=[];

% Check for previous molecule picking on same file
[p,fname] = fileparts( filename );
inds_fname = [p filesep fname '_picked_inds.txt'];

if exist(inds_fname,'file'),
    text = 'These traces have been binned before.  ';
    text = [text 'Do you want to reload your selections?'];
    answer = questdlg(text, 'Load picking selections?','Yes','No','No');
    
    if strcmp(answer,'Yes'),
        fid = fopen(inds_fname,'r');
        
        handles.NoFRETs_indexes = sscanf( fgetl(fid), '%f' )';
        handles.FRETs_indexes   = sscanf( fgetl(fid), '%f' )';
        handles.Best_indexes    = sscanf( fgetl(fid), '%f' )';
        
        set(handles.btnSave,'Enable','on');
    end
end

% Initialize picking boxes
set(handles.editBin1,'String', num2str(numel(handles.NoFRETs_indexes)) );
set(handles.editBin2,'String', num2str(numel(handles.FRETs_indexes))   );
set(handles.editBin3,'String', num2str(numel(handles.Best_indexes))    );
set(handles.chkBin1,'Enable','on', 'Value', 0);
set(handles.chkBin2,'Enable','on', 'Value', 0);
set(handles.chkBin3,'Enable','on', 'Value', 0);
set(handles.btnSelAll1,'Enable','on');
set(handles.btnSelAll2,'Enable','on');
set(handles.btnSelAll3,'Enable','on');

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

set(handles.btnPrint,    'Enable','on' );
set(handles.btnLoadDWT,  'Enable','on' );

if isChannel(handles.data,'acceptor2') && ~isChannel(handles.data,'donor2')
    isThreeColor = 'on';
else
    isThreeColor = 'off';
end
set(handles.edCrosstalk2, 'Enable',isThreeColor,'String','0' );
set(handles.sldCrosstalk2,'Enable',isThreeColor,'Value',0 );


% Initialize array for tracking FRET donor-blinking threshold value.
% The default value of zero is a marker that the value hasn't been
% calculated yet (but should be using traceStat).
% We don't calculate it here because then there would a long delay loading
% the file; small delays for each trace are not perceptible. 
handles.fretThreshold = zeros( handles.data.nTraces, 1  );
set( handles.sldThreshold, 'min', 0, 'max', 200, 'sliderstep', [0.01 0.1] );

% Set data correction starting values.
% The crosstalk value here reflects *the correction that has already been
% made* -- the actual data are modified each time
handles.crosstalk = zeros( handles.data.nTraces, 2  );


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
dwt_fname = [p filesep fname '.qub.dwt'];
dwt_fname2 = [p filesep fname '.dwt'];

if exist( dwt_fname, 'file' ),
    handles = loadDWT_ex( handles, dwt_fname );
elseif exist( dwt_fname2, 'file' ),
    handles = loadDWT_ex( handles, dwt_fname2 );
else
    handles.idl = [];
end


% Got to and plot first molecule.
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
    ylabel(handles.axFret, 'FRET');
    ylim(handles.axFret, [-0.1 1]);
else
    ylabel(handles.axFret, '');
    ylim(handles.axFret, 'auto');
end


% END FUNCTION OpenTracesFile




%=========================================================================%
%============================   NAVIGATION   =============================%


%----------GO TO MOLECULE----------%
function handles = editGoTo_Callback(hObject, eventdata, handles)
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
    set(handles.btnNextTop,'Enable','off');
    set(handles.btnNextBottom,'Enable','off');
else
    set(handles.btnNextTop,'Enable','on');
    set(handles.btnNextBottom,'Enable','on');
end
if mol-1>=1
    set(handles.btnPrevTop,'Enable','on');
    set(handles.btnPrevBottom,'Enable','on');
else
    set(handles.btnPrevTop,'Enable','off');
    set(handles.btnPrevBottom,'Enable','off');
end

% Reset these values for the new trace.
fluorNames = handles.data.channelNames( handles.data.idxFluor );
handles.backgrounds = zeros( 1,numel(fluorNames) );

% If no value has been calculated for FRET threshold, do it now.
trace = handles.data.getSubset(mol);
handles.stats = traceStat(trace);
    
if handles.data.isChannel('fret') && handles.fretThreshold(mol) == 0,
    total = zeros( size(trace) );
    
    for i=1:numel(fluorNames),
        total = total + trace.(fluorNames{i});
    end
    
    constants = cascadeConstants;
    s = handles.stats.lifetime + 5;
    range = s:min(s+constants.NBK,handles.data.nFrames);
    
    if numel(range)<10,
        handles.fretThreshold(mol) = 100; %arbitrary
    else
        handles.fretThreshold(mol) = constants.blink_nstd * std(total(range));
    end
    
    % Adjust scroll bar range if the new value falls outside of it.
    sldMax = get( handles.sldThreshold, 'max' );
    sldMax = max( sldMax, 2*handles.fretThreshold(mol) );
    set( handles.sldThreshold, 'max', sldMax );
end

% Set bin checkboxes
set(handles.chkBin1,'Value', any(handles.NoFRETs_indexes==mol) );
set(handles.chkBin2,'Value', any(handles.FRETs_indexes==mol) );
set(handles.chkBin3,'Value', any(handles.Best_indexes==mol) );

% Re-initialize figure objects.
set( handles.edCrosstalk1,  'String', sprintf('%.3f',handles.crosstalk(mol,1)) );
set( handles.sldCrosstalk1, 'Value',  handles.crosstalk(mol,1) );

set( handles.edCrosstalk2,  'String', sprintf('%.3f',handles.crosstalk(mol,2)) );
set( handles.sldCrosstalk2, 'Value',  handles.crosstalk(mol,2) );

set( handles.edThreshold,  'String', sprintf('%.2f',handles.fretThreshold(mol)) );
set( handles.sldThreshold, 'Value',  handles.fretThreshold(mol) );

set(handles.btnSubDonor,    'Enable','on' );
set(handles.btnSubBoth,     'Enable','on' );
set(handles.btnSubAcceptor, 'Enable','on' );
set(handles.btnSubUndo,     'Enable','off');

guidata(hObject,handles);
plotter(handles);



%----------GO BACK TO PREVIOUS MOLECULE----------%
% --- Executes on button press in btnPrevTop.
function btnPrevTop_Callback(hObject, eventdata, handles)
% User clicked "previous molecule" button.
set( handles.editGoTo,'String',num2str(handles.molecule_no-1) );
editGoTo_Callback( handles.editGoTo, [], handles );



%----------GO TO NEXT MOLECULE----------%
% --- Executes on button press in btnNextTop - 'Next Molecule'.
function btnNextTop_Callback(hObject, eventdata, handles)
% User clicked "next molecule" button.
set( handles.editGoTo,'String',num2str(handles.molecule_no+1) );
editGoTo_Callback( handles.editGoTo, [], handles );




%=========================================================================%
%=========================   MOLECULE BINNING   ==========================%

% --- Executes on button press in chkBin1.
function addToBin_Callback(hObject, eventdata, handles, index)
% User clicked on one of the check boxes associated with each bin.
% The last parameter determines which bin was indicated.

mol = handles.molecule_no;
val = get(hObject,'Value');

if val==1,  %checking
    if index==1,
        handles.NoFRETs_indexes = [handles.NoFRETs_indexes mol];
    elseif index==2,
        handles.FRETs_indexes = [handles.FRETs_indexes mol];
    elseif index==3,
        handles.Best_indexes = [handles.Best_indexes mol];
    end
    
else  %unchecking
    if index==1,
        handles.NoFRETs_indexes = handles.NoFRETs_indexes( ...
                                  handles.NoFRETs_indexes~=mol );
    elseif index==2,
        handles.FRETs_indexes = handles.FRETs_indexes( ...
                                handles.FRETs_indexes~=mol );       
    elseif index==3,
        handles.Best_indexes = handles.Best_indexes( ...
                               handles.Best_indexes~=mol );
    end
end


NoFRET_no = numel( handles.NoFRETs_indexes );
FRET_no = numel( handles.FRETs_indexes );
Best_no = numel( handles.Best_indexes );

set(handles.editBin1,'String',num2str(NoFRET_no));
set(handles.editBin2,'String',num2str(FRET_no));  
set(handles.editBin3,'String',num2str(Best_no));

set(handles.btnSave,'Enable','on');
guidata(hObject,handles);


% --- Executes on button press in chkBin1.
function toggleBin(hObject, eventdata, handles, index)
% User clicked on one of the check boxes associated with each bin.
% The last parameter determines which bin was indicated.

mol = handles.molecule_no;
chk_name = ['chkBin' num2str(index)];
val = ~get( handles.(chk_name), 'Value' );
set( handles.(chk_name), 'Value', val );

if val==1,  %if unchecked, check
    if index==1,
        handles.NoFRETs_indexes = [handles.NoFRETs_indexes mol];
    elseif index==2,
        handles.FRETs_indexes = [handles.FRETs_indexes mol];
    elseif index==3,
        handles.Best_indexes = [handles.Best_indexes mol];
    end
    
else  %if checked, uncheck
    if index==1,
        handles.NoFRETs_indexes = handles.NoFRETs_indexes( ...
                                  handles.NoFRETs_indexes~=mol );
    elseif index==2,
        handles.FRETs_indexes = handles.FRETs_indexes( ...
                                handles.FRETs_indexes~=mol );       
    elseif index==3,
        handles.Best_indexes = handles.Best_indexes( ...
                               handles.Best_indexes~=mol );
    end
end


NoFRET_no = numel( handles.NoFRETs_indexes );
FRET_no = numel( handles.FRETs_indexes );
Best_no = numel( handles.Best_indexes );

set(handles.editBin1,'String',num2str(NoFRET_no));
set(handles.editBin2,'String',num2str(FRET_no));  
set(handles.editBin3,'String',num2str(Best_no));

set(handles.btnSave,'Enable','on');
guidata(hObject,handles);




%----------SAVE TRACES----------%
% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% User clicked "Save Traces".
% Traces files are saved for each bin in which there are picked molecules.


%--- Save indexes of picked molecules to file
[p,f]=fileparts(handles.filename);
baseFilename = [p filesep f];

fid = fopen( [baseFilename '_picked_inds.txt'], 'w' );
fprintf(fid, '%d ', handles.NoFRETs_indexes); fprintf(fid,'\n');
fprintf(fid, '%d ', handles.FRETs_indexes);   fprintf(fid,'\n');
fprintf(fid, '%d ', handles.Best_indexes);    fprintf(fid,'\n');
fclose(fid);


%--- Save files
if ~isempty(handles.Best_indexes),
    filename = [baseFilename '_best_fret.traces'];
    savePickedTraces( handles, filename, handles.Best_indexes );
end

if ~isempty(handles.FRETs_indexes),
    filename = [baseFilename '_all_fret.traces'];
    savePickedTraces( handles, filename, handles.FRETs_indexes );
end

if ~isempty(handles.NoFRETs_indexes),
    filename = [baseFilename '_no_fret.traces'];
    savePickedTraces( handles, filename, handles.NoFRETs_indexes );
end


% Finish up
set(hObject,'Enable','off');



function savePickedTraces( handles, filename, indexes )
% Save picked traces and idealizations to file.

% Sort indexes so they are in the same order as in the file, rather than in
% the order selected.
indexes = sort(indexes);

% Put together the subset of selected traces for saving.
data = handles.data.getSubset(indexes);  %create a copy

[p,f,e] = fileparts(filename);
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
    dwtFilename = [p filesep f '.qub.dwt'];
    
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

% end function savePickedTraces





%=========================================================================%
%========================   TRACE CORRECTIONS   ==========================%


%----------HANDLE BACKGROUND SUBSTRACTION BUTTONS----------%
function btnSubBoth_Callback(hObject, eventdata, handles, mode)
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
plotter(handles);




%----------ADJUST CROSSTALK WITH SLIDER----------%
% --- Executes on slider movement.
function sldCrosstalk_Callback(hObject, eventdata, handles, ch )
% Called when user changes the scroll bar for specifying FRET threshold.
%

assert( ch==1 || ch==2 );
mol = handles.molecule_no;

oldCrosstalk = handles.crosstalk(mol,ch);
handles.crosstalk(mol,ch) = get(hObject,'Value');
delta = oldCrosstalk - handles.crosstalk(mol,ch);

% Adjust the acceptor fluorescence to subtract donor->acceptor crosstalk
% according to the new value. The trace has already been adjusted according
% to the old value, so that has to be "undone" first.
chNames = handles.data.channelNames;  %we assume these are in order of wavelength!!
ch1 = chNames{ch};
ch2 = chNames{ch+1};

handles.data.(ch2)(mol,:) = handles.data.(ch2)(mol,:) + ...
                                          delta * handles.data.(ch1)(mol,:);
                         
% Save and display the result
name = sprintf('edCrosstalk%d',ch);
set( handles.(name), 'String',sprintf('%.3f',handles.crosstalk(mol,ch)) );

handles = updateTraceData( handles );
guidata(hObject,handles);
plotter(handles);



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
function sldThreshold_Callback(hObject, eventdata, handles)
% Called when user changes the scroll bar for specifying FRET threshold.
%

mol = handles.molecule_no;

% Update edit box with new value
handles.fretThreshold(mol) = get(hObject,'Value');
set(handles.edThreshold,'String', sprintf('%.2f',handles.fretThreshold(mol)) );

% Plot data with new threshold value
handles = updateTraceData( handles );
guidata(hObject,handles);
plotter(handles);



function edThreshold_Callback(hObject, eventdata, handles)
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
plotter(handles);





function handles = updateTraceData( handles )
% Recalculate FRET and stats. This is called following any changes made to
% the fluorescence traces (bg subtraction, crosstalk, etc) or the FRET
% threshold.

if ~isChannel(handles.data,'fret'),
    return;
end

m        = handles.molecule_no;
donor    = handles.data.donor(m,:);
acceptor = handles.data.acceptor(m,:);

% Calculate total intensity and donor lifetime.
if isChannel(handles.data,'acceptor2') && ~isChannel(handles.data,'donor2'),
    isThreeColor = true;
    acceptor2 = handles.data.acceptor2(m,:);
    total = donor+acceptor+acceptor2;
else
    isThreeColor = false;
    total = donor+acceptor;
end

lt = max(1, calcLifetime(total) );

% Calculate FRET
fret = acceptor./total;
fret( total<handles.fretThreshold(m) ) = 0;
fret( lt:end ) = 0;
handles.data.fret(m,:) = fret;

if isThreeColor,
    fret2 = acceptor2./total;
    fret2( total<handles.fretThreshold(m) ) = 0;
    fret2( lt:end ) = 0;
    handles.data.fret2(m,:) = fret2;
end

% Recalculate stats.
handles.stats = traceStat( handles.data.getSubset(m) );


% END FUNCTION updateTraceData



%=========================================================================%
%===========================   PLOT & PRINT   ============================%


%----------PRINT TRACE----------%
% --- Executes on button press in btnPrint.
function btnPrint_Callback(hObject, eventdata, handles)
% Opens the print system dialog to print the current trace.
printdlg(handles.figure1);



%----------PLOT TRACES----------%
function plotter(handles)
% Draw traces for current molecule and update displayed stats.

m     = handles.molecule_no;
chNames = handles.data.channelNames;
fluorCh = chNames(handles.data.idxFluor);
nCh = numel(fluorCh);

% Determine colors to user for plotting fluorescence.
if isfield(handles.data.fileMetadata,'wavelengths'),
    chColors = Wavelength_to_RGB( handles.data.fileMetadata.wavelengths );
elseif ismember('fret',handles.data.channelNames),
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

% Get trace properties and reset GUI values with these results
stats = handles.stats;
lt = stats.lifetime;
FRETlifetime = stats.acclife;
snr = stats.snr;
CC = stats.corr;

set(handles.editCorrelation,'String', sprintf('%.2f',CC) );
set(handles.editSNR,'String', sprintf('%.2f',snr) );

if ismember('acceptor2',handles.data.channelNames),  %isThreeColor,
    fret2Lifetime = stats.fret2Lifetime;
    set(handles.editLifetime,'String',  sprintf('%d, %d, %d', [FRETlifetime fret2Lifetime lt]));
else
    set(handles.editLifetime,'String',  sprintf('%d, %d', [FRETlifetime lt]));
end


[p,name,ext] = fileparts( handles.filename );
data_fname = [name ext];

% Plot fluorophore traces
time = handles.data.time;
if time(1)~=1, %first time is 1 if in frame number (not ms)
    time = time/1000; %in seconds
end

cla( handles.axFluor );

for c=1:numel(fluorCh),
    trace = handles.data.(fluorCh{c})(m,:);
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
    plot( handles.axTotal, time, repmat(handles.fretThreshold(m),1,handles.data.nFrames), 'b-');
    
    mean_on_signal  = mean( total(1:lt) );
    mean_off_signal = mean( total(lt+5:end) );

    simplified_cy3=[mean_on_signal*ones(1,lt)...
            mean_off_signal*ones(1,handles.data.nFrames-lt)];
    simplified_cy5=[mean_on_signal*ones(1,FRETlifetime)...
            mean_off_signal*ones(1,handles.data.nFrames-FRETlifetime)];

    plot( handles.axTotal, time,simplified_cy5,'r' );
    plot( handles.axTotal, time,simplified_cy3,'g' );
end



% Plot FRET efficiency
cla( handles.axFret );
if ismember('fret',chNames),
    plot( handles.axFret, time,handles.data.fret(m,:), 'b-');
end

if ismember('fret2',chNames),
    plot( handles.axFret, time,handles.data.fret2(m,:), 'm-');
end

if isfield(handles,'idl') && ~isempty(handles.idl),
    stairs( handles.axFret, time, handles.idl(m,:), 'r-', 'LineWidth',1 );
end


xlim([time(1) time(end)]);



% end function plotter.







% --- Executes on button press in btnLoadDWT.
function btnLoadDWT_Callback(hObject, eventdata, handles)
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
    [idl,handles.dwtIDs] = dwtToIdl( dwt, traceLen, offsets );
    
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








% --- Executes on button press in btnSelAll3.
function btnSelAll_Callback(hObject, eventdata, handles, index)
% User clicked the "select all" button above one of the bins.
% This is dangerous because all existing selections in that bin could be
% lost if this was accidental, so a warning dialog was added.

result = questdlg('Are you sure? All existing selections in this bin will be lost', ...
                            'Select all traces','OK','Cancel','Cancel');
if ~strcmp(result,'OK'),
    return;
end

if index==1,
    handles.NoFRETs_indexes = 1:handles.data.nTraces;
    set(handles.chkBin1,'Value',1);
elseif index==2,
    handles.FRETs_indexes = 1:handles.data.nTraces;
    set(handles.chkBin2,'Value',1);
elseif index==3,
    handles.Best_indexes = 1:handles.data.nTraces;
    set(handles.chkBin3,'Value',1);
end

NoFRET_no = numel( handles.NoFRETs_indexes );
FRET_no = numel( handles.FRETs_indexes );
Best_no = numel( handles.Best_indexes );

set(handles.editBin1,'String',num2str(NoFRET_no));
set(handles.editBin2,'String',num2str(FRET_no));  
set(handles.editBin3,'String',num2str(Best_no));

set(handles.btnSave,'Enable','on');
guidata(hObject,handles);




% --- Executes on button press in btnSelAll3.
function navKeyPress_Callback(hObject, eventdata, handles)
% Handles keyboard shortcut commands for moving through traces and putting
% them into bins. Called when keys are pressed when one of the navigation
% buttons has active focus.
% FIXME: not working for top buttons yet.

ch = get(gcf,'CurrentCharacter');

switch ch
    case 28, %left arrow key
    btnPrevTop_Callback( hObject, eventdata, handles ); 
    
    case 29, %right arrow key
    btnNextTop_Callback( hObject, eventdata, handles );
    
    case 'a',
    toggleBin( hObject, eventdata, handles, 1 );
    
    case 's',
    toggleBin( hObject, eventdata, handles, 2 );
    
    case 'd',
    toggleBin( hObject, eventdata, handles, 3 );
end


%end function navKeyPress_Callback




% --- Executes on button press in btnGettraces.
function btnGettraces_Callback(hObject, eventdata, handles)
%

% Get the filename and movie coordinates of the selected trace.
% FIXME: this assumes new format traces!
m = handles.molecule_no;
id = handles.data.traceMetadata(m).ids;

if any( id=='#' ),
    output = split('#',id);
    [movieFilename,traceID] = deal( output{:} );
else
    warning('No gettraces metadata available. re-run gettraces!');
    return;
end


% Verify file actually exists in the specified location. If not, given a
% warning and try to find it in the current location.
if ~exist( movieFilename, 'file' ),
    warning('Movie file specified in trace metadata doesn''t exist!');
        
    [p,f,e] = fileparts(movieFilename);
    altFilename = [pwd filesep f e];
    if exist( altFilename, 'file' ),
        movieFilename = altFilename;
    else        
        disp( ['Unable to find associated movie file: ' movieFilename] );
        disp( 'Please find the associated movie file manually.' )
        [f,p] = uigetfile( '*.stk', 'Manually find associated movie file', [pwd filesep f e] );
        movieFilename = [p f];
        
        % Verify the selected file exists.
        if ~ischar(f) || ~exist(movieFilename,'file'),
            return;
        end
    end
end

gettraces_gui( movieFilename, handles.data.traceMetadata(m) );


% end function btnGettraces_Callback
