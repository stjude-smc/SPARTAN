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

%   Copyright 2007-2016 Cornell University All Rights Reserved.

% Last Modified by GUIDE v2.5 05-May-2016 12:10:15


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

updateSpartan; %check for updates

disp('Keyboard shortcuts:');
disp('   Left arrow key  = Previous molecule');
disp('   Right arrow key = Next molecule');
disp('   a, s, d = Put molecule in No, All, or Best FRET, respectively');
disp('   z       = Zoom in to trace. Press again to zoom out');    

% Choose default command line output for sorttraces
handles.output = hObject;
handles.vals = [];

constants = cascadeConstants();
set( handles.figure1, 'Name', [mfilename ' - ' constants.software] );

% Setup FRET threshold line dragging callbacks
set(handles.figure1, 'WindowButtonMotionFcn',@(a,b)sorttraces('threshDrag_Callback',a,b,guidata(a)), ...
                     'WindowButtonUpFcn',@(a,b)sorttraces('threshUp_Callback',a,b,guidata(a)));
handles.threshDrag = false;

% Link x-axes - zooming on one plot will automatically zoom on the other
ax = [handles.axFluor handles.axTotal handles.axFret];
linkaxes(ax,'x');

% SETUP AXES labels and settings.
% Hold is needed so we don't lose the grid/zoom/etc settings, which are
% lost when a new plot is generated on the axes...
for i=1:numel(ax),
    zoom( ax(i), 'on' );
    hold( ax(i), 'on' );
    set( zoom(ax(i)), 'ActionPostCallback',@zoom_callback);
end

ylabel( handles.axFluor, 'Fluorescence' );
ylabel( handles.axTotal, 'Total Fluorescence' );
[handles.axFOV, handles.idl, handles.idlFret] = deal([]);
guidata(hObject, handles);

% If called by autotrace, 2nd argument is filename of traces output file.
if numel(varargin) > 1
    handles = btnOpen_Callback( [], [], handles, varargin{2} );
        
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
varargout{1} = handles.output;
%--- END FUNCTION sorttraces_OutputFcn




%=========================================================================%
%=========================   LOAD TRACES FILES   =========================%

%----------"OPEN TRACES FILE" Button----------%
function handles = btnOpen_Callback(~, ~, handles, filename)
% Open a user-selected traces file.

if nargin<4,
    filename = getFile;
end
if isempty(filename), return; end
set(handles.figure1,'pointer','watch'); drawnow;

% Load the file
data = loadTraces( filename );

if isempty(data),
    warndlg('File is empty, so it cannot be loaded');
    set(handles.figure1,'pointer','arrow'); drawnow;
    return;
end

if ~isfield(data.fileMetadata,'zeroMethod'),
    data.fileMetadata.zeroMethod = 'threshold';
end
useThresh = strcmpi(data.fileMetadata.zeroMethod,'threshold');


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
handles.adjusted = false( data.nTraces, 1 ); %bool if trace was changed

if isa(data,'TracesFret4') && isChannel(data,'acceptor2'),
    choiceChanged = data.verifyFretGeometry();
    handles.adjusted(:) = choiceChanged;
end


handles.filename = filename;
handles.data = data;

% Set default correction values, which correspond to no change.
if useThresh,
    handles.fretThreshold = NaN( data.nTraces, 1  );
else
    handles.fretThreshold = zeros( data.nTraces, 1  );
end
    

nFluor = 3;  %numel(data.idxFluor);
handles.background = zeros( data.nTraces, nFluor ); %not used (yet)
handles.crosstalk  = zeros( data.nTraces, nFluor, nFluor );
handles.gamma      = ones(  data.nTraces, nFluor );


% Trace indexes of binned molecules
handles.binNames = {'No FRET', 'All FRET', 'Best FRET'};
handles.bins = cell(3,1);

% Check saved state file, including corrections and molecule selections.
% If no .mat file is found, look for .txt file (old version with only selections)
[p,fname] = fileparts( filename );
inds_fname = fullfile(p, [fname '_savedState.mat']);
if ~exist(inds_fname,'file'),
    inds_fname = fullfile(p, [fname '_picked_inds.txt']);  %legacy
end

if exist(inds_fname,'file'),
    answer = questdlg('Load previously saved state?', mfilename, 'Yes','No', 'Yes');
    
    if strcmp(answer,'Yes'),
        handles = loadSavedState(handles, inds_fname);
    end
end

% Initialize picking boxes
set(handles.mnuBin, 'Enable','on');
for i=1:numel(handles.binNames),
    set( handles.(sprintf('editBin%d',i)), 'String',num2str(numel(handles.bins{i})) );
    set( handles.(sprintf('chkBin%d',i)), 'String',handles.binNames{i}, ...
         'Enable','on', 'Value', 0);
end

% Turn on other controls that can now be used now that a file is loaded.
set( [handles.edCrosstalk1 handles.sldCrosstalk1 handles.edGamma1 handles.sldGamma1], ...
      'Enable', onoff(ismember('fret',data.channelNames)) );

set( [handles.edCrosstalk2  handles.sldCrosstalk2 handles.edCrosstalk3 ...
      handles.sldCrosstalk3 handles.edGamma2      handles.sldGamma2], ...
      'Enable', onoff(ismember('fret2',data.channelNames)) );

set( [handles.tbLoadIdl handles.mnuLoadIdl handles.tbGettraces ...
      handles.mnuGettraces handles.mnuZeroMethod], 'Enable','on');

% Reset x-axis label to reflect time or frame-based.
time = data.time;
if time(1)==1,
    xlabel( handles.axFret, 'Frame Number' );
else
    time = time/1000; %convert from ms to seconds
    xlabel( handles.axFret, 'Time (s)' );
end
xlim( handles.axFret, [time(1) time(end)] );

% Look for an corresponding idealization file and load it if found.
try
    dwt_fname = findDwt({filename});
    handles = loadDWT_ex( handles, dwt_fname{1} );
catch
    handles.idl = [];
    handles.idlFret = [];
end

% Molecule location
% FIXME: assumes donor is at the origin (upper left).
if isfield(data.traceMetadata,'donor_x'),
    set(handles.axLocation, 'Visible','on');
    axis(handles.axLocation,'square');  %FIXME
    if isfield(data.fileMetadata,'fieldSize')
        szMovie = data.fileMetadata.movieSize;
        axis(handles.axLocation, [0 szMovie(1) 0 szMovie(2)]);
    else
        xlim(handles.axLocation, [0 max([data.traceMetadata.donor_x])] );
        ylim(handles.axLocation, [0 max([data.traceMetadata.donor_y])] );
    end
else
    cla(handles.axLocation);
    set(handles.axLocation, 'Visible','off');
end

% Got to and plot first molecule.
% The "GoTo" callback ends up calling guidata() to save handles.
handles.molecule_no = 1;
set(handles.editGoTo,'Enable','on','String','1');
editGoTo_Callback(handles.editGoTo, [], handles);


% Add legends to the plotted traces.
% We want to do this once here because legend() is slow.
if numel(data.idxFluor)>1,
    legend( handles.axFluor, data.channelNames{data.idxFluor} );
else
    legend( handles.axFluor, 'off' );
end

if numel(data.idxFret)>1,
    legend( handles.axFret, data.channelNames{data.idxFret} );
else
    legend( handles.axFret, 'off' );
end


% Adjust bottom axis, depending on the type of data.
if ismember('fret',data.channelNames),
    ylabel(handles.axFret, data.fretAxisLabel);
else
    ylabel(handles.axFret, '');
    ylim(handles.axFret, 'auto');
end
zoom reset;  %remember new axis limits when zooming out.

% Enable buttons
set([handles.mnuSaveAs handles.mnuExportText handles.mnuSubDonor ...
     handles.mnuSubAcceptor handles.mnuSubBoth handles.mnuResetBG ...
     handles.mnuCorrResetAll handles.btnSubBoth], 'Enable','on');
set(handles.figure1,'pointer','arrow'); drawnow;


% END FUNCTION OpenTracesFile




%----------LOAD SAVED STATE FROM FILE----------%
function handles = loadSavedState(handles, inds_fname)
% Load trace selections and adjustments from file.

[~,~,ext] = fileparts(inds_fname);

if strcmp('.mat',ext),
    savedState = load( inds_fname );
    
    % Verify that the savedState.mat file is valid.
    requiredFields = {'bins','adjusted','crosstalk','gamma','background','fretThreshold'};
    if ~all( isfield(savedState,requiredFields) ),
        warning('Invalid sorttraces saved state file. Ignoring.');
        return;
    end
    
    nTraces = numel(savedState.adjusted);
    if nTraces~=handles.data.nTraces,
        warning('Number of traces in saved state does not match loaded .traces file');
        return;
    end
    
    % Handle older files that do not have the full crosstalk matrix.
    if ndims(savedState.crosstalk)<3, %#ok<ISMAT>
        nFluor = 3; %numel(handles.data.idxFluor);
        crosstalk = zeros(nTraces, nFluor, nFluor);
        crosstalk(:,1,2) = savedState.crosstalk(:,1);
        crosstalk(:,2,3) = savedState.crosstalk(:,2);
        savedState.crosstalk = crosstalk;
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

set([handles.tbSave handles.mnuSave], 'Enable','on');

%END FUNCTION loadSavedState




%=========================================================================%
%============================   NAVIGATION   =============================%


%----------GO TO MOLECULE----------%
function handles = editGoTo_Callback(hObject, ~, handles)
% Called when user changes the molecule number textbox. Jump to and plot
% the indicated trace.
constants = cascadeConstants;

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

% Update molecule number display. FIXME: use two separate controls.
[~,name,ext] = fileparts(handles.filename);
data_fname = strrep([name ext],'_',' '); %avoid LaTeX interpretation
if numel(data_fname)>35,
    data_fname = ['...' data_fname(numel(data_fname)-34:end)];
end
titleTxt = sprintf('Molecule %d of %d in %s',mol, handles.data.nTraces, data_fname);
title(handles.axFluor,titleTxt);

% Make sure that the molecule selected actually exists.
set(handles.btnNextBottom, 'Enable', onoff(mol+1<=handles.data.nTraces) );
set(handles.btnPrevBottom, 'Enable', onoff(mol-1>=1) );

% Reset these values for the new trace.
fluorNames = handles.data.channelNames( handles.data.idxFluor );
handles.backgrounds = zeros( 1,numel(fluorNames) );

% If no value has been calculated for FRET threshold, do it now.
% FIXME: this calculation is duplicated from TracesFret*.
trace = handles.data.getSubset(mol);
handles.stats = traceStat(trace);
handles.trace = trace;

if strcmpi(handles.data.fileMetadata.zeroMethod,'threshold'),
    if trace.isChannel('fret') && isnan(handles.fretThreshold(mol)),
        total = zeros( size(trace) );

        % FIXME: this should only consider channels contributing to FRET (donor,
        % acceptor, acceptor2).
        for i=1:numel(fluorNames),
            total = total + trace.(fluorNames{i});
        end

        s = handles.stats.lifetime + 5;
        range = s:min(s+constants.NBK,trace.nFrames);

        if numel(range)<10,
            handles.fretThreshold(mol) = 150; %arbitrary
        else
            handles.fretThreshold(mol) = constants.blink_nstd * std(total(range));
        end
    end
end

% Set bin checkboxes
for i=1:numel(handles.bins),
    chkName = sprintf('chkBin%d',i);
    set(handles.(chkName),'Value', any(handles.bins{i}==mol) );
end

% Set correction controls to saved values for this trace.
from = [1 2 1];
to   = [2 3 3];

for i=1:numel(from),
    crosstalk = handles.crosstalk(mol,from(i),to(i));
    name = sprintf('edCrosstalk%d',i);
    set( handles.(name),  'String', sprintf('%.3f',crosstalk) );
    name = sprintf('sldCrosstalk%d',i);
    set( handles.(name), 'Value', crosstalk );
end

set( handles.edGamma1, 'String', sprintf('%.2f',handles.gamma(mol,2)) );
set( handles.sldGamma1, 'Value', handles.gamma(mol,2) );

set( handles.edGamma2, 'String', sprintf('%.2f',handles.gamma(mol,3)) );
set( handles.sldGamma2,'Value',  handles.gamma(mol,3) );

handles = plotter(handles);
guidata(hObject,handles);

% END FUNCTION editGoTo_Callback


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

set([handles.tbSave handles.mnuSave], 'Enable','on');
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
function btnSave_Callback(~, ~, handles) %#ok<DEFNU>
% Save selected traces in a new file for each bin.

[p,f]=fileparts(handles.filename);
baseFilename = fullfile(p,f);

% Save the current state of the GUI (selections, adjustments, etc)
saveState(handles);

% Save selected traces in one file per bin
binNames = lower( strrep(handles.binNames,' ','_') );  %FIXME: may need to remove other special characters.

for i=1:numel(handles.bins),
    filename = [baseFilename '_' binNames{i} '.traces'];
    savePickedTraces( handles, filename, handles.bins{i} );
end

set([handles.tbSave handles.mnuSave], 'Enable','off');

% END FUNCTION btnSave_Callback




% --- Executes on button press in btnSaveInPlace.
function btnSaveInPlace_Callback(~, ~, handles) %#ok<DEFNU>
% Save all traces to a new file with all corrections applied.

% Save the current state of the GUI (selections, adjustments, etc)
saveState(handles);

% Ask the user for a target filename to save as
[p,f] = fileparts(handles.filename);
filename = fullfile(p,[f '_adjusted.traces']);

[f,p] = uiputfile('.traces','Save current file as:',filename);
if f==0, return; end

% Apply all corrections and save traces and idealization to file.
savePickedTraces( handles, fullfile(p,f), 1:handles.data.nTraces );

% END FUNCTION savePickedTraces




function savePickedTraces( handles, filename, indexes )
% Save picked traces and idealizations to file.

[p,f] = fileparts(filename);
baseFilename = fullfile(p,f);
dwtfname = [baseFilename '.qub.dwt'];

% If no traces remain, there is nothing to save. Rather than save an empty file,
% save nothing and delete any previously saved files.
if isempty(indexes),
    if exist(filename,'file'),  delete(filename);  end
    if exist(dwtfname,'file'),  delete(dwtfname);  end
    return;
end

set(handles.figure1,'pointer','watch'); drawnow;

% Select traces, apply all corrections, and save to file.
indexes = sort(indexes);  %ensure traces are in original order
saveTraces( filename, adjustTraces(handles,indexes) );

% Save idealizations of selected traces, if available.
if ~isempty(handles.idl),
    [dwt,offsets] = idlToDwt( handles.idl(indexes,:) );
    saveDWT( dwtfname, dwt, offsets, handles.dwtModel, handles.data.sampling );
end

set(handles.figure1,'pointer','arrow');

% END FUNCTION savePickedTraces




function saveState( handles )
% Save a .mat containing the current GUI state to be recovered later.

constants = cascadeConstants;
savedState.version = constants.version;

fields = {'binNames','bins','adjusted','crosstalk','gamma','background','fretThreshold'};
for i=1:numel(fields),
    savedState.(fields{i}) = handles.(fields{i});
end

[p,f] = fileparts(handles.filename);
filename = fullfile(p,[f '_savedState.mat']);
save(filename, '-struct', 'savedState');

% END FUNCTION saveState




%=========================================================================%
%========================   TRACE CORRECTIONS   ==========================%


%----------HANDLE BACKGROUND SUBSTRACTION BUTTONS----------%
function btnSubBoth_Callback(~, ~, handles, mode) %#ok<DEFNU>
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

updateTraceData( handles );

% END FUNCTION sldCrosstalk_Callback




%----------ADJUST CROSSTALK WITH SLIDER----------%
% --- Executes on slider movement.
function sldCrosstalk_Callback(hObject, ~, handles )
% Called when user changes the scroll bar for specifying FRET threshold.
%
tag = get(hObject,'Tag');
ch = tag(end)-'0';

from = [1 2 1];
to   = [2 3 3];

mol = handles.molecule_no;
crosstalk = get(hObject,'Value');
handles.crosstalk(mol, from(ch), to(ch)) = crosstalk;
                         
% Save and display the result
name = sprintf('edCrosstalk%d',ch);
set( handles.(name), 'String',sprintf('%.3f',crosstalk) );

updateTraceData( handles );

% END FUNCTION sldCrosstalk_Callback



function edCrosstalk_Callback(hObject, ~, handles ) %#ok<DEFNU>
% Called when user changes the text box for specifying FRET threshold.
%
tag = get(hObject,'Tag');
ch = tag(end)-'0';

crosstalk = str2double( get(hObject,'String') );
name = sprintf('sldCrosstalk%d',ch);  %slider control name

% Adjust slider range to contain the new value.
sldMax = get( handles.(name), 'max' );
sldMin = get( handles.(name), 'min' );

if crosstalk<sldMin,
    set( handles.(name), 'min', crosstalk );
elseif crosstalk>sldMax
    set( handles.(name), 'max', crosstalk );
end

set( handles.(name), 'Value',crosstalk );

% Make the corrections and update plots.
sldCrosstalk_Callback( handles.(name), [], handles );



%---------------------------------------------------%
function edGamma_Callback(hObject, ~, handles) %#ok<DEFNU>
% Text box for adjusting apparent gamma was changed. Scale acceptor1.
% A value of 1 means no adjustment. A value of 5 will multiply the acceptor
% by a factor of 5.

% Determine the acceptor channel number from control name
tag = get(hObject,'Tag');
ch = tag(end)-'0';
assert( ch==1 | ch==2, 'Invalid acceptor channel number' );
name = sprintf('sldGamma%d',ch);  %slider control name

gamma = str2double( get(hObject,'String') );

% Adjust slider range to contain the new value.
sldMax = get( handles.(name), 'max' );
sldMin = get( handles.(name), 'min' );
gamma = max(gamma,sldMin);
if gamma>sldMax,
    set( handles.(name), 'max', gamma );
end

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

updateTraceData( handles );

% END FUNCTION sldGamma_Callback



% --- Executes on button press in btnResetAllCorrections.
function btnResetAllCorrections_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Clear trace corrections for ALL TRACES.
a = questdlg( 'This will reset corrections on ALL TRACES. Are you sure?', ...
              'Reset all corrections', 'OK','Cancel', 'OK' );

% Reset corrections variables to their defaults.
if strcmp(a,'OK'),
    handles.fretThreshold = NaN( size(handles.fretThreshold)  );
    handles.adjusted   = false( size(handles.adjusted)   );
    handles.background = zeros( size(handles.background) );
    handles.crosstalk  = zeros( size(handles.crosstalk)  );
    handles.gamma      = ones(  size(handles.gamma)      );
end

% Update GUI controls and redraw the trace.
editGoTo_Callback( hObject, [], handles );

% END FUNCTION btnResetAllCorrections_Callback




function handles = updateTraceData( handles )
% Recalculate FRET, update plots, and save GUI data (handles).
% This is called following any changes made to the fluorescence traces
% (bg subtraction, crosstalk, etc) or the FRET threshold.

handles.adjusted(handles.molecule_no) = true;

set([handles.mnuSave handles.tbSave], 'Enable','on');
handles = plotter(handles);
guidata(handles.figure1, handles);

% END FUNCTION updateTraceData




function output = adjustTraces( handles, indexes )
% DATA = adjustTraces( handles, INDEXES ) returns a new Traces object DATA
% applying all adjustments (background, gamma, crosstalk, and FRET threshold)
% to the data loaded from file (handles.data).
% Traces that were not adjusted (handles.adjusted is false) are unaltered.
% If trace INDEXES are given, DATA only includes those traces.

% Extract traces of interest
if nargin<2,
    indexes = 1:handles.data.nTraces;
end
output = handles.data.getSubset(indexes);  %creates a copy

% Determine which traces have adjustments that must be applied.
idxAdjusted = to_row( find(handles.adjusted(indexes)) );  %ones to adjust
if numel(idxAdjusted)==0, return; end


% Make all adjustments to raw data
crosstalk   = handles.crosstalk(indexes,:,:);
gamma       = handles.gamma(indexes,:);
thresholds  = handles.fretThreshold(indexes);

chNames = output.channelNames(output.idxFluor);  %increasing wavelength order.
nFlour = numel(output.idxFluor);

for i=idxAdjusted,
    % Subtract forward crosstalk (blue to red). The order matters.
    for src=1:nFlour,
        for dst=1:nFlour,
            if src>=dst, continue; end  %only consider forward crosstalk
            
            ch1 = chNames{src};
            ch2 = chNames{dst};
            output.(ch2)(i,:) = output.(ch2)(i,:) - crosstalk(i,src,dst)*output.(ch1)(i,:);
        end
    end
    
    % Scale channels so their relative brightness is the same.
    for ch=1:nFlour,
        output.(chNames{ch})(i,:) = output.(chNames{ch})(i,:)*gamma(i,ch);
    end

end %for each selected, adjusted trace


% Recalculate FRET from the corrected data, removing donor dark states.
if isfield(output.fileMetadata,'zeroMethod') && strcmpi(output.fileMetadata.zeroMethod,'threshold'),
    output.recalculateFret( idxAdjusted, thresholds(idxAdjusted) );
else
    output.recalculateFret( idxAdjusted );
end


% END FUNCTION adjustTraces





%=========================================================================%
%===========================   PLOT & PRINT   ============================%


%----------PRINT TRACE----------%
% --- Executes on button press in btnPrint.
function btnPrint_Callback(~, ~, handles) %#ok<DEFNU>
% Opens the print system dialog to print the current trace.
printdlg(handles.figure1);



%----------PLOT TRACES----------%
function handles = plotter(handles)
% Draw traces for current molecule and update displayed stats.


% Adjust data if some setting was changed.
m = handles.molecule_no;

if handles.adjusted(m),
    data = adjustTraces(handles,m);
else
    data = handles.trace;
end

chNames = data.channelNames;
fluorCh = chNames(data.idxFluor);
nCh = numel(fluorCh);


% Show the molecule location within field of view, if the window is open.
if ishandle(handles.axFOV),
    output = strsplit(data.traceMetadata.ids,'#');
    
    if strcmp(handles.movieFilename, output{1}),
        % Get x/y location of current molecule in each fluorescence channel.
        x = nan(nCh,1);  y = nan(nCh,1);
        for i=1:nCh,
            if isfield(data.traceMetadata, [fluorCh{i} '_x']) && ...
                    isfield(data.traceMetadata, [fluorCh{i} '_y']),
                x(i) = data.traceMetadata.([fluorCh{i} '_x']);
                y(i) = data.traceMetadata.([fluorCh{i} '_y']);
            end
        end
        
        % Draw markers on selection points.
        viscircles( handles.axFOV, [x y], repmat(3,numel(x),1), 'EdgeColor','w' );
    else
        % Close if current molecule is not from the loaded movie.
        % FIXME: should instead load the new movie and continue seamlessly.
        close( get(handles.axFOV,'Parent') );
    end
end

% Draw molecule location in small panel
if strcmp( get(handles.axLocation, 'Visible'), 'on'),
    cla(handles.axLocation);
    loc = [data.traceMetadata.donor_x data.traceMetadata.donor_y];
    line( loc(1),loc(2), 'MarkerSize',4, 'LineStyle','none', 'Marker','o',...
                         'Color','r', 'Parent',handles.axLocation );
end


% Determine colors for plotting fluorescence.
if isfield(data.fileMetadata,'wavelengths'),
    chColors = Wavelength_to_RGB( data.fileMetadata.wavelengths );
else
    chColors = [0 1 0; 1 0 0; 0.5 0 0];
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


% Plot fluorophore traces
time = data.time;
if time(1)~=1, %first time is 1 if in frame number (not ms)
    time = time/1000; %in seconds
end

cla( handles.axFluor );

for c=1:numel(fluorCh),
    trace = data.(fluorCh{c});
    plot( handles.axFluor, time,trace, 'Color',chColors(c,:) );
end
total = data.total; %Excludes 'factor' and misc channels

% Plot total fluorescence
cla( handles.axTotal );
plot( handles.axTotal, time,total,'k' );
axis([handles.axTotal handles.axFluor],'auto y');

% Draw lines representing donor (green) and acceptor (red) alive times
if ismember('fret',chNames),
    mean_on_signal = mean( total(data.fret~=0) );
    mean_off_signal = mean( total(lt+5:end) );

    simplified_cy3 = mean_on_signal*(data.fret~=0);
    simplified_cy3(data.fret==0) = mean_off_signal;
        
    simplified_cy5=[mean_on_signal*ones(1,FRETlifetime)...
            mean_off_signal*ones(1,data.nFrames-FRETlifetime)];
    
    plot( handles.axTotal, time,simplified_cy5,'r' );
    plot( handles.axTotal, time,simplified_cy3,'g' );
    
    % Show adjustable FRET threshold line, if applicable.
    if strcmpi(data.fileMetadata.zeroMethod,'threshold')
        handles.hThresh = plot( handles.axTotal, time([1 data.nFrames]), ...
                repmat(handles.fretThreshold(m),1,2), 'b-', ...
                'ButtonDownFcn',@(a,b)sorttraces('threshDown_Callback',a,b,guidata(a)));
    end
end


% Plot FRET efficiency
cla( handles.axFret );
if ismember('fret',chNames),
    plot( handles.axFret, time,data.fret, 'b-');
end

if ismember('fret2',chNames),
    plot( handles.axFret, time,data.fret2, 'm-');
end

if ~isempty(handles.idlFret),
    stairs( handles.axFret, time, handles.idlFret(m,:), 'r-', 'LineWidth',1 );
end

xlim(handles.axFret, [time(1) time(end)]);
ylim(handles.axFret, [-0.1 1]);

drawnow;

% end function plotter.


function zoom_callback(hObject, ~)
% Called as ActionPostCallback for any of the axes objects upon zooming.
% Updates any stats (correlation) calculated over the zoomed region.

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
    
end

%END FUNCTION



%=========================================================================%
%===================   FRET THRESHOLD LINE DRAGGING   ====================%
function threshDown_Callback(hObject, ~, handles) %#ok<DEFNU>
% Called when the user STARTS dragging the FRET threshold line in axTotal.
% This allows the user to the manually set the threshold.
handles.threshDrag = true;
guidata(hObject,handles);

function threshUp_Callback(hObject, ~, handles) %#ok<DEFNU>
% Called when the user STOPS dragging the FRET threshold line in axTotal.
if handles.threshDrag,
    cp = get(handles.axTotal,'CurrentPoint');
    handles.fretThreshold(handles.molecule_no) = cp(1,2);
    handles.threshDrag = false;
    guidata(hObject,handles);
    
    handles = updateTraceData(handles);
    guidata(hObject,handles);
end

function threshDrag_Callback(~, ~, handles) %#ok<DEFNU>
% Called repeatedly while the FRET threshold line in axTotal is dragged.
% Allows the user to the manually set the threshold.
if handles.threshDrag,
    cp = get(handles.axTotal,'CurrentPoint');
    set( handles.hThresh, 'YData', [cp(1,2) cp(1,2)] );
end



%=========================================================================%
%==================   IDEALIZATION BUTTON CALLBACKS   ====================%
function btnLoadDWT_Callback(hObject, ~, handles) %#ok<DEFNU>
% Loads an idealization (.dwt file) for later plotting in plotter().

handles = loadDWT_ex( handles );
handles = plotter(handles);
guidata(hObject,handles);

% END FUNCTION btnLoadDWT_Callback


function handles = loadDWT_ex( handles, varargin)
% This function actually loads the dwell-time information.

time = handles.data.time;
sampling = time(2)-time(1);
nTraces  = handles.data.nTraces;
traceLen = handles.data.nFrames;

% Get filename for .dwt file from user and load it.
% FIXME: suggest a file name (or at least location) automatically).
[dwt,dwtSampling,offsets,model] = loadDWT( varargin{:} );
handles.dwtModel = model;
handles.idl = [];

% Verify that the DWT matches the trace data.
if isempty(dwt), return; end

if sampling~=double(dwtSampling), 
    msgbox('Data and idealization sampling intervals do not match!','Error loading idealization','Error')

elseif offsets(end)>(nTraces*traceLen)
    msgbox('Idealization size does not match trace data.','Error loading idealization','Error');

else
    % Convert DWT to idealization (state sequence).
    idl = dwtToIdl( dwt, offsets, traceLen, nTraces );
    handles.idl = idl;
    
    % Convert state sequence to idealized FRET trace.
    for i=1:nTraces
        if iscell(model),
            m = model{i};
        else
            m = model;
        end
        
        fretValues = [NaN; m(:,1)];
        idl(i,:) = fretValues( idl(i,:)+1 );
    end
    handles.idlFret = idl;
    
end %if errors

set(handles.mnuClearIdl,'Enable','on');

% END FUNCTION loadDWT_ex



% --- Executes on button press in btnClearIdl.
function btnClearIdl_Callback(hObject, ~, handles) %#ok<DEFNU>
% Clear the currently loaded idealization if any.

handles.idl = [];
handles.idlFret = [];
set(handles.mnuClearIdl,'Enable','off');

handles = plotter(handles);
guidata(hObject,handles);

% END FUNCTION btnClearIdl_Callback
% ----------------------------------------------------------------------




% --- Executes on button press in btnSelAll.
function btnSelAll_Callback(hObject, ~, handles, index) %#ok<DEFNU>
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

set([handles.tbSave handles.mnuSave],'Enable','on');
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

    set([handles.tbSave handles.mnuSave], 'Enable','on');
    guidata(hObject,handles);
end

%end function btnSelClear_Callback





% --- Executes on button press in btnSelAll3.
function navKeyPress_Callback(hObject, eventdata, handles) %#ok<DEFNU>
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
function btnGettraces_Callback(hObject, ~, handles) %#ok<DEFNU>
% Display an image of the field-of-view from the movie that the current trace
% came from and its physical location in each fluorescence channel. Iterating
% over traces will then update the molecule location.


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
    output = strsplit(id,'#');
    [movieFilename,~] = deal( output{:} );
    handles.movieFilename = movieFilename; %base name from IDs.
else
    disp('No gettraces metadata available. re-run gettraces!');
    return;
end

% Handle IDs that use the trace file name instead of the movie.
[p,f,e] = fileparts(movieFilename);
if ~strcmp(e,'.stk') && ~strcmp(e,',tiff') && ~strcmp(e,'.tif'),
    e = '.tif';
    movieFilename = fullfile(p,[f e]);
end


% Find the movie data given in metadata, if it exists, plot an image from
% the first few frames, and indicate the location of the molecule.
if isempty(handles.axFOV) || ~ishandle(handles.axFOV),    
    
    % If the movie file doesn't exist, allow the user to look for it.
    if ~exist( movieFilename, 'file' ),
        % First look in the current directory
        movieFilename = fullfile(pwd, [f e]);
        
        if ~exist( movieFilename, 'file' ),
            [f,p] = uigetfile( '*.tif;*.tiff;*.stk', 'Manually find associated movie file', ...
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

    % Load an image from the first 10 frames of the movie.
    stkData = gettraces( movieFilename );
    image_t = stkData.stk_top-stkData.background;

    % Display the field of view
    sort_px = sort(stkData.stk_top(:));
    val = sort_px( floor(0.98*numel(sort_px)) );
    imshow( image_t, [0 val], 'Parent',handles.axFOV );
    colormap(handles.axFOV, gettraces_colormap);
    
    zoom(handles.axFOV,'on');
end


% Show molecule location
handles = plotter(handles);
guidata(hObject,handles);

% end function btnGettraces_Callback




% --- Executes when user attempts to close figure1.
function sorttraces_CloseRequestFcn(hObject, ~, handles) %#ok<DEFNU>
% Give the user a chance to save current state before closing.

if strcmpi(get(handles.mnuSave,'Enable'),'on') && ...
   strcmpi(get(handles.mnuSaveAs,'Enable'),'on'),

    a = questdlg('Save current state before exiting?', 'Close sorttraces', ...
                 'Yes','No','Cancel', 'Yes' );
    switch a
        case 'Yes'
            saveState(handles);
        case 'Cancel'
            return; %do not close yet
    end
end

% Close the molecule location window if open
if ~isempty(handles.axFOV) && ishandle(handles.axFOV),
    close( get(handles.axFOV,'Parent') );
end

% Close the sorttraces window
delete(hObject);

% END FUNCTION sorttraces_CloseRequestFcn



% --------------------------------------------------------------------
function mnuZeroMethod_Callback(hObject, ~, handles) %#ok<DEFNU>
% Set method for detecting donor blinks (setting FRET to zero).

zeroMethod = handles.data.fileMetadata.zeroMethod;
current = find( strcmpi(zeroMethod,TracesFret.zeroMethodNames) );

[sel,ok] = listdlg('PromptString','Method to detect donor blinking and set FRET to zero:', ...
                   'SelectionMode','single','ListSize',[300 120], ...
                   'ListString',TracesFret.zeroMethodDesc, ...
                   'InitialValue',current );

if ok && sel~=current,
    set(handles.figure1,'pointer','watch'); drawnow;
    
    % Recalculate FRET with the new mode.
    handles.data.fileMetadata.zeroMethod = TracesFret.zeroMethodNames{sel};
    handles.data.recalculateFret();
    
    % Update GUI controls and redraw current trace.
    if sel==2,
        handles.fretThreshold(:) = NaN;
    else
        handles.fretThreshold(:) = 0;
    end
    guidata(hObject,handles);
    
    editGoTo_Callback(hObject, [], handles);
    set(handles.figure1,'pointer','arrow'); drawnow;
end

% END FUNCTION mnuZeroMethod_Callback
