function varargout = sorttraces(varargin)

% SORTTRACES M-file for sorttraces.fig
%      SORTTRACES takes *.traces files as input, and allows the user to
%      correct the background and crosstalk, and to place each trace into
%      different bins. Three text files are outputed: "_all_fret",
%      "_no_fret", and "_best_fret", all with extension .txt, and all space
%      delimited.
%
%      The figure window was created using GUIDE. Sorttraces.fig must be
%      present in the same directory as the m-file. 3/2006 -JBM.
%
%      Total fluorescence axis, threshold, and correlation coefficient
%      calculator added. Correlation coefficient is printed to the
%      same file as lifetime and SNR. 7/2006 -JBM.
%
%      Ability to open files from gettraces, and text files from autotrace 
%      with unique identifiers was added. A bug in the background
%      correction was fixed too. 12/2006 -JBM.

% Depends on: sorttraces.fig, LoadTraces.m, CorrectTraces, cascadeConstants,
%    trace_stat (which requires: RLE_filter, CalcLifetime)

% Last Modified by GUIDE v2.5 20-Aug-2012 16:10:29


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sorttraces (see VARARGIN)

firstRun = 0;

% Initialize GUI if sorrtraces is being launched for the first time.
if ~isfield(handles,'constants')
    firstRun = 1;
    
    % Choose default command line output for sorttraces
    handles.output = hObject;
    handles.vals = [];

    % Initialize some variables, utilizing the handles structure
    handles.constants = cascadeConstants();
    handles.defaultFretThreshold=3000;
    handles.default_crosstalk=0;
    set(handles.sldCrosstalk, 'Value',  handles.default_crosstalk);
    set(handles.sldThreshold, 'Value',  handles.defaultFretThreshold);
    set(handles.edCrosstalk,  'String', '0');

    % Link x-axes - zooming on one plot will automatically zoom on the other
    linkaxes([handles.axFluor handles.axTotal handles.axFret],'x');

    warning off MATLAB:divideByZero;
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
    
elseif ~firstRun
    disp('Sorttraces is already running!');
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
%======================   LOAD TRACES FILES   ======================%


%----------"OPEN TRACES FILE" Button----------%
function btnOpen_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get traces filename by menu driven input
[handles.datafile,handles.datapath]=...
    uigetfile({'*.traces;*.txt'},'Choose a traces file');
if handles.datafile==0, return; end
handles.filename=[handles.datapath handles.datafile];

% Load the file and initialize sorttraces
handles = OpenTracesFile( handles.filename, handles );


guidata(hObject, handles);
% END FUNCTION btnOpen_Callback




%----------OPEN TRACES FROM FILE----------%
% Called from btnOpen_Callback and sorttraces_OpeningFcn.
% Setting guidata is done there.
function handles = OpenTracesFile( filename, handles )

handles.filename = filename;

% Load the file
data = loadTraces( filename );
handles.donor    = data.donor;
handles.acceptor = data.acceptor;
handles.fret     = data.fret;
handles.time     = data.time;
handles.traceMetadata = data.traceMetadata;

[handles.Ntraces,handles.len] = size(handles.donor);

if size(handles.donor,1)<1,
    error('File is empty');
end

% Make sure time axis is in seconds (not frames)
if handles.time(1)==1,
    f = inputdlg('What is the sampling interval (in ms) for this data?');
    sampling = str2double(f)
    handles.time = sampling.*(0:handles.len-1);
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
        
        handles.NoFRETs_indexes = str2num( fgetl(fid) );
        handles.FRETs_indexes   = str2num( fgetl(fid) );
        handles.Best_indexes    = str2num( fgetl(fid) );
        
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


% Set data correction starting values
handles.fretThreshold = handles.defaultFretThreshold;
handles.crosstalk = handles.default_crosstalk;

% Initialize arrays for background subtraction.
handles.donor_background = 0;
handles.acceptor_background = 0;

% Initialize arrays for tracking FRET donor-blinking threshold value
handles.fretThresholdArray = repmat( handles.defaultFretThreshold, ...
                                    handles.Ntraces,1  );

% Re-initialize figure objects
handles.molecule_no = 1;
set(handles.editGoTo,'Enable','on','String',num2str(handles.molecule_no));

set(handles.edThreshold, 'Enable','on', 'String',num2str(handles.fretThreshold));
set(handles.sldThreshold,'Enable','on', 'Value', handles.fretThreshold);

set(handles.edCrosstalk, 'Enable','on', 'String',num2str(handles.crosstalk));
set(handles.sldCrosstalk,'Enable','on', 'Value', handles.crosstalk);


set(handles.btnNextTop,'Enable','on');
set(handles.btnSubDonor,'Enable','on');
set(handles.btnSubBoth,'Enable','on');
set(handles.btnSubAcceptor,'Enable','on');
set(handles.btnPrevTop,'Enable','off');
set(handles.btnNextBottom,'Enable','on');
set(handles.btnPrevBottom,'Enable','off');
set(handles.btnPrint, 'Enable','on');
set(handles.btnLoadDWT, 'Enable','on');

handles.idl = [];

plotter(handles);

% END FUNCTION OpenTracesFile




%=========================================================================%
%======================   NAVIGATION   ======================%


%----------GO TO MOLECULE----------%
function handles = editGoTo_Callback(hObject, eventdata, handles)
% hObject    handle to editGoTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get trace ID from GUI
mol=str2double(get(hObject,'String'));

% If trace ID is invalid, reset it to what it was before.
if isnan(mol) || mol>handles.Ntraces || mol<1,
    disp('WARNING in sorttraces: Invalid trace number. Resetting.');
    set( hObject,'String',num2str(handles.molecule_no) );
    return;
else
    handles.molecule_no = mol;
end

% Make sure that the molecule selected actually exists.
if handles.molecule_no+1>handles.Ntraces
    set(handles.btnNextTop,'Enable','off');
    set(handles.btnNextBottom,'Enable','off');
else
    set(handles.btnNextTop,'Enable','on');
    set(handles.btnNextBottom,'Enable','on');
end
if handles.molecule_no-1>=1
    set(handles.btnPrevTop,'Enable','on');
    set(handles.btnPrevBottom,'Enable','on');
else
    set(handles.btnPrevTop,'Enable','off');
    set(handles.btnPrevBottom,'Enable','off');
end

% Reset these values for the new trace.
% handles.fretThreshold = handles.defaultFretThreshold;
handles.fretThreshold = handles.fretThresholdArray(mol);
handles.crosstalk     = handles.default_crosstalk;

handles.donor_background = 0;
handles.acceptor_background = 0;

% Set bin checkboxes
set(handles.chkBin1,'Value', any(handles.NoFRETs_indexes==mol) );
set(handles.chkBin2,'Value', any(handles.FRETs_indexes==mol) );
set(handles.chkBin3,'Value', any(handles.Best_indexes==mol) );

% Re-initialize figure objects.
set(handles.edCrosstalk,'String',num2str(handles.crosstalk));
set(handles.sldCrosstalk,'Value',handles.crosstalk);

set(handles.edThreshold,'String',num2str(handles.fretThreshold));
set(handles.sldThreshold,'Value',handles.fretThreshold);

set(handles.btnSubDonor,'Enable','on');
set(handles.btnSubBoth,'Enable','on');
set(handles.btnSubAcceptor,'Enable','on');
set(handles.btnSubUndo,'Enable','off');

guidata(hObject,handles);
plotter(handles);


%----------GO BACK TO PREVIOUS MOLECULE----------%
% --- Executes on button press in btnPrevTop.
function btnPrevTop_Callback(hObject, eventdata, handles)
% hObject    handle to btnPrevTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set( handles.editGoTo,'String',num2str(handles.molecule_no-1) );
editGoTo_Callback( handles.editGoTo, [], handles );



%----------GO TO NEXT MOLECULE----------%
% --- Executes on button press in btnNextTop - 'Next Molecule'.
function btnNextTop_Callback(hObject, eventdata, handles)
% hObject    handle to btnNextTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set( handles.editGoTo,'String',num2str(handles.molecule_no+1) );
editGoTo_Callback( handles.editGoTo, [], handles );




%=========================================================================%
%======================   MOLECULE BINNING   ======================%

% --- Executes on button press in chkBin1.
function addToBin_Callback(hObject, eventdata, handles, index)
% hObject    handle to chkBin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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




%----------SAVE TRACES----------%
% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

if ~isempty(handles.Best_indexes),
    filename = [baseFilename '_no_fret.traces'];
    savePickedTraces( handles, filename, handles.NoFRETs_indexes );
end


% Finish up
set(hObject,'Enable','off');



function savePickedTraces( handles, filename, indexes )
% Save picked traces and idealizations to file.

data.time     = handles.time;
data.donor    = handles.donor(indexes,:);
data.acceptor = handles.acceptor(indexes,:);
data.fret     = handles.fret(indexes,:);
data.traceMetadata = handles.traceMetadata(indexes);

saveTraces( filename, 'traces', data );

% Save idealizations of selected traces, if available.
if isfield(handles,'idl') && ~isempty(handles.idl),
    traceLen = numel(data.time);
    
    [p,f] = fileparts(filename);
    dwtFilename = [p filesep f '.qub.dwt'];
    
    saveDWT( dwtFilename, handles.dwt(indexes), ...
             (0:numel(indexes)-1)*traceLen, handles.dwtModel, handles.dwtSampling );
end

% end function savePickedTraces





%=========================================================================%
%======================   TRACE CORRECTIONS   ======================%


%----------SUBTRACT BOTH BACKGROUNDS----------%
% --- Executes on button press in btnSubBoth.
function btnSubBoth_Callback(hObject, eventdata, handles, mode)
% hObject    handle to btnSubBoth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

m = handles.molecule_no;

% Same as above, but only for both donor and acceptor.
xlim=get(handles.axFluor,'XLim');

% Convert axis limits to frames
if handles.time(1)~=1,
    dt = handles.time(2)-handles.time(1);
    xlim = floor(xlim./(dt/1000));
end

if xlim(1)<1, xlim(1)=1; end
if xlim(2)>handles.len, xlim(2)=handles.len; end

% Subtract donor background
if mode==1 || mode==3
    handles.donor_background=...
        mean(handles.donor(m,xlim(1):xlim(2)));
    handles.donor(m,:)=...
        handles.donor(m,:)-...
        handles.donor_background;
end

% Subtract acceptor background
if mode==2 || mode==3
    handles.acceptor_background=...
        mean(handles.acceptor(m,xlim(1):xlim(2)));
    handles.acceptor(m,:)=...
        handles.acceptor(m,:)-...
        handles.acceptor_background;
end

% Undo background subtraction
if mode==4
    handles.donor(m,:) = handles.donor(m,:) + handles.donor_background;
    handles.acceptor(m,:)=...
        handles.acceptor(m,:)+...
        handles.acceptor_background;
end

if mode<4
    set(handles.btnSubUndo,'Enable','on');    %undo
else  %undo
    set(handles.btnSubUndo,'Enable','off');    %undo
end

handles = updateTraceData( handles );
guidata(hObject,handles);
plotter(handles);




%----------ADJUST CROSSTALK WITH SLIDER----------%
% --- Executes on slider movement.
function sldCrosstalk_Callback(hObject, eventdata, handles)
% hObject    handle to sldCrosstalk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mol = handles.molecule_no;

% First, restore trace to it's original state w/ no crosstalk correction
oldCrosstalk = handles.crosstalk;
handles.acceptor(mol,:) = ...
    handles.acceptor(mol,:) + oldCrosstalk*handles.donor(mol,:);

% Second, make crosstalk correction again using new value
handles.crosstalk=get(hObject,'Value');
set(handles.edCrosstalk,'String',num2str(handles.crosstalk));

handles.acceptor(mol,:) = ...
    handles.acceptor(mol,:) - handles.crosstalk*handles.donor(mol,:);

% Save and display the result
handles = updateTraceData( handles );
guidata(hObject,handles);
plotter(handles);



%---------------------------------------------------%

% --- Executes on slider movement.
function sldThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to sldThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mol = handles.molecule_no;

% Update edit box with new value
handles.fretThreshold = get(hObject,'Value');
set(handles.edThreshold,'String', num2str(handles.fretThreshold) );

handles.fretThresholdArray(mol) = handles.fretThreshold;

% Plot data with new threshold value
handles = updateTraceData( handles );
guidata(hObject,handles);
plotter(handles);



function edThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edThreshold as text
%        str2double(get(hObject,'String')) returns contents of edThreshold as a double

% Update slider with new value
handles.fretThreshold = str2double( get(hObject,'String') );
set( handles.sldThreshold, 'Value', handles.fretThreshold );

handles.fretThresholdArray(mol) = handles.fretThreshold;

% Plot data with new threshold value
handles = updateTraceData( handles );
guidata(hObject,handles);
plotter(handles);





function handles = updateTraceData( handles )
% 

% Create a new FRET thresholded FRET signal
m        = handles.molecule_no;
donor    = handles.donor(m,:);
acceptor = handles.acceptor(m,:);
fret     = handles.fret(m,:);
total    = donor+acceptor;

stats = traceStat( donor,acceptor,fret );
lt = stats.lifetime;
fret = acceptor./total;
fret( total<handles.fretThreshold ) = 0;
if lt>0,
    fret( lt:end ) = 0;
end

handles.fret(m,:) = fret;



% END FUNCTION updateTraceData




%=========================================================================%
%======================   PLOT & PRINT   ======================%


%----------PRINT TRACE----------%
% --- Executes on button press in btnPrint.
function btnPrint_Callback(hObject, eventdata, handles)
% hObject    handle to btnPrint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

printdlg(handles.figure1);



%----------PLOT TRACES----------%
function plotter(handles)

constants = handles.constants;
% gamma = constants.gamma;
gamma = 1;

m        = handles.molecule_no;
donor    = handles.donor(m,:);
acceptor = handles.acceptor(m,:);
fret     = handles.fret(m,:);
total    = donor+acceptor;

try
    stats = traceStat( donor,acceptor,fret, constants );

    % Save calculated values
    %handles.minlt = minlt;
    lt = stats.lifetime;
    FRETlifetime = stats.acclife;
    snr = stats.snr;
    CC = stats.corr;
catch err
    lt = 0; FRETlifetime=0; snr=0; CC=0;
    disp('Could not calculate values: error in traceStat');
end


% Reset GUI values with these results
set(handles.editCorrelation,'String', sprintf('%.2f',CC) );
set(handles.editLifetime,'String',  sprintf('%d, %d', [FRETlifetime lt]));
set(handles.editSNR,'String', sprintf('%.2f',snr) );


[p,name] = fileparts( handles.filename );
data_fname = strrep(name, '_', '\_');

% Plot fluorophore traces
signal = donor+acceptor;
% time = 1:handles.len;

time = handles.time;
inFrames = (handles.time(1)==1);
if ~inFrames
    time = time/1000; %in seconds
end

axes(handles.axFluor);
cla;
hold off;
plot( time,gamma*donor,'g', time,acceptor,'r' );
title(['Molecule ' num2str(m) ' of '...
    num2str(handles.Ntraces) ' of ' data_fname]);
ylabel('Fluorescence');
grid on;
zoom on;
hold on;


% Plot total fluorescence
axes(handles.axTotal);
hold off;
cla;
plot( time,gamma*donor+acceptor,'k',...
    time,handles.fretThreshold*ones(1,handles.len),'m');
ylabel('Total Fluorescence');
grid on;
zoom on;
hold on;


% Draw lines representing Cy3 (green) and Cy5 (red) alive times
mean_on_signal  = mean( signal(1:lt) );
mean_off_signal = mean( signal(lt+5:end) );

simplified_cy3=[mean_on_signal*ones(1,lt)...
        mean_off_signal*ones(1,handles.len-lt)];
simplified_cy5=[mean_on_signal*ones(1,FRETlifetime)...
        mean_off_signal*ones(1,handles.len-FRETlifetime)];

axes(handles.axTotal);
plot(time,simplified_cy5,'r');
plot(time,simplified_cy3,'g');
zoom on;


% Plot FRET efficiency
axes(handles.axFret);
cla;
plot(time,fret,'b');  hold on;
if isfield(handles,'idl') && ~isempty(handles.idl),
    stairs( time, handles.idl(m,:), 'r-', 'LineWidth',1 );
end

xlabel('Frame Number');
ylabel('FRET Efficiency');
ylim([-0.1 1]);
grid on;
zoom on;
hold off;

if inFrames,
    xlabel('Frame Number');
else
    xlabel('Time (sec)');
end







% --- Executes on button press in btnLoadDWT.
function btnLoadDWT_Callback(hObject, eventdata, handles)
% Loads an idealization (.dwt file) for later plotting in plotter().

time = handles.time;
sampling = time(2)-time(1);
[nTraces,traceLen] = size(handles.fret);

% Get filename for .dwt file from user and load it.
[dwt,dwtSampling,offsets,model] = loadDWT;
if isempty(dwt), return; end

handles.dwt = dwt;
handles.dwtSampling = dwtSampling;
handles.dwtModel = model;

dwtSampling = double(dwtSampling);

% Verify that the DWT matches the trace data.
handles.idl = [];

if sampling~=dwtSampling, 
    msgbox('Data and idealization sampling intervals do not match!','Error loading idealization','Error')

elseif offsets(end)>(nTraces*traceLen)
    msgbox('Idealization size does not match trace data.','Error loading idealization','Error');

else
    % Convert DWT to idealization (state sequence).
    idl = dwtToIdl( dwt, traceLen, offsets );

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
    
%     handles.idl = reshape( idl, size(handles.fret) );
    handles.idl = idl;
    
end %if errors


% Save data to GUI
guidata(hObject,handles);

% Update GUI
plotter(handles);







% --- Executes on button press in btnSelAll3.
function btnSelAll_Callback(hObject, eventdata, handles, index)

if index==1,
    handles.NoFRETs_indexes = 1:handles.Ntraces;
    set(handles.chkBin1,'Value',1);
elseif index==2,
    handles.FRETs_indexes = 1:handles.Ntraces;
    set(handles.chkBin2,'Value',1);
elseif index==3,
    handles.Best_indexes = 1:handles.Ntraces;
    set(handles.chkBin3,'Value',1);
end
    
%unchecking
% if index==1,
%     handles.NoFRETs_indexes = [];
% elseif index==2,
%     handles.FRETs_indexes = [];       
% elseif index==3,
%     handles.Best_indexes = [];
% end


NoFRET_no = numel( handles.NoFRETs_indexes );
FRET_no = numel( handles.FRETs_indexes );
Best_no = numel( handles.Best_indexes );

set(handles.editBin1,'String',num2str(NoFRET_no));
set(handles.editBin2,'String',num2str(FRET_no));  
set(handles.editBin3,'String',num2str(Best_no));

set(handles.btnSave,'Enable','on');
guidata(hObject,handles);




% --- Executes on button press in btnGettraces.
function btnGettraces_Callback(hObject, eventdata, handles)
%

% Get the filename and movie coordinates of the selected trace.
% FIXME: this assumes new format traces!
m = handles.molecule_no;
id = handles.traceMetadata(m).ids;

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

gettraces_gui( movieFilename, handles.traceMetadata(m) );


% end function btnGettraces_Callback





