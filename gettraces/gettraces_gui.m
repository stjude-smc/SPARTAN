function varargout = gettraces_gui(varargin)
% GETTRACES  Extract fluorescence traces from wide-field fluorescence movies
%
%   This GUI provides an interface to view wide-field fluoerscence movies
%   (.tif and .stk files), detect peaks of fluorescence intensity within the
%   field of view, align peaks from separate spectral channels, and save
%   single molecule fluorescence traces by summing the intensity over a fixed
%   window around each peak in each frame of the movie.
%
%   Imaging profiles describe how the movie is divided into fluorescence
%   channels (side-by-side tiles, sequential frames, etc), descriptions of these
%   channels, and various analysis settings. Default profiles are defined in
%   cascadeConstants.m and can be modified for your lab's setup. Custom profiles
%   can also be saved within the GUI that are persisent across sessions, but
%   have less flexibility.
%
%   For algorithm details, see the MovieParser class and associated methods.
%
%   See also: MovieParser, Movie_TIFF, subfield, tirfProfile.

%   Copyright 2007-2016 Cornell University All Rights Reserved.

% Last Modified by GUIDE v2.5 03-May-2018 15:09:13


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gettraces_OpeningFcn, ...
                   'gui_OutputFcn',  @gettraces_OutputFcn, ...
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


% --------------------- GUI INITIALIZATION --------------------- %
% --- Executes just before gettraces is made visible.
function gettraces_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.

updateSpartan; %check for updates
constants = cascadeConstants();
set( handles.figure1, 'Name', ['gettraces - ' constants.software] );
handles.output = hObject;

% Load colormap for image viewer
set( handles.figure1, 'Colormap',gettraces_colormap() );

% Load built-in and saved user-custom.
profiles = constants.gettraces_profiles;
handles.nStandard = numel(profiles);

if ispref('SPARTAN','gettraces_customProfiles')
    try
        profiles = [profiles getpref('SPARTAN','gettraces_customProfiles')];
    catch
        warndlg('Previously saved custom profiles could not be loaded due to a version incompatibility');
        rmpref('SPARTAN','gettraces_customProfiles');
    end
end

% Add profiles from cascadeConstants to settings menu.
for i=1:numel(profiles),
    if i<=numel(constants.gettraces_profiles), pos=i; else, pos=i+1; end
    
    hMenu(i) = uimenu(handles.mnuProfiles, 'Label',profiles(i).name, ...
                      'Position',pos, 'Callback',@mnuProfiles_Callback); %#ok<AGROW>
end

% Put customization menu items in the correct spots.
set(handles.mnuSettingsCustom,'Position',handles.nStandard+1);
handles.profile = constants.gettraces_defaultProfile;  %index to current profile (FIXME: rename)
set( hMenu(handles.profile), 'Checked','on' );

% Context menus for field-specific settings (names, wavelength, etc).
hZoom = zoom(handles.figure1);
set(hZoom, 'UIContextMenu', handles.mnuField);
zoom(handles.figure1,'on');

% Setup default values for parameter values -- 2-color FRET.
handles.hProfileMenu = hMenu;
handles.profiles = profiles;
handles.params = profiles(handles.profile);
guidata(hObject, handles);

% Set up GUI elements to reflect the internal parameter values.
mnuProfiles_Callback( hMenu(handles.profile), [], handles);

% END FUNCTION gettraces_OpeningFcn



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, handles) %#ok<DEFNU>
% Save custom profiles for future sessions before exit
% FIXME: setpref can interfere with older versions. enable with caution.
% setpref('SPARTAN','gettraces_customProfiles',handles.profiles(handles.nStandard+1:end));
delete(hObject);
% END FUNCTION figure1_CloseRequestFcn



% --- Outputs from this function are returned to the command line.
function varargout = gettraces_OutputFcn(~, ~, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;




%==========================================================================
%=========================  OPEN and DISPLAY MOVIE  =======================
%==========================================================================

% --- Executes on button press in openstk.
function openstk_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Get filename of input data from user. 
% Multi-select is for multi-part movies (ordinary TIFFs limited to 2GB).
[datafile,datapath] = uigetfile( '*.stk;*.tif;*.tiff', 'Choose a movie file', ...
                                 'MultiSelect','on' );
if ~iscell(datafile),
    if datafile==0, return; end  %user hit cancel
end
handles.stkfile = strcat(datapath,datafile);

% Load the movie
OpenStk( handles.stkfile, handles, hObject );

% END FUNCTION openstk_Callback



% --------------------- OPEN SINGLE MOVIE --------------------- %
function handles = OpenStk(filename, handles, hObject)
% Load movie data from file. filename is a cell array of files.

if isempty(filename), return; end
if ~iscell(filename), filename = {filename}; end

% Remove "-file00x" extension for multi-file TIFF stacks.
[p,f,e] = fileparts( filename{1} );
f = regexprep(f,'-file[0-9]*$','');

% If a single file is selected, look for the others in multi-file TIFFs.
% If the user selected multiple files, we assume that they got all of them.
% FIXME: not used anymore and will fail if the filename has any special
% characters in it (e.g., plus sign). Uncomment with caution.
% This essentially forces the user to select all files in a multi-part TIFF.
% if numel(filename)==1,
%     d = dir( [f '*.tif*'] );
%     if numel(d)>1,
%         d = regexpi( {d.name}, [f '(-file[0-9]*)?\.tiff?$'], 'match' );
%         filename = [d{:}];
%     end
% end

% Trancate name if too long ot fit into window without wrapping.
fnameText = fullfile(p, [f e]);
if numel(fnameText)>80,
    fnameText = ['...' fnameText(end-80:end)];
end

set( handles.txtFilename, 'String',fnameText);
set( handles.txtAlignWarning, 'Visible','off' );

set([handles.txtOverlapStatus handles.txtIntegrationStatus, ...
     handles.txtWindowOverlap handles.txtPSFWidth handles.nummoles], ...
     'String', '');
set(handles.panAlignment, 'ForegroundColor',[0 0 0]);
set(handles.tblAlignment, 'Data',{});


% Load movie data, clearing original to save memory
set(handles.figure1,'pointer','watch'); drawnow;

stkData = MovieParser( filename, handles.params );  %does not draw from gettraces params!
handles.stkData = stkData;

set( handles.sldScrub, 'Min',1, 'Max',stkData.nFrames, 'Value',1, ...
     'SliderStep',[1/stkData.nFrames,0.02] );

% Setup slider bar (adjusting maximum value in image, initially 2x max)
idxFields = handles.params.idxFields;
stk_top = cat(3, stkData.stk_top{idxFields} );
sort_px = sort(stk_top(:));
val = 2*sort_px( floor(0.99*numel(sort_px)) );
high = min( ceil(val*10), 32000 );  %uint16 maxmimum value

set(handles.scaleSlider,'min',0, 'max',high, 'value', val);
set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));

% Create axes for sub-fields (listed in column-major order, like stk_top)
delete( findall(handles.figure1,'type','axes') );  %remvoe old axes
axopt = {'Visible','off'};

switch numel(idxFields)
case 1
    ax = [];
    handles.axTotal = axes( handles.panView, 'Position',[0.02 0 0.95 0.95], axopt{:} );
    
case 2
    ax(1)           = axes( handles.panView, 'Position',[0     0    0.325 0.95], axopt{:} );  %L
    ax(2)           = axes( handles.panView, 'Position',[0.335 0    0.325 0.95], axopt{:} );  %R
    handles.axTotal = axes( handles.panView, 'Position',[0.67  0    0.325 0.95], axopt{:} );
    
case {3,4}
    ax(1)           = axes( handles.panView, 'Position',[0.0   0.5  0.325 0.47], axopt{:} );  %TL
    ax(2)           = axes( handles.panView, 'Position',[0     0    0.325 0.47], axopt{:} );  %BL
    ax(3)           = axes( handles.panView, 'Position',[0.335 0.5  0.325 0.47], axopt{:} );  %TR
    ax(4)           = axes( handles.panView, 'Position',[0.335 0    0.325 0.47], axopt{:} );  %BR
    handles.axTotal = axes( handles.panView, 'Position',[0.67  0.25 0.325 0.47], axopt{:} );
    
otherwise
    error('Invalid field geometry');
end

% Show fluorescence fields for all channels
axopt = {'YDir','reverse', 'Color',get(handles.figure1,'Color'), 'Visible','off'};
handles.himshow = [];
handles.ax = ax;

if ~isempty(ax)
    for i=1:numel(stkData.stk_top),
        handles.himshow(i) = image( ax(i), stkData.stk_top{i}, 'CDataMapping','scaled' );
        set(ax(i), 'UserData',i, 'CLim',[0 val], axopt{:});
    end
    linkaxes( [to_row(ax) handles.axTotal] );
end
set( handles.himshow, 'UIContextMenu',handles.mnuField );

% Show total fluorescence channel
total = sum( cat(3,stkData.stk_top{idxFields}), 3);
handles.himshow(end+1) = image( handles.axTotal, total, 'CDataMapping','scaled' );
set(handles.axTotal, 'CLim',[0 val*2], axopt{:} );

if abs(log2( size(total,2)/size(total,1) )) < 1
    % For roughly symmetric movies, allows for better window resizing.
    axis( [ax handles.axTotal], 'image' );  %adjust axes size to match image
else
    % Allows for better zooming of narrow (high time resolution) movie.
    axis( [ax handles.axTotal], 'equal' );  %keep the axes size fixed.
end
setAxTitles(handles);

guidata(hObject,handles);
set( handles.figure1, 'pointer','arrow');
set([handles.btnPickPeaks handles.mnuPick handles.mnuViewMetadata handles.mnuTirfProfile], 'Enable','on');
set([handles.btnSave handles.mnuFileSave handles.mnuFileSaveAs],'Enable','off');

%end function OpenStk



% --- Executes on slider movement.
function sldScrub_Callback(hObject, ~, handles) %#ok<DEFNU>
% Allows user to scroll through the movie

% Stop any currently playing movies.
set(handles.btnPlay,'String','Play');

% Read frame of the movie
idxFrame = round( get(hObject,'Value') );
allgeo = true( size(handles.params.geometry) );
fields = subfield( handles.stkData.movie, allgeo, idxFrame );

for f=1:numel(fields)
    field = single(fields{f}) - handles.stkData.background{f};
    set( handles.himshow(f), 'CData',field );
end

% END FUNCTION sldScrub_Callback



% --- Executes on button press in btnPlay.
function btnPlay_Callback(~, ~, handles) %#ok<DEFNU>
% Play movie

% Clicking when the button is labeled 'stop' causes the loop below to terminate.
if strcmpi( get(handles.btnPlay,'String') ,'Play' )
    set(handles.btnPlay,'String','Stop');
else
    set(handles.btnPlay,'String','Play');
    return;
end

stkData = handles.stkData;
startFrame = round( get(handles.sldScrub,'Value') );
allgeo = true( size(handles.params.geometry) );

for i=startFrame:stkData.movie.nFrames
    fields = subfield( stkData.movie, allgeo, i );
    
    for f=1:numel(fields)
        field = single(fields{f}) - stkData.background{f};
        set( handles.himshow(f), 'CData',field );
    end
    
    set(handles.sldScrub,'Value',i);
    drawnow;
    
    % Terminate early if the user clicks the 'Stop' button.
    state = get(handles.btnPlay,'String');
    if strcmpi(state,'Play'), return; end
end

set(handles.btnPlay,'String','Play');

% END FUNCTION btnPlay_Callback




% --------------- OPEN ALL STK'S IN DIRECTORY (BATCH) --------------- %

% --- Executes on button press in batchmode.
function batchmode_Callback(hObject, ~, handles,direct)
% direct     target directory location to look for new files

% Get input parameter values
skipExisting = strcmpi( get(handles.mnuBatchOverwrite,'Checked'), 'on');
recursive = strcmpi( get(handles.mnuBatchRecursive,'Checked'), 'on');
auto = strcmpi( get(handles.mnuBatchAuto,'Checked'), 'on');
skipExisting = auto || skipExisting;  %always skip if auto-detect new files.

% Get location of files for gettraces to process
if nargin>=4 && exist(direct,'dir'),
    % Get a fresh copy of handles. The one passed in arguments is an old
    % copy made when the timer was created. Kind of ugly...
    handles = guidata(handles.mnuBatchRecursive);
else
    direct=uigetdir('','Choose directory:');
    if direct==0, return; end
    disp(direct);
end

% Get list of files in current directory (option: and all subdirectories)
movieFiles = regexpdir(direct,'^.*\.(tiff?|stk)$',recursive);
nFiles = length(movieFiles);

% Wait for 100ms to give sufficient time for polling file sizes in the
% main loop below.
pause(0.1);


% ---- Run the ordinary gettraces procedure on each file
nTraces  = zeros(nFiles,1); % number of peaks found in file X
existing = zeros(nFiles,1); % true if file was skipped (traces file exists)

for i=1:nFiles
    text = sprintf('Batch Processing: %d/%d (%.0f%% complete)', i,nFiles, 100*((i-1)/nFiles) );
    set(handles.txtProgress,'String',text);
    
    stk_fname = movieFiles(i).name;
    handles.stkfile = stk_fname;
    
    % Skip if previously processed (.traces file exists)
    [p,name] = fileparts(stk_fname);
    traceFname = fullfile(p, [name '.rawtraces']);
    
    if skipExisting && exist(traceFname,'file'),
        existing(i) = 1;
        continue;
    end
    
    % Poll the file to make sure it isn't changing.
    % This could happen when a file is being saved during acquisition.
    d = dir(stk_fname);
    if movieFiles(i).datenum ~= d(1).datenum,
        disp( ['Skipping (save in process?): ' stk_fname] );
        existing(i) = 1;
        continue;
    end
    
    % Load movie file
    handles = OpenStk(handles.stkfile,handles, hObject);
    
    % Pick molecules using default parameter values
    handles = getTraces_Callback(hObject, [], handles);
    
    % Save the traces to file
    mnuFileSave_Callback(hObject, [], handles);
    nTraces(i) = handles.num;
end


% ----- Create log file with results
log_fid = fopen( fullfile(direct,'gettraces.log'), 'wt' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');
fprintf(log_fid, '%s', evalc('disp(handles.params)'));
%FIXME: structure parameters are not displayed here (alignment!)

% Log list of files processed by gettraces
fprintf(log_fid,'\n%s\n\n%s\n%s\n\n%s\n',date,'DIRECTORY',direct,'FILES');

for i=1:nFiles
    if existing(i),
        fprintf(log_fid, 'SKIP %s\n', movieFiles(i).name);
    else
        fprintf(log_fid, '%.0f %s\n', nTraces(i), movieFiles(i).name);
    end
end

fclose(log_fid);
set(handles.txtProgress,'String','Batch processing: finished.');

% END FUNCTION batchmode_Callback





%==========================================================================
%=======================  PICK PEAKS and SAVE TRACES  =====================
%==========================================================================

function handles = getTraces_Callback(hObject, ~, handles)
% Find peak locations from total intensity

% Do nothing if no file has been loaded. This function may be triggered by
% changing the settings fields before a file is open.
if ~isfield(handles, 'stkfile')
    return;
end
set(handles.tblAlignment, 'Data',{}, 'RowName',{});
set(handles.figure1,'pointer','watch'); drawnow;
stkData = handles.stkData;

% Clear any existing selection markers from previous calls.
delete(findobj(handles.figure1,'type','line')); drawnow;

% Locate single molecules
try
    getPeaks(stkData, handles.params);
catch e
    if ~strcmpi(e.identifier,'parfor_progressbar:cancelled')
        errordlg( ['Error: ' e.message], 'Gettraces' );
    end
    set(handles.figure1,'pointer','arrow'); drawnow;
    return;
end

% The alignment may involve shifting (or distorting) the fields to get a
% registered donor+acceptor field. Show this distorted imaged so the user
% can see what the algorithm is doing.
val = get(handles.scaleSlider,'value');
minimum = get(handles.scaleSlider,'min');
val = max(val,minimum+1);

set( handles.himshow(end), 'CData',stkData.total_t );
set( handles.axTotal, 'CLim',[minimum*2 val*2] );


% If no alignment data given (for example in single-channel recordings),
% don't display any status messages.
set( handles.txtAlignWarning, 'Visible','off' );

if isempty(stkData.alignStatus),
    set(handles.panAlignment, 'Visible','off', 'ForegroundColor', [0 0 0]);
    
% Display alignment status to inform user if realignment may be needed.
% Format: translation deviation (x, y), absolute deviation (x, y)
else
    a = stkData.alignStatus;
    handles.params.alignment = a;
    
    tableData = cell(numel(a)-1,6);
    fmt = {'% 0.2f','% 0.2f','% 0.2f','% 0.2f %%','%0.2f','%0.2f'};  %sprintf formats for each field
    idxShow = find(  ~arrayfun(@(x)isempty(x.theta), a)  );  %field numbers
    
    for idxRow=1:numel(idxShow),  %row in displayed table
        %FIXME: quality is not defined for all alignments (!)
        i = idxShow(idxRow);  %field number
        row = [a(i).dx a(i).dy a(i).theta 100*(a(i).sx-1) a(i).quality a(i).abs_dev];
        tableData(idxRow,:) = arrayfun( @(x)sprintf(fmt{x},row(x)), 1:numel(fmt), 'Unif',false )';

        if a(i).quality==0,
            tableData{idxRow,5}='';
        end
    end
    
    set( handles.tblAlignment, 'Data',tableData(1:numel(a)-1,:) );
    set( handles.tblAlignment, 'RowName',handles.params.chDesc(idxShow) );
    
    % If the alignment quality (confidence) is low, warn the user.
    methods = {'Alignment Disabled','Aligned from File', ...
               'Auto Aligned (ICP)','Alignment Memorized'};
    text = methods{handles.params.alignMethod};
    
    if isfield(a,'quality') &&  any( [a.quality]<1.1 & [a.quality]>0 ),
        text = [text sprintf(': LOW CONFIDENCE!')];
    end
    set(handles.panAlignment, 'Title',text);

    % Color the text to draw attention to it if the alignment is bad.
    % FIXME: this should depend on the nhood/window size. 1 px may be small.
    d = max(0,min(1,  1.75*(max([a.abs_dev])-0.25)  ));
    set( handles.tblAlignment, 'ForegroundColor', d*[1 0 0] );
    set( handles.panAlignment, 'ForegroundColor', d*[1 0 0] );
    
    % Show a big warning for total misalignment (no corresponding peaks).
    if any( [a.abs_dev] >=0.7 ),
        set( handles.txtAlignWarning, 'Visible','on' );
    end
end


% Fraction of molecules close enough for PSFs to overlap (overcrowding).
percentOverlap = stkData.fractionOverlapped*100;
c = max(0,min(1, (percentOverlap-40)/10 ));

set(  handles.txtOverlapStatus, 'ForegroundColor', c*[0.9 0 0], ...
       'String', sprintf('Molecules rejected: %0.0f%%', percentOverlap)  );


% Fraction of overlapping integration windows (also overcrowding).
percentWinOverlap = mean(stkData.fractionWinOverlap*100);
set(  handles.txtWindowOverlap, 'String', ...
      sprintf('Residual win. overlap: %0.1f%%', percentWinOverlap)  );
set( handles.txtWindowOverlap, 'ForegroundColor', (percentWinOverlap>10)*[0.9 0 0] );


% Get (approximate) average fraction of fluorescence collected within the
% integration window of each molecule. Set the text color to red where the
% intensity is not well collected at the current integration window size.
eff = stkData.integrationEfficiency;
set(  handles.txtIntegrationStatus, 'String', ...
      sprintf('Intensity collected: %0.0f%% ', eff)  );
set( handles.txtIntegrationStatus, 'ForegroundColor', (eff<70)*[0.9 0 0] );


% Estimate the peak width from pixel intensity distribution.
set( handles.txtPSFWidth, 'String', sprintf('PSF size: %0.1f px',stkData.psfWidth) );
set( handles.txtPSFWidth, 'ForegroundColor', (stkData.psfWidth>handles.params.nPixelsToSum-1)*[0.9 0 0] );


% Graphically show peak centers
highlightPeaks( handles );

% Update GUI controls
handles.num = size(stkData.peaks,1);
set( handles.nummoles, 'String', sprintf('%d (of %d)',handles.num, ...
                      size(stkData.rejectedPicks,1)+handles.num) );

set( [handles.btnSave handles.mnuFileSave handles.mnuFileSaveAs ...
                                     handles.mnuHidePeaks], 'Enable','on');
set( [handles.mnuAlignSave handles.mnuAlignKeep], ...
                   'Enable',onoff(handles.params.alignMethod>1) );

set(handles.figure1,'pointer','arrow');
guidata(hObject,handles);

% end function



function highlightPeaks(handles)
% Draw circles around each selected fluorescence spot.

style = {'LineStyle','none','marker','o'};
stkData = handles.stkData;
idxField = find(handles.params.geometry);

if ~isscalar(handles.params.geometry)
    for i=1:size(stkData.peaks,3)
        ax = handles.ax( idxField(i) );

        line( stkData.peaks(:,1,i), stkData.peaks(:,2,i), ...
                style{:}, 'color','w', 'Parent',ax );
        line( stkData.rejectedPicks(:,1,i), stkData.rejectedPicks(:,2,i), ...
                style{:}, 'color',[0.4,0.4,0.4], 'Parent',ax );
    end
end

% Draw markers on selection points (total intensity composite image).
line( stkData.total_peaks(:,1), stkData.total_peaks(:,2), ...
        style{:}, 'color','y',  'Parent',handles.axTotal );
line( stkData.rejectedTotalPicks(:,1), stkData.rejectedTotalPicks(:,2), ...
        style{:}, 'color',[0.4,0.4,0.0], 'Parent',handles.axTotal );

% end function highlightPeaks


% --- Executes on button press in btnHidePicks.
function btnHidePicks_Callback(~, ~, handles)  %#ok<DEFNU>
% Hide the circles drawn to indicate molecule locations.
delete(findobj(handles.figure1,'type','line'));
% END FUNCTION btnHidePicks_Callback



% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

% --- Executes on button press in mnuFileSave.
function mnuFileSave_Callback(~, ~, handles, prompt)

% Create output file name
filename = handles.stkfile;
if iscell(filename), filename=filename{1}; end
[p,f] = fileparts(filename);
f = regexprep(f,'-file[0-9]*$',''); %remove multi-part TIFF extension
filename = fullfile(p,[f '.rawtraces']);

% Prompt user for target filename, if called from the "File->Save As" menu.
if nargin>=4 && prompt
    [f,p] = uiputfile('*.traces', 'Save traces as:', filename);
    if f==0, return; end  %user hit cancel
    filename = fullfile(p,f);
end

set(handles.figure1,'pointer','watch'); drawnow;

% Integrate fluorophore PSFs into fluorescence traces and save to file.
try
    integrateAndSave(handles.stkData, filename, handles.params);
catch e
    if ~strcmpi(e.identifier,'parfor_progressbar:cancelled')
        errordlg( ['Error: ' e.message], 'Gettraces' );
    end
end
set(handles.figure1,'pointer','arrow'); drawnow;

% END FUNCTION mnuFileSave_Callback





%========================================================================
%======================  VIEW/SETTINGS CALLBACKS  =======================
%========================================================================

% --- Executes on slider movement.
function scaleSlider_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Update axes color limits from new slider value
val = get(hObject,'value');
set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));
txtMaxIntensity_Callback(handles.txtMaxIntensity, [], handles);
% END FUNCTION scaleSlider_Callback


function txtMaxIntensity_Callback(hObject, ~, handles)
% Update axes color limits from new slider value
val = str2double( get(hObject,'String') );
minimum = get(handles.scaleSlider,'min');
maximum = get(handles.scaleSlider,'max');

val = max(val,minimum+1); %prevent errors in GUI
maximum = max(val,maximum);

set( handles.scaleSlider, 'Value',val );
set( handles.scaleSlider, 'max',maximum );
set( handles.ax, 'CLim',[minimum val] );

if isscalar(handles.params.geometry)  %Single-channel recordings
    set( handles.axTotal, 'CLim',[minimum val] );
else
    set( handles.axTotal, 'CLim',[minimum*2 val*2] );
end

set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));
guidata(hObject,handles);

% END FUNCTION txtMaxIntensity_Callback



% --------------------------------------------------------------------
function mnuProfiles_Callback(hObject, ~, handles)
% Called when any imaging profile is selected from the "Settings" menu.
% Apply imaging profile settings, including dividing the movie into quadrants.

if nargin<3,  handles = guidata(hObject);  end

% Save any changes to previous profile
handles.profiles(handles.profile) = handles.params;

% Mark the new profile as current
set(handles.hProfileMenu,'Checked','off');
set(hObject, 'Checked','on');
handles.profile = find(hObject==handles.hProfileMenu);

% If running, stop the "auto detect" timer. Otherwise, it may be triggered by
% the change in settings.
fileTimer = timerfind('Name','gettraces_fileTimer');
if ~isempty(fileTimer),
    stop(fileTimer);
    delete(fileTimer);
    set(handles.mnuBatchAuto,'Checked','off');
end

% Get parameter values associated with the selected profile.
params = handles.profiles(handles.profile);
handles.params = params;
guidata(hObject,handles);

% Set all GUI to defaults of currently selected profile.
set( handles.mnuBatchRecursive,   'Checked', onoff(params.recursive)    );
set( handles.mnuBatchOverwrite,   'Checked', onoff(params.skipExisting) );

set( findobj('Parent',handles.mnuAlign), 'Checked','off' );
set( findobj('Parent',handles.mnuAlign,'Position',params.alignMethod), ...
     'Checked','on' );

% Enable alignment, crosstalk, and scale controls only in multi-color.
isMultiColor = ~isscalar(handles.params.geometry);
set( [handles.mnuAlign handles.btnCrosstalk handles.btnScaleAcceptor], ...
                       'Enable',onoff(isMultiColor) );
set( handles.tblAlignment, 'Visible',onoff(isMultiColor) );

set( handles.panAlignment, 'Title','Software Alignment', ...
                           'Visible',onoff(isMultiColor) );

% If a movie has already been loaded, reload movie with new setup.
if isfield(handles,'stkfile')
    handles = OpenStk( handles.stkfile, handles, hObject );
    set([handles.mnuAlignSave handles.mnuAlignKeep], 'Enable',onoff(params.alignMethod>1));
end

% END FUNCTION cboGeometry_Callback



% --------------------------------------------------------------------
function mnuSettingsCustom_Callback(hObject, ~, handles) %#ok<DEFNU>
% Called when Settings->Customize... menu clicked.
% Allows the user to temporarily alter settings for the current profile.

params = handles.params;
oldParams = params;

% Create options for bgTraceField, which can be empty or a number.
% (settingdlg doesn't like this, so we have to convert it to string. FIXME)
nCh = numel(handles.params.chNames);
fopt = cellfun(@num2str, num2cell(1:nCh), 'uniform',false);
fopt = [{''} fopt];
params.bgTraceField = num2str(params.bgTraceField);

% Create dialog for changing imaging settings
prompt = {'Name:', 'Threshold (0 for auto):', 'Auto picking sensitivity', ...
          'Integration window size (px):', 'Integration neighbhorhood (px):', ...
          'Minimum separation (px):', 'ADU/photon conversion:', ...
          'Donor blink detection method:', 'Background trace field:', ...
          'Frames to average for picking:'};
fields = {'name', 'don_thresh', 'thresh_std', 'nPixelsToSum', 'nhoodSize', ...
          'overlap_thresh', 'photonConversion', 'zeroMethod', 'bgTraceField', ...
          'nAvgFrames'};
isInt = @(x)~isnan(x) && isreal(x) && isscalar(x) && x==floor(x);
isNum = @(x)~isnan(x) && isreal(x) && isscalar(x);
types = {[],isNum,isNum,isInt,isInt,isNum,isNum,{'off','threshold','skm'},fopt,isInt};

if handles.profile > handles.nStandard
    prompt{1} = 'Name (clear to remove profile):';
end

params = settingdlg(params, fields, prompt, types);
if isempty(params), return; end  %user hit cancel

if ~isempty(params.bgTraceField),
    params.bgTraceField = str2double(params.bgTraceField);
end

% If a standard profile is renamed, make a custom setting instead.
if handles.profile < handles.nStandard  && ~strcmpi(params.name,handles.params.name)
    if isempty(params.name),
        errordlg('Built-in profiles cannot be deleted. Modify cascadeConstants instead');
        return;
    end
    
    % Overwrite if name matches an existing profile, except built-in.
    handles.params = params;
    handles.profiles(end+1) = params;
    handles.profile = numel(handles.profiles);

    % Update settings menu with new item.
    set( handles.hProfileMenu, 'Checked','off' );  %uncheck all
    handles.hProfileMenu(end+1) = uimenu( handles.mnuProfiles, 'Checked','on', ...
              'Label',handles.params.name, 'Callback',@mnuProfiles_Callback );
    guidata(hObject,handles);

% Modify existing custom profile
elseif ~isempty(params.name)
    set( handles.hProfileMenu(handles.profile), 'Label',params.name );
    handles.params = params;
    guidata(hObject,handles);
    
% Delete profile if name was cleared
else
    delete( handles.hProfileMenu(handles.profile) );
    handles.hProfileMenu(handles.profile) = [];
    handles.profiles(handles.profile) = [];
    mnuProfiles_Callback( handles.hProfileMenu(1), [], handles );
    return;
end

if isfield(handles,'stkData')  %if a movie has been loaded
    
    % Reload movie if changing number of frames to average
    if params.nAvgFrames~=oldParams.nAvgFrames
        OpenStk( handles.stkfile, handles, hObject );

    % If molecules were already picked and settings have changed, re-pick.
    elseif isfield(handles,'stkData') && ~isempty(handles.stkData.peaks)
        getTraces_Callback(hObject, [], handles);
    end
    
end

% END FUNCTION mnuSettingsCustom_Callback



% --- Executes on button press in btnMetadata.
function btnMetadata_Callback(~, ~, handles)  %#ok<DEFNU>
% Display a simple diaglog with the MetaMorph metadata for the first frame.

stkData = handles.stkData;
if isempty(stkData), return; end  %disable button?

% MetaMorph specific metadata
if isfield(stkData.movie.header,'MM')
    metadata = stkData.movie.header.MM(1);
    if isfield( stkData.movie.header, 'MM_wavelength' ),
        wv = stkData.movie.header.MM_wavelength;
        metadata.MM_wavelength = wv(wv>100);
    end
    
% Generic TIFF format metadata
else
    metadata = stkData.movie.header;
end

% Display metadata fields as a list in a message box.
fields = fieldnames(metadata);
output = {};

for i=1:numel(fields),
    fname = fields{i};
    data = metadata.(fname);
    if isempty(data), continue; end
    if strcmpi(fname,'DateTime'), continue; end  % Ignore frame-specific fields
    
    % Truncate long text fields
    if isnumeric(data)
        data = num2str(data);
    else
        idxline = find(data==newline,1,'first');
        if ~isempty(idxline),
            data = data(1:idxline-1);
        end
    end
    if numel(data)>70
        data = [data(1:70) ' ...'];
    end
    
    output{end+1} = sprintf( '%s:  %s', fname, data ); %#ok
    disp( output{end} );
end

msgbox( output, 'Movie metadata' );

% END FUNCTION btnMetadata_Callback



function mnuTirfProfile_Callback(~, ~, handles) %#ok<DEFNU>
% Executes when the "View->Illumination Profile" menu is clicked.
% FIXME: consider directly integration only from stk_top instead.
% This would not require saving traces, which takes a long time.

[p,f] = fileparts(handles.stkfile);
outname = fullfile(p,[f '.rawtraces']);
if exist(outname,'file'),
    tirfProfile(outname);
else
    disp('No .rawtraces file. Save one first to get a profile');
end

% END FUNCTION mnuTirfProfile_Callback



% --------------------------------------------------------------------
function mnuViewMontage_Callback(varargin) %#ok<DEFNU>
% Display multiple movies simultaneously for direct comparison.

% Get movie paths from user
filter = {'*.tif;*.tiff;*.stk','Movie Files (*.tif;*.tiff)'; ...
          '*.*','All Files (*.*)'};
f = getFiles(filter,'Movie Montage: Select Files');
if isempty(f), return; end  %user hit cancel.

% Reshape input for common sizes. FIXME
if numel(f)==4, reshape(f,2,2); end

% Create a new window to view all movies simultaneously
m = MovieMontage(f);
m.show();

% END FUNCTION mnuViewMontage_Callback





%========================================================================
%======================  FIELD SETTINGS CALLBACKS  ======================
%========================================================================

function mnuFieldSettings_Callback(hObject, ~, handles) %#ok<DEFNU>
% Context menu to alter field-specific settings (name, wavelength, etc).
% FIXME: alter settingdlg to make this work somehow?

fieldID = get(gca,'UserData');  %quadrant
chID = find( find(handles.params.geometry)==fieldID );  %index in parameter list
input = [];

% Prompt for new values and verify validity.
if isequal(hObject,handles.mnuFieldSettings)
    prompt = {'Role:', 'Description:', 'Wavelength (nm):', 'Scale intensity by:'};

    if ~isempty(chID)
        currentopt = {handles.params.chNames{chID} handles.params.chDesc{chID} ...
                      num2str(handles.params.wavelengths(chID)) ...
                      num2str(handles.params.scaleFluor(chID)) };
    else
        currentopt = {'','','','1'};
    end

    answer = inputdlg(prompt, 'Change settings', 1, currentopt);
    if isempty(answer), return; end   %user hit cancel

    % Validate input
    if ~isempty(answer{1})
        if ~ismember(answer{1}, properties(TracesFret4)),
            errordlg( ['Invalid channel name ' answer{1}] );
            return;
        elseif isnan(str2double(answer{3}))
            errordlg('Invalid wavelength');
            return;
        end
        input = struct( 'chNames',answer{1}, 'chDesc',answer{2}, ...
           'wavelengths',str2double(answer{3}), 'scaleFluor',str2double(answer{4}) );
    end
end

% Save the new parameters
handles.params = gettraces_setch(handles.params, fieldID, input);
guidata(hObject,handles);
setAxTitles(handles);

% Reset picking
if ~isempty( findobj(handles.figure1,'type','line') )
    getTraces_Callback(hObject,[],handles);
end

% END FUNCTION mnuFieldSettings_Callback



function setAxTitles(handles)
% Set axes titles from imaging profile settings, including colors.

% Clear any existing titles
for i=1:numel(handles.ax),
    title(handles.ax(i),'');
end

% Create new titles
if numel(handles.ax)>0
    p = handles.params;
    
    for i=1:numel(p.idxFields)  %i is channel index
        chColor = Wavelength_to_RGB( p.wavelengths(i) );

        title( handles.ax( p.idxFields(i) ), ...
               sprintf('%s (%s) #%d',p.chNames{i},p.chDesc{i},i), ...
               'BackgroundColor',chColor, 'FontSize',10, 'Visible','on', ...
               'Color',(sum(chColor)<1)*[1 1 1] ); % White text for dark backgrounds.
    end
end
title(handles.axTotal,'Total Intensity', 'FontSize',10, 'Visible','on');

% END FUNCTION setTitles





%========================================================================
%=================  SOFTWARE ALIGNMENT MENU CALLBACKS  ==================
%========================================================================

function cboAlignMethod_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Change software alignment mode and re-pick molecules.
% 1=off, 2=load from file, 3=Auto (ICP), 4=memorize (keep using).

assert( ~isscalar(handles.params.geometry) );
sel = get(hObject,'Position');  %position within the menu, 1=top.

% Load alignment from file, if requested.
if sel==2
    % Ask user for a filename.
    [f,p] = uigetfile('*.mat','Select an alignment settings file');
    if f==0, return; end
    input = load( fullfile(p,f) );

    % Verify the alignment file is valid. FIXME: check for matching geometry!
    idx = find(  ~arrayfun(@(x)isempty(x.tform), input.alignment), 1 );  %first non-empty entry
    if isfield(input.alignment(idx).tform,'tdata'),
        errordlg('Alignment files from version 2.8 and earlier are not supported', 'Gettraces', 'modal');
        return;
    elseif ~isa(input.alignment(idx).tform,'affine2d'),
        errordlg('Invalid alignment: unrecognized tform class', 'Gettraces', 'modal');
        return;
    end

    handles.params.alignment = input.alignment;
    
elseif sel==4
    % Can only align memorize a valid alignment.
    if isempty(handles.params.alignment), return; end
end

% Re-pick molecules and update GUI with new mode.
handles.params.alignMethod = sel;
guidata(hObject,handles);
getTraces_Callback(hObject,[],handles);

set(findobj('Parent',handles.mnuAlign), 'Checked','off');
set(hObject, 'Checked','on');

% END FUNCTION cboAlignMethod_Callback



function btnSaveAlignment_Callback(~, ~, handles)  %#ok<DEFNU>
% Save current software alignment settings (which may be set to do nothing
% at all) to file so they can be reloaded later.

if isempty(handles.params.alignment) || isscalar(handles.params.geometry)
    return;  %can't save an invalid alignment
end

[p,f] = fileparts(handles.stkfile);
alignfile = fullfile( p, [f '_align.mat'] );

[f,p] = uiputfile('*.mat','Save software alignment settings',alignfile);

if f,
    % FIXME: should also remove abs_dev.
    alignment = rmfield( handles.params.alignment, {'quality'} );   %#ok<NASGU>
    save( fullfile(p,f), 'alignment' );
end

%end function btnSaveAlignment_Callback




%========================================================================
%===================  SPECTRAL CORRECTION CALLBACKS  ====================
%========================================================================

function btnCrosstalk_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Callback for button to set crosstalk correction settings.

params = handles.params;
nCh = numel(params.chNames);
if nCh<2, return; end  %no crosstalk

% Enumerate all possible crosstalk pairs
[src,dst] = find( triu(true(nCh),1) );

% Handle two-color special case giving a single value.
if numel(params.crosstalk)==1,
    c = zeros(2);
    c(1,2) = params.crosstalk;
    params.crosstalk = c;
end

% Prompt the user crosstalk parameters.
% Channel names MUST be in order of wavelength!
prompts  = cell(numel(src), 1);
defaults = cell(numel(src), 1);

for i=1:numel(src),
    prompts{i} = sprintf('%s to %s', params.chDesc{src(i)}, params.chDesc{dst(i)} );
    defaults{i} = num2str( params.crosstalk(src(i),dst(i)) );
end

result = inputdlg(prompts, 'gettraces', 1, defaults);
if isempty(result), return; end  %user hit cancel
result = cellfun(@str2double, result);

% Verify the inputs are valid and save them.
if any( isnan(result) | result>1 | result<0 ),
    warndlg('Invalid crosstalk values');
    return;
end

for i=1:numel(src),
    handles.params.crosstalk(src(i), dst(i)) = result(i);
end

guidata(hObject,handles);

%end function btnCrosstalk_Callback



function btnScaleAcceptor_Callback(hObject, ~, handles) %#ok<DEFNU>
% Set values for scaling the fluorescence intensity of acceptor channels so
% that they all have the same apparent brightness (gamma=1).

% Prompt the user for new multipliers for gamma correction.
params = handles.params;
prompts = cellfun( @(a,b)sprintf('%s (%s):',a,b), params.chNames, ...
                                       params.chDesc, 'UniformOutput',false );
defaults = cellfun(@num2str, num2cell(params.scaleFluor), 'UniformOutput',false);

result = inputdlg(prompts, 'gettraces', 1, defaults);
if isempty(result), return; end

% Verify the inputs are valid and save them.
result = cellfun(@str2double, result);
if any(isnan(result)),
    warndlg('Invalid scaling values');
    return;
end

handles.params.scaleFluor = to_row(result);
guidata(hObject,handles);

%end function btnCrosstalk_Callback





%========================================================================
%======================  AUTO-DETECT NEW FILES  =========================
%========================================================================

% --- Executes on button press in chkAutoBatch.
function chkAutoBatch_Callback(hObject, ~, ~)  %#ok<DEFNU>
% 
    
% If another timer is running, stop it.
fileTimer = timerfind('Name','gettraces_fileTimer');
if ~isempty(fileTimer),
    stop(fileTimer);
    delete(fileTimer);
end

% Start a new timer if requested
if strcmpi(get(hObject,'Checked'), 'off'),
    % Ask the user for a directory location
    targetDir = uigetdir('','Choose directory:');
    if targetDir==0, return; end
    disp(targetDir);
    
    % Start a thread that will periodically check for new movies every 5
    % seconds and process them automatically.
    fileTimer = timer('ExecutionMode','fixedSpacing','StartDelay',1, ...
                              'Name','gettraces_fileTimer', 'TimerFcn', ...
                              {@updateFileTimer,hObject,targetDir}, ...
                              'StopFcn',{@stopFileTimer,hObject}, ...
                              'Period',2.0,'BusyMode','drop');
    start(fileTimer);
    set(hObject, 'Checked','on');
else
    set(hObject, 'Checked','off');
end

% END FUNCTION chkAutoBatch_Callback


function stopFileTimer(~,~,hObject)
% Called on error during timer callback or when the timer is stopped.
if ishandle(hObject)
    handles = guidata(hObject);
    set(handles.mnuBatchOverwrite,'Checked','off');
end
% END FUNCTION stopFileTimer


function updateFileTimer(~,~,hObject,targetDir)
% This function runs each time the timer is fired, looking for any new
% movies that may have appeared on the path.

if ~ishandle(hObject),
fileTimer = timerfind('Name','gettraces_fileTimer');
    stop(fileTimer);
    delete(fileTimer);
    return;
end

batchmode_Callback( hObject, [], guidata(hObject), targetDir );
% END FUNCTION updateFileTimer
