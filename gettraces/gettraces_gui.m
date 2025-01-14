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

% Last Modified by GUIDE v2.5 24-Sep-2024 08:53:58


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
handles.stkData = MovieParser();

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
handles.stkData.params = profiles(handles.profile);
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

filter = {'*.tif;*.tiff;*.pma','All supported movie formats (*.tif,*.tiff,*.pma)'; ...
          '*.tif;*.tiff','TIFF image stacks (*.tif, *.tiff)'; ...
          '*.pma','Legacy raw frame data (*.pma)'
          '*.*','All Files (*.*)'};
[datafile,datapath] = uigetfile( filter, 'Open movie file', 'MultiSelect','on' );
if ~iscell(datafile)
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


% Load movie data
set(handles.figure1,'pointer','watch'); drawnow;
stkData = handles.stkData;
stkData.openStk(filename);


% Launch or update viewer
if ~isfield(handles,'viewer')
    handles.viewer = MovieViewer( stkData.chExtractor );
else
    handles.viewer.chExtractor = stkData.chExtractor;
end
handles.viewer.subtractBGImage = stkData.params.subtractBGImage;
handles.viewer.show(handles.panView);


guidata(hObject,handles);
set( handles.figure1, 'pointer','arrow');
set([handles.btnPickPeaks handles.mnuPick handles.mnuViewMetadata handles.mnuTirfProfile], 'Enable','on');
set([handles.btnSave handles.mnuFileSave handles.mnuFileSaveAs],'Enable','off');

%end function OpenStk



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
movieFiles = regexpdir(direct,'^.*\.(tiff?|stk|pma)$',recursive);
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
    nTraces(i) = size(handles.stkData.peaks,1);
end


% ----- Create log file with results
log_fid = fopen( fullfile(direct,'gettraces.log'), 'wt' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');
fprintf(log_fid, '%s', evalc('disp(handles.stkData.params)'));
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

function handles = btnHidePicks_Callback(~, ~, handles) %#ok<*DEFNU>
handles.viewer.hidePeaks();
%end function btnHidePicks_Callback



function handles = getTraces_Callback(~, ~, handles)
% Find peak locations from total intensity

% Do nothing if no file has been loaded. This function may be triggered by
% changing the settings fields before a file is open.
if isempty(handles.stkData.chExtractor),  return;  end

set(handles.tblAlignment, 'Data',{}, 'RowName',{});
set(handles.figure1,'pointer','watch'); drawnow;
stkData = handles.stkData;

% Clear any existing selection markers from previous calls.
handles.viewer.hidePeaks();

% Locate single molecules
try
    stkData.getPeaks();
catch e
    if ~strcmpi(e.identifier,'parfor_progressbar:cancelled')
        errordlg( ['Error: ' e.message], 'Gettraces' );
    end
    set(handles.figure1,'pointer','arrow'); drawnow;
    return;
end

set(handles.figure1,'pointer','arrow');
nmol = size(stkData.peaks,1);

% The alignment may involve shifting (or distorting) the fields to get a
% registered donor+acceptor field. Show this distorted imaged so the user
% can see what the algorithm is doing.
% val = get(handles.scaleSlider,'value');
% minimum = get(handles.scaleSlider,'min');
% val = max(val,minimum+1);

% set( handles.himshow(end), 'CData',stkData.total_t );
% set( handles.axTotal, 'CLim',[minimum*2 val*2] );
% 
set( handles.viewer.hImg(end), 'CData',stkData.total_t );
% set( handles.viewer.ax(end), 'CLim',[minimum*2 val*2] );  %incorrect for single color!

% If no alignment data given (for example in single-channel recordings),
% don't display any status messages.
set( handles.txtAlignWarning, 'Visible','off' );

if isempty(stkData.alignStatus),
    set(handles.panAlignment, 'ForegroundColor',[0 0 0], 'Title','No Alignment');
    
% Display alignment status to inform user if realignment may be needed.
% Format: translation deviation (x, y), absolute deviation (x, y)
else
    a = stkData.alignStatus;
    stkData.params.alignment = a;
    
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
    set( handles.tblAlignment, 'RowName',{stkData.chExtractor.channels(idxShow).name} );
    
    % If the alignment quality (confidence) is low, warn the user.
    methods = {'Alignment Disabled','Aligned from File', ...
               'Auto Aligned (ICP)','Alignment Memorized'};
    text = methods{stkData.params.alignMethod};
    
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

% Update GUI controls
set( handles.nummoles, 'String', sprintf('%d (of %d)',nmol, size(stkData.rejectedPicks,1)+nmol) );
set( [handles.btnSave handles.mnuFileSave handles.mnuFileSaveAs handles.mnuHidePeaks], 'Enable',onoff(nmol>0) );
set( [handles.mnuAlignSave handles.mnuAlignKeep], 'Enable',onoff(stkData.params.alignMethod>1) );
set( [handles.txtPSFWidth handles.txtOverlapStatus handles.txtWindowOverlap handles.txtIntegrationStatus], 'String','' );
if nmol<1, return; end


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
set(handles.txtSaturation,'String',sprintf('Saturated pixels: %0.1f%%',100*stkData.fractionSaturated), ...
    'ForegroundColor', (stkData.fractionSaturated>0.01)*[0.9 0 0]);

% Estimate the peak width from pixel intensity distribution.
set( handles.txtPSFWidth, 'String', sprintf('PSF size: %0.1f px',stkData.psfWidth) );
set( handles.txtPSFWidth, 'ForegroundColor', (stkData.psfWidth>stkData.params.nPixelsToSum-1)*[0.9 0 0] );

% Graphically show peak centers
handles.viewer.highlightPeaks( stkData.peaks, stkData.rejectedPicks, stkData.total_peaks, stkData.rejectedTotalPicks );


% end function



% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

% --- Executes on button press in mnuFileSave.
function mnuFileSave_Callback(~, ~, handles, prompt)

if isempty(handles.stkData.peaks)
    disp('No molecules selected! Skipping file');
    return;
end

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
    integrateAndSave(handles.stkData, filename);
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



% --------------------------------------------------------------------
function mnuProfiles_Callback(hObject, ~, handles)
% Called when any imaging profile is selected from the "Settings" menu.
% Apply imaging profile settings, including dividing the movie into quadrants.

if nargin<3,  handles = guidata(hObject);  end

% Save any changes to previous profile
%handles.profiles(handles.profile) = handles.stkData.params;

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
handles.stkData.params = params;
guidata(hObject,handles);

% Set all GUI to defaults of currently selected profile.
set( handles.mnuBatchRecursive,   'Checked', onoff(params.recursive)    );
set( handles.mnuBatchOverwrite,   'Checked', onoff(params.skipExisting) );

set( findobj('Parent',handles.mnuAlign), 'Checked','off' );
set( findobj('Parent',handles.mnuAlign,'Position',params.alignMethod), ...
     'Checked','on' );

% Enable alignment, crosstalk, and scale controls only in multi-color.
isMultiColor = ~isscalar(handles.stkData.params.geometry);
set( [handles.mnuAlign handles.mnuCrosstalk handles.mnuScaleAcceptor], ...
                       'Enable',onoff(isMultiColor) );
set( handles.tblAlignment, 'Visible',onoff(isMultiColor) );

set( handles.panAlignment, 'Title','Software Alignment', ...
                           'Visible',onoff(isMultiColor) );

% If a movie has already been loaded, reload movie with new setup.
if ~isempty(handles.stkData.chExtractor)
    handles = OpenStk( handles.stkfile, handles, hObject );
    set([handles.mnuAlignSave handles.mnuAlignKeep], 'Enable',onoff(params.alignMethod>1));
end

% END FUNCTION cboGeometry_Callback



% --------------------------------------------------------------------
function mnuSettingsCustom_Callback(hObject, ~, handles) %#ok<DEFNU>
% Called when Settings->Customize... menu clicked.
% Allows the user to temporarily alter settings for the current profile.

params = handles.stkData.params;
oldParams = params;

% Create options for bgTraceField
if handles.stkData.nChannels>0
    fopt = [{'None'} handles.stkData.chExtractor.channels.name];
else
    fopt = {'None'};
end

% Create dialog for changing imaging settings
prompt = {'Name:', 'Picking Threshold Value:', 'Use Automatic Threshold:', 'Auto picking sensitivity:', ...
          'Integration window size (px):', 'Integration neighbhorhood (px):', ...
          'Minimum separation (px):', ...
          'Donor blink detection method:', 'Background trace field:', ...
          'Frames to average for picking:', 'Subtract background image:'};
fields = {'name', 'don_thresh', 'autoThresh', 'thresh_std', 'nPixelsToSum', 'nhoodSize', ...
          'overlap_thresh', 'zeroMethod', 'bgTraceField', 'nAvgFrames','subtractBGImage'};
isInt = @(x)~isnan(x) && isreal(x) && isscalar(x) && x==floor(x);
isNum = @(x)~isnan(x) && isreal(x) && isscalar(x);
isBool = @(x)islogical(x);
types = {[],[],isBool,isNum,isInt,isInt,isNum,{'off','threshold','skm'},fopt,isInt,isBool};

if handles.profile > handles.nStandard
    prompt{1} = 'Name (clear to remove profile):';
end

if ~isempty(params.bgTraceField)
    params.bgTraceField = fopt( params.bgTraceField+1 );
else
    params.bgTraceField = 'None';
end

params = settingdlg(params, fields, prompt, types);
if isempty(params), return; end  %user hit cancel

if ~isempty(params.bgTraceField)
    params.bgTraceField = find( strcmp(params.bgTraceField,fopt) )-1;
    if params.bgTraceField==0, params.bgTraceField=''; end
end

% If a standard profile is renamed, make a custom setting instead.
if handles.profile < handles.nStandard  && ~strcmpi(params.name,handles.stkData.params.name)
    if isempty(params.name),
        errordlg('Built-in profiles cannot be deleted. Modify cascadeConstants instead');
        return;
    end
    
    % Overwrite if name matches an existing profile, except built-in.
    handles.stkData.params = params;
    handles.profiles(end+1) = params;
    handles.profile = numel(handles.profiles);

    % Update settings menu with new item.
    set( handles.hProfileMenu, 'Checked','off' );  %uncheck all
    handles.hProfileMenu(end+1) = uimenu( handles.mnuProfiles, 'Checked','on', ...
              'Label',handles.stkData.params.name, 'Callback',@mnuProfiles_Callback );
    guidata(hObject,handles);

% Modify existing custom profile
elseif ~isempty(params.name)
    set( handles.hProfileMenu(handles.profile), 'Label',params.name );
    handles.stkData.params = params;
    guidata(hObject,handles);
    
% Delete profile if name was cleared
else
    delete( handles.hProfileMenu(handles.profile) );
    handles.hProfileMenu(handles.profile) = [];
    handles.profiles(handles.profile) = [];
    mnuProfiles_Callback( handles.hProfileMenu(1), [], handles );
    return;
end

if ~isempty(handles.stkData.chExtractor)  %if a movie has been loaded
    
    % Reload movie if required
    if params.nAvgFrames~=oldParams.nAvgFrames || params.subtractBGImage~=oldParams.subtractBGImage
        OpenStk( handles.stkfile, handles, hObject );

    % If molecules were already picked and settings have changed, re-pick.
    elseif ~isempty(handles.stkData.peaks)
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
if isfield(stkData.chExtractor.movie.header,'MM')
    metadata = stkData.chExtractor.movie.header.MM(1);
    if isfield( stkData.chExtractor.movie.header, 'MM_wavelength' ),
        wv = stkData.chExtractor.movie.header.MM_wavelength;
        metadata.MM_wavelength = wv(wv>100);
    end
    
% Generic TIFF format metadata
else
    metadata = stkData.chExtractor.movie.header;
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
        data = mat2str(data);
    elseif ~ischar(data)
        %FIXME: need a recursive routine to draw out nested struct array.
        continue;
    end
    
    output{end+1} = sprintf( '%s:  %s', fname, data ); %#ok
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
filter = {'*.tif;*.tiff;*.pma','Movie Files (*.tif,*.tiff,*.pma)'; ...
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

chID = get(gca,'UserData');
if isempty(chID), return; end  %total intensity field
try
    handles.stkData.updateChannel( chID );
catch e
    errordlg(e.message);
    return;
end

% Update GUI, including picked molecules (if any).
handles.viewer.setAxTitles();
if ~isempty( findobj(handles.figure1,'type','line') )
    getTraces_Callback(hObject,[],handles);
end

% END FUNCTION mnuFieldSettings_Callback



function mnuFieldArrangement_Callback(~, ~, handles) %#ok<DEFNU>
% Context menu to reset which fields are in the current movie.

try
    success = handles.stkData.updateFieldArrangement();
catch e
    errordlg(e.message);
    return;
end

% Reload interface to handle changed number of axes.
if success
    % Reset GUI to state w/o picked molecules.
    % This should be a function!
    set( handles.txtAlignWarning, 'Visible','off' );
    set([handles.txtOverlapStatus handles.txtIntegrationStatus, ...
         handles.txtWindowOverlap handles.txtPSFWidth handles.nummoles], ...
         'String', '');
    set(handles.panAlignment, 'ForegroundColor',[0 0 0]);
    set(handles.tblAlignment, 'Data',{});
    set([handles.btnSave handles.mnuFileSave handles.mnuFileSaveAs],'Enable','off');

    % Update viewer
    handles.viewer.show(handles.panView);
end

% END FUNCTION mnuFieldSettings_Callback


function mnuUpdateAlex_Callback(~, ~, handles)
% Context menu to reset which fields are in the current movie.
% try
    handles.stkData.updateAlex();
% catch e
%     errordlg(e.message);
%     return;
% end
% END FUNCTION mnuUpdateAlex_Callback




%========================================================================
%=================  SOFTWARE ALIGNMENT MENU CALLBACKS  ==================
%========================================================================

function cboAlignMethod_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Change software alignment mode and re-pick molecules.
% 1=off, 2=load from file, 3=Auto (ICP), 4=memorize (keep using).

assert( handles.stkData.nChannels>1 );
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

    handles.stkData.params.alignment = input.alignment;
    
elseif sel==4
    % Can only align memorize a valid alignment.
    if isempty(handles.stkData.params.alignment), return; end
end

% Re-pick molecules and update GUI with new mode.
handles.stkData.params.alignMethod = sel;
guidata(hObject,handles);
getTraces_Callback(hObject,[],handles);

set(findobj('Parent',handles.mnuAlign), 'Checked','off');
set(hObject, 'Checked','on');

% END FUNCTION cboAlignMethod_Callback



function btnSaveAlignment_Callback(~, ~, handles)  %#ok<DEFNU>
% Save current software alignment settings (which may be set to do nothing
% at all) to file so they can be reloaded later.

if isempty(handles.stkData.params.alignment) || isscalar(handles.stkData.params.geometry)
    return;  %can't save an invalid alignment
end

[p,f] = fileparts(handles.stkfile);
alignfile = fullfile( p, [f '_align.mat'] );

[f,p] = uiputfile('*.mat','Save software alignment settings',alignfile);

if f,
    % FIXME: should also remove abs_dev.
    alignment = rmfield( handles.stkData.params.alignment, {'quality'} );   %#ok<NASGU>
    save( fullfile(p,f), 'alignment' );
end

%end function btnSaveAlignment_Callback




%========================================================================
%===================  SPECTRAL CORRECTION CALLBACKS  ====================
%========================================================================

% --------------------------------------------------------------------
function mnuCrosstalk_Callback(~, ~, handles)
% Callback for button to set crosstalk correction settings.
handles.stkData.updateCrosstalk();


% --------------------------------------------------------------------
function mnuScaleAcceptor_Callback(~, ~, handles)
% Set values for scaling the fluorescence intensity of acceptor channels so
% that they all have the same apparent brightness (gamma=1).
handles.stkData.updateScaling();





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



% --------------------------------------------------------------------
function mnuLoadROI_Callback(hObject, ~, handles)
% Load region of interest descriptor from a file.

[f,p] = uigetfile('*.txt','gettraces:load ROI');
if ~isequal(f,0)
    data = load( fullfile(p,f) );
    nX = handles.stkData.chExtractor.nX;
    nY = handles.stkData.chExtractor.nY;
    
    assert( size(data,2)==2, 'Invalid ROI; must be list x-y coordinates' );
    assert( all(data(:)>0) & all(data(:,1)<=nX) & all(data(:,2)<=nY), 'Invalid ROI: outside of image bounds' );
    handles.stkData.roi = data;
    getTraces_Callback(hObject,[],handles);
end

%END FUNCTION mnuLoadROI


% --------------------------------------------------------------------
function mnuDrawROI_Callback(hObject, ~, handles)
% Manually draw an ROI and select spots within it.

ROI = drawpolygon(handles.viewer.ax(end), 'FaceAlpha',0, 'Color','red');
handles.stkData.roi = ROI.Position;
delete(ROI);
getTraces_Callback(hObject,[],handles);

% END FUNCTION mnuDrawROI_Callback


% --------------------------------------------------------------------
function mnuClearROI_Callback(hObject, ~, handles)
% Clear any existing ROI

handles.stkData.roi = [];
getTraces_Callback(hObject,[],handles);

% END FUNCTION mnuClearROI_Callback


% --------------------------------------------------------------------
function mnuCenterQuadROI_Callback(hObject, ~, handles)

nX = handles.stkData.chExtractor.nX;
nY = handles.stkData.chExtractor.nY;

handles.stkData.roi = [nX/4 nY/4; nX/4 3*nY/4; 3*nX/4 3*nY/4; 3*nX/4 nY/4];
getTraces_Callback(hObject,[],handles);

% END FUNCTION mnuCenterQuadROI_Callback


% --------------------------------------------------------------------
function mnuCircleROI_Callback(hObject, ~, handles)

nX = handles.stkData.chExtractor.nX;
nY = handles.stkData.chExtractor.nY;

roi = images.roi.Circle('Center',[nX/2,nY/2], 'Radius',nX/2);
handles.stkData.roi = roi.Vertices;
getTraces_Callback(hObject,[],handles);

% END FUNCTION mnuCircleROI_Callback



% --------------------------------------------------------------------
function mnuAutoMeta_Callback(hObject, ~, handles)
% If checked, automatically adjust analysis settings to the experimental
% parameters (field arrangement, channel assignments, etc.) as described
% in movie metadata. If not selected, a series of movies can be analyzed
% with the same settings.

status = strcmpi( get(hObject,'Checked'), 'on' );
handles.stkData.autoFromMetadata = ~status;
set(hObject,'Checked',~status);
msgbox('Reload movie for setting to take effect.');

% END FUNCTION mnuAutoMeta_Callback


% --------------------------------------------------------------------
function mnuCountOverTime_Callback(~, ~, handles)

handles.stkData.getParticleCountTrace();

% END FUNCTION mnuCountOverTime_Callback

