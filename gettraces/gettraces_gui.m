function varargout = gettraces_gui(varargin)
% GETTRACES M-file for gettraces.fig
%      
%      Converts Metamorph .stk files into single-molecule fluorescence and
%      FRET traces. The picking algorithm was taken from gui3.m, written by
%      Harold Kim. I fixed problems with the batch mode, and reconfigured
%      the startup to open .stks instead of .pmas. I also reconfigured the
%      output, so that now the result is a new .traces file with a unique 
%      id for each molecule donor fluorescence, acceptor fluorescence, and
%      fret for each molecule.
%
%      To run in batch mode, first load a single .stk file and adjust the
%      threshold and txtOverlap rejection. Then select batch mode and choose a
%      directory. All the .stk files in that directory will be converted to
%      separate .traces files using the established threshold and txtOverlap
%      rejection. Batch mode is not recursive.

%   Copyright 2007-2016 Cornell University All Rights Reserved.

% Last Modified by GUIDE v2.5 01-Nov-2017 13:23:16


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

% Load colormap for image viewer
handles.colortable = gettraces_colormap();

% Choose default command line output for gettraces
handles.output = hObject;

% Load built-in and saved user-custom.
profiles = constants.gettraces_profiles;
nStandard = numel(profiles);

if ispref('SPARTAN','gettraces_customProfiles')
    profiles = [profiles getpref('SPARTAN','gettraces_customProfiles')];
% else
%     error('STUB');
end

% Add profiles from cascadeConstants to settings menu.
for i=1:numel(profiles),
    if i<=numel(constants.gettraces_profiles), pos=i; else, pos=i+2; end
    
    hMenu(i) = uimenu(handles.mnuProfiles, 'Label',profiles(i).name, ...
                      'Position',pos, 'Callback',@mnuProfiles_Callback); %#ok<AGROW>
end

% Put customization menu items in the correct spots.
set(handles.mnuSettingsCustom,'Position',nStandard+1);
set(handles.mnuSettingsSave,  'Position',nStandard+2);
if numel(hMenu)>nStandard
    set( hMenu(nStandard+1), 'Separator','on' );
end

handles.profile = constants.gettraces_defaultProfile;  %index to current profile (FIXME: rename)
set( hMenu(handles.profile), 'Checked','on' );



% Context menus for field-specific settings (names, wavelength, etc).
hZoom = zoom(handles.figure1);
set(hZoom, 'UIContextMenu', handles.mnuField);
zoom(handles.figure1,'on');

% Setup default values for parameter values -- 2-color FRET.
handles.hProfileMenu = hMenu;
handles.profiles = profiles;
handles.alignment = [];  %current alignment parameters (status)
guidata(hObject, handles);

% Set up GUI elements to reflect the internal parameter values.
cboGeometry_Callback(hObject, [], handles);

% END FUNCTION gettraces_OpeningFcn




% --- Outputs from this function are returned to the command line.
function varargout = gettraces_OutputFcn(~, ~, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------- OPEN SINGLE MOVIE (CALLBACK) ---------------- %

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
if isappdata(handles.figure1,'stkData')
    rmappdata(handles.figure1,'stkData');
end
set(handles.figure1,'pointer','watch'); drawnow;

[stkData] = gettraces( filename, handles.params );
setappdata(handles.figure1,'stkData',stkData);

set( handles.sldScrub, 'Min',1, 'Max',stkData.movie.nFrames, 'Value',1, ...
     'SliderStep',[1/stkData.movie.nFrames,0.02] );

% Setup slider bar (adjusting maximum value in image, initially 2x max)
% low = min(min(stkData.stk_top));
sort_px = sort(stkData.stk_top(:));
val = sort_px( floor(0.98*numel(sort_px)) );
high = min( ceil(val*10), 32000 );  %uint16 maxmimum value

set(handles.scaleSlider,'min',0, 'max',high, 'value', val);
set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));

% Create axes for sub-fields
delete( findall(handles.figure1,'type','axes') );  %remvoe old axes
ax = [];
fields = subfield(stkData.stk_top-stkData.background, handles.params.geometry);
spopt = {handles.panView, 'XTickLabel','', 'YTickLabel','', 'Units','normalized'};

switch handles.params.geometry
case 1,  %---- Single-Channel (full-chip)
    handles.axTotal = axes( spopt{:}, 'Position',[0.02 0 0.95 0.95] );
    
case 2,  %---- Dual-Channel (L/R)
    ax(1)           = axes( spopt{:}, 'Position',[0      0     0.325  0.95] );
    ax(2)           = axes( spopt{:}, 'Position',[0.335  0     0.325  0.95] );
    handles.axTotal = axes( spopt{:}, 'Position',[0.67   0     0.325  0.95] );
    
case 3,  %---- Dual-Channel (T/B)
    ax(1)           = axes( spopt{:}, 'Position',[0.02   0.51  0.48   0.47] );
    ax(2)           = axes( spopt{:}, 'Position',[0.02   0     0.48   0.47] );
    handles.axTotal = axes( spopt{:}, 'Position',[0.51   0.25  0.48   0.47] );
    
case 4,  %---- Quad-Channel (TL/TR/BL/BR)
    ax(1)           = axes( spopt{:}, 'Position',[0.0    0.5   0.325  0.47] );
    ax(2)           = axes( spopt{:}, 'Position',[0.335  0.5   0.325  0.47] );
    ax(4)           = axes( spopt{:}, 'Position',[0.335  0     0.325  0.47] );
    ax(3)           = axes( spopt{:}, 'Position',[0      0     0.325  0.47] );
    handles.axTotal = axes( spopt{:}, 'Position',[0.67   0.25  0.325  0.47] );
otherwise
    error('Invalid imaging geometry');
end

% Show fluorescence fields for all channels
handles.ax = ax;
handles.himshow = [];

if handles.params.geometry>1,
    for i=1:numel(fields),
        h = imshow( fields{i}, [0 val], 'Parent',ax(i) );
        set(h,'UIContextMenu',handles.mnuImage);
        colormap(ax(i),handles.colortable);
        set(ax(i),'UserData',i);
        handles.himshow(i) = h;
    end
    setAxTitles(handles);
    linkaxes( [to_row(ax) handles.axTotal] );
end

guidata(hObject,handles);

% Show total fluorescence channel
total = sum( cat(3,fields{:}), 3);
h = imshow( total, [0 val*2], 'Parent',handles.axTotal );
set(h,'UIContextMenu',handles.mnuImage);
colormap(handles.axTotal,handles.colortable);
title(handles.axTotal,'Total Intensity', 'FontSize',10);

set(handles.figure1,'pointer','arrow');
set([handles.btnPickPeaks handles.mnuPick handles.mnuViewMetadata], 'Enable','on');
set([handles.btnSave handles.mnuFileSave handles.mnuFileSaveAs],'Enable','off');

%end function OpenStk



% --- Executes on slider movement.
function sldScrub_Callback(hObject, ~, handles) %#ok<DEFNU>
% Allows user to scroll through the movie

% Stop any currently playing movies.
set(handles.btnPlay,'String','Play');

% Read frame of the movie
stkData  = getappdata(handles.figure1,'stkData');
idxFrame = round( get(hObject,'Value') );
frame  = double( stkData.movie.readFrame(idxFrame) );
fields = subfield(frame-stkData.background, handles.params.geometry);

for i=1:numel(fields)
    set( handles.himshow(i), 'CData',fields{i} );
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

stkData = getappdata(handles.figure1,'stkData');
startFrame = round( get(handles.sldScrub,'Value') );

for i=startFrame:stkData.movie.nFrames
    frame  = double( stkData.movie.readFrame(i) );
    fields = subfield(frame-stkData.background, handles.params.geometry);
    
    for f=1:numel(fields)
        set( handles.himshow(f), 'CData',fields{f} );
    end
    
    set(handles.sldScrub,'Value',i);
    drawnow;
    
    % Terminate early if the user clicks the 'Stop' button.
    state = get(handles.btnPlay,'String');
    if strcmpi(state,'Play'), return; end
end

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






% ------------------------ PICK INTENSITY PEAKS ------------------------ %

function handles = getTraces_Callback(hObject, ~, handles)
%

%----- Find peak locations from total intensity

% Do nothing if no file has been loaded. This function may be triggered by
% changing the settings fields before a file is open.
if ~isfield(handles, 'stkfile')
    return;
end
set(handles.tblAlignment, 'Data',{}, 'RowName',{});

set(handles.figure1,'pointer','watch'); drawnow;

% Locate single molecules
stkData = getappdata(handles.figure1,'stkData');
stkData = getPeaks(stkData, handles.params);
stkData = getIntegrationWindows(stkData, handles.params);
setappdata(handles.figure1,'stkData',stkData);

% The alignment may involve shifting (or distorting) the fields to get a
% registered donor+acceptor field. Show this distorted imaged so the user
% can see what the algorithm is doing.
val = get(handles.scaleSlider,'value');
minimum = get(handles.scaleSlider,'min');
val = max(val,minimum+1);

imshow( stkData.total_t, [minimum*2 val*2], 'Parent',handles.axTotal );
colormap( handles.axTotal, handles.colortable );
title( handles.axTotal, 'Total Intensity' );


% If no alignment data given (for example in single-channel recordings),
% don't display any status messages.
set( handles.txtAlignWarning, 'Visible','off' );

if ~isfield(stkData,'alignStatus') || isempty(stkData.alignStatus),
    set(handles.panAlignment, 'Visible','off', 'ForegroundColor', [0 0 0]);
    handles.alignment = [];
    
% Display alignment status to inform user if realignment may be needed.
% Format: translation deviation (x, y), absolute deviation (x, y)
else
    a = stkData.alignStatus;
    handles.alignment = a;
    
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
    % FIXME: try to color individual rows according to degree of misalignment,
    % for example using HTML tags in text (rg, <html><b><font color='red',
    % or color="#FF00FF"). <center> tag might also be useful.
    if any( [a.abs_dev] > 0.25 ),
        d = max( [a.abs_dev] );
        set( handles.tblAlignment, 'ForegroundColor', [(3/2)*min(2/3,d) 0 0] );
        set( handles.panAlignment, 'ForegroundColor', min(1,d*[1 0 0]) );
    else
        set( handles.tblAlignment, 'ForegroundColor', [0 0 0] );
        set( handles.panAlignment, 'ForegroundColor', [0 0 0] );
    end
    
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
% percentTracesWinOverlap = 100*sum(stkData.fractionWinOverlap>0)/numel(stkData.fractionWinOverlap);

set(  handles.txtWindowOverlap, 'String', ...
      sprintf('Residual win. overlap: %0.1f%%', percentWinOverlap)  );
set( handles.txtWindowOverlap, 'ForegroundColor', (percentWinOverlap>10)*[0.9 0 0] );


% Get (approximate) average fraction of fluorescence collected within the
% integration window of each molecule. Set the text color to red where the
% intensity is not well collected at the current integration window size.
eff = 100*stkData.integrationEfficiency(:,handles.params.nPixelsToSum);
eff = nanmean(eff);
set(  handles.txtIntegrationStatus, 'String', ...
      sprintf('Intensity collected: %0.0f%% ', eff)  );
set( handles.txtIntegrationStatus, 'ForegroundColor', (eff<70)*[0.9 0 0] );


% Estimate the peak width from pixel intensity distribution.
eff = stkData.integrationEfficiency;
eff = eff( ~any(isnan(eff')), : );  %ignore NaN values, which can happen in with empty fields.
decay = zeros( size(eff,1), 1 ); %number pixels to integrate to get 70% intensity integrated.

for i=1:size(eff,1),
    decay(i) = find( eff(i,:)>=0.7, 1, 'first' );%default 3-color channel assignments.
end

set( handles.txtPSFWidth, 'String', sprintf('PSF size: %0.1f px',mean(decay)) );
set( handles.txtPSFWidth, 'ForegroundColor', (mean(decay)>handles.params.nPixelsToSum)*[0.9 0 0] );



%----- Graphically show peak centers
handles.x = stkData.peaks(:,1);
handles.y = stkData.peaks(:,2);
handles.total_x = stkData.total_peaks(:,1);
handles.total_y = stkData.total_peaks(:,2);

handles.rx = stkData.rejectedPicks(:,1);
handles.ry = stkData.rejectedPicks(:,2);
handles.rtotal_x = stkData.rejectedTotalPicks(:,1);
handles.rtotal_y = stkData.rejectedTotalPicks(:,2);

highlightPeaks( handles );

% Update GUI controls
handles.num = numel(handles.total_x);
set( handles.nummoles, 'String', sprintf('%d (of %d)',handles.num, ...
                                    numel(handles.rtotal_x)+handles.num) );

set( [handles.btnSave handles.mnuFileSave handles.mnuFileSaveAs ...
                                     handles.mnuHidePeaks], 'Enable','on');

set(handles.figure1,'pointer','arrow');
guidata(hObject,handles);

% end function



function highlightPeaks(handles)
% Draw circles around each selected fluorescence spot.
    
% Parameters for drawing circles around each detected molecule.
style1 = {'LineStyle','none','marker','o','color','w'};
style2 = {'LineStyle','none','marker','o','color','y'};

% Parameters for overlapped molecules that won't actually be recorded.
% These are darker to de-emphasize them, but yet still be visible.
style1b = {'LineStyle','none','marker','o','color',[0.4,0.4,0.4]};
style2b = {'LineStyle','none','marker','o','color',[0.4,0.4,0.0]};

stkData = getappdata(handles.figure1,'stkData');
[nrow,ncol] = size(stkData.stk_top);

% Clear any existing selection markers from previous calls.
delete(findobj(handles.figure1,'type','line'));
    
if handles.params.geometry==2, %dual-channel, L/R
    % Assign each spot to the corresponding field image.
    indL  = find( handles.x<=(ncol/2) );
    indR  = find( handles.x> (ncol/2) );
    rindL = find( handles.rx<=(ncol/2) );
    rindR = find( handles.rx> (ncol/2) );
    
    % Draw markers on selection points (donor side)
    axes(handles.ax(1));
    line( handles.x(indL),   handles.y(indL),   style1{:},  'Parent',handles.ax(1) );
    line( handles.rx(rindL), handles.ry(rindL), style1b{:}, 'Parent',handles.ax(1) );

    % Draw markers on selection points (acceptor side)
    axes(handles.ax(2));
    line( handles.x(indR)-(ncol/2),   handles.y(indR),   style1{:},  'Parent',handles.ax(2) );
    line( handles.rx(rindR)-(ncol/2), handles.ry(rindR), style1b{:}, 'Parent',handles.ax(2) );
    
elseif handles.params.geometry==3, %dual-channel, T/B
    error('STUB');
    
elseif handles.params.geometry==4, %quad-channel
    % Assign each spot to a quadrant. FIXME: can this be a loop?
    indUL = find( handles.x<=(ncol/2) & handles.y<=(nrow/2) );
    indLL = find( handles.x<=(ncol/2) & handles.y> (nrow/2) );
    indLR = find( handles.x> (ncol/2) & handles.y> (nrow/2) );
    indUR = find( handles.x> (ncol/2) & handles.y<=(nrow/2) );
    
    rindUL = find( handles.rx<=(ncol/2) & handles.ry<=(nrow/2) );
    rindLL = find( handles.rx<=(ncol/2) & handles.ry> (nrow/2) );
    rindLR = find( handles.rx> (ncol/2) & handles.ry> (nrow/2) );
    rindUR = find( handles.rx> (ncol/2) & handles.ry<=(nrow/2) );
    
    axes(handles.ax(1));
    line( handles.x(indUL),   handles.y(indUL),   style1{:},  'Parent',handles.ax(1)  );
    line( handles.rx(rindUL), handles.ry(rindUL), style1b{:}, 'Parent',handles.ax(1) );
    
    axes(handles.ax(3));
    line( handles.x(indLL),   handles.y(indLL)-(nrow/2),   style1{:},  'Parent',handles.ax(3)  );
    line( handles.rx(rindLL), handles.ry(rindLL)-(nrow/2), style1b{:}, 'Parent',handles.ax(3) );
    
    axes(handles.ax(4));
    line( handles.x(indLR)-(ncol/2),   handles.y(indLR)-(nrow/2),   style1{:},  'Parent',handles.ax(4)  );
    line( handles.rx(rindLR)-(ncol/2), handles.ry(rindLR)-(nrow/2), style1b{:}, 'Parent',handles.ax(4) );
    
    axes(handles.ax(2));
    line( handles.x(indUR)-(ncol/2),   handles.y(indUR),   style1{:},  'Parent',handles.ax(2)  );
    line( handles.rx(rindUR)-(ncol/2), handles.ry(rindUR), style1b{:}, 'Parent',handles.ax(2) );
end

% Draw markers on selection points (total intensity composite image).
axes(handles.axTotal);
line( handles.total_x,  handles.total_y,  style2{:},  'Parent',handles.axTotal  );
line( handles.rtotal_x, handles.rtotal_y, style2b{:}, 'Parent',handles.axTotal );

% end function highlightPeaks


% --- Executes on button press in btnHidePicks.
function btnHidePicks_Callback(~, ~, handles)  %#ok<DEFNU>
% Hide the circles drawn to indicate molecule locations.
delete(findobj(handles.figure1,'type','line'));
% END FUNCTION btnHidePicks_Callback



% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

% --- Executes on button press in mnuFileSave.
function mnuFileSave_Callback(~, ~, handles, prompt)

if nargin<4, prompt=false; end

if iscell(handles.stkfile),
    filename = handles.stkfile{1};
else
    filename = handles.stkfile;
end

% Remove file extension for multi-part movies.
filename = regexprep(filename,'-file[0-9]*.[A-Za-z0-9]*$','');

% Prompt user for target filename, if called from the "File->Save As" menu.
if prompt,
    [p,f] = fileparts(filename);
    [f,p] = uiputfile('*.traces', 'Save traces as:', fullfile(p,[f '.rawtraces']));
    if f==0, return; end
    handles.params.outFilename = fullfile(p,f);
end

% Integrate fluorophore point-spread functions, generate fluorescence
% traces, and save to file.
set(handles.figure1,'pointer','watch'); drawnow;

try
    stkData = getappdata(handles.figure1,'stkData');
    if isfield(stkData,'peaks')
        integrateAndSave(stkData, filename, handles.params);
    else
        gettraces( stkData, handles.params, filename );
    end
catch e
    if strcmpi(e.identifier,'parfor_progressbar:cancelled')
        disp('Gettraces: Operation cancelled by user.');
    elseif ~cascadeConstants('debug')
        errordlg( ['Error: ' e.message], 'Gettraces' );
    else
        rethrow(e);
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

if handles.params.geometry==1, %Single-channel recordings
    set( handles.axTotal, 'CLim',[minimum val] );
else
    set( handles.axTotal, 'CLim',[minimum*2 val*2] );
end

set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));
guidata(hObject,handles);

% END FUNCTION txtMaxIntensity_Callback


function mnuProfiles_Callback(hObject, ~)
% Called when any imaging profile is selected from the "Settings" menu.
% Marks only the active profile and breaks up the movie.

handles = guidata(hObject);
set(handles.hProfileMenu,'Checked','off');
set(hObject, 'Checked','on');
handles.profile = find(hObject==handles.hProfileMenu);

cboGeometry_Callback(hObject, [], handles);

% END FUNCTION mnuProfiles_Callback



% --- Executes on selection change in cboGeometry.
function handles = cboGeometry_Callback(hObject, ~, handles)
%

% If running, stop the "auto detect" timer. Otherwise, it may be triggered by
% the change in settings.
fileTimer = timerfind('Name','gettraces_fileTimer');
if ~isempty(fileTimer),
    stop(fileTimer);
    delete(fileTimer);
    set(handles.mnuBatchAuto,'Checked','off');
end

% Get parameter values associated with the selected profile.
% Warning: if cascadeConstants is changed to add a new profile or rearrange
% profiles, this can have unpredictable effects...
params = handles.profiles(handles.profile);
handles.params = params;

% Set all GUI to defaults of currently selected profile.
set( handles.mnuBatchRecursive,   'Checked', onoff(params.recursive)    );
set( handles.mnuBatchOverwrite,   'Checked', onoff(params.skipExisting) );

set( findobj('Parent',handles.mnuAlign), 'Checked','off' );
set( findobj('Parent',handles.mnuAlign,'Position',params.alignMethod), ...
     'Checked','on' );
set([handles.mnuAlignSave handles.mnuAlignKeep], 'Enable',onoff(params.alignMethod>1));


% Enable alignment, crosstalk, and scale controls only in multi-color.
isMultiColor = numel(handles.params.idxFields)>1;
set( [handles.mnuAlign handles.btnCrosstalk handles.btnScaleAcceptor], ...
                       'Enable',onoff(isMultiColor) );
set( handles.tblAlignment, 'Visible',onoff(isMultiColor) );

set( handles.panAlignment, 'Title','Software Alignment', ...
                           'Visible',onoff(isMultiColor) );

% If a movie has already been loaded, reload movie with new setup.
if isfield(handles,'stkfile'),
    handles = OpenStk( handles.stkfile, handles, hObject );
end

guidata(hObject,handles);



% --------------------------------------------------------------------
function mnuSettingsCustom_Callback(hObject, ~, handles) %#ok<DEFNU>
% Called when Settings->Customize... menu clicked.
% Allows the user to temporarily alter settings for the current profile.

% All options for fluorescence channel identifiers. See subfield_mask.m.
% FIXME: ideally this should be the natural field name (e.g., 'acceptor').
fopt = { {''}, {'','L','R'}, {'','T','B'}, {'','TL','TR','BL','BR'} };
fopt = fopt{ handles.params.geometry };

prompt = {'Name:','Threshold (0 for auto):', 'Integration window size (px):', ...
          'Minimum separation (px):', 'ADU/photon conversion:', ...
          'Donor blink detection method:', 'Integration neighbhorhood (px):', ...
          'Background trace field'};
fields = {'name','don_thresh', 'nPixelsToSum', 'overlap_thresh', ...
          'photonConversion', 'zeroMethod', 'nhoodSize','bgTraceField'};
types = {[],[],[],[],[],{'off','threshold','skm'},[],fopt};
params = settingdlg(handles.params,fields,prompt,types);
if isempty(params), return; end

% Update profile name
constants = cascadeConstants;
nStandard = numel(constants.gettraces_profiles);

if handles.profile < nStandard  && ~strcmpi(params.name,handles.params.name)
    % Standard profile's name was modified, so we assume they want to create a
    % new custom-named profile, rather than giving an error.
    % FIXME
    errordlg('Can''t rename built-in profiles');
    return;
end

% Save modified profile list.
set( handles.hProfileMenu(handles.profile), 'Label',params.name );
handles.params = params;
handles.profiles(handles.profile) = params;
guidata(hObject,handles);

setpref('SPARTAN','gettraces_customProfiles',handles.profiles(nStandard+1:end));

% END FUNCTION mnuSettingsCustom_Callback



% --------------------------------------------------------------------
function mnuSettingsSave_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Save current settings as a new imaging profile in the Settings menu.

constants = cascadeConstants;
nStandard = numel(constants.gettraces_profiles);

% Ask for the new name
newName = inputdlg('New profile name:', mfilename, 1, {handles.params.name});
if isempty(newName), return; end  %user hit cancel
newName = newName{1};

% Overwrite if name matches an existing profile, except built-in.
% FIXME: make this a shortcut to adding a new custom profile?
if strcmpi( newName, handles.params.name )
    if handles.profile <= nStandard,
        errordlg('Built-in profiles cannot be modified. Alter cascadeConstants.m instead.');
        return;
    else
        handles.profiles(handles.profile) = handles.params;
    end

% Save as a new profile if the name is new.
else
    handles.params.name = newName;
    handles.profiles(end+1) = handles.params;
    handles.profile = numel(handles.profiles);
    
    % Update settings menu with new item.
    set( handles.hProfileMenu, 'Checked','off' );  %uncheck all
    handles.hProfileMenu(end+1) = uimenu( handles.mnuProfiles, 'Checked','on', ...
              'Label',handles.params.name, 'Callback',@mnuProfiles_Callback );
    
    % If this is the first in the list, add separator.
    if handles.profile == nStandard+1
        set( handles.hProfileMenu(end), 'Separator','on' );
    end
end

% Save the updated profile list.
setpref('SPARTAN','gettraces_customProfiles',handles.profiles(nStandard+1:end));
guidata(hObject,handles);

% END FUNCTION mnuSettingsSave_Callback



% --- Executes on button press in btnMetadata.
function btnMetadata_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Display a simple diaglog with the MetaMorph metadata for the first frame.

stkData = getappdata(handles.figure1,'stkData');

if isempty(stkData) || ~isfield(stkData.movie.header,'MM')
    if isfield( stkData.movie.header,'ImageDescription' ),
        msgbox( stkData.movie.header.ImageDescription, 'ImageDescription' );
    else
        % This movie doesn't have any metadata
        set( hObject, 'Enable','off' );
    end
    return;
end

% Grab extra MetaMorph fields outside the MM struct.
metadata = stkData.movie.header.MM(1);

if isfield( stkData.movie.header, 'MM_wavelength' ),
    wv = stkData.movie.header.MM_wavelength;
    metadata.MM_wavelength = wv(wv>100);
end

% Convert structure data into text to display.
fields = fieldnames(metadata);
flen = max( cellfun(@numel,fields) );
nFields = numel(fields);

output = cell( nFields, 1 );

for i=1:nFields, %for each field
    fname = fields{i};
    data  = metadata.(fname);
    
    if isnumeric(data),
        data = num2str(data);
    end
    
    output{i} = sprintf( ['%' num2str(flen) 's:  %s'], fname, data );
end

% Display the dialog.
msgbox( output, 'MetaMorph metadata' );

% END FUNCTION btnMetadata_Callback



function mnuTirfProfile_Callback(~, ~, handles) %#ok<DEFNU>
% Executes when the "View->Illumination Profile" menu is clicked.
% FIXME: no direct way to get the output filename, or to know whether it has
% been saved? For now just guess.

[p,f] = fileparts(handles.stkfile);
outname = fullfile(p,[f '.rawtraces']);
if exist(outname,'file'),
    tirfProfile(outname);
else
    disp('No .rawtraces file. Save one first to get a profile');
end

% END FUNCTION mnuTirfProfile_Callback




%========================================================================
%======================  FIELD SETTINGS CALLBACKS  ======================
%========================================================================

function mnuFieldSettings_Callback(hObject, ~, handles) %#ok<DEFNU>
% Context menu to alter field-specific settings (name, wavelength, etc).
% FIXME: alter settingdlg to make this work somehow?
% FIXME: some redundant code here.

idxField = get(gca,'UserData');  %quadrant
if isempty(idxField), return; end  %total intensity field
chID = find(handles.params.idxFields==idxField);  %index in parameter list

% Prompt for new values and verify validity.
prompt = {'Role (ex: donor):', 'Description (ex: Cy3):', 'Wavelength (nm):', ...
          'Scale intensity by:'};
if ~isempty(chID)
    currentopt = {handles.params.chNames{chID} handles.params.chDesc{chID} ...
                  num2str(handles.params.wavelengths(chID)) ...
                  num2str(handles.params.scaleFluor(chID)) };
else
    currentopt = {'','','','1'};
end

answer = inputdlg(prompt, 'Change settings', 1, currentopt);
if isempty(answer), return; end   %user hit cancel

if ~ismember(answer{1}, properties(TracesFret4)),
    errordlg( ['Invalid channel name ' answer{1}] );
    return;
elseif isnan(str2double(answer{3}))
    errordlg('Invalid wavelength');
    return;
end

% Save the new parameters
input = struct( 'chNames',answer{1}, 'chDesc',answer{2}, ...
                'wavelengths',str2double(answer{3}), 'scaleFluor',str2double(answer{4}) );
handles.params = gettraces_setch(handles.params, idxField, input);
guidata(hObject,handles);
setAxTitles(handles);

% Reset picking
if ~isempty( findobj(handles.figure1,'type','line') )
    getTraces_Callback(hObject,[],handles);
end

% END FUNCTION mnuFieldSettings_Callback



function mnuFieldRemove_Callback(hObject, ~, handles) %#ok<DEFNU>
% Context menu to remove a field from consideration.

idxField = get(gca,'UserData');  %quadrant
if isempty(idxField), return; end

fieldID = find(handles.params.idxFields==idxField,1);  %index in parameter list
if isempty(fieldID), return; end

if numel(handles.params.chNames)<2,
    warndlg('At least one channel must be defined');
    return;
end

handles.params = gettraces_setch(handles.params, idxField, []);
guidata(hObject,handles);
setAxTitles(handles);

% Reset picking
if ~isempty( findobj(handles.figure1,'type','line') )
    getTraces_Callback(hObject,[],handles);
end

% END FUNCTION mnuFieldRemove_Callback



function setAxTitles(handles)
% Set axes titles from imaging profile settings, including colors.

% Clear any existing titles
for i=1:numel(handles.ax),
    title(handles.ax(i),'');
end

% Create new titles
p = handles.params;
fields = p.idxFields;

for i=1:numel(fields)
    chColor = Wavelength_to_RGB(p.wavelengths(i));

    title(handles.ax(fields(i)), ...
          sprintf('%s (%s) #%d',p.chNames{i},p.chDesc{i},fields(i)), ...
          'BackgroundColor',chColor, 'FontSize',10, ...
          'Color',(sum(chColor)<1)*[1 1 1] ); % White text for dark backgrounds.
end

% END FUNCTION setTitles





%========================================================================
%=================  SOFTWARE ALIGNMENT MENU CALLBACKS  ==================
%========================================================================

function cboAlignMethod_Callback(hObject, ~, handles)  %#ok<DEFNU>
% Change software alignment mode and re-pick molecules.
% 1=off, 2=load from file, 3=Auto (ICP), 4=memorize (keep using).

assert(handles.params.geometry>1);
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

    handles.alignment = input.alignment; %GUI state
    handles.params.alignment = input.alignment; %gettraces() input
    
elseif sel==4
    handles.params.alignment = rmfield( handles.alignment, {'quality'} );
    %FIMXE: only run getTraces_Callback if molecules not picked yet.
end

% Re-pick molecules and update GUI with new mode.
handles.params.alignMethod = sel;
getTraces_Callback(hObject,[],handles);

set(findobj('Parent',handles.mnuAlign), 'Checked','off');
set(hObject, 'Checked','on');
set([handles.mnuAlignSave handles.mnuAlignKeep], 'Enable',onoff(sel>1));

% END FUNCTION cboAlignMethod_Callback



function btnSaveAlignment_Callback(~, ~, handles)  %#ok<DEFNU>
% Save current software alignment settings (which may be set to do nothing
% at all) to file so they can be reloaded later.

assert( isfield(handles,'alignment') && ~isempty(handles.alignment) && handles.params.geometry>1 );

[p,f] = fileparts(handles.stkfile);
alignfile = fullfile( p, [f '_align.mat'] );

[f,p] = uiputfile('*.mat','Save software alignment settings',alignfile);

if f,
    % FIXME: should also remove abs_dev.
    alignment = rmfield( handles.alignment, {'quality'} );   %#ok<NASGU>
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




