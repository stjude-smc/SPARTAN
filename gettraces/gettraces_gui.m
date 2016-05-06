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

% Last Modified by GUIDE v2.5 06-May-2016 17:12:37


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

% Add profiles from cascadeConstants to settings menu.
profiles = {constants.gettraces_profiles.name};
for i=1:numel(profiles),
    h = uimenu(handles.mnuProfiles, 'Label',profiles{i}, 'Callback',@mnuProfiles_Callback);
    if i>1 && lower(profiles{i}(1))~=lower(profiles{i-1}(1)),
        set(h, 'Separator','on');
    end
    if i==constants.gettraces_defaultProfile,
        set(h, 'Checked','on');
    end
end
set(handles.mnuSettingsCustom,'Position',numel(profiles)+1);

% Setup default values for parameter values -- 2-color FRET.
handles.alignment = [];  %current alignment parameters (status)
handles.profile = constants.gettraces_defaultProfile;
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
if numel(filename)==1,
    d = dir( [f '*.tif*'] );
    if numel(d)>1,
        d = regexpi( {d.name}, [f '(-file[0-9]*)?\.tiff?$'], 'match' );
        filename = [d{:}];
    end
end

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
set(handles.panAlignment, 'Title','Software Alignment', 'ForegroundColor',[0 0 0]);
set(handles.tblAlignment, 'Data',{});


set(handles.figure1,'pointer','watch'); drawnow;

% Clear the original stack to save memory
if isappdata(handles.figure1,'stkData')
    rmappdata(handles.figure1,'stkData');
end

% Load movie data
[stkData] = gettraces( filename, handles.params );
handles.stk_top = stkData.stk_top;

% Since the image stack is very large, it is stored in ApplicationData
% instead of GUIData for memory efficiency
setappdata(handles.figure1,'stkData', stkData);


% Setup slider bar (adjusting maximum value in image, initially 2x max)
% low = min(min(handles.stk_top));
sort_px = sort(handles.stk_top(:));
low=0;
val = sort_px( floor(0.98*numel(sort_px)) );
high = min( ceil(val*10), 32000 );  %uint16 maxmimum value

set(handles.scaleSlider,'min',low);
set(handles.scaleSlider,'max',high);
set(handles.scaleSlider,'value', val);
set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));

%
image_t = handles.stk_top-stkData.background;
[nrow,ncol] = size(image_t);
fields = {};  ax = [];

%---- Show fields for Single-Channel (full-chip) recordings.
if handles.params.geometry==1,
    handles.total_t = image_t;
    
    % Show full field of view.
    handles.axTotal = subplot( 1,1,1, 'Parent',handles.panView, 'Position',[0.2 0 0.6 0.95] );
    
    
%---- Show fields for Dual-Channel (half-chip, L/R) recordings.
elseif handles.params.geometry==2,    
    fields{1} = image_t(:,1:ncol/2);  %left, donor
    fields{2} = image_t(:,(ncol/2)+1:end);  %right, acceptor
    handles.total_t = fields{1}+fields{2};

    % Setup axes
    handles.axDonor    = subplot( 1,3,1, 'Parent',handles.panView, 'Position',[0.025 0 0.3 0.95] );
    handles.axAcceptor = subplot( 1,3,2, 'Parent',handles.panView, 'Position',[0.350 0 0.3 0.95] );
    ax = [handles.axDonor,handles.axAcceptor];
        
    % Show total intensity image
    handles.axTotal = subplot( 1,3,3, 'Parent',handles.panView, 'Position',[0.675 0 0.3 0.95] );
    
    
%---- Show fields for Quad-Channel (quarter-chip, L/R/U/D) recordings.
elseif handles.params.geometry>2,
    % Split up channels and combine, assuming perfect alignment.
    fields{1} = image_t( 1:nrow/2, 1:ncol/2 );              %upperLeft
    fields{2} = image_t( 1:nrow/2, (ncol/2)+1:end );        %upperRight
    fields{3} = image_t( (nrow/2)+1:end, 1:ncol/2 );        %lowerLeft
    fields{4} = image_t( (nrow/2)+1:end, (ncol/2)+1:end );  %lowerRight
    handles.total_t = fields{1} + fields{2} + fields{3} + fields{4};
    
    % Setup axes
    handles.axUL = subplot( 2,3,1, 'Parent',handles.panView, 'Position',[0.025 0.5   0.3 0.45] );
    handles.axUR = subplot( 2,3,2, 'Parent',handles.panView, 'Position',[0.35  0.5   0.3 0.45] );
    handles.axLL = subplot( 2,3,4, 'Parent',handles.panView, 'Position',[0.025 0.025 0.3 0.45] );
    handles.axLR = subplot( 2,3,5, 'Parent',handles.panView, 'Position',[0.35  0.025 0.3 0.45] );
    ax = [handles.axUL handles.axUR handles.axLL handles.axLR];

    % Also plot the total intensity image.
    handles.axTotal = subplot( 2,3,3, 'Parent',handles.panView, 'Position',[0.7 0.25 0.3 0.45] );
end

handles.ax = ax;

% Show fluorescence fields for all channels
chColors = Wavelength_to_RGB(handles.params.wavelengths);
chNames = handles.params.chNames;

for i=1:numel(fields),
    % i is the index of physical CCD chip locations.
    % idxCh is the corresponding index into list of channels (there may be none).
    idxCh = find( handles.params.idxFields==i ); 
        
    imshow( fields{i}, [low val], 'Parent',ax(i) );
    colormap(ax(i),handles.colortable);

    if ~isempty(idxCh) && ~isempty(chNames{idxCh}),
        % Give each field a title with the background color matching the
        % wavelength of that channel.        
        h = title( ax(i), [chNames{idxCh} ' (' handles.params.chDesc{idxCh} ') #' num2str(i)], ...
                   'BackgroundColor',chColors(idxCh,:), 'FontSize',10 );

        % Use white text for very dark background colors.
        if sum(chColors(idxCh,:)) < 1,
            set(h,'Color',[1 1 1]);
        end
        
        set(ax(i), 'UserData',idxCh);
    end
end

% Link axes so zooming one zooms all.
linkaxes( [ax handles.axTotal] );
    
% Context menus for field-specific settings (names, wavelength, etc).
% FIXME: this also adds one to "total intensity"...
zoom(handles.figure1,'off');
hZoom = zoom(handles.figure1);
set(hZoom, 'UIContextMenu', handles.mnuField);
zoom(handles.figure1,'on');

% Show total fluorescence channel
imshow( handles.total_t, [low*2 val*2], 'Parent',handles.axTotal );
colormap(handles.axTotal,handles.colortable);
zoom(handles.axTotal,'on');
title(handles.axTotal,'Total Intensity', 'FontSize',10);

% Finish up
set(handles.figure1,'pointer','arrow');
set([handles.btnPickPeaks handles.mnuPick handles.mnuViewMetadata], 'Enable','on');

guidata(hObject,handles);

%end function OpenStk




% --------------- OPEN ALL STK'S IN DIRECTORY (BATCH) --------------- %

% --- Executes on button press in batchmode.
function batchmode_Callback(hObject, ~, handles,direct)
% direct     target directory location to look for new files

% Get input parameter values
skipExisting = strcmpi( get(handles.mnuBatchOverwrite,'Checked'), 'on');
recursive = strcmpi( get(handles.mnuBatchRecursive,'Checked'), 'on');

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

set(handles.figure1,'pointer','watch'); drawnow;

% Locate single molecules
stkData = getappdata(handles.figure1,'stkData');
[stkData,peaks] = gettraces( stkData, handles.params );

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
    set(handles.panAlignment, 'Title','Software Alignment');
    handles.alignment = [];
    
% Display alignment status to inform user if realignment may be needed.
% Format: translation deviation (x, y), absolute deviation (x, y)
else
    a = stkData.alignStatus;
    handles.alignment = a;
    
    if handles.params.alignMethod==1,
        text = 'Alignment Deviation:';
    else
        text = 'Alignment Applied:';
    end
    
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
    if isfield(a,'quality'),
        if any( [a.quality]<1.1 & [a.quality]>0 ),
            text = [text sprintf(' (low confidence!)')];
        end
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

set(  handles.txtOverlapStatus, 'String', ...
      sprintf('Molecules rejected: %0.0f%%', percentOverlap)  );
set( handles.txtOverlapStatus, 'ForegroundColor', (percentOverlap>30)*[0.9 0 0] );


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
handles.x = peaks(:,1);
handles.y = peaks(:,2);
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

set( [handles.btnSave handles.mnuFileSave handles.mnuHidePeaks], 'Enable','on');

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

[nrow,ncol] = size(handles.stk_top);

% Clear any existing selection markers from previous calls.
delete(findobj(handles.figure1,'type','line'));
    
if handles.params.geometry==2, %dual-channel
    % Assign each spot to the corresponding field image.
    indL  = find( handles.x<=(ncol/2) );
    indR  = find( handles.x> (ncol/2) );
    rindL = find( handles.rx<=(ncol/2) );
    rindR = find( handles.rx> (ncol/2) );
    
    % Draw markers on selection points (donor side)
    axes(handles.axDonor);
    line( handles.x(indL),   handles.y(indL),   style1{:},  'Parent',handles.axDonor );
    line( handles.rx(rindL), handles.ry(rindL), style1b{:}, 'Parent',handles.axDonor );

    % Draw markers on selection points (acceptor side)
    axes(handles.axAcceptor);
    line( handles.x(indR)-(ncol/2),   handles.y(indR),   style1{:},  'Parent',handles.axAcceptor );
    line( handles.rx(rindR)-(ncol/2), handles.ry(rindR), style1b{:}, 'Parent',handles.axAcceptor );
    
elseif handles.params.geometry>2, %quad-channel
    % Assign each spot to a quadrant. FIXME: can this be a loop?
    indUL = find( handles.x<=(ncol/2) & handles.y<=(nrow/2) );
    indLL = find( handles.x<=(ncol/2) & handles.y> (nrow/2) );
    indLR = find( handles.x> (ncol/2) & handles.y> (nrow/2) );
    indUR = find( handles.x> (ncol/2) & handles.y<=(nrow/2) );
    
    rindUL = find( handles.rx<=(ncol/2) & handles.ry<=(nrow/2) );
    rindLL = find( handles.rx<=(ncol/2) & handles.ry> (nrow/2) );
    rindLR = find( handles.rx> (ncol/2) & handles.ry> (nrow/2) );
    rindUR = find( handles.rx> (ncol/2) & handles.ry<=(nrow/2) );
    
    axes(handles.axUL);
    line( handles.x(indUL),   handles.y(indUL),   style1{:},  'Parent',handles.axUL  );
    line( handles.rx(rindUL), handles.ry(rindUL), style1b{:}, 'Parent',handles.axUL );
    
    axes(handles.axLL);
    line( handles.x(indLL),   handles.y(indLL)-(nrow/2),   style1{:},  'Parent',handles.axLL  );
    line( handles.rx(rindLL), handles.ry(rindLL)-(nrow/2), style1b{:}, 'Parent',handles.axLL );
    
    axes(handles.axLR);
    line( handles.x(indLR)-(ncol/2),   handles.y(indLR)-(nrow/2),   style1{:},  'Parent',handles.axLR  );
    line( handles.rx(rindLR)-(ncol/2), handles.ry(rindLR)-(nrow/2), style1b{:}, 'Parent',handles.axLR );
    
    axes(handles.axUR);
    line( handles.x(indUR)-(ncol/2),   handles.y(indUR),   style1{:},  'Parent',handles.axUR  );
    line( handles.rx(rindUR)-(ncol/2), handles.ry(rindUR), style1b{:}, 'Parent',handles.axUR );
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
function mnuFileSave_Callback(~, ~, handles)

set(handles.figure1,'pointer','watch'); drawnow;

if iscell(handles.stkfile),
    filename = handles.stkfile{1};
else
    filename = handles.stkfile;
end

% Remove file extension for multi-part movies.
filename = regexprep(filename,'-file[0-9]*.[A-Za-z0-9]*$','');

% Integrate fluorophore point-spread functions, generate fluorescence
% traces, and save to file.
stkData = getappdata(handles.figure1,'stkData');
gettraces( stkData, handles.params, filename );

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

if handles.params.geometry==1, %Single-channel recordings
    set( handles.axTotal, 'CLim',[minimum val] );
    
elseif handles.params.geometry==2, %Dual-channel recordings
    set( [handles.axDonor handles.axAcceptor], 'CLim',[minimum val] );
    set( handles.axTotal, 'CLim',[minimum*2 val*2] );
    
elseif handles.params.geometry>2,
    set( [handles.axUL handles.axUR handles.axLL handles.axLR], ...
                                                   'CLim',[minimum val] );
    set( handles.axTotal, 'CLim',[minimum*2 val*2] );
end

set(handles.txtMaxIntensity,'String', sprintf('%.0f',val));
guidata(hObject,handles);

% END FUNCTION txtMaxIntensity_Callback


function mnuProfiles_Callback(hObject, ~)
% Called when any imaging profile is selected from the "Settings" menu.
% Marks only the active profile and breaks up the movie.

handles = guidata(hObject);
set(findobj('Parent',get(hObject,'Parent')), 'Checked','off');
set(hObject, 'Checked','on');
handles.profile = get(hObject,'Position');

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
constants = cascadeConstants;
params = constants.gettraces_profiles(handles.profile);
handles.params = params;


% Set all GUI to defaults of currently selected profile.
set( handles.mnuBatchRecursive,   'Checked', onoff(params.recursive)    );
set( handles.mnuBatchOverwrite,   'Checked', onoff(params.skipExisting) );

set( findobj('Parent',handles.mnuAlign), 'Checked','off' );
set( findobj('Parent',handles.mnuAlign,'Position',params.alignMethod), ...
     'Checked','on' );
set([handles.mnuAlignSave handles.mnuAlignKeep], 'Enable',onoff(params.alignMethod>1));

set( handles.edScaleAcceptor, 'String', num2str(params.scaleAcceptor) );
set( handles.txtDACrosstalk,  'String', num2str(params.crosstalk) );

% Enable alignment, crosstalk, and scale controls only in multi-color.
nCh = numel(handles.params.idxFields);
set( [handles.mnuAlign handles.txtDACrosstalk ...
      handles.btnCrosstalk handles.edScaleAcceptor  ...
      handles.btnScaleAcceptor],  'Enable',onoff(nCh>1) );

set(handles.txtDACrosstalk,   'Visible', onoff(numel(params.crosstalk)<2) );
set(handles.btnCrosstalk,     'Visible', onoff(numel(params.crosstalk)>1)  );
set(handles.edScaleAcceptor,  'Visible', onoff(numel(params.scaleAcceptor)<2) );
set(handles.btnScaleAcceptor, 'Visible', onoff(numel(params.scaleAcceptor)>1)  );


% If a movie has already been loaded, reload movie with new setup.
if isfield(handles,'stkfile'),
    handles = OpenStk( handles.stkfile, handles, hObject );
end

guidata(hObject,handles);



% --------------------------------------------------------------------
function mnuSettingsCustom_Callback(hObject, ~, handles) %#ok<DEFNU>
% Called when Settings->Customize... menu clicked.
% Allows the user to temporarily alter settings for the current profile.

prompt = {'Threshold (0 for auto):', 'Integration window size (px):', ...
          'Minimum separation (px):', 'ADU/photon conversion:', ...
          'Donor blink detection method:', 'Integration neighbhorhood (px):'};
fields = {'don_thresh', 'nPixelsToSum', 'overlap_thresh', ...
          'photonConversion', 'zeroMethod', 'nhoodSize'};
types = {[],[],[],[],{'off','threshold','skm'}};
handles.params = settingsDialog(handles.params,fields,prompt,types);
guidata(hObject,handles);

% END FUNCTION mnuSettingsCustom_Callback



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



% --------------------------------------------------------------------
function mnuFieldSettings_Callback(hObject, ~, handles) %#ok<DEFNU>
% Context menu to alter field-specific settings (name, wavelength, etc).
% FIXME: alter settingsDialog to make this work somehow?

fieldID = get(gca,'UserData');
if isempty(fieldID), return; end  %total intensity field

% Prompt for new values and verify validity.
prompt = {'Role (ex: donor):', 'Description (ex: Cy3):', 'Excitation wavelength (nm):'};
currentopt = {handles.params.chNames{fieldID} handles.params.chDesc{fieldID} ...
              num2str(handles.params.wavelengths(fieldID)) };

answer = inputdlg(prompt, 'Change settings', 1, currentopt);
if isempty(answer), return; end   %user hit cancel

if ~ismember(answer{1}, properties(TracesFret4)),
    errordlg( ['Invalid channel name ' answer{1}] );
    return;
end

% Save the new parameters
handles.params.chNames{fieldID}     = answer{1};
handles.params.chDesc{fieldID}      = answer{2};
handles.params.wavelengths(fieldID) = str2double(answer{3});
guidata(hObject,handles);

% Update the GUI. FIXME should call a function.
axID = handles.params.idxFields(fieldID);
chColor = Wavelength_to_RGB(handles.params.wavelengths(fieldID));

h = title( handles.ax(axID), [answer{1} ' (' answer{2} ') #' num2str(axID)], ...
           'BackgroundColor',chColor, 'FontSize',10 );
       
% White text for dark backgrounds.
set(h,'Color',(sum(chColor)<1)*[1 1 1]); 

% END FUNCTION mnuFieldSettings_Callback




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
    sel = 2; %basically the same thing
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
% When there are more than 2 channels, the crosstalk is more than just a
% scalar and can't be represented in the text box easily, so this button
% will launch a dialog to show all the possible parameter values and allow
% the user to change them.

params = handles.params;
assert( params.geometry>1 && numel(params.crosstalk)>1 );


% Prompt the user crosstalk parameters for 3-color.
% We assume channel names are in order of wavelength.
src = [1 1 2]; %Cy3->Cy5, Cy3->Cy7, Cy5->Cy7
dst = [2 3 3];

prompts  = cell( numel(src), 1 );
defaults = cell( numel(src), 1 );

for i=1:numel(src)
    prompts{i} = sprintf('%s (%s) -> %s (%s)', ...
                         params.chNames{src(i)}, params.chDesc{src(i)}, ...
                         params.chNames{dst(i)}, params.chDesc{dst(i)} );
    defaults{i} = num2str( params.crosstalk(src(i),dst(i)) );
end

result = inputdlg( prompts, 'Enter crosstalk values', 1, defaults );
if isempty(result), return; end


% Process the new crosstalk values
crosstalk = zeros( numel(params.chNames) );

for i=1:numel(src),
    c = str2double( result{i} );
    
    if isnan(c) || c>1 || c<0,
        fprintf( 'Error: invalid crosstalk value %s)n', result{i} );
        return;
    end
    
    crosstalk( src(i), dst(i) ) = c;
end

handles.params.crosstalk = crosstalk;
guidata(hObject,handles);

%end function btnCrosstalk_Callback



function btnScaleAcceptor_Callback(hObject, ~, handles) %#ok<DEFNU>
% Set values for scaling the fluorescence intensity of acceptor channels so
% that they all have the same apparent brightness (gamma=1). This button
% should only be visit for multi-color FRET where there are multiple
% channels that should be scaled.
% FIXME: should be extended to factor and donor2 channels??

params = handles.params;

idx = ~cellfun( @isempty, strfind(params.chNames, 'acceptor') );
if ~isfield(params,'scaleAcceptor') || isempty(params.scaleAcceptor),
    params.scaleAcceptor = ones(numel(idx), 1);
end

% Prompt the user for new multipliers for gamma correction.
prompts = cellfun( @(a,b)sprintf('%s (%s):',a,b), params.chNames(idx), ...
                                 params.chDesc(idx), 'UniformOutput',false );
defaults = cellfun(@num2str, num2cell(params.scaleAcceptor), 'UniformOutput',false);

result = inputdlg(prompts, 'Gettraces: scale acceptor', 1, defaults);
if isempty(result), return; end

% Verify the inputs are valid and save them.
result = cellfun(@str2double, result);
if any(isnan(result)),
    disp('Gettraces: ignoring invalid scale acceptor values.');
    return;
end

handles.params.scaleAcceptor = result;
guidata(hObject,handles);

%end function btnCrosstalk_Callback





%========================================================================
%======================  AUTO-DETECT NEW FILES  =========================
%========================================================================

% --- Executes on button press in chkAutoBatch.
function chkAutoBatch_Callback(hObject, ~, handles)  %#ok<DEFNU>
% 
    
% If another timer is running, stop it.
fileTimer = timerfind('Name','gettraces_fileTimer');

if ~isempty(fileTimer),
    stop(fileTimer);
    delete(fileTimer);
    %disp('Timer deleted.');
end

% Start a new timer if requested
if strcmpi(get(hObject,'Checked'), 'off'),
    % Ask the user for a directory location
    targetDir = uigetdir('','Choose directory:');
    if targetDir==0, return; end
    disp(targetDir);
    
    % Force "skip existing"; process only new movies in each iteration.
    set(handles.mnuBatchOverwrite,'Checked','on');
    
    % Start a thread that will periodically check for new movies every 5
    % seconds and process them automatically.
    %disp('Timer started.');
    fileTimer = timer('ExecutionMode','fixedSpacing','StartDelay',1, ...
                              'Name','gettraces_fileTimer', 'TimerFcn', ...
                              {@updateFileTimer,hObject,targetDir}, ...
                              'StopFcn',{@stopFileTimer,hObject}, ...
                              'Period',2.0,'BusyMode','drop');
    start(fileTimer);
    set(hObject, 'Checked','on');
    %FIXME: add an error/stop function to clear the checkbox.
else
    set(hObject, 'Checked','off');
end

% END FUNCTION chkAutoBatch_Callback


function stopFileTimer(~,~,hObject)
% Called on error during timer callback or when the timer is stopped.
handles = guidata(hObject);
set(handles.mnuBatchOverwrite,'Checked','off');
% END FUNCTION stopFileTimer


function updateFileTimer(~,~,hObject,targetDir)
% This function runs each time the timer is fired, looking for any new
% movies that may have appeared on the path.
% disp('Timer fired');
batchmode_Callback( hObject, [], guidata(hObject), targetDir );
% END FUNCTION updateFileTimer




