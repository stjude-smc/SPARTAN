function varargout = gettraces(varargin)
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
%      threshold and overlap rejection. Then select batch mode and choose a
%      directory. All the .stk files in that directory will be converted to
%      separate .traces files using the established threshold and overlap
%      rejection. Batch mode is not recursive.
%
%               -JBM 12/06
%
% Changes:
% 
% 0. Much faster run time, less memory usage (AppData), comments
% 1. Overlap distances now calculated for centroid of fluor peak.
% 2. Full bg correction used for peak picking
% 3. Peak search uses plus-shaped window instead of 3x3
% 
%           -DT 10/22

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help gettraces

% Last Modified by GUIDE v2.5 20-Aug-2008 15:45:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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
function gettraces_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gettraces (see VARARGIN)

% Choose default command line output for gettraces
handles.output = hObject;

% Initialize global variables
% handles.CODE_PATH='/home/dsterry/Documents/cornell/blanchard/code/';
handles.CODE_PATH = '';

% Update handles structure
guidata(hObject, handles);

set(handles.saveTraces,'Enable','off');
set(handles.getTraces,'Enable','off');
% set(handles.getTracesCy5,'Enable','off');
% set(handles.batchmode,'Enable','off');


% --- Outputs from this function are returned to the command line.
function varargout = gettraces_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------------- OPEN STK MOVIE (CALLBACK) --------------------- %

% --- Executes on button press in openstk.
function openstk_Callback(hObject, eventdata, handles)
% hObject    handle to openstk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load stk file
[datafile,datapath]=uigetfile( ...
    '*.stk;*.stk.bz2','Choose a stk file');
if datafile==0, return; end

handles.stkfile = strcat(datapath,datafile);

% Load the movie
handles = OpenStk( handles.stkfile, handles, hObject );


set(handles.getTraces,'Enable','on');
% set(handles.getTracesCy5,'Enable','on');
% set(handles.batchmode,'Enable','on');

guidata(hObject,handles);



% --------------------- OPEN STK MOVIE --------------------- %
function handles = OpenStk(filename, handles, hObject)


set(handles.txtFilename,'String',filename);

% Clear the original stack to save memory
setappdata(handles.figure1,'stk',[]);

% Load colormap for image viewer
fid=fopen('colortable.txt','r');
colortable=fscanf(fid,'%d',[3 256]);
colortable=colortable'/255;
fclose(fid);
handles.colortable = colortable;

% Load movie data
[stk,handles.stk_top,handles.background] = OpenStk2(filename);

 % Since the image stack is very large, it is stored in ApplicationData
 % instead of GUIData for memory efficiency
setappdata(handles.figure1,'stk', stk); 
clear stk;


% Setup slider bar (adjusting maximum value in image, initially 2x max)
low = min(min(handles.stk_top));
high = max(max(handles.stk_top));
high = ceil(high*2);

set(handles.scaleSlider,'min',low);
set(handles.scaleSlider,'max',high);
set(handles.scaleSlider,'value', (low+high)/2);

%
image_t    = handles.stk_top-handles.background+mean2(handles.background);
[nrow,ncol] = size(image_t);

donor_t    = image_t(:,1:ncol/2);
acceptor_t = image_t(:,(ncol/2)+1:end);
total_t    = donor_t+acceptor_t;

% Show donor image
axes(handles.axDonor);
imshow( donor_t, [low (high+low)/2] );
colormap(colortable);  zoom on;

% Show acceptor image
axes(handles.axAcceptor);
imshow( acceptor_t, [low (high+low)/2] );
colormap(colortable);  zoom on;

% Show total intensity image
axes(handles.axTotal);
imshow( total_t, [low*2 (high+low)] );
colormap(colortable);  zoom on;

linkaxes( [handles.axDonor handles.axAcceptor handles.axTotal] );


% Finish up
guidata(hObject,handles);




function [stk,stk_top,background] = OpenStk2(filename)

% If the movie is compressed, deflate it to a temporary location first
if strfind(filename,'.stk.bz2'),
    
    % decompress
    z_fname = filename;
    status = system( ['bunzip2 --keep "' z_fname '"'] );
    if status ~= 0, 
        error('Error uncompressing STK: %s',z_fname);
    end
    
    % Load stk movie -- stk(X,Y,frame number)
    filename = strrep(z_fname,'.stk.bz2','.stk');
    [stk, stkX, stkY, Nframes] = tiffread(filename);
    
    % Delete the compressed movie -- we no longer need it!
    delete( z_fname );

else
    % Load stk movie -- stk(X,Y,frame number)
    [stk, stkX, stkY, Nframes] = tiffread(filename);
end


% Create an average image of the first 10 frames (stk_top)
averagesize = min([10 Nframes]);
movie = stk(:,:,1:averagesize);
stk_top = mean(movie,3);



% Create an estimated background image by:
% 1. Divide the image into den*den squares
% 2. For each square, find the fifth lowest number
% 3. Rescaling these values back to the original image size
if stkX == 128
    den=4;  %not optimized
elseif stkX == 170
    den=5;
%     set(handles.overlap,'String', 2.1);
elseif stkX == 256
    den=8;
%     set(handles.overlap,'String', 2.5);
elseif stkX == 512
    den=16;
%     set(handles.overlap,'String', 4.5);
end

background = stk_top;  %**
temp = zeros( stkY/den, stkX/den );

for i=1:stkY/den
    for j=1:stkX/den
        sort_temp = background(den*(i-1)+1:den*i,den*(j-1)+1:den*j);
        sort_temp = reshape(sort_temp,1,den*den);  % make into a vector
        sort_temp = sort(sort_temp);

        temp(i,j) = sort_temp(den);  % get the 1/den % smallest value
    end
end

% Rescale the image back to the actual size
background=imresize(temp,[stkX stkY],'bilinear');
% handles.background = mean( stk(:,:,end-4:end), 3 );


% END FUNCTION OpenStk




% --------------- OPEN ALL STK'S IN DIRECTORY (BATCH) --------------- %

% --- Executes on button press in batchmode.
function batchmode_Callback(hObject, eventdata, handles)
% hObject    handle to batchmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

skipExisting = get(handles.chkOverwrite,'Value');
don_thresh = str2double(get(handles.donthresh,'String'));
overlap_thresh = str2double(get(handles.overlap,'String'));

%iptsetpref('imshowInitialMagnification','fit');

direct=uigetdir('','Choose directory:');
if direct==0, return; end

% Create header for log file
log_fid = fopen( [direct filesep 'gettraces.log'], 'w' );
fprintf(log_fid, 'Donor Thresh = %.1f\n',don_thresh);
fprintf(log_fid, 'Overlap = %.1f\n',overlap_thresh);
% fprintf(log_fid, 'Acceptor Thresh = %.1f\n',acc_thresh);

fprintf(log_fid,'\n%s\n\n%s\n%s\n\n%s\n',date,'DIRECTORY',direct,'FILES');

disp(direct);

% Create progress bar at bottom of window
set(handles.txtProgress,'String','Creating traces, please wait...');

% Get list of files in current directory (option: and all subdirectories)
recursive = get(handles.chkRecursive,'Value');

if recursive
    stk_files  = rdir([direct filesep '**' filesep '*.stk*']);
else
    stk_files  = rdir([direct filesep '*.stk*']);
end

tic;
% For each file in the user-selected directory
i = 1;
for file = stk_files'
    handles.stkfile = file.name;
    
    % Skip if previously processed (.traces file exists)
    stk_fname = strrep(file.name,'.bz2','');
    [p,name]=fileparts(stk_fname);
    traceFname = [p filesep name '.traces'];
    
    if skipExisting && exist(traceFname,'file'),
        disp( ['Skipping (already processed): ' stk_fname] );
        continue;
    end
    
    % Load STK file
    handles = OpenStk(handles.stkfile,handles, hObject);
    
    % Pick molecules using default parameter values
    handles = getTraces_Callback(hObject, [], handles);
    
    % Save the traces to file
    saveTraces_Callback(hObject, [], handles);
    
    % Save entry in log file
    fprintf(log_fid, '%.0f %s\n', handles.num, file.name);
    
    
%     text = sprintf('Creating traces: %.0f%%', 100*(i/size(stk_files,1)) );
%     set(handles.txtProgress,'String',text);
%     guidata(hObject,handles);
end

set(handles.txtProgress,'String','Finished.');
fclose(log_fid);
toc
guidata(hObject,handles);







% --------------- PICK MOLECULES CALLBACKS --------------- %


%------------- Pick Cy3 spots (CALLBACK) ----------------- 
% --- Executes on button press in getTraces.
function handles = getTraces_Callback(hObject, eventdata, handles)
% hObject    handle to getTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% NOTE: ignores 3 pixels from all edges, only looks at left side (cy3)


%----- Load picking parameters
overlap_thresh = str2double(get(handles.overlap,'String'));

% Donor threshold is relative to background level.
% If not chosen by user, is 7x standard dev. of background
don_thresh = str2double(get(handles.donthresh,'String'));

image_t=handles.stk_top-handles.background;
[nrow ncol] = size(image_t);

if don_thresh==0
   % would be better to use std over time than space?
   don_thresh = 10*std2( handles.background(1+3:nrow-3,1+3:ncol/2-3) );
else
    don_thresh = don_thresh-mean2(handles.background);
end


%----- Find peak locations from total intensity

[handles.x,handles.y,total_t] = getPeaks( image_t, don_thresh, overlap_thresh,handles );
handles.num = numel(handles.x)/2;

% axes(handles.axTotal);
% % low = 1236;  high = 28966;
% imshow( total_t );
% colormap(handles.colortable);  zoom on;

%----- GUI stuff

% Clear selection markers
delete(findobj(gcf,'type','line'));

% Draw markers on selection points (donor side)
l = 1:2:numel(handles.x);

axes(handles.axDonor);
line(handles.x(l),handles.y(l),'LineStyle','none','marker','o','color','y','EraseMode','background');

% Draw markers on selection points (acceptor side)
ll = 2:2:numel(handles.x);

axes(handles.axAcceptor);
line(handles.x(ll)-(ncol/2),handles.y(ll),'LineStyle','none','marker','o','color','w','EraseMode','background');

% Draw markers on selection points (acceptor side)
axes(handles.axTotal);
line(handles.x(l),handles.y(l),'LineStyle','none','marker','o','color','w','EraseMode','background');

% Update GUI controls
set(handles.nummoles,'String',num2str(handles.num));
set(handles.saveTraces,'Enable','on');

guidata(hObject,handles);





%------------- Pick fluorescence intensity peaks ----------------- 
function [picksX,picksY,total_t] = getPeaks( image_t, threshold, overlap_thresh, handles )
% Localizes the peaks of molecules in the Cy3 channel and infers the
% positions of the Cy5 peaks by applying a transofmration to the Cy3 peak
% positions.
%
% picksX - X-coords of all molcules, as Cy3,Cy5,Cy3,Cy5 in order
% picksY - Y-coords ...

[nrow ncol] = size(image_t);


%---- 1. Create translation mapping function to correct for misalignment.
% Translation is rounded to the nearest integer
[dx,dy] = mapPoints( image_t, threshold, overlap_thresh );
% dx = round(dx);
% dy = round(dy);
dx = 0;
dy = -1;


%---- 2. Generate combined (donor+acceptor) image after transform
donor_t    = image_t(:,1:ncol/2);
acceptor_t = image_t(:,(ncol/2)+1:end); %untranslated image

% Generate acceptor image that is padded
dedge = 2; %padding size, maximum transform size
assert( abs(dx)<=dedge && abs(dy)<=dedge, ...
        'Gettraces-getPeaks: transform is too large' );

padded = zeros( size(donor_t)+2*dedge );

% Copy acceptor channel image into center, after translation
padded( dedge+dy+(1:nrow), dedge+dx+(1:ncol/2) ) = acceptor_t;

% Copy out center region, which is the new transformed image
% and add it to donor image to create adjusted total intensity image
transAcceptor_t = padded( 1+dedge:end-dedge, 1+dedge:end-dedge );
total_t = donor_t + transAcceptor_t;

% Find intensity peaks above background -- DONOR CHANNEL VERSION
[don_x,don_y] = pickPeaks( donor_t, threshold, overlap_thresh );
nPicked = numel(don_x);

acc_x = don_x + (ncol/2) +dx;
acc_y = don_y            +dy;

% Find intensity peaks above background -- TOTAL INTENSITY VERSION
% [don_x,don_y] = pickPeaks( total_t, threshold, overlap_thresh );
% nPicked = numel(don_x);
% 
% acc_x = don_x + (ncol/2) -dx;
% acc_y = don_y            -dy;

%---- 3. Output resulting peak selection coordinates
picksX = zeros(nPicked*2,1);
picksY = zeros(nPicked*2,1);

for i=1:nPicked,
    picksX(2*i-1) = don_x(i);
    picksY(2*i-1) = don_y(i);
    picksX(2*i)   = acc_x(i);
    picksY(2*i)   = acc_y(i);
end    


% end function getTraces



function wavg = WeightedAvg( vals, weights )
% Produces a weighted average of <vals>.
% Used by getTraces to find fluor centroid position

weights = weights / sum(weights);  % normalize

wavg = sum( vals.*weights );

% end function WEIGHTEDAVG


function [dx,dy, f] = mapPoints( image_t, threshold, overlap_thresh )
% mapPoints  Develop mapping function

dx = 0; dy = 0; %default to no correction

[nrow ncol] = size(image_t);
donor_t = image_t(:,1:ncol/2);
acceptor_t = image_t(:,(ncol/2)+1:end);

[don_x,don_y] = pickPeaks( donor_t, threshold, overlap_thresh );
nPicked = numel(don_x);

%---- 2. Estimate coordinates of acceptor-side peaks

% Initial guess: simple translation to right half (no rotation, scale)
acc_x = don_x + (ncol/2);
acc_y = don_y;

for j=1:nPicked,

    % Refine acceptor peak positions by finding local maxima
    % within the 3x3 grid around the initial guess.
    temp = image_t( acc_y(j)-1:acc_y(j)+1, acc_x(j)-1:acc_x(j)+1 );
    [maxy, maxx] = find(temp==max(max(temp)),1);
    acc_x(j) = acc_x(j) +maxx-2;
    acc_y(j) = acc_y(j) +maxy-2;
    
end


%---- 3. Create mapping function: donor to acceptor side
% ...to account for misalignment of DualView.  The transform is then 
% applied to get predicted Cy5 peak positions.
% 
%   Tranform matrix format:
%   scale X | shear Y | 0
%   shear X | scale Y | 0
%   trans X | trans Y | 1
% 
% NOTE that the image is a grid, so transforms of non-integer amounts
% have stochastic effects going from continuous to discrete cooredinates.
% 
% AFFINE works best
% 

if nPicked<5, 
    % If less than 5 molecules, don't attempt to derive a transform,
    % just use translation
    warning( 'Too few molecules to create transformation' );
%     f = [1 0 0 ; 0 1 0 ; ncol/2 0 1];
    dx = 0; dy = 0; %default to no correction
    return;
else
    donor_points    = [don_x ; don_y]';
    acceptor_points = [acc_x ; acc_y]';

    tform = cp2tform(acceptor_points,donor_points,'linear conformal');  %image proc TK
    f = tform.tdata.Tinv;
    
    u = [0 1]; v=[0 0];
    [tx,ty] = tformfwd(tform, u,v);
    dx = tx(2)-tx(1)
    dy = ty(2)-ty(1)
    angle = (180/pi) * atan2(dy, dx) 
    scale = 1 / sqrt(dx^2 + dy^2)
end

if abs(dx)>1 || abs(dy)>1
    disp('Warning: Dual-View alignment off by more than 1 pixel');
end
    

%---- 4. Use the transform to get acceptor-side coordinates
% for i=1:nPicked,
%     
%     % Apply transform to get acceptor-side coordinates
%     points=[don_x(i) don_y(i) 1]*f;
%     points=round(points);
% 
%     acc_x(i) = points(1);
%     acc_y(i) = points(2);
% end    


% end function mapPoints



function [don_x,don_y] = pickPeaks( image_t, threshold, overlap_thresh )
% Localizes the peaks of fluorescence, removing any too close together.
% NOTE: only the donor side is passed to this function

[nrow ncol] = size(image_t);


% 1. For each pixel (excluding edges), pick those above threshold that have
% greater intensity than their neighbors (3x3,local maxima) -- donor only.
% tempxy is peak position, centroidxy is estimated true molecule position.
nMols=0;
for i=1+3:nrow-3,
    for j=1+3:ncol-3,
        block = image_t(i-1:i+1,j-1:j+1);
        cross = block( [2,4,5,6,8] );
        
        if image_t(i,j)>threshold && image_t(i,j)==max(cross)
            
            % Calc centroid position using intensity-weighted average of
            % position
            off_x = WeightedAvg( 1:3, sum(block,1)  ) -1.5;
            off_y = WeightedAvg( 1:3, sum(block,2)' ) -1.5;
%                   
%             line(j+off_x-0.5,i+off_y-0.5,'marker','o','color','y','EraseMode','background')
            
            % Save the position of this molecule
            nMols=nMols+1;
            tempx(nMols)=j;
            tempy(nMols)=i;
            
            centroidx(nMols) = j+ off_x;
            centroidy(nMols) = i+ off_y;
        end
    end
end

if nMols<1,
    don_x = [];
    don_y = [];
    return;
end


% 2. Remove pick pairs that are closer than the cutoff radius:
% Signals may have cross-contamination. -- donor only
% THIS STEP NEEDS IMPROVEMENT - if molecules are too close, they will not
% be picked as seperate peaks.
overlap=zeros(1,nMols);

for i=1:nMols,
    for j=i+1:nMols,
        % quickly ignore molecules way beyond the threshold
        if centroidy(j)-centroidy(i) > overlap_thresh
            break;
        end
        % otherwise, we have to check more precisely
        if (centroidx(i)-centroidx(j))^2 + (centroidy(i)-centroidy(j))^2 ...
           < overlap_thresh^2
            overlap(i)=1;
            overlap(j)=1;
        end
    end
end


indexes = find( overlap==0 ); %find peaks with overlap less than threshold
don_x = tempx(indexes);
don_y = tempy(indexes);


% END FUNCTION pickPeaks








% --------------------- SAVE PICKED TRACES TO FILE --------------------- %
% --- Executes on button press in saveTraces.
function saveTraces_Callback(hObject, eventdata, handles)
% hObject    handle to saveTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

peaks = [handles.x handles.y]';
assert( size(peaks,1)==2 );

stk_top = handles.stk_top;
stk = getappdata(handles.figure1,'stk');

integrateAndSave( stk, stk_top, peaks, handles.stkfile );

clear stk;

%guidata(hObject,handles);



% function saveTraces( handles, hObject )
function integrateAndSave( stk, stk_top, peaks, stk_fname )
% NOTE: can find which pixels to use by correlation (original code did so?)

[stkX,stkY,Nframes] = size(stk);

if stkX == 128
    squarewidth=3;
    NumPixels=5;  %unoptimized
elseif stkX == 170
    squarewidth=3;
    NumPixels=4;
elseif stkX == 256
    squarewidth=3;
    NumPixels=4;
elseif stkX == 512
    squarewidth=5;
    NumPixels=16;
end

% Get x,y coordinates of picked peaks
Npeaks = size(peaks,2)/2;
x = peaks(1,:);
y = peaks(2,:);

regions=zeros(NumPixels,2,2*Npeaks);  %pixel#, dimension(x,y), peak#

% Define regions over which to integrate each peak --
% Done separately for each channel!
for m=1:2*Npeaks
    
    hw = floor(squarewidth/2);
    
    % Get pixels around picked point (squarewidth^2 pixels)
    peak = stk_top( ...
            y(m)-hw:y(m)+hw, ...
            x(m)-hw:x(m)+hw  );
    center = sort( peak(:) );  %vector of sorted intensities
    
    % Get pixels whose intensity is greater than the median (max=NumPixels).
    % We just want the centroid to avoid adding noise ...
    % A is x-coord, B is Y-coord of the top <NumPixels> pixels
    [A,B]=find( peak>=center(squarewidth*squarewidth-NumPixels+1), NumPixels );
    
    % Define a region over which to integrate each peak
    regions(:,:,m) = [ A+y(m)-hw-1, B+x(m)-hw-1  ];
end


% Create a trace for each molecule across the entire movie
traces = zeros(2*Npeaks,Nframes);

s = size(stk);
idx = sub2ind( s(1:2), regions(:,1,:), regions(:,2,:) );
% bg = handles.background;

for k=1:Nframes,
    frame = double(stk(:,:,k));
    traces(:,k) = sum( frame(idx) );
end

clear stk;


% Subtract background from last point of all traces so data
% can fit into an int16 matrix
donor    = traces(1:2:end,:);
acceptor = traces(2:2:end,:);

donor    = donor    - mean( donor(:,end)    );
acceptor = acceptor - mean( acceptor(:,end) );


% Save data to file.
stk_fname = strrep(stk_fname,'.bz2','');
[p,name]=fileparts(stk_fname);
save_fname = [p filesep name '.traces'];

saveTraces( save_fname, 'traces', donor,acceptor );


% Save the locations of the picked peaks for later lookup.
% FORMAT:  mol_name, don x, don y, acc x, acc y
% filename=strrep(handles.stkfile,'.stk','.loc.txt');
% fid = fopen(filename,'w');
% 
% for j=1:nTrances/2;
%     don_x = x(2*j-1);
%     don_y = y(2*j-1);
%     acc_x = x(2*j);
%     acc_y = y(2*j);
%     
%     fprintf(fid, '%s_%d- %d %d %d %d\n', name, j, don_x,don_y,acc_x,acc_y);
% end
% 
% fclose(fid);


% END FUNCTION saveTraces_Callback



% --------------------- MISC. GUI CALLBACK FUNCTIONS --------------------- %

% --- Executes on slider movement.
function scaleSlider_Callback(hObject, eventdata, handles)
% hObject    handle to scaleSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set( handles.axDonor,    'CLim',[get(hObject,'min') get(hObject,'value')] );
set( handles.axAcceptor, 'CLim',[get(hObject,'min') get(hObject,'value')] );
set( handles.axTotal,    'CLim',[get(hObject,'min')*2 get(hObject,'value')*2] );




