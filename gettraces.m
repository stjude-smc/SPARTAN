function [stkData,peaks,image_t] = gettraces(varargin)
% GETTRACES  Extract smFluorescence traces from movies
%
%    GETTRACES()
%    Launches the gettraces graphical user interface.
%
%    [STK,PEAKS,IMG] = GETTRACES( FILENAME, PARAMS )
%    Loads a movie from FILENAME, finds fluorescence peaks, and returns
%    their locations (PEAKS) as a Nx2 matrix (x,y).  IMG is the image used
%    for selecting intensity peaks.
%
%    [STK,PEAKS,IMG] = gettraces( FILENAME, PARAMS, OUTFILE )
%    As above, but also extracts traces from the PEAKS locations and saves
%    the resulting traces to OUTFILE.
%
%    STK can also be passed instead of FILENAME if the movie has already
%    been loaded by gettraces. This saves on load time.
%
%    GETTRACES( DIRECTORY, PARAMS )
%    For each filename in DIRECTORY, the movie is loaded, processed, and
%    a .traces file is saved automatically. No data is returned.
%    This is "batch mode".  Additional (optional) fields in PARAMS include
%    option to check all child folders (recursive) and option to skip
%    movies that have already been processed (skipExisting).
%
%    PARAMS is a struct with any of the following options:
%     - don_thresh:     total (D+A) intensity threshold for pick selection
%                       (if 0, threshold is automatically selected).
%     - overlap_thresh: peaks are rejected if they are closer than this radius.
%     - nPixelsIntegrated: number of pixels proximal to each peak to sum
%          to produce fluorescence traces. Higher values capture more
%          intensity, but also capture more noise...
%     - saveLocations: write a text file with the locations of each peak to
%          file.
%     - geometry: imaging geometry can be single-channel/full-chip (1), 
%                 dual-channel/left-right (2), or quad-channel (3).
%                 Default: dual-channel (2).
%

% TODO: also return unfiltered peak list (overlap=0) so overlap statistics
% can be displayed in scripts that call this function.





%% Parse input arguments

% If calling directly from command line,
% launch the GUI.
if nargin==0,
    gettraces_gui;
    return;
end


%------ Load parameter values (always second parameter!)
global params;

if nargin<2, %provide defaults if no parameters given.
    params.don_thresh = 0;
    params.overlap_thresh = 2.1;
else
    params = varargin{2};
end

% PARAMTERES FOR SOFTWARE ALIGNMENT!
% The numbers are the displacements of the right (acceptor) field relative to
% the left (donor) field that are required to get the two fields aligned.
params.align_dx =  0;
params.align_dy =  0;


% If any parameters are not specified, give a default value of 0.
% This way, we don't have to constantly check if a parameter exists.
paramNames = {'skipExisting','recursive','don_thresh','overlap_thresh','saveLocations','quiet'};
for i=1:numel(paramNames),
    if ~isfield(params,paramNames{i})
        params.(paramNames{i}) = 0;
    end
end

if ~isfield(params,'geometry'),
    params.geometry=2; %dual-channel (L/R).
end


%------ If a structure is specified, load as STK data
if nargin>=1 && isstruct(varargin{1}),
    stkData = varargin{1};

%------ Otherwise, load an individual file
elseif nargin>=1 && ischar(varargin{1}) && ~isdir(varargin{1}),
    stkData = OpenStk( varargin{1} );

%------ If a directory name is given, run in batch mode and then terminate.
elseif nargin>=1 && ischar(varargin{1}) && isdir(varargin{1})
    batchmode( varargin{1} );
    return;

else
    error('gettraces: Invalid param 1');
end




%------ Find peaks of intensity corrosponding to isolated single molecules

% Generate image that will be used to select peaks
image_t = stkData.stk_top - stkData.background;

if ~params.don_thresh
    if ~isfield(params,'thresh_std')
        constants = cascadeConstants;
        thresh_std = constants.gettracesThresholdStd;
    else
        thresh_std = params.thresh_std;
    end

    % Automatically threshold setting based on variance of
    % background intensity at end of movie.
    endBG = sort( stkData.endBackground(:) );
    endBG_lowerHalf = endBG( 1:floor(numel(endBG)*0.75) );
    params.don_thresh = thresh_std*std( endBG_lowerHalf );
else
    params.don_thresh = params.don_thresh-mean2(stkData.background);
end

% Find peak locations from total intensity
[peaksX,peaksY] = getPeaks( image_t, params );
peaks = [peaksX peaksY];



%------ Extract traces using peak locations
if nargin>=3 && ischar(varargin{3}),
    if ~isfield(params,'nPixelsToSum'),
        params.nPixelsToSum = 4;
    end
    
    outputFilename = varargin{3};
    integrateAndSave( stkData.movie, stkData.stk_top, peaks', outputFilename );
end
    
return;


%% --------------------- OPEN STK MOVIE (CALLBACK) --------------------- %



function [stkData] = OpenStk(filename)

global params;

if nargin<2,
    params.geometry=2;
end

% If the movie is compressed, deflate it to a temporary location first
if strfind(filename,'.bz2'),
    
    % decompress
    z_fname = filename;
    status = system( ['bunzip2 --keep "' z_fname '"'] );
    if status ~= 0, 
        error('Error uncompressing STK: %s',z_fname);
    end
    
    filename = strrep(z_fname,'.bz2','');
    
    % Delete the compressed movie -- we no longer need it!
    delete( z_fname );

elseif strfind(filename,'.movie'),
    error('No simulated movie support (for now)');
%     fid = fopen( filename, 'r' );
%     stkSize = fread( fid, 3, 'int16' );
%     stkY = stkSize(1); stkX = stkSize(2); nFrames = stkSize(3);
% %     time = fread( fid, nFrames, 'float' );
%     stk = fread( fid, stkY*stkX*nFrames, 'int16' );
%     stk = reshape(stk,stkSize');
%     fclose(fid);
%     time = 1:nFrames;
end

if strfind(filename,'.stk') || strfind(filename,'.tif') || strfind(filename,'.lsm'),
    % Load stk movie -- stk(X,Y,frame number)
    movie = Movie_TIFF( filename );
    time = movie.timeAxis;
    stkX = movie.nX;
    stkY = movie.nY;
    nFrames = movie.nFrames;
end


% Create an average image of the first 10 frames (stk_top)
averagesize = min([10 nFrames]);
stk_top = movie.readFrames(1:averagesize);
stk_top = mean(stk_top,3);



% Create an estimated background image by:
% 1. Divide the image into den*den squares
% 2. For each square, find the fifth lowest number
% 3. Rescaling these values back to the original image size
if stkX <= 128
    den=4;  %not optimized
elseif stkX <= 170
    den=5;
elseif stkX <= 256
    den=8;
elseif stkX <= 512
    den=16;
end

background = stk_top;  %**
temp = zeros( floor(stkY/den), floor(stkX/den) );

for i=1:size(temp,1),
    for j=1:size(temp,2),
        sort_temp = background(den*(i-1)+1:den*i,den*(j-1)+1:den*j);
        sort_temp = reshape(sort_temp,1,den*den);  % make into a vector
        sort_temp = sort(sort_temp);

        temp(i,j) = sort_temp(den);  % get the 1/den % smallest value
    end
end

% Rescale the image back to the actual size
background=imresize(temp,[stkY stkX],'bilinear');
% handles.background = mean( stk(:,:,end-4:end), 3 );

% Also create a background image from the last few frames;
% useful for defining a threshold.
endBackground = movie.readFrames(nFrames-11:nFrames-1);

if params.geometry==1,
    % No combining needed for single-color imaging.
elseif params.geometry==2,
    % Combine left and right side of field for two-color imaging.
    endBackground = endBackground(:,1:stkX/2,:) + endBackground(:,(1+stkX/2):end,:);
elseif params.geometry>2,
    % Combine fluorescence from the four quadrants. For 3-color imaging, one
    % should be left out. FIXME
    endBackground = endBackground(1:stkY/2,1:stkX/2,:) + endBackground(1:stkY/2,(1+stkX/2):end,:) + ...
                    endBackground((1+stkY/2):end,1:stkX/2,:) + endBackground((1+stkY/2):end,(1+stkX/2):end,:);
end

% Combine the image stack data with the extras just calculated for later
% processing and return them.
% All of this could be combined into the Movie class! TODO
stkData.movie = movie;
stkData.stk_top = stk_top;
stkData.background = background;
stkData.endBackground = double(endBackground);
stkData.time = time;
stkData.stkY=stkY;
stkData.stkX=stkX;
stkData.nFrames=nFrames;

% END FUNCTION OpenStk




% --------------- OPEN ALL STK'S IN DIRECTORY (BATCH) --------------- %

function batchmode(direct)

global params;

% Get list of files in current directory (option: and all subdirectories)
if params.recursive
    movieFilenames  = rdir([direct filesep '**' filesep '*.stk']);
    movieFilenames  = [movieFilenames; rdir([direct filesep '**' filesep '*.stk.bz2'])];
    movieFilenames  = [movieFilenames; rdir([direct filesep '**' filesep '*.movie'])];
else
    movieFilenames  = rdir([direct filesep '*.stk']);
    movieFilenames  = [movieFilenames; rdir([direct filesep '*.stk.bz2'])];
    movieFilenames  = [movieFilenames; rdir([direct filesep '*.movie'])];
end

nFiles = length(movieFilenames);

% Wait for 100ms to give sufficient time for polling file sizes in the
% main loop below.
pause(0.1);


% Process each file in the user selected directory.
nTraces  = zeros(nFiles,1); % number of peaks found in file X
existing = zeros(nFiles,1); % true if file was skipped (traces file exists)

for i=1:nFiles,
    
    % Skip if previously processed (.traces file exists)
    stk_fname = strrep(movieFilenames(i).name,'.bz2','');
    [p,name]=fileparts(stk_fname);
    traceFname = [p filesep name '.traces'];
    
    if params.skipExisting && exist(traceFname,'file'),
        if ~params.quiet,
            disp( ['Skipping (already processed): ' stk_fname] );
        end
        existing(i) = 1;
        continue;
    end
    
    % Poll the file size to make sure it isn't changing.^M
    % This could happen when a file is being saved during acquisition.^M
    d = dir(movieFilenames(i).name);
    if movieFilenames(i).bytes ~= d(1).bytes,
        disp( ['Skipping (save in process?): ' movieFilenames(i).name] );
        existing(i) = 1;
        continue;
    end
    
    % Show waitbar only when new data must be loaded.
    if ~exist('h','var'),
        h = waitbar( (i-1)/nFiles,'Extracting traces from movies...');
    end
    
    % Load STK file
    try
        stkData = OpenStk( movieFilenames(i).name );
    catch e
        disp(e);
        disp('Skipping file: corrupted, missing, or not completely saved.');
        existing(i) = 1;
        continue;
    end
    
    % Pick molecules using default parameter values
    image_t = stkData.stk_top - stkData.background;

    if ~params.don_thresh
        if ~isfield(params,'thresh_std')
            constants = cascadeConstants;
            thresh_std = constants.gettracesThresholdStd;
        else
            thresh_std = params.thresh_std;
        end
        
        % Automatically threshold setting based on variance of
        % background intensity at end of movie.
        endBG = sort( stkData.endBackground(:) );
        endBG_lowerHalf = endBG( 1:floor(numel(endBG)*0.75) );
        don_thresh = thresh_std*std( endBG_lowerHalf );
    else
        don_thresh = params.don_thresh-mean2(stkData.background);
    end

    % Find peak locations from total intensity
    params2 = params;
    params2.don_thresh = don_thresh;
    
    [peaksX,peaksY] = getPeaks( image_t, params2 );
    peaks = [peaksX peaksY];
    nTraces(i) = numel(peaksX)/2;
    
    % Save the traces to file
    integrateAndSave( stkData.movie, stkData.stk_top, peaks', traceFname );
    
    waitbar(i/nFiles, h); drawnow;
end

if exist('h','var'),  close(h);  end

% If no new data was loaded, nothing more to do.
if all(existing),  return;  end


% ----- Create log file with results
log_fid = fopen( [direct filesep 'gettraces.log'], 'w' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');

names = fieldnames(  params );
vals  = struct2cell( params );

for i=1:numel(names),
    fprintf(log_fid, '  %15s:  %.2f\n', names{i}, vals{i});
end

% Log list of files processed by gettraces
fprintf(log_fid,'\n%s\n\n%s\n%s\n\n%s\n',date,'DIRECTORY',direct,'FILES');

for i=1:nFiles
    if existing(i),
        fprintf(log_fid, 'SKIP %s\n', movieFilenames(i).name);
    else
        fprintf(log_fid, '%.0f %s\n', nTraces(i), movieFilenames(i).name);
    end
end


% Clean up
fclose(log_fid);







% --------------- PICK MOLECULES CALLBACKS --------------- %

%------------- Pick Cy3 spots ----------------- 
function [picksX,picksY,total_t] = getPeaks( image_t, params )
% Localizes the peaks of molecules in the Cy3 channel and infers the
% positions of the Cy5 peaks by applying a transofmration to the Cy3 peak
% positions.
%
% picksX - X-coords of all molcules, as Cy3,Cy5,Cy3,Cy5 in order
% picksY - Y-coords ...


% If using a single-channel setup, use a simplified proceedure.
if params.geometry==1,
    total_t = image_t;
    [picksX,picksY] = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
    picksX = reshape( picksX, numel(picksX),1 );
    picksY = reshape( picksY, numel(picksY),1 );
    return;
end
% Otherwise, assume dual-channel. FIXME


[nrow ncol] = size(image_t);
ncol = ncol/2;

%---- 1. Transform acceptor side so that the two channels align properly.
donor_t    = image_t( :, 1:ncol );
acceptor_t = image_t( :, (ncol+1):end );

dx = round(params.align_dx);
dy = round(params.align_dy);

total_t = zeros( size(donor_t) );
total_t( max(1,1+dy):min(nrow,nrow+dy), max(1,1+dx):min(ncol,ncol+dx) ) = ...
    acceptor_t( max(1,1-dy):min(nrow,nrow-dy), max(1,1-dx):min(ncol,ncol-dx) );

total_t = total_t + donor_t; %sum the two fluorescence channels.
%FIXME: regions w/o acceptor intensity should be erased!


%---- 2. Get coordinates of intensity peaks using the summed image.
[don_x,don_y] = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
nPicked = numel(don_x);

acc_x = don_x + ncol +dx;
acc_y = don_y        +dy;


%---- 3. Save results for output
picksX = zeros(nPicked*2,1);
picksY = zeros(nPicked*2,1);

for i=1:nPicked,
    picksX(2*i-1) = don_x(i);
    picksY(2*i-1) = don_y(i);
    picksX(2*i)   = acc_x(i);
    picksY(2*i)   = acc_y(i);
end    


%---- 4. Re-estimate coordinates of acceptor-side peaks to verify alignment.
align_acc_x = acc_x;
align_acc_y = acc_y;

for j=1:nPicked,
    % Refine acceptor peak positions by finding local maxima
    % within the 3x3 grid around the initial guess.
    temp = image_t( acc_y(j)-1:acc_y(j)+1, acc_x(j)-1:acc_x(j)+1 );
    [maxy, maxx] = find(temp==max(max(temp)),1);
    align_acc_x(j) = align_acc_x(j) +maxx-2;
    align_acc_y(j) = align_acc_y(j) +maxy-2;
end

% Verify the alignment
x_align = mean( align_acc_y-acc_y );
y_align = mean( align_acc_y-acc_y );
if abs(x_align)+abs(y_align)>0.5,
    warning('gettraces:badAlignment','Fluorescence fields may be out of alignment.');
end


% end function getTraces



function wavg = WeightedAvg( vals, weights )
% Produces a weighted average of <vals>.
% Used by getTraces to find fluor centroid position

weights = weights / sum(weights);  % normalize

wavg = sum( vals.*weights );

% end function WEIGHTEDAVG



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
        if abs(centroidy(j)-centroidy(i)) > overlap_thresh
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


% function saveTraces( handles, hObject )
function integrateAndSave( movie, stk_top, peaks, stk_fname )
% NOTE: can find which pixels to use by correlation

global params;


nFrames = movie.nFrames;
time = movie.timeAxis;
wbh = waitbar(0,'Extracting traces from movie data');

% Specify the number of most intense proximal pixels to sum when generating
% fluorescence traces (depends on experimental point-spread function).
assert( isfield(params,'nPixelsToSum'), 'Missing nPixelsToSum' );

swChoices = (1:2:19);
idx = find( params.nPixelsToSum<=(swChoices.^2), 1,'first' );
squarewidth = swChoices(idx);

% Get x,y coordinates of picked peaks
Npeaks = size(peaks,2);
x = peaks(1,:);
y = peaks(2,:);

regions=zeros(params.nPixelsToSum,2,Npeaks);  %pixel#, dimension(x,y), peak#

% Define regions over which to integrate each peak --
% Done separately for each channel!
for m=1:Npeaks
    
    hw = floor(squarewidth/2);
    
    % Get pixels around picked point (squarewidth^2 pixels)
    peak = stk_top( ...
            y(m)-hw:y(m)+hw, ...
            x(m)-hw:x(m)+hw  );
    center = sort( peak(:) );  %vector of sorted intensities
    
    % Get pixels whose intensity is greater than the median (max=NumPixels).
    % We just want the centroid to avoid adding noise ...
    % A is x-coord, B is Y-coord of the top <NumPixels> pixels
    [A,B]=find( peak>=center(squarewidth*squarewidth-params.nPixelsToSum+1), ...
                params.nPixelsToSum );
    
    % Define a region over which to integrate each peak
    regions(:,:,m) = [ A+y(m)-hw-1, B+x(m)-hw-1  ];
end


% Create a trace for each molecule across the entire movie
traces = zeros(Npeaks,nFrames);

idx = sub2ind( [movie.nY movie.nX], regions(:,1,:), regions(:,2,:) );

for k=1:nFrames,
    frame = double( movie.readFrame(k) );
    if params.nPixelsToSum>1
        traces(:,k) = sum( frame(idx) );
    else
        traces(:,k) = diag( frame(y,x) );
    end
    
    if mod(k,100)==0,
        waitbar( k/nFrames, wbh );
    end
end


if params.geometry==1, %single-channel
    donor    = traces;
    acceptor = zeros( size(traces) );
elseif params.geometry==2, %dual-channel
    donor    = traces(1:2:end,:);
    acceptor = traces(2:2:end,:);
else  %quad-channel
    %TODO
end


% Subtract background from last point of all traces so data
% can fit into an int16 matrix
donor    = donor    - mean( donor(:,end)    );
acceptor = acceptor - mean( acceptor(:,end) );


% Save data to file.
stk_fname = strrep(stk_fname,'.bz2','');
[p,name]=fileparts(stk_fname);
save_fname = [p filesep name '.traces'];

saveTraces( save_fname, 'traces', donor,acceptor, [], time );


% Save the locations of the picked peaks for later lookup.
% FORMAT:  mol_name, don x, don y, acc x, acc y
if params.saveLocations,
    filename=strrep(stk_fname,'.stk','.loc.txt');
    fid = fopen(filename,'w');

    for j=1:Npeaks/2,
        don_x = x(2*j-1);
        don_y = y(2*j-1);
        acc_x = x(2*j);
        acc_y = y(2*j);

        fprintf(fid, '%s_%d %d %d %d %d\n', name, j, don_x,don_y,acc_x,acc_y);
    end

    fclose(fid);
end

close( wbh );




