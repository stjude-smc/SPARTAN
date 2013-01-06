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
%    a .rawtraces file is saved automatically. No data is returned.
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
%     - crosstalk: donor-to-acceptor channel fluorescence channel
%                  bleedthrough (crosstalk) as a fraction. For correction.
%     - photonConversion: fluorescence AU/photon conversion factor.
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
    params.overlap_thresh = 2.3;
    params.nPixelsToSum = 4;
else
    params = varargin{2};
end


% If any parameters are not specified, give a default value of 0.
% This way, we don't have to constantly check if a parameter exists.
paramNames = {'skipExisting','recursive','don_thresh','overlap_thresh','saveLocations','quiet','alignTranlate','refineAlign'};
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

% If all we want is the movie data, don't pick yet.
if nargout==1,
    return;
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
peaks = getPeaks( image_t, params );

% Generate integration windows for later extracting traces.
[stkData.regions,stkData.integrationEfficiency] = ...
                            getIntegrationWindows(image_t, peaks', params);


%------ Extract traces using peak locations
if nargin>=3 && ischar(varargin{3}),
    if ~isfield(params,'nPixelsToSum'),
        params.nPixelsToSum = 4;
    end
    
    outputFilename = varargin{3};
    integrateAndSave( stkData, peaks', outputFilename );
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

params.bgBlurSize = den;


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
    
    % Skip if previously processed (.rawtraces file exists)
    stk_fname = strrep(movieFilenames(i).name,'.bz2','');
    [p,name]=fileparts(stk_fname);
    traceFname = [p filesep name '.rawtraces'];
    
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
    
    peaks = getPeaks( image_t, params2 );
    nTraces(i) = size(peaks,1);
    
    % Generate integration windows for later extracting traces.
    [stkData.regions,stkData.integrationEfficiency] = ...
                                getIntegrationWindows(image_t, peaks', params);
    
    % Save the traces to file
    integrateAndSave( stkData, peaks', traceFname );
    
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
        fprintf(log_fid, '%.0f %s\n', nTraces(i)/2, movieFilenames(i).name);
    end
end


% Clean up
fclose(log_fid);







% --------------- PICK MOLECULES CALLBACKS --------------- %

%------------- Pick single molecule spots ----------------- 
function [picks,total_t, align] = getPeaks( image_t, params )
% Localizes the peaks of molecules from a summed image of the two channels
% (in FRET experiments).
%
% picksX - X-coords of all molcules, as Cy3,Cy5,Cy3,Cy5 in order
% picksY - Y-coords ...


% If using a single-channel setup, use a simplified proceedure
% since we don't have to worry about alignment, etc.
if params.geometry==1,
    total_t = image_t;
    picks = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
    assert( size(picks,2)==2 ); %x and y columns, molecules in rows.
    return;
end
% Otherwise, assume dual-channel. FIXME


[nrow ncol] = size(image_t);
ncol = ncol/2;


%---- 1. Pick molecules as peaks of intensity from summed (D+A) image)
donor_t    = image_t( :, 1:ncol );
acceptor_t = image_t( :, (ncol+1):end );

total_t = acceptor_t + donor_t; %sum the two fluorescence channels.

donor_picks = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
nPicked = size(donor_picks,1);

% Define acceptor size peaks as a simple translation across the chip.
picks = zeros(nPicked*2,2); %donor, acceptor alternating; 2 columns = x,y

picks(1:2:end,:) = donor_picks; %donor
picks(2:2:end,1) = donor_picks(:,1) + ncol; %acceptor
picks(2:2:end,2) = donor_picks(:,2);        %acceptor

% Refine peak locations to determine if realignment is needed.
refinedPicks = refinePeaks( image_t, picks );

x_align = mean( refinedPicks(2:2:end,1)-refinedPicks(1:2:end,1)-ncol );
y_align = mean( refinedPicks(2:2:end,2)-refinedPicks(1:2:end,2)      );
x_align_abs = mean(abs( refinedPicks(2:2:end,1)-refinedPicks(1:2:end,2)-ncol )); %this will detect rotation as well
y_align_abs = mean(abs( refinedPicks(2:2:end,2)-refinedPicks(1:2:end,2)      ));
align = [x_align y_align x_align_abs y_align_abs];

% If specified, use the refined peak positions to handle slight misalignment.
if params.refineAlign,
    picks = refinedPicks;
end

% If the alignment is close (by translation), we're done.
if abs(x_align)<0.5 && abs(y_align)<0.5,
    return;
end


%---- 2. Transform acceptor side so that the two channels align properly.
% If the alignment is off by a significant margin, the fields are
% realigned in software. This will only handle translation. Rotation and
% other complex distortions are harder. FIXME.
warning('gettraces:badAlignment','Fluorescence fields may be out of alignment.');
fprintf( 'X Translation (Absolute): %.1f  (%.1f)\n', x_align, x_align_abs );
fprintf( 'Y Translation (Absolute): %.1f  (%.1f)\n', y_align, y_align_abs );

% Just give a warning unless asked to do software alignment in settings.
if ~params.alignTranslate,  return;  end 


dx = round( x_align );
dy = round( y_align );

% Sum the donor and acceptor fields, after translating the acceptor size to
% deal with the misalignment.
align_t = zeros( size(donor_t) );
align_t( max(1,1+dy):min(nrow,nrow+dy), max(1,1+dx):min(ncol,ncol+dx) ) = ...
    acceptor_t( max(1,1-dy):min(nrow,nrow-dy), max(1,1-dx):min(ncol,ncol-dx) );

total_t_old = total_t;
total_t = align_t + donor_t; %sum the two fluorescence channels.

% Pick peaks from the aligned image.
[don_x,don_y] = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
nPicked = numel(don_x);

% Since the alignment image has been shifted to compensate for misalignment
% already (above), adjust the output coordinates so they are relative to
% the actual fields, not the adjusted fields.
picks = zeros(nPicked*2,2); %donor, acceptor alternating; 2 columns = x,y

picks(1:2:end,1) = don_x;            %donor
picks(2:2:end,1) = don_x + ncol -dx; %acceptor
picks(1:2:end,2) = don_y;     %donor
picks(2:2:end,2) = don_y -dy; %acceptor


%---- 4. Re-estimate coordinates of acceptor-side peaks to verify alignment.
refinedPicks = refinePeaks( image_t, picks );

% If specified, use the refined peak positions to handle slight misalignment.
if params.refineAlign,
    picks = refinedPicks;
end

% Verify the alignment
x_align = mean( refinedPicks(2:2:end,1)-refinedPicks(1:2:end,1)-ncol );
y_align = mean( refinedPicks(2:2:end,2)-refinedPicks(1:2:end,2)      );
x_align_abs = mean(abs( refinedPicks(2:2:end,1)-refinedPicks(1:2:end,2)-ncol )); %this will detect rotation as well
y_align_abs = mean(abs( refinedPicks(2:2:end,2)-refinedPicks(1:2:end,2)      ));
align = [x_align y_align x_align_abs y_align_abs];

if x_align_abs>0.5 || y_align_abs>0.5,
    warning('gettraces:badAlignment','Fluorescence fields are STILL out of alignment. Rotation is off?');
    fprintf( 'X Translation (Absolute): %.1f  (%.1f)\n', x_align, x_align_abs );
    fprintf( 'Y Translation (Absolute): %.1f  (%.1f)\n', y_align, y_align_abs );
end


% end function getTraces



function wavg = WeightedAvg( vals, weights )
% Produces a weighted average of <vals>.
% Used by getTraces to find fluor centroid position

weights = weights / sum(weights);  % normalize

wavg = sum( vals.*weights );

% end function WEIGHTEDAVG



function [clean_picks,all_picks] = pickPeaks( image_t, threshold, overlap_thresh )
% Localizes the peaks of fluorescence, removing any too close together.

[nrow ncol] = size(image_t);


% 1. For each pixel (excluding edges), pick those above threshold that have
% greater intensity than their neighbors (3x3,local maxima).
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
    donor_picks = zeros(0,2);
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

clean_picks = [don_x' don_y'];
all_picks = [tempx' tempy'];


% END FUNCTION pickPeaks


function peaks = refinePeaks( image_t, peaks )
% pickPeaks simply finds peaks of intensity in the total (D+A) image. Here
% the peak locations are refined to account for slight differences due to
% misalignment, where the donor and acceptor peaks may be in different
% relative positions. This can be used to re-align the images in software.
% The input image and peak locations are listed as:
%     donor/acceptor/donor/acceptor/etc.

for j=1:size(peaks,1),
    % Refine acceptor peak positions by finding local maxima
    % within the 3x3 grid around the initial guess.
    temp = image_t( peaks(j,2)-1:peaks(j,2)+1, peaks(j,1)-1:peaks(j,1)+1 );
    [maxy, maxx] = find(temp==max(temp(:)),1);
    peaks(j,1) = peaks(j,1) +maxx-2;  %X
    peaks(j,2) = peaks(j,2) +maxy-2;  %Y
end

% END FUNCTION refinePeaks



% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

function [regions,integrationEfficiency] = getIntegrationWindows( stk_top, peaks, params )

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
integrationEfficiency = zeros(Npeaks,squarewidth^2);

for m=1:Npeaks
    
    hw = floor(squarewidth/2);
    
    % Get pixels around picked point (squarewidth^2 pixels)
    peak = stk_top( ...
            y(m)-hw:y(m)+hw, ...
            x(m)-hw:x(m)+hw  );
    center = sort( peak(:) );  %vector of sorted intensities in peak.
    
    % Estimate the fraction of intensity in each pixel.
    % This is used in the GUI to show the efficiency of collecting
    % intensity at a given integratin window size and to estimate the size
    % of the point-spread function.
    integrationEfficiency(m,:) = integrationEfficiency(m,:) + cumsum( center(end:-1:1)/sum(center) )';
    
    % Get pixels whose intensity is greater than the median (max=NumPixels).
    % We just want the centroid to avoid adding noise ...
    % A is x-coord, B is Y-coord of the top <NumPixels> pixels
    [A,B]=find( peak>=center(squarewidth*squarewidth-params.nPixelsToSum+1), ...
                params.nPixelsToSum );
    
    % Define a region over which to integrate each peak
    regions(:,:,m) = [ A+y(m)-hw-1, B+x(m)-hw-1  ];
end


% end function getIntegrationWindows




% function saveTraces( handles, hObject )
function integrateAndSave( stkData, peaks, stk_fname )
% NOTE: can find which pixels to use by correlation

global params;

movie = stkData.movie;
nFrames = movie.nFrames;
data.time = movie.timeAxis;
wbh = waitbar(0,'Extracting traces from movie data');

% Get x,y coordinates of picked peaks
Npeaks = size(peaks,2);
x = peaks(1,:);
y = peaks(2,:);

regions = stkData.regions;  % pixel#, dimension(x,y), peak#


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


% Convert fluorescence to arbitrary units to photon counts.
constants = cascadeConstants;
if isfield(params,'photonConversion'),
    donor    = donor./params.photonConversion;
    acceptor = acceptor./params.photonConversion;
end

% Make an adjustment for crosstalk on the camera.
if ~isfield(params,'crosstalk'),
    params.crosstalk = constants.crosstalk;
end
acceptor = acceptor - params.crosstalk*donor;

% Subtract background and calculate FRET
[data.donor,data.acceptor,data.fret] = correctTraces(donor,acceptor,constants);


% Correct for variable sensitivity across acceptor-channel camera,
% according to measurements with DNA oligos -- QZ
% This is very specific to our equipment, so it should realistically be put
% somewhere outside (cascadeConstants) as an option.
% if params.acceptorFieldRescale,
    % creat a Look up table for intensity correction
    for j=1:Npeaks/2,
        acc_y = y(2*j);
        yCorrection = 0.87854+acc_y*9.45332*10^(-4);
        acceptor(j,:) = acceptor(j,:)/yCorrection;
    end
% end


% ---- Metadata: save various metadata parameters from movie here.

% -- Fields specific to this movie.
% data.movieMetadata.filename  = strrep(stk_fname,'.bz2','');
% data.movieMetadata.crosstalk = params.crosstalk;

% -- Fields specific to each trace:
% Save the locations of the picked peaks for later lookup.
% Also save indexes to map traces to movie metadata (here, everything is 1
% because there is only one movie).
nTraces = size(data.donor,1);

if params.geometry==1
    data.traceMetadata = struct( 'donor_x',num2cell(x), 'donor_y',num2cell(y) );
elseif params.geometry==2,
    data.traceMetadata = struct( ...
        'donor_x',    num2cell( x(1:2:end) ), 'donor_y',    num2cell( y(1:2:end) ), ...
        'acceptor_x', num2cell( x(2:2:end) ), 'acceptor_y', num2cell( y(2:2:end) ) ...
    );
elseif params.geometry>2,
    % TODO
end
% z = num2cell( ones(1,nTraces) );
% [data.traceMetadata.movieIndex] = deal( z{:} );

% -- Create trace identifiers.
% This is just the full path to the movie plus a trace number. This can be
% used to later find the corresponding original movie data for each
% individual trace, even after many rounds of processing.
stk_fname = strrep(stk_fname,'.bz2','');

for i=1:nTraces,
    data.traceMetadata(i).ids = sprintf( '%s#%d', stk_fname, i );
end

% ---- Save data to file.
[p,name]=fileparts(stk_fname);
save_fname = [p filesep name '.rawtraces'];

saveTraces( save_fname, 'traces', data );


close( wbh );

% end function integrateAndSave


