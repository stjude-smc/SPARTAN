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
%                 dual-channel/left-right (2), or quad-channel (3/4).
%                 Default: dual-channel (2).
%     - crosstalk: donor-to-acceptor channel fluorescence channel
%                  bleedthrough (crosstalk) as a fraction. For correction.
%     - photonConversion: fluorescence AU/photon conversion factor.
%     - fieldNames: cell array of assignments of fluorescence field names
%                     (e.g., donor, acceptor, factor). If one is left blank
%                     that field is disregarded. OPTIONAL.
%


constants = cascadeConstants;


% If calling directly from command line, launch the GUI.
if nargin==0,
    gettraces_gui;
    return;
end


%------ Load parameter values (always second parameter!)
global params;


% Set default parameter values and merge with user-specified values.
% The user's options will override any existing defaults.
params = constants.gettracesDefaultParams;

if nargin>=2,
    userParams = varargin{2};
    % FIXME: verify all required parameter settings are given!    
    params = catstruct( params, userParams );
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
[peaks,stkData.total_t,stkData.alignStatus,stkData.total_peaks,stkData.fractionOverlapped, ...
    stkData.rejectedPicks,stkData.rejectedTotalPicks] = getPeaks( image_t, params );

% Generate integration windows for later extracting traces.
[stkData.regions,stkData.integrationEfficiency] = ...
                            getIntegrationWindows(image_t, peaks, params);


%------ Extract traces using peak locations
if nargin>=3 && ischar(varargin{3}),
    outputFilename = varargin{3};
    integrateAndSave( stkData, peaks, outputFilename );
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
end

[p,f,ext] = fileparts(filename);

if ~isempty( strfind(ext,'.stk') ),
    % Load stk movie -- stk(X,Y,frame number)
    movie = Movie_STK( filename );
    time = movie.timeAxis;
    stkX = movie.nX;
    stkY = movie.nY;
    nFrames = movie.nFrames;
end

if ~isempty( strfind(ext,'.tif') ),
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
% 
% *** den should be much larger than the PSF size. A good rule of thumb is
% at least 3x the size of the width (not std) of the PSF. This is difficult
% to generalize just from the size of the frame.
den = 6;
params.bgBlurSize = den;


background = stk_top;  %**
temp = zeros( floor(stkY/den), floor(stkX/den) );

for i=1:size(temp,1),
    for j=1:size(temp,2),
        sort_temp = background(den*(i-1)+1:den*i,den*(j-1)+1:den*j);
        sort_temp = sort( sort_temp(:) );

        temp(i,j) = sort_temp( den );  % get the 1/den % smallest value
    end
end

% Rescale the image back to the actual size
background=imresize(temp,[stkY stkX],'bicubic');

% Also create a background image from the last few frames;
% useful for defining a threshold.
% FIXME: this causes problems in experiments where the background at the
% end of the movies is higher than at the beginning: tRNA selection.
endBackground = movie.readFrames(nFrames-11:nFrames-1);

if params.geometry==1,
    % No combining needed for single-color imaging.
    stkData.nChannels = 1;
elseif params.geometry==2,
    % Combine left and right side of field for two-color imaging.
    endBackground = endBackground(:,1:stkX/2,:) + endBackground(:,(1+stkX/2):end,:);
    stkData.nChannels = 2;
elseif params.geometry>2,
    % Combine fluorescence from the four quadrants. For 3-color imaging, one
    % should be left out. FIXME
    endBackground = endBackground(1:stkY/2,1:stkX/2,:) + endBackground(1:stkY/2,(1+stkX/2):end,:) + ...
                    endBackground((1+stkY/2):end,1:stkX/2,:) + endBackground((1+stkY/2):end,(1+stkX/2):end,:);
    stkData.nChannels = 4;
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
        endBG = sort( stkData.background(:) );
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
                                getIntegrationWindows(image_t, peaks, params);
    
    % Save the traces to file
    integrateAndSave( stkData, peaks, traceFname );
    
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
    if iscell( vals{i} )
        f = repmat( '%s, ', 1,numel(vals{i}));
        f = f(1:end-2);
        fprintf(log_fid, ['  %15s:  ' f '\n'], names{i}, vals{i}{:});
    else
        fprintf(log_fid, '  %15s:  %.2f\n', names{i}, vals{i});
    end
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
function [picks,total_t,align,total_picks,fractionOverlapped,rejectedPicks,rejectedTotalPicks] ...
             = getPeaks( image_t, params )
% Localizes the peaks of molecules from a summed image of the two channels
% (in FRET experiments). The selection is made on the total fluorescence
% intensity image (summing all channels into a single image) to minimize
% bias (e.g., molecules with a dim donor and bright acceptor would be
% missed if we selected just on the donor side). The error in alignment is
% estimated by looking in the immediate neighborhood of each spot for the
% actual intensity maximum, which could be different on each side if the
% two fields are not aligned or if there are optical distortions.
%
% (If params.alignTranslate and alignRotate are set:) If the fields are not
% closely aligned, alignSearch() will try many possible alignments over a
% range of values, returning the best one. This process is pretty slow and
% the rotation adds noise to peak locations (picks), so rotation is turned
% off by default.
%
%    picks   = peak center positions for each channel listed in order:
%            For dual-color (FRET), it is Cy3,Cy5,Cy3,Cy5,...
%            For quad-color, UL,LL,LR,UR is the order.
%    total_t = total intensity image that was used for selecting peaks.
%    align   = alignment information as struct(dx,dy,theta,mag,abs_dev).
%            If software alignment is used, these numbers correspond
%            to the alignment applied.
%    total_picks = peak locations in the total intensity "channel" used for
%            picking in the first place.
%    fractionOverlapped = fraction of traces removed because they were too
%            close to a neighboring peak.
%    rejectedPicks = locations of peaks that will not be considered because
%            they are too close to a neighbor.
%    rejectedTotalPicks = rejected peaks in total intensity image.
% 

% nChannels = params.geometry;
align = [];
abs_dev = 0;


% Define each channel's dimensions and sum fields together.
% For dual-channel, we assume they are arranged L/R.
% For quad-channel, it is UL/LL/LR/UR
[nrow ncol] = size(image_t);

% Indexes into the list of all possible fields.
indD = find( strcmp(params.chNames,'donor') );
indA = find( strcmp(params.chNames,'acceptor') );

if params.geometry==1,
    total_t = image_t;

elseif params.geometry==2,
    left  = image_t( :, 1:ncol/2 );
    right = image_t( :, (ncol/2+1):end );
    allFields = {left,right};
    
    donor_t = allFields{indD};
    acceptor_t = allFields{indA};
    
    total_t = acceptor_t + donor_t; %sum the two fluorescence channels.

elseif params.geometry>2,
    upperLeft  = image_t( 1:nrow/2, 1:ncol/2 );
    lowerLeft  = image_t( (nrow/2)+1:end, 1:ncol/2 );
    lowerRight = image_t( (nrow/2)+1:end, (ncol/2)+1:end );
    upperRight = image_t( 1:nrow/2, (ncol/2)+1:end );
    
    % Extract donor and acceptor peaks.
    allFields = {upperLeft, lowerLeft, lowerRight, upperRight};
    donor_t = allFields{indD};
    acceptor_t = allFields{indA};
    
    total_t = upperLeft + lowerLeft + lowerRight + upperRight;
end


%---- 1. Pick molecules as peaks of intensity from summed (D+A) image)
[total_picks,rejected] = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
nPicked = size(total_picks,1);

assert( all(total_picks(:))>0, 'bad peak locations' );


% Predict the locations of peaks in all channels.
quadrants = find( ~cellfun(@isempty,params.chNames) );  %graduant number for each channel.
channelNames = params.chNames(quadrants);
nChannels = numel(channelNames);  %# of channels TO USE.
nCh = nChannels;

if params.geometry==1,
    picks = total_picks;
    %nChannels = 1; nCh=1;

elseif params.geometry>1,
    % For each channel, predict peak locations as straight translations 
    % from the peak locations in the total intensity image.
    picks = zeros(nCh*nPicked, 2);
    
    for i=1:nChannels,  
        picks(i:nCh:end,:) = translatePeaks( total_picks, size(total_t), quadrants(i) );
    end
end



%%%%% Estimation of misalignment (for dual or multi-color)
% FIXME: this only works for donor/acceptor (FRET) channels. The factor
% channel may not have any intensity to do the alignment calculation...
% There is also no crosstalk that could be used to measure it...
if params.geometry>1,
    
    % Get indexes of the donor and acceptor channels into the final picks
    % array, which only includes the channels IN USE -- channels without a
    % name are ignored.
    indD = find( strcmp(channelNames,'donor') );
    indA = find( strcmp(channelNames,'acceptor') );
    
    % Refine peak locations to determine if realignment is needed.
    refinedPicks = refinePeaks( image_t, picks );
    
    % Calculate mean deviations from donor->acceptor fields.
    if params.geometry==2,
        r_mod = [ mod(refinedPicks(:,1),ncol/2) refinedPicks(:,2) ]; %2-color
    else
        r_mod = [ mod(refinedPicks(:,1),ncol/2) mod(refinedPicks(:,2),nrow/2) ]; %4-color
    end
    x_diff = r_mod(indA:nCh:end,1)-r_mod(indD:nCh:end,1);
    y_diff = r_mod(indA:nCh:end,2)-r_mod(indD:nCh:end,2);
    
    % Use the picked peak locations to create a simple transformation
    % (including translation, rotation, and scaling) from donor to acceptor.
    tform = cp2tform( r_mod(indD:nCh:end,:), r_mod(indA:nCh:end,:), ...
                      'nonreflective similarity');
    T = tform.tdata.T;
    dx = T(3,1);  dy = T(3,2);  sx = T(1,1);  sy = T(2,2);
    theta = asin(T(1,2)) *180/pi;
    mag = mean([sx sy]);
    
    abs_dev = mean( sqrt( x_diff.^2 + y_diff.^2 )  ); % mean displacement.
    align = struct('dx',dx,'dy',dy,'theta',theta,'mag',mag,'abs_dev',abs_dev);
    
end


%%%%% Optional software alignment algorithm (for dual-color only!)
% If the alignment is close, no need to adjust.
% Just give a warning unless asked to do software alignment in settings.
if params.geometry>1 && abs_dev>0.5 && (params.alignTranslate || params.alignRotate),
    % 
    %---- 2. Transform acceptor side so that the two channels align properly.
    
    % Try out all possible alignments within a range and find the one with
    % the best donor-acceptor intensity overlap. The quality score is the
    % mean aligned peak magnitude vs random alignment. If the score is low,
    % reject it and just say "we don't know".
    [align_reg,total_reg,quality] = alignSearch( [donor_t acceptor_t], params );
    
    if ~( align_reg.dx==0 && align_reg.dy==0 && align_reg.theta==0 ),
        % The alignment is clearly off, but we may not have high confidence
        % in the predicted alignment. Just give a warning if it is below
        % some arbitrary, but carefully chosen, threshold.
        if quality<1.12,
            warning('Low confidence alignment. Parameters are out of range or data quality is poor.');
            disp(quality);
        end

        align = align_reg;
        total_t = total_reg;
        align.quality=quality;

        % Pick peaks from the aligned, total intensity image.
        [total_picks,rejected] = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );

        % Transform peak location in the total intensity image to where
        % they fall in the original images.
        T = [ sx*cos(align.theta*pi/180)     sin(align.theta*pi/180)  0 ; ...
                -sin(align.theta*pi/180)  sy*cos(align.theta*pi/180)  0 ; ...
                     align.dx                      align.dy           1 ];
        tform_a = maketform('affine',T);
        
        % Predict peak locations of the other channels assuming simple
        % translation -- assuming fields are aligned.
        [acceptor_picks,remove] = translatePeaks( total_picks, size(total_t), ...
                                               quadrants(indA), tform_a );
                                           
        % Remove peaks that were moved outside the field boundries in at
        % least one of the fields due to the software alignment.
        total_picks  = total_picks(~remove,:);
        rejected     = rejected(~remove);
        nPicked = size(total_picks,1);
        
        % For each channel, predict peak locations as straight translations 
        % from the peak locations in the total intensity image.
        picks = zeros( nCh*nPicked, 2 );

        for i=1:nCh,  
            picks(i:nCh:end,:) = translatePeaks( total_picks, size(total_t), quadrants(i) );
        end

        % Replace acceptor locations with aligned ones.
        picks(indA:nCh:end,:) = acceptor_picks(~remove,:);


        %---- 4. Re-estimate coordinates of acceptor-side peaks to verify alignment.
        refinedPicks = refinePeaks( image_t, picks );

        if params.geometry==2,
            r_mod = [ mod(refinedPicks(:,1),ncol/2) refinedPicks(:,2) ]; %2-color
        else
            r_mod = [ mod(refinedPicks(:,1),ncol/2) mod(refinedPicks(:,2),nrow/2) ]; %4-color
        end

        x_diff = r_mod(indA:nCh:end,1)-r_mod(indD:nCh:end,1);
        y_diff = r_mod(indA:nCh:end,2)-r_mod(indD:nCh:end,2);
        align.abs_dev = mean(  sqrt( x_diff.^2 + y_diff.^2 )  );

        % Calculate residual alignment deviation after software alignment.
        % It is easy to calculate deviation in each peak from where it
        % should be (assuming they are now aligned), but this is not the
        % same as the distance between the two peaks (D,A) in the aligned
        % image. FIXME!!!
        %residuals = picks-refinedPicks;
        %align.residual_dev = 2*mean(  sqrt( residuals(:,1).^2 + residuals(:,2).^2 )  );
    end
end


% Save output
rejectedTotalPicks = total_picks(rejected,:);
total_picks = total_picks(~rejected,:);

good = repmat( ~rejected, [nChannels,1] ); good=logical(good(:));
rejectedPicks = picks( ~good,: );
picks = picks( good,: );

fractionOverlapped = size(rejectedTotalPicks,1) / nPicked;


% end function getTraces



function wavg = WeightedAvg( vals, weights )
% Produces a weighted average of <vals>.
% Used by getTraces to find fluor centroid position

weights = weights / sum(weights);  % normalize

wavg = sum( vals.*weights );

% end function WEIGHTEDAVG



function [picks,boolRejected] = pickPeaks( image_t, threshold, overlap_thresh )
% Localizes the peaks of fluorescence.
%   picks = locations (x,y) of all molecules selected.
%   rejectedPicks = locations of molecules that are too close to a neighbor
%                        and should be ignored in analysis.
%
[nrow ncol] = size(image_t);


% 1. For each pixel (excluding edges), pick those above threshold that have
% greater intensity than their neighbors (3x3,local maxima).
% tempxy is peak position, centroidxy is estimated true molecule position.
bf = 3; %pixels around the edges to ignore.

[rows,cols] = find( image_t(1+bf:nrow-bf,1+bf:ncol-bf)>threshold );
rows = rows+bf;
cols = cols+bf;
off_x = 0; off_y = 0;

nMols=0;
for n=1:numel(rows),
    i = rows(n);
    j = cols(n);
    
    block = image_t(i-1:i+1,j-1:j+1);
    cross = block( [2,4,5,6,8] );

    if image_t(i,j)==max(cross)
        % Calc centroid position using intensity-weighted average of
        % position. This gives some additional information that is useful
        % for alignment...I think.
        if overlap_thresh~=0,
            off_x = WeightedAvg( 1:3, sum(block,1)  ) -2;
            off_y = WeightedAvg( 1:3, sum(block,2)' ) -2;
        end

        % Save the position of this molecule
        nMols=nMols+1;
        tempx(nMols)=j;
        tempy(nMols)=i;

        centroidx(nMols) = j+ off_x;
        centroidy(nMols) = i+ off_y;
    end
end

if nMols<1,
    picks = zeros(0,2);
    boolRejected = false(0,2);
    return;
end


% 2. Remove pick pairs that are closer than the cutoff radius, but still
% easily distinguishable as distinct molecules. The loop iterates in a way
% to avoid searching each pair twice and once the peaks are pretty far away
% (in terms of y-axis/rows), it stops searching for overlap for that
% molecule.
overlap=zeros(1,nMols);

if overlap_thresh~=0,
    for i=1:nMols,
        for j=i+1:nMols,
            % quickly ignore molecules way beyond the threshold
            if abs(centroidx(j)-centroidx(i)) > 1+ceil(overlap_thresh),
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
end


picks = [tempx' tempy'];
boolRejected = overlap==1; %find peaks that overlap (are rejected)


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




function [bestAlign,bestReg,quality] = alignSearch( image_t, params )
% ALIGNSEARCH    Aligns two wide-field images (donor+acceptor in FRET)
%
%    [bestAlign,bestReg] = alignSearch( image_t, params )
%
% Full search for a transformation (rotation+translation) of the acceptor
% side relative to the donor. For each combination of dx, dy, and dtheta
% (across a range in each), transform the acceptor side image, sum the
% donor and acceptor images, and calculate peak intensities (above a
% threshold). Find the one with maximum intensity, which should correspond
% to the well-aligned image.
%
% This method is slow and should only be used when the fields are far out
% of alignment.
%
% Note that this method is biased toward zero rotation because the imrotate
% method blurs out the acceptor image, so the peak intensities are
% necessarily lower. The method may also be biased toward low-FRET:
% The imrotate method blurs the acceptor field, so that high-FRET molecules
% are somewhat less likely to be picked at a set threshold. If the
% threshold is way below all the molecules, it will have no biasing effect,
% but if it is close and some molecules are not discovered, it will cause
% problems. This might be fixed by rotating BOTH fields...
%
% WARNING: This method does not handle difference in magnification across
% the fields. FIXME.
%

% TODO: displacement ranges should be stored in params/cascadeConstants.
% 
tic;

% Seperate fluorescence channels
[nrow,ncol] = size(image_t);
donor_t    = image_t( 1:nrow, 1:ncol/2 );
acceptor_t = image_t( 1:nrow, ((ncol/2)+1):end );

[nrow,ncol] = size(donor_t);

don_thresh = params.don_thresh;


% Determine search range based on parameters.
if isfield(params,'alignRotate') && params.alignRotate,
    theta_range = -4:0.25:4;
else
    theta_range = 0;
end

if isfield(params,'alignTranslate') && params.alignTranslate,
    trans_range = -5:1:5;
else
    trans_range = 0;
end

ntheta = numel(theta_range);
ntrans = numel(trans_range);


if ~params.quiet && ntheta>1,
    h = waitbar(0,'Searching for the optimal alignment');
%     tic;
end


% Space to store the best alignment for each possible rotation value.
% This is required for parfor to work correctly because all threads must be
% totally independent.
bestScores = zeros(1,ntheta);
avgScores = zeros(1,ntheta);
bestAligns = cell(1,ntheta);
bestRegs = cell(1,ntheta);


% for t=1:ntheta,  %use this for single-threaded
parfor t=1:ntheta,
    theta = theta_range(t);
    
    totalScore = 0;
    
    % Rotate the image, removing excess around the edges ('crop').
    % Rotation is done first for speed since imrotate is pretty slow and
    % the rotated image can be reused for the translations.
    % imrotate blurs out the image, which makes detecting peaks harder at a
    % set threshold. Rotating only the acceptor field, as we do here, may
    % introduce some bias, but in practice it is small. Rotating both
    % fields improves the bias, but the method doens't work very well.
    rot_a = imrotate( acceptor_t, theta, 'bicubic', 'crop' );
            
    for dx=trans_range,
        for dy=trans_range,
            % Translate the image, also removing the excess.
            registered = zeros( size(acceptor_t) );

            registered( max(1,1-dy):min(nrow,nrow-dy), max(1,1-dx):min(ncol,ncol-dx) ) = ...
                rot_a( max(1,1+dy):min(nrow,nrow+dy), max(1,1+dx):min(ncol,ncol+dx) );

            % Calculate the median intensity of pixels above the picking
            % threshold in the (aligned) total intensity image.
            % This will be higher when intensity is gathered from both
            % donor and acceptor side and low when they are seperated.
            total = donor_t+registered;
            picks = total>don_thresh;
            I = ( mean( total(picks) )-mean( total(~picks) ) ) / std(total(~picks));
            
            %scores(dx==trans_range) = I;
            totalScore = totalScore+I;
            
            % If this is the best alignment so far, save it.
            if I>bestScores(t),
                bestScores(t) = I;
                bestAligns{t} = struct('dx',dx,'dy',dy,'theta',theta);
                bestRegs{t} = total;
            end
        end
    end
    
    % Collect averages, etc for this parfor run
    avgScores(t) = totalScore/(ntrans*ntrans);
     
%     if ~params.quiet && ntheta>1,  %breaks parfor
%         waitbar( t/ntheta, h );
%     end
end
% disp(toc);


% Over all possible rotations, find the best one.
[bestScore,bestIdx] = max(bestScores);
bestAlign = bestAligns{bestIdx};
bestReg   = bestRegs{bestIdx};


% If the correct alignment is out of range or the data are just random, we
% will always get a result, but it will be meaningless. To indicate the
% quality of (or confidence in) the optimal alignment, pass along the score
% magnitude relative to other choices. Ideally, this should be more than
% 1.1-1.15 (10-15% higher intensity than a random alignment).
quality = bestScore / mean(avgScores);


if ~params.quiet && ntheta>1,
    close(h);
end





% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

function [regions,integrationEfficiency] = getIntegrationWindows( stk_top, peaks, params )

% Specify the number of most intense proximal pixels to sum when generating
% fluorescence traces (depends on experimental point-spread function).
swChoices = (1:2:19);
idx = find( params.nPixelsToSum<=(swChoices.^2), 1,'first' );
squarewidth = swChoices(idx);

% Get x,y coordinates of picked peaks
Npeaks = size(peaks,1);
x = peaks(:,1);
y = peaks(:,2);

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
% NOTE: can find which pixels to use by correlation.
% Up to this point the code is basically generic about what the channels
% mean (donor/acceptor/etc). Here with multi-color measurements, we need to
% be able to handle many possibly configurations gracefully. Donor and
% acceptor properties in the data structure should remain. Then additional
% channels can augment that main signal (single FRET pair). This is most
% often in the form of factor binding. Also have to handle each specific
% quadrant may have a different assignment (donor in one, factor or
% acceptor in another) for each experiment. GOOD TIMES! --d

global params;

movie = stkData.movie;
nFrames = movie.nFrames;
data.time = movie.timeAxis;
wbh = waitbar(0,'Extracting traces from movie data');

% Get x,y coordinates of picked peaks
Npeaks = size(peaks,1);
x = peaks(:,1);
y = peaks(:,2);

regions = stkData.regions;  % pixel#, dimension(x,y), peak#


% Create a trace for each molecule across the entire movie.
% The estimated background image is also subtracted to help with molecules
% that do not photobleaching during the movie.
traces = zeros(Npeaks,nFrames);

idx = sub2ind( [movie.nY movie.nX], regions(:,1,:), regions(:,2,:) );

for k=1:nFrames,
    frame = double( movie.readFrame(k) )  -stkData.background;
    if params.nPixelsToSum>1
        traces(:,k) = sum( frame(idx) );
    else
        traces(:,k) = diag( frame(y,x) );
    end
    
    if mod(k,100)==0,
        waitbar( k/nFrames, wbh );
    end
end


% Convert fluorescence to arbitrary units to photon counts.
constants = cascadeConstants;
if isfield(params,'photonConversion') && ~isempty(params.photonConversion) && params.photonConversion~=0,
    traces = traces./params.photonConversion;
else
    warning('Conversion from ADU to photons was not performed!');
end

% Extract individual channels from the traces matrix.
% For channel names, ignore empty strings that are placeholders for
% unused channels.
data.channelNames = params.chNames( ~cellfun( @isempty, params.chNames ) );
nCh = numel(data.channelNames);

if params.geometry==1, %single-channel    
    data.donor    = traces;
    data.acceptor = zeros( size(traces) );
    
elseif params.geometry>1, %two- or three- or four-color
    % Add each fluorescence channel to the data structure.
    for i=1:nCh,
        ch    = data.channelNames{i};
        indCh = find( strcmp(data.channelNames,ch) );
        data.(ch) = traces(indCh:nCh:end,:);
    end
    
    % Make an adjustment for crosstalk on the camera.
    % For now we assume this is the same regardless of geometry.
    if ~isfield(params,'crosstalk'),
        params.crosstalk = constants.crosstalk;
    end
    data.acceptor = data.acceptor - params.crosstalk*data.donor;
end

% Save color assignments for display in sorttraces.
emptyCh = cellfun(@isempty,data.channelNames); %unused channels
data.fileMetadata.wavelengths = params.wavelengths(~emptyCh);


% Correct for variable sensitivity across acceptor-channel camera,
% according to measurements with DNA oligos -- QZ
% This is very specific to our equipment, so it should realistically be put
% somewhere outside (cascadeConstants) as an option.
if params.geometry==2,
    % creat a Look up table for intensity correction
    for j=1:Npeaks/2,
        acc_y = y(2*j);
        yCorrection = 0.87854+acc_y*9.45332*10^(-4);
        data.acceptor(j,:) = data.acceptor(j,:)/yCorrection;
    end
end

% Subtract background and calculate FRET
[data.donor,data.acceptor,data.fret] = correctTraces(data.donor,data.acceptor);
if params.geometry>1,
    data.channelNames = [data.channelNames 'fret'];
end

% Also calculate FRET for multi-color FRET.
% FIXME: only donor->acceptor crosstalk is handled!
if isfield(data,'donor2') && isfield(data,'acceptor2'),
    [data.donor2,data.acceptor2,data.fret2] = correctTraces( ...
                                              data.donor2, data.acceptor2);
elseif isfield(data,'acceptor2')
    [~,accs,frets] = correctTraces( data.donor, {data.acceptor,data.acceptor2} );
    data.acceptor  = accs{1};
    data.acceptor2 = accs{2};
    data.fret  = frets{1};
    data.fret2 = frets{2};
end

if isfield(data,'fret2'),
    data.channelNames = [data.channelNames 'fret2'];
end



% ---- Metadata: save various metadata parameters from movie here.

data.fileMetadata.crosstalk = params.crosstalk;
% data.fileMetadata.chDesc = params.chDesc; %fixme: cells not supported!


% -- Fields specific to each trace:
% Save the locations of the picked peaks for later lookup.
nTraces = size(data.donor,1);

if params.geometry==1
    data.traceMetadata = struct( 'donor_x',num2cell(x), 'donor_y',num2cell(y) );

elseif params.geometry>1,
    % Add each fluorescence channel to the data structure.
    for i=1:nCh,
        ch    = data.channelNames{i};
        indCh = find( strcmp(data.channelNames,ch) );
        data.traceMetadata.([ch '_x']) = x(indCh:nCh:end);
        data.traceMetadata.([ch '_y']) = y(indCh:nCh:end);
    end
end


% -- Create unique trace identifiers.
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


