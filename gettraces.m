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

% TODO: also return unfiltered peak list (overlap=0) so overlap statistics
% can be displayed in scripts that call this function.


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
    
    % Specify default channel assignments if none given.
    % Setting it here is better in case geometry but not names are set in
    % userParams because the default includes names for 2-channel.
    if ~isfield(userParams,'chNames') && params.geometry>2
        params.chNames = constants.gettraces_chNames4;
        params.chDesc  = constants.gettraces_chDesc4;
    end
    
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
[peaks,stkData.total_t,stkData.alignStatus] = getPeaks( image_t, params );

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
den = min( floor(stkX/32), 32 );
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
function [picks,total_t, align] = getPeaks( image_t, params )
% Localizes the peaks of molecules from a summed image of the two channels
% (in FRET experiments). The selection is made on the total fluorescence
% intensity image (summing all channels into a single image) to minimize
% bias (e.g., molecules with a dim donor and bright acceptor would be
% missed if we selected just on the donor side). This inital selection is
% refined by looking in the immediate neighborhood of each spot for the
% actual intensity maximum, which could be different on each side if the
% two fields are not aligned or if there are optical distortions.
% 
% If translation is detected, the acceptor channel image is translated, and
% the channels are summed again for molecule selection. When the molecules
% are selected a final time, the deviations are used once again to create a
% transformation describing any deviations. This map is then used to choose
% the final acceptor-side locations (optional -- "refine align"). "Refine
% align" may not work very well since it uses a transformation of donor to
% acceptor side, when what we want is total->donor and total->acceptor and
% move BOTH sets of selections... FIXME
% 
%
% picksX - X-coords of all molcules, as Cy3,Cy5,Cy3,Cy5 in order
% picksY - Y-coords ...
% For quad-color, UL,LL,LR,UR is the order.

nChannels = params.geometry;
align = [];

% If using a single-channel setup, use a simplified proceedure
% since we don't have to worry about alignment, etc.
if params.geometry==1,
    total_t = image_t;
    picks = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
    assert( size(picks,2)==2 ); %x and y columns, molecules in rows.
    return;
end
% Otherwise, assume dual-channel. FIXME


% Define each channel's dimensions and sum fields together.
% For dual-channel, we assume they are arranged L/R.
% For quad-channel, it is UL/LL/LR/UR
[nrow ncol] = size(image_t);

if params.geometry==2,
    donor_t    = image_t( :, 1:ncol/2 );
    acceptor_t = image_t( :, (ncol/2+1):end );

    total_t = acceptor_t + donor_t; %sum the two fluorescence channels.

elseif params.geometry>2,
    upperLeft  = image_t( 1:nrow/2, 1:ncol/2 );
    lowerLeft  = image_t( (nrow/2)+1:end, 1:ncol/2 );
    lowerRight = image_t( (nrow/2)+1:end, (ncol/2)+1:end );
    upperRight = image_t( 1:nrow/2, (ncol/2)+1:end );
    
    total_t = upperLeft + lowerLeft + lowerRight + upperRight;
end



%---- 1. Pick molecules as peaks of intensity from summed (D+A) image)
total_picks = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
nPicked = size(total_picks,1);

assert( all(total_picks(:))>0, 'bad peak locations' );

% Define acceptor size peaks as a simple translation across the chip.
% Picks is a cell array of peak locations, one per channel. Some channels
% (cell elements) may be empty if that channel isn't used.
picks = ones(nChannels*nPicked, 2); %donor, acceptor alternating; 2 columns = x,y

if params.geometry==2
    picks(1:2:end,:) = total_picks; % L (donor)

    picks(2:2:end,1) = total_picks(:,1) + ncol/2; % R (acceptor) x
    picks(2:2:end,2) = total_picks(:,2);          % R (acceptor) y

elseif params.geometry>2,
    chNames = params.chNames;
    k = 1;
    
    % FIXME: choose which fields to pick based on naming!!
    if ~isempty(chNames{1}),
        picks(k:nChannels:end,:) = total_picks; % UL
        k = k+1;
    end
    
    if ~isempty(chNames{2}),
        picks(k:nChannels:end,1) = total_picks(:,1);          % LL x
        picks(k:nChannels:end,2) = total_picks(:,2) + nrow/2; % LL y
        k = k+1;
    end
    
    if ~isempty(chNames{3}),
        picks(k:nChannels:end,1) = total_picks(:,1) + ncol/2; % LR x
        picks(k:nChannels:end,2) = total_picks(:,2) + nrow/2; % LR y
        k = k+1;
    end
    
    if ~isempty(chNames{4}),
        picks(k:nChannels:end,1) = total_picks(:,1) + ncol/2; % UR x
        picks(k:nChannels:end,2) = total_picks(:,2);          % UR y
    end
end


align = [];

% FIXME: this only works for donor/acceptor (FRET) channels. The factor
% channel may not have any intensity to do the alignment calculation...
% There is also no crosstalk that could be used to measure it...
if params.geometry>1,
    
    % Get indexes for donor and acceptor channels. This is really only
    % needed for three-color.
    channelNames = params.chNames( ~cellfun( @isempty, params.chNames ) );
    nCh = numel(channelNames);
    assert( nCh==nChannels );
    indD = find( strcmp(channelNames,'donor') );
    indA = find( strcmp(channelNames,'acceptor') );
    
    % Refine peak locations to determine if realignment is needed.
    % The _abs alignment calculation will detect rotation and other distortions.
    % Using rem() deals with the offsets for each channel in the field-of-view.
    refinedPicks = refinePeaks( image_t, picks );
    
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
    rot = asin(T(1,2)) *180/pi;
    
    align_abs = mean( sqrt( x_diff.^2 + y_diff.^2 )  );    
    align = [dx dy align_abs rot mean([sx,sy])]
    
    
    save('images.mat','donor_t','acceptor_t','total_t','tform');
    

    % If the alignment is close (by translation), we're done.
    if all( abs(align(1:2))<=0.5 ),  return;  end

    warning('gettraces:badAlignment','Fluorescence fields may be out of alignment.');
    fprintf( 'Alignment deviation: %.1f (X), %0.1f (Y), %.1f (absolute)\n', [dx dy align_abs] );
   
    
    % If specified, use the refined peak positions to handle slight misalignment.
    % This generally causes more harm than good and is not recommended!!
    if params.refineAlign,
       picks = refinedPicks;
    end
    
    
    % Just give a warning unless asked to do software alignment in settings.
    % FIXME: code below only works for two-channel FRET.
    if ~params.alignTranslate || params.geometry>2,  return;  end 
         

    %---- 2. Transform acceptor side so that the two channels align properly.
    % If the alignment is off by a significant margin, the fields are
    % realigned in software. This will only handle translation by 1 pixel.
    % Rotation and other complex distortions are harder. FIXME.
    
    dx = round(dx); dy = round(dy); %can only shift fields by integer amounts.
    
    % Sum the donor and acceptor fields, after translating the acceptor size to
    % deal with the misalignment.
    align_t = zeros( size(donor_t) );
    align_t( max(1,1-dy):min(nrow,nrow-dy), max(1,1-dx):min(ncol/2,ncol/2-dx) ) = ...
        acceptor_t( max(1,1+dy):min(nrow,nrow+dy), max(1,1+dx):min(ncol/2,ncol/2+dx) );

    total_t = align_t + donor_t; %sum the two fluorescence channels.

    % Pick peaks from the aligned image.
    total_picks = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
    nPicked = numel(total_picks)/2;
    
    % Since the alignment image has been shifted to compensate for misalignment
    % already (above), adjust the output coordinates so they are relative to
    % the actual fields, not the adjusted fields.
    picks = zeros(nPicked*2,2); %donor, acceptor alternating; 2 columns = x,y

    picks(1:2:end,1) = total_picks(:,1);              %donor x
    picks(2:2:end,1) = total_picks(:,1) + ncol/2 +dx; %acceptor x
    picks(1:2:end,2) = total_picks(:,2);     %donor y
    picks(2:2:end,2) = total_picks(:,2) +dy; %acceptor y


    %---- 4. Re-estimate coordinates of acceptor-side peaks to verify alignment.
    refinedPicks = refinePeaks( image_t, picks );

    % Use the picked peak locations to create a simple transformation
    % (including translation, rotation, and scaling) from donor to acceptor.
    r_mod = [ mod(refinedPicks(:,1),ncol/2) refinedPicks(:,2) ]; %2-color
    tform = cp2tform( r_mod(indD:nCh:end,:), r_mod(indA:nCh:end,:), ...
                      'nonreflective similarity');

    T = tform.tdata.T;
    dx = T(3,1);  dy = T(3,2);  sx = T(1,1);  sy = T(2,2);
    rot = asin(T(1,2)) *180/pi;

    % Verify the alignment
    x_diff = r_mod(indA:nCh:end,1)-r_mod(indD:nCh:end,1);
    y_diff = r_mod(indA:nCh:end,2)-r_mod(indD:nCh:end,2);
    align_abs = mean(  sqrt( x_diff.^2 + y_diff.^2 )  );
    align = [dx dy align_abs rot mean([sx sy])]
    

    % Use actual peak locations if option is specified.
    if params.refineAlign,
       picks = refinedPicks;
    end

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

disp(threshold);


% 1. For each pixel (excluding edges), pick those above threshold that have
% greater intensity than their neighbors (3x3,local maxima).
% tempxy is peak position, centroidxy is estimated true molecule position.
[rows,cols] = find( image_t(1+3:nrow-3,1+3:ncol-3)>threshold );
rows = rows+3;
cols = cols+3;

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

if nMols<1,
    clean_picks = zeros(0,2);
    all_picks = clean_picks;
    return;
end


% 2. Remove pick pairs that are closer than the cutoff radius, but still
% easily distinguishable as distinct molecules. The loop iterates in a way
% to avoid searching each pair twice and once the peaks are pretty far away
% (in terms of y-axis/rows), it stops searching for overlap for that
% molecule.
overlap=zeros(1,nMols);

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

indDonor    = find( strcmp(data.channelNames,'donor') );
indAcceptor = find( strcmp(data.channelNames,'acceptor') );
indFactor   = find( strcmp(data.channelNames,'factor') );

if params.geometry==1, %single-channel    
    donor    = traces;
    acceptor = zeros( size(traces) );
    
elseif params.geometry==2, %dual-channel
    donor    = traces(indDonor:2:end,:);
    acceptor = traces(indAcceptor:2:end,:);

elseif params.geometry==3, %three-color
    donor    = traces(indDonor:3:end,:);
    acceptor = traces(indAcceptor:3:end,:);
    factor   = traces(indFactor:3:end,:);
end

% Make an adjustment for crosstalk on the camera.
% For now we assume this is the same regardless of geometry.
% TODO: in extreme cases where the fields are not aligned, there will be no
% crosstalk on the other side, giving negative fluorescence values. If
% these were detected, it could be a good indication for a warning.
if ~isfield(params,'crosstalk'),
    params.crosstalk = constants.crosstalk;
end
acceptor = acceptor - params.crosstalk*donor;

% Background substraction for third channel for factor binding.
% This is tricky because it isn't dependant on the donor and very often may
% not have a distinct "photobleaching" event.
% So I just naively substract from the last 10 frames...
% There are better ways to do this. FIXME
if params.geometry==3,
    data.factor = factor;
    for i=1:size(factor,1),
        data.factor(i,:) = factor(i,:) - mean( factor(i,end-9:end) );
    end
end

% Correct for variable sensitivity across acceptor-channel camera,
% according to measurements with DNA oligos -- QZ
% This is very specific to our equipment, so it should realistically be put
% somewhere outside (cascadeConstants) as an option.
if params.geometry==2,
    % creat a Look up table for intensity correction
    for j=1:Npeaks/2,
        acc_y = y(2*j);
        yCorrection = 0.87854+acc_y*9.45332*10^(-4);
        acceptor(j,:) = acceptor(j,:)/yCorrection;
    end
end

% Subtract background and calculate FRET
[data.donor,data.acceptor,data.fret] = correctTraces(donor,acceptor,constants);
data.channelNames = [data.channelNames 'fret'];


% ---- Metadata: save various metadata parameters from movie here.

% -- Fields specific to this movie.
% data.movieMetadata.filename  = strrep(stk_fname,'.bz2','');
% data.movieMetadata.crosstalk = params.crosstalk;

% -- Fields specific to each trace:
% Save the locations of the picked peaks for later lookup.
% Also save indexes to map traces to movie metadata (here, everything is 1
% because there is only one movie).
% FIXME: don't assume channel orders! This is wrong! Should make this
% generic by using channel names here!
nTraces = size(data.donor,1);

if params.geometry==1
    data.traceMetadata = struct( 'donor_x',num2cell(x), 'donor_y',num2cell(y) );
    
elseif params.geometry==2,
    data.traceMetadata = struct( ...
        'donor_x',    num2cell( x(indDonor:2:end) ),    'donor_y',    num2cell( y(indDonor:2:end) ), ...
        'acceptor_x', num2cell( x(indAcceptor:2:end) ), 'acceptor_y', num2cell( y(indAcceptor:2:end) )  ...
    );

elseif params.geometry==3,
    data.traceMetadata = struct( ...
        'donor_x',    num2cell( x(indDonor:3:end) ),    'donor_y',    num2cell( y(indDonor:3:end)    ), ...
        'acceptor_x', num2cell( x(indAcceptor:3:end) ), 'acceptor_y', num2cell( y(indAcceptor:3:end) ), ...
        'factor_x',   num2cell( x(indFactor:3:end) ),   'factor_y',   num2cell( y(indFactor:3:end)   )  ...
    );
end
% z = num2cell( ones(1,nTraces) );
% [data.fileMetadata.movieIndex] = deal( z{:} );


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


