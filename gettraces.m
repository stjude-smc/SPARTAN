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

%------ Otherwise, load a list of files (cell array) or single file.
elseif nargin>=1 && iscell(varargin{1}) || (ischar(varargin{1}) && ~isdir(varargin{1})),
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

if ~iscell(filename),  filename = {filename};  end

[~,~,ext] = fileparts(filename{1});

% Load movie data from file.
if ~isempty( strfind(ext,'.stk') ),
    movie = Movie_STK( filename );
elseif ~isempty( strfind(ext,'.tif') ),
    movie = Movie_TIFF( filename );    
else
    error('Unrecognized file type');
end

time = movie.timeAxis;
stkX = movie.nX;
stkY = movie.nY;
nFrames = movie.nFrames;


% Average the first 10 frames of the movie to create an image (stk_top)
% to use for finding molecules.
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
    movieFilenames  = [ movieFilenames ; rdir([direct filesep '**' filesep '*.tif*']) ];
else
    movieFilenames  = rdir([direct filesep '*.stk']);
    movieFilenames  = [ movieFilenames ; rdir([direct filesep '*.tif*']) ];
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
    stk_fname = movieFilenames(i).name;
    [p,name]=fileparts(stk_fname);
    traceFname = fullfile(p, [name '.rawtraces']);
    
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
log_fid = fopen( fullfile(direct,'gettraces.log'), 'w' );

% Log parameter values used in gettraces
fprintf(log_fid,'GETTRACES PARAMETERS:\n');

output = evalc('disp(params)');
fprintf(log_fid, '%s', output);
%FIXME: structure parameters are not displayed here (alignment!)

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

% A note on notation: i, indD, indA, etc are indexes into the list of channels
% as they will appear in the output data (donor,acceptor). params.idxFields and
% quadrants identify the physical position of each channel on the camera chip.
% When looking into the image, use idxFields.
quadrants = params.idxFields;
channelNames = params.chNames;
nCh = numel(channelNames);  %# of channels TO USE.


% Define each channel's dimensions and sum fields together.
[nrow ncol] = size(image_t);  %full-chip size.

if params.geometry==1,
    allFields = image_t;

elseif params.geometry==2,
    left  = image_t( :, 1:ncol/2 );
    right = image_t( :, (ncol/2+1):end );
    allFields = cat( 3, left, right );

elseif params.geometry>2,
    upperLeft  = image_t( 1:nrow/2,       1:ncol/2       );
    upperRight = image_t( 1:nrow/2,       (ncol/2)+1:end );
    lowerLeft  = image_t( (nrow/2)+1:end, 1:ncol/2       );
    lowerRight = image_t( (nrow/2)+1:end, (ncol/2)+1:end );
    
    allFields = cat( 3, upperLeft, upperRight, lowerLeft, lowerRight );
end

% Sum fields to get a total intensity image for finding molecules.
total_t = sum( allFields(:,:,quadrants), 3 );
[nrow,ncol] = size(total_t); %from now on, this is the size of subfields.



%---- 1. Pick molecules as peaks of intensity from summed (D+A) image)
[total_picks,rejected] = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );
nPicked = size(total_picks,1);

assert( all(total_picks(:))>0, 'bad peak locations' );


% Predict the locations of peaks in all channels.
if params.geometry==1,
    picks = total_picks;

elseif params.geometry>1,
    % For each channel, predict peak locations as straight translations 
    % from the peak locations in the total intensity image.
    picks = zeros(nCh*nPicked, 2);
    
    for i=1:nCh,  
        picks(i:nCh:end,:) = translatePeaks( total_picks, size(total_t), quadrants(i) );
    end
end



%%%%% Estimation of misalignment (for dual or multi-color)
align = struct('dx',{},'dy',{},'theta',{},'sx',{},'sy',{},'abs_dev',{});
indD = find( strcmp(channelNames,'donor') ); %donor channel to align to.

if params.geometry>1,
    % Refine peak locations and how much they deviate. This helps determine
    % if realignment is needed.
    refinedPicks = refinePeaks( image_t, picks );
    r_mod = [ mod(refinedPicks(:,1),ncol) mod(refinedPicks(:,2),nrow) ];
    residuals = refinedPicks-picks;
    
    % For each channel, find a crude alignment using control points. This
    % helps determine if software alignment is needed.
    d = r_mod(indD:nCh:end,:);  %donor (reference) points
    
    for i=1:nCh,
        if i==indD, continue; end %don't try to align donor to itself.
        
        % Calculate mean deviations from donor->acceptor fields.
        dev = residuals(i:nCh:end,:)-residuals(indD:nCh:end,:);
        abs_dev = mean(  sqrt( dev(:,1).^2 + dev(:,2).^2 )  );

        % Use the picked peak locations to create a simple transformation
        % (including translation, rotation, and scaling) from donor to
        % each of the other fields.
        a = r_mod(i:nCh:end,:);

        tform = cp2tform( d(~rejected,:), a(~rejected,:), 'nonreflective similarity' );
        T = tform.tdata.T;
        theta = asind( T(1,2) );
        align(i) = struct( 'dx',T(3,1), 'dy',T(3,2), 'theta',theta, ...
                           'sx',T(1,1)/cosd(theta), 'sy',T(2,2)/cosd(theta), ...
                           'abs_dev',abs_dev );
        
        % Measure the "quality" of the alignment as the magnitude increase in
        % score compared to a "random" alignment.
        %quality =  weberQuality(ref_img,reg_img,params.don_thresh)
    end
end


%%%%% Optional software alignment algorithm (for dual-color only!)
% If the alignment is close, no need to adjust.
% Just give a warning unless asked to do software alignment in settings.
if params.geometry>1 && params.alignMethod>1,
    % FIXME (?): the user may expect the alignment to be applied even if
    % the deviation is small when a specific alignment is loaded!
    
    % 
    %---- 2. Transform acceptor side so that the two channels align properly.
    
    % Try out all possible alignments within a range and find the one with
    % the best donor-acceptor intensity overlap. The quality score is the
    % mean aligned peak magnitude vs random alignment. If the score is low,
    % reject it and just say "we don't know".
    quality = zeros(nCh,1);
    registered_t = cell(nCh,1);
    newAlign = struct('dx',{},'dy',{},'theta',{},'sx',{},'sy',{},'abs_dev',{},'tform',{});
    tform = cell(nCh,1);
    
    donor_t = allFields( :,:, params.idxFields(indD) ); %target field to align to
    total_t = donor_t;
    
    for i=1:nCh,
        if i==indD, continue; end %don't try to align donor to itself.
        target_t = allFields(:,:,params.idxFields(i));
        
        % Search for an optimal alignment of the selected field vs donor.        
        if params.alignMethod==2,
            % Nothing to search, just apply the alignment.
            newAlign(i) = params.alignment(i);
            registered_t{i} = imtransform( target_t, newAlign(i).tform, ...
                                'bicubic', 'XData',[1 ncol], 'YData',[1 nrow]);
            
        elseif params.alignMethod==3,
            % Use peaks of fluorescence as control points.
            [newAlign(i),registered_t{i}] = alignSearch_cpt( ...
                                                   donor_t, target_t, params );
                                 
        elseif params.alignMethod==4,
            % Old, slow brute force method.
            [newAlign(i),registered_t{i}] = alignSearch_weber( ...
                                                   donor_t, target_t, params );
        end
        
        total_t = total_t + registered_t{i};
        tform{i} = newAlign(i).tform;
        
        % Measure the "quality" of the alignment as the magnitude increase in
        % score compared to a "random" alignment.
        quality(i) =  weberQuality(donor_t,registered_t{i},params.don_thresh);
    end
    
    % If the optimal alignment is not trivial, re-pick molecule locations and
    % derive alignment deviation score.
    if ~( all([newAlign.dx]==0) && all([newAlign.dy]==0) && all([newAlign.theta]==0) ),
        % Give a warning for poor quality alignment.
        if any( quality<1.1 & quality>0 ),
            warning('gettraces:lowConfidenceAlignment', ...
                    'Low confidence alignment. Parameters are out of range or data quality is poor.');
            disp(quality);
        end
        
        % Pick peaks from the aligned, total intensity image.
        [total_picks,rejected] = pickPeaks( total_t, params.don_thresh, params.overlap_thresh );

        % For each channel, predict peak locations as straight translations 
        % from the peak locations in the total intensity image.
        nPicked = size(total_picks,1);
        picks = zeros( nCh*nPicked, 2 );
        remove = zeros( nPicked,1 );

        for i=1:nCh,
            [picks(i:nCh:end,:),r] = translatePeaks( total_picks, ...
                                 size(total_t), quadrants(i), tform{i} );
            remove = remove | r;
        end
        
        % Remove peaks that were moved outside the field boundries in at
        % least one of the fields due to the software alignment.
        total_picks  = total_picks(~remove,:);
        rejected     = rejected(~remove);
        nPicked = size(total_picks,1);

        % Remove molecules if the peaks in any field are out of range.
        z = repmat( remove', [nCh,1] ); z = z(:);
        picks = picks(~z,:);


        %---- 4. Re-estimate coordinates of acceptor-side peaks to verify alignment.
        % Then normalize so that the "expected" location of each peak is (0,0).
        % The rmsd is then the distance between each channel and the donor
        refinedPicks = refinePeaks( image_t, picks );
        residuals = refinedPicks-picks; 
        
        for i=1:nCh,
            dev = residuals(i:nCh:end,:)-residuals(indD:nCh:end,:);
            newAlign(i).abs_dev = mean(  sqrt( dev(:,1).^2 + dev(:,2).^2 )  );
            newAlign(i).quality = quality(i);
        end
        align = newAlign;
    end
end


% Save output
rejectedTotalPicks = total_picks(rejected,:);
total_picks = total_picks(~rejected,:);

good = repmat( ~rejected, [nCh,1] ); good=logical(good(:));
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

assert( ~any(tempx(:)<bf|tempy(:)<bf) ); 


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




function [align,reg_img] = alignSearch_cpt( ref_img, target_img, params )
% ALIGNSEARCH    Align two images using a simple control point algorithm.
%
%    [bestAlign,bestReg,quality] = alignSearch( REF, TARGET, params)
%
% Finds all bright spots in the field-of-view, picking on each field separately,
% then uses nearest neighbors in each field as control points for creating a
% transformation matrix. This allows translation, rotation, scaling, and some
% non-linear distortions. Generally this is only possible with beads samples,
% where there is a high SNR signal on all channels and no background junk to
% confuse the algorithm. This will not work well for single fluorophores.
%
% TODO: allow multiple arguments for aligning multiple fields. The first image
% will be the reference.
%

[nrow,ncol] = size(ref_img);
assert( all(size(ref_img)==size(target_img)) );

% 1) Detect beads as brights spots above background.
ref_peaks    = pickPeaks( ref_img,    0.7*params.don_thresh, params.overlap_thresh );
target_peaks = pickPeaks( target_img, 0.7*params.don_thresh, params.overlap_thresh );


% 2) For each peak in the reference, find the corresponding peak in the target.
% K-nearest neighbor search (Statistics toolbox), where K=1. Using K>1 could be
% used to verify there is no uncertainty (other pints are far away).
% This will not work if the beads are too dense or if the distortion is too
% severe -- the "correct" choice has to be closer than nearby beads.
[idx,dist] = knnsearch( target_peaks, ref_peaks );
target_peaks = target_peaks(idx,:);
sel = dist<params.nPixelsToSum;

idx = idx(sel);
target_peaks = target_peaks(sel,:);
ref_peaks    = ref_peaks(sel,:);

% Find any reference peaks that use the same target peak multiple times. Such
% points clearly have some uncertainty and should be thrown out.
% Apparently this does not help. 
if ~params.quiet,
    u = unique(idx);
    n = histc(idx,u); %count how many times each index is used
    fprintf('alignSearch_cpt: %.0f (%0.0f%%) were used; of those %.0f (%.0f%%) are ambiguous\n', ...
             numel(idx), 100*numel(idx)/numel(dist), sum(n>1), 100*sum(n>1)/numel(n) );
end

% idx = u(n==1);
% ref_peaks    = ref_peaks(idx,:);
% target_peaks = target_peaks(idx,:);


% 3) Generate the transformation matrix from control points
assert( all(size(ref_peaks)==size(target_peaks)) );
tform = cp2tform( target_peaks, ref_peaks, 'nonreflective similarity' );


% 4) Transform the target field image so it is aligned with the reference.
% X/YData are needed so imtransform crops the output to the size of the input.
reg_img = imtransform( target_img, tform, 'bicubic', ...
                          'XData',[1 ncol], 'YData',[1 nrow]);
                      
                      
% 5) Calculate distortion parameters from tform.
ss = tform.tdata.Tinv(2,1);
sc = tform.tdata.Tinv(1,1);
tx = tform.tdata.Tinv(3,1);
ty = tform.tdata.Tinv(3,2);
scale = sqrt(ss*ss + sc*sc);
theta = atan2(ss,sc)*180/pi;

align = struct( 'dx',tx,'dy',ty, 'theta',theta, 'sx',scale,'sy',scale, ...
                'abs_dev',0, 'tform',tform );

% Measure the "quality" of the alignment as the magnitude increase in
% score compared to a "random" alignment.
%quality =  weberQuality(ref_img,reg_img,params.don_thresh);



% END FUNCTION %function alignSearch_cpt




function [bestAlign,bestReg] = alignSearch_weber( donor_t, acceptor_t, params )
% ALIGNSEARCH    Aligns two wide-field images (donor+acceptor in FRET)
%
%    [bestAlign,bestReg] = alignSearch( image_t, params, [forceAlign] )
%
% Full search for a transformation (rotation+translation) of the acceptor
% side relative to the donor. For each combination of dx, dy, and dtheta
% (across a range in each), transform the acceptor side image, sum the
% donor and acceptor images, and calculate a contrast score. The one with
% the brightest, sharpest peaks wins.
%


% These parameters used to be in cascadeConstants, but since they are very
% specific to this algorithm, I moved them here. This makes some of the code
% below redundant!
params.alignment.theta = -4:0.1:4;
params.alignment.dx    = -4:1:4;
params.alignment.dy    = -4:1:4;  %all dx and dy must be integers
params.alignment.sx    = 1;  %magnification (x)
params.alignment.sy    = 1;  %magnification (y)
params.alignRotate    = true;
params.alignTranslate = true;

[nrow,ncol] = size(donor_t);


% Get parameter search ranges from input arguments, if available.
if isfield(params,'alignment') && ~isempty(params.alignment),
    dx_range = params.alignment.dx;
    dy_range = params.alignment.dy;
    theta_range = params.alignment.theta;
else
    % If no parameters given, disable search.
    theta_range = 0;
    dx_range = 0;
    dy_range = 0;
end

% Disable search for parameters if the settings say not to search.
if ~isfield(params,'alignRotate') || ~params.alignRotate,
    theta_range = 0;
end

if ~isfield(params,'alignTranslate') || ~params.alignTranslate,
    dx_range = 0;
    dy_range = 0;
end


% Round translation to integers; the algorithm do subpixel shifts.
dx_range = round(dx_range);
dy_range = round(dy_range);

ntheta = numel(theta_range);
ndx = numel(dx_range);
ndy = numel(dy_range);

if ~params.quiet && ntheta>1,
    h = waitbar(0,'Searching for the optimal alignment');
end


% Reserve space for the best alignment for each possible rotation value.
% This is required for parfor to work correctly because all threads must be
% totally independent.
%scores = zeros(ntheta,ndx,ndy);
bestScores = zeros(ntheta,1);
bestAligns = cell(ntheta,1);
bestRegs = cell(ntheta,1);

don_thresh = params.don_thresh;

for t=1:ntheta,  %use this for single-threaded
% parfor t=1:ntheta,
    theta = theta_range(t);
    scoreTemp = zeros(ndx,ndy); %contrast scores for all translations.
    
    % Rotate the image, removing excess around the edges ('crop').
    % Rotation is done first for speed since the rotated image can be
    % reused for the translations. imrotate blurs out the image, which
    % lowers the contrast ratio, so there is a slight bias toward no
    % theta=0. Rotating both fields by some large amount outside the search
    % range can remove the bias, but in practice this isn't necessary.
    % There may also be a minor bias toward low FRET using this since the
    % acceptor side will have weaker intensity.
    rot_a = imrotate( acceptor_t, theta, 'bicubic', 'crop' );
    
    for i=1:ndx,
        dx = dx_range(i);
            
        for j=1:ndy,
            dy = dy_range(j);
            
            % Translate the image, also removing the excess.
            registered = zeros( size(acceptor_t) );

            registered( max(1,1-dy):min(nrow,nrow-dy), max(1,1-dx):min(ncol,ncol-dx) ) = ...
                rot_a( max(1,1+dy):min(nrow,nrow+dy), max(1,1+dx):min(ncol,ncol+dx) );

            % Calculate Weber contrast as a score. When the fields are not
            % aligned, the peak intensities go down and some peaks fall
            % below the threshold of background, lowering contrast.
            total = donor_t+registered;
            picks = total>don_thresh;    %pixels above background
            Ib = mean( total(~picks) );  %background intensity
            S = ( mean(total(picks)) -Ib ) / Ib;  %Weber contrast score
            
            % If this is the best alignment so far, save it.
            if all( S>scoreTemp(:) ),
                bestAligns{t} = struct('dx',dx,'dy',dy,'theta',theta,'sx',1,'sy',1,'abs_dev',0);
                bestRegs{t} = registered;
                bestScores(t) = S;
            end
            
            scoreTemp(i,j) = S;  %fixme for direct indexing.
        end
    end
    
    %scores(t,:,:) = scoreTemp;
    
    if ~params.quiet && ntheta>1,  %do not use with parfor
        waitbar( t/ntheta, h );
    end
end


% Over all possible rotations, find the best one.
[~,bestIdx] = max(bestScores);
bestAlign = bestAligns{bestIdx};
bestReg   = bestRegs{bestIdx};
        

% Create a transformation matrix from the alignment parameters.
Ts = ones(3);
Ts(1,1) = 1;  % x scaling
Ts(2,2) = 1;  % y scaling

T = [  cosd(bestAlign.theta)  sind(bestAlign.theta)  0 ; ...
      -sind(bestAlign.theta)  cosd(bestAlign.theta)  0 ; ...
            bestAlign.dx           bestAlign.dy      1 ];
T = T.*Ts;
bestAlign.tform = fliptform( maketform('affine',T) );


% Measure the "quality" of the alignment as the magnitude increase in
% score compared to a "random" alignment.
% quality = max(bestScores) / weberRandom(donor_t,acceptor_t,don_thresh);
%quality =  weberQuality(donor_t,bestReg,don_thresh);


if ~params.quiet && ntheta>1,
    close(h);
end


% END FUNCTION alignSearch



function [quality,randomScore] = weberQuality( base, registered, thresh )
% Calculates a Weber contrast score for random scampling of the target image.
% This provides a way to calculate the score expected for a random (poor)
% alignment. quality = weberScore/weberRandom(X,Y).
% Measure the "quality" of the alignment as the magnitude increase in
% score compared to a "random" alignment, which is approximated by
% scrambling the donor and acceptor field data. In practice, this gives
% similar scores to the mean of all scores from an alignment search across
% a wide parameter range. We do not just average the scores in the search
% because there may only be one parameter value or a small range. Quality
% scores above 1.1-1.15 (10-15% higher contrast than a random) are
% generally acceptable; any lower and the data may have other problems.
%

N = 6; %number of repititions for averaging.
S = zeros(N,1);

for i=1:N,
    if i==1, %base case, no scrambling
        total = base(:)+registered(:);
    else
        total = base( randperm(numel(base)) ) + registered( randperm(numel(registered)) );
    end
    
    Ib = mean( total(total<thresh) );  %background intensity
    S(i) = ( mean(total(total>thresh)) -Ib ) / Ib;  %Weber contrast score
end

randomScore = mean(S(2:end));
quality = S(1) / randomScore;


% END FUNCTION weberRandom




% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

function [regions,integrationEfficiency] = getIntegrationWindows( stk_top, peaks, params )
% For each molecule location in "peaks", find the most intense pixels in
% its immediate neighborhood (defined by params.nPixelsToSum). These
% regions are used by integrateAndSave() to sum most of the intensity for
% each peak and generate fluorescence-time traces.
%

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




function integrateAndSave( stkData, peaks, stk_fname )
% For each location in "peaks", sum a region that includes most of the
% intensity for that spot and repeat for each time point to get a
% fluorescence-time trace for each peak. Background subtraction, crosstalk
% correction, and calculation of derived signals (FRET traces) is all done
% here. Then the result is saved as a .rawtraces file with metadata.
%

global params;

movie = stkData.movie;
nFrames = movie.nFrames;
wbh = waitbar(0,'Extracting traces from movie data');

% Get x,y coordinates of picked peaks
Npeaks = size(peaks,1);
x = peaks(:,1);
y = peaks(:,2);

regions = stkData.regions;  % pixel#, dimension(x,y), peak#



% Create channel name list for the final data file. This includes FRET channels,
% which are not in the movie. chNames includes only fluorescence fields.
chNames = params.chNames( ~cellfun(@isempty,params.chNames) );
nCh = numel(chNames);
nTraces = Npeaks/nCh;

dataNames = chNames;

if ismember('acceptor',dataNames),
    dataNames = [dataNames 'fret'];
end

if ismember('acceptor2',dataNames),
    dataNames = [dataNames 'fret2'];
end

% Create traces object, where the data will be stored.
if params.geometry>2 && ismember('donor',dataNames)
    data = TracesFret4(nTraces,nFrames,dataNames);
elseif ismember('donor',dataNames)
    data = TracesFret(nTraces,nFrames,dataNames);
else
    data = TracesFluor(nTraces,nFrames,dataNames);
end

data.time = movie.timeAxis;

% Ask the user
if ~params.quiet && data.time(1)==1,
    disp('Time axis information is not present in movie!');
    a = inputdlg('Time resolution (ms):','No time axis in movie!');
    a = str2double(a);
    if ~isempty(a), %empty if user hit Cancel or entered an invalid number
        data.time = a*( 0:nFrames-1 );
    else
        disp('Using frames as the time axis');
    end
end


% Create a trace for each molecule across the entire movie.
% The estimated background image is also subtracted to help with molecules
% that do not photobleaching during the movie.
traces = zeros(Npeaks,nFrames,'single');

idx = sub2ind( [movie.nY movie.nX], regions(:,1,:), regions(:,2,:) );

for k=1:nFrames,
    frame = single( movie.readFrame(k) )  -stkData.background;
    
    if params.nPixelsToSum>1
        traces(:,k) = sum( frame(idx) );
    else
        traces(:,k) = diag( frame(y,x) );
    end
    
    if mod(k,100)==0,
        waitbar( 0.9*k/nFrames, wbh );
    end
end


% Convert fluorescence to arbitrary units to photon counts.
if isfield(params,'photonConversion') && ~isempty(params.photonConversion) && params.photonConversion~=0,
    traces = traces./params.photonConversion;
else
    warning( 'gettraces:noPhotonConversion', ...
             'Conversion from ADU to photons was not performed!' );
end


% Extract individual channels from the traces matrix.
% For channel names, ignore empty strings that are placeholders for
% unused channels.
if params.geometry==1, %single-channel    
    data.donor    = traces;
    data.acceptor = zeros( size(traces), 'single' );
    
elseif params.geometry>1, %two- or three- or four-color
    % Add each fluorescence channel to the data structure.
    for i=1:nCh,
        data.( data.channelNames{i} ) = traces(i:nCh:end,:);
    end
    
    % Make an adjustment for crosstalk on the camera. If there are 3 or 4
    % colors, params.crosstalk is a matrix of (src,dst) values, where most
    % are 0 (no crosstalk).
    if numel(params.crosstalk)==1,
        data.acceptor = data.acceptor - params.crosstalk*data.donor;
    else
        [src,dst] = find( params.crosstalk~=0 );
        
        for i=1:numel(src),
            chSrc = data.channelNames{ src(i) };
            chDst = data.channelNames{ dst(i) };
            val = params.crosstalk(src(i),dst(i));
            data.(chDst) = data.(chDst) - val*data.(chSrc);
        end
    end
end


% Correct for non-uniform sensitivity across the field-of-view in the donor
% and acceptor fields. This is generally fixed and determined by the
% optical properties in the light path and sometimes the CCD chips. The
% parameters in cascadeConstants (.biasCorrection) are lamda functions that
% give a scaling factor for each point in the field-of-view. To generate
% this, find functions that give a flat response profile. The elements in
% the lamba function cell array are in the same order as the channels.
if isfield(params,'biasCorrection') && ~isempty(params.biasCorrection),
    
    for i=1:nCh,
        chName = data.channelNames{i};
        corr = params.biasCorrection{i}( x(i:nCh:end), y(i:nCh:end) );
        data.(chName) = data.(chName)  ./  repmat( corr, [1 nFrames] );
    end
    
end


% Subtract background, correct for crosstalk, and calculate FRET
data = correctTraces(data);



% ---- Metadata: save various metadata parameters from movie here.

% Keep only parameters for channels that are being analyzed, not all
% possible channels in the configuration.
chToKeep = ~cellfun(@isempty,params.chNames);

data.fileMetadata(1).wavelengths = params.wavelengths(chToKeep);

if numel(params.crosstalk)>1
    data.fileMetadata.crosstalk = params.crosstalk(chToKeep,chToKeep);
end
% data.fileMetadata.chDesc = params.chDesc(chToKeep); %fixme: cells not supported!


% -- Fields specific to each trace:
% Save the locations of the picked peaks for later lookup.
nTraces = size(data.donor,1);

if params.geometry==1
    data.traceMetadata = struct( 'donor_x',num2cell(x), 'donor_y',num2cell(y) );

elseif params.geometry>1,
    % Add each fluorescence channel to the data structure.
    for i=1:nCh,
        ch = data.channelNames{i};
        data.traceMetadata().([ch '_x']) = x(i:nCh:end);
        data.traceMetadata().([ch '_y']) = y(i:nCh:end);
    end
end


% -- Create unique trace identifiers.
% This is just the full path to the movie plus a trace number. This can be
% used to later find the corresponding original movie data for each
% individual trace, even after many rounds of processing.
for i=1:nTraces,
    data.traceMetadata(i).ids = sprintf( '%s#%d', stk_fname, i );
end

% ---- Save data to file.
[p,name]=fileparts(stk_fname);
save_fname = fullfile(p, [name '.rawtraces']);

saveTraces( save_fname, 'traces', data );


close( wbh );

% end function integrateAndSave


