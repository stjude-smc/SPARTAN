function [stkData,peaks,image_t] = gettraces(varargin)
% GETTRACES  Extract smFluorescence traces from movies
%
%    GETTRACES()
%    Launches the gettraces graphical user interface.
%
%    [STK,PEAKS,IMG] = GETTRACES( FILENAME, PARAMS )
%    Loads a movie from FILENAME, finds fluorescence peaks, and returns
%    their locations (PEAKS) as a Nx2 matrix (x,y).  IMG is the image used
%    for selecting intensity peaks (sum of aligned fluorescence channels).
%    STK is the stkData structure containing the internal state.
%
%    [STK,PEAKS,IMG] = gettraces( FILENAME, PARAMS, OUTFILE )
%    As above, but also extracts traces from the PEAKS locations and saves
%    the resulting traces to OUTFILE.
%
%    STK can also be passed instead of FILENAME if the movie has already
%    been loaded by gettraces. This saves on load time.
%
%    GETTRACES( DIRECTORY, PARAMS )
%    This is "batch mode".  For each filename in DIRECTORY, the movie is loaded,
%    processed, and a .rawtraces file is saved automatically. Additional
%    (optional) fields in PARAMS include option to check all child folders
%    (recursive) and option to skip processed movies (skipExisting).
%
%
%  PARAMS is a struct with any of the following options:
% 
% - don_thresh:     total (D+A) intensity threshold for pick selection
%                   (if 0, threshold is automatically selected).
% 
% - overlap_thresh: peaks are rejected if they are closer than this radius.
% 
% - nPixelsIntegrated: number of pixels proximal to each peak to sum
%      to produce fluorescence traces. Higher values capture more
%      intensity, but also capture more noise...
% 
% - saveLocations: write a text file with the locations of each peak to
%      file.
% 
% - geometry: imaging geometry can be single-channel/full-chip (1), 
%             dual-channel/left-right (2), or quad-channel (3/4).
%             Default: dual-channel (2).
% 
% - crosstalk: donor-to-acceptor channel fluorescence channel
%              bleedthrough (crosstalk) as a fraction. For correction.
% 
% - photonConversion: fluorescence ADU/photon conversion factor.
% 
% - fieldNames: cell array of assignments of fluorescence field names
%                 (e.g., donor, acceptor, factor). If one is left blank
%                 that field is disregarded. OPTIONAL.
% 

% NOTE: because of the structure of this file, all variables defined in this
% initial section are global to the other sub-functions in the file. This way,
% common parameters (params, stkData, constants) do not have to be passed
% around. FIXME: avoid using common names where possible to pervent confusion,
% including image_t, peaks, and picks.


% If calling directly from command line, launch the GUI.
if nargin==0,
    gettraces_gui;
    return;
end


%------ Load parameter values (always second parameter!)
% Set default parameter values and merge with user-specified values.
% The user's options will override any existing defaults.
constants = cascadeConstants;
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

% Find peak locations from total intensity.
% outputs: total_t, alignStatus, total_peaks, fractionOverlapped, 
%    rejectedPicks, rejectedTotalPicks
peaks = getPeaks( image_t );

% Generate integration windows for later extracting traces.
% outputs: regions, integrationEfficiency, fractionOverlap
getIntegrationWindows(image_t, peaks);


%------ Extract traces using peak locations
if nargin>=3 && ischar(varargin{3}),
    outputFilename = varargin{3};
    integrateAndSave( peaks, outputFilename );
end
    
return;









%% --------------------- OPEN STK MOVIE (CALLBACK) --------------------- %

function [stkData] = OpenStk(filename)


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
endBackground = movie.readFrames(nFrames-11:nFrames-1);

if params.geometry==1,
    % No combining needed for single-color imaging.
    stkData.nChannels = 1;
elseif params.geometry==2,
    % Combine left and right side of field for two-color imaging.
    endBackground = endBackground(:,1:stkX/2,:) + endBackground(:,(1+stkX/2):end,:);
    stkData.nChannels = 2;
elseif params.geometry>2,
    % Combine fluorescence from the four quadrants.
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

end %FUNCTION OpenStk




% --------------- OPEN ALL STK'S IN DIRECTORY (BATCH) --------------- %

function batchmode(direct)

% Get list of files in current directory (option: and all subdirectories)
movieFiles = regexpdir(direct,'^.*\.(tiff?|stk)$',params.recursive);
nFiles = length(movieFiles);

% Wait for 100ms to give sufficient time for polling file sizes in the
% main loop below.
pause(0.1);


% Process each file in the user selected directory.
nTraces  = zeros(nFiles,1); % number of peaks found in file X
existing = zeros(nFiles,1); % true if file was skipped (traces file exists)

for i=1:nFiles,
    
    % Skip if previously processed (.rawtraces file exists)
    stk_fname = movieFiles(i).name;
    [p,name]=fileparts(stk_fname);
    traceFname = fullfile(p, [name '.rawtraces']);
    
    if params.skipExisting && exist(traceFname,'file'),
        if ~params.quiet,
            disp( ['Skipping (already processed): ' stk_fname] );
        end
        existing(i) = 1;
        continue;
    end
    
    % Poll the file to make sure it isn't changing.
    % This could happen when a file is being saved during acquisition.
    d = dir(movieFiles(i).name);
    if movieFiles(i).datenum ~= d(1).datenum,
        disp( ['Skipping (save in process?): ' movieFiles(i).name] );
        existing(i) = 1;
        continue;
    end
    
    % Show waitbar only when new data must be loaded.
    if ~exist('h','var') && ~params.quiet,
        h = waitbar( (i-1)/nFiles,'Extracting traces from movies...');
    end
    
    % Load STK file
    try
        [stkData,peaks,image_t] = gettraces( movieFiles(i).name, params, traceFname );
    catch e
        disp('Skipping file: corrupted, missing, or not completely saved.');
        disp(movieFiles(i).name);
        disp(e);
        existing(i) = 1;
        continue;
    end
    
    if exist('h','var'),
        waitbar(i/nFiles, h); drawnow;
    end
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
        fprintf(log_fid, 'SKIP %s\n', movieFiles(i).name);
    else
        fprintf(log_fid, '%.0f %s\n', nTraces(i)/2, movieFiles(i).name);
    end
end


% Clean up
fclose(log_fid);

end %FUNCTION batchmode





% --------------- PICK MOLECULES CALLBACKS --------------- %

%------------- Pick single molecule spots ----------------- 
function picks = getPeaks( image_t )
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


% If the threshold for detecting intensity peaks is not given, calculate it
% automatically from the std of background regions at the end of the movie.
% FIXME: This must be recalcualted after software alignment, if applicable.
if ~params.don_thresh
    if ~isfield(params,'thresh_std')
        thresh_std = constants.gettracesThresholdStd;
    else
        thresh_std = params.thresh_std;
    end

    % Calculate threshold from variance in background intensity at end of movie.
    % FIXME: this does not work well when background levels change during the
    % movie, for example with injection of Cy5-labeled tRNA.
    endBG = sort( stkData.endBackground(:) );
    endBG_lowerHalf = endBG( 1:floor(numel(endBG)*0.75) );
    params.don_thresh = thresh_std*std( endBG_lowerHalf );
else
    params.don_thresh = params.don_thresh-mean2(stkData.background);
end


% A note on notation: i, indD, indA, etc are indexes into the list of channels
% as they will appear in the output data (donor,acceptor). params.idxFields and
% quadrants identify the physical position of each channel on the camera chip.
% When looking into the image, use idxFields.
quadrants = params.idxFields;
channelNames = params.chNames;
nCh = numel(channelNames);  %# of channels TO USE.


% Define each channel's dimensions and sum fields together.
[nrow,ncol] = size(image_t);  %full-chip size.

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
% For now, we assume everything is aligned. FIXME: if an alignment file is
% loaded, there's no need to check first.
total_t = sum( allFields(:,:,quadrants), 3 );
[nrow,ncol] = size(total_t); %from now on, this is the size of subfields.



%---- 1. Pick molecules as peaks of intensity from summed (D+A) image)
[total_picks,rejected] = pickPeaks( total_t, params.don_thresh, ...
                                     params.nhoodSize, params.overlap_thresh );
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
quality = zeros(nCh,1);

if params.geometry>1 && numel(picks)>0,
    % Refine peak locations and how much they deviate. This helps determine
    % if realignment is needed. FIXME: this could be improved (?) by using
    % centroid locations for translatePeaks above?
    refinedPicks = getCentroids( image_t, picks, params.nhoodSize );
    r_mod = [ mod(refinedPicks(:,1),ncol) mod(refinedPicks(:,2),nrow) ];
    residuals = refinedPicks-picks;
    
    % For each channel, find a crude alignment using control points. This
    % helps determine if software alignment is needed.
    d = r_mod(indD:nCh:end,:);  %donor (reference) points
    donor_t = allFields( :,:, params.idxFields(indD) ); %target field to align to
    
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
        target_t = allFields(:,:,params.idxFields(i));
        quality(i) =  weberQuality(donor_t,target_t,0.7*params.don_thresh);
    end
end


%%%%% Optional software alignment algorithm (for dual-color only!)
% If the alignment is close, no need to adjust.
% Just give a warning unless asked to do software alignment in settings.
if params.geometry>1 && params.alignMethod>1 && numel(picks)>0,
    % FIXME (?): the user may expect the alignment to be applied even if
    % the deviation is small when a specific alignment is loaded!
    
    % 
    %---- 2. Transform acceptor side so that the two channels align properly.
    
    % Try out all possible alignments within a range and find the one with
    % the best donor-acceptor intensity overlap. The quality score is the
    % mean aligned peak magnitude vs random alignment. If the score is low,
    % reject it and just say "we don't know".
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
        
        elseif params.alignMethod==3,
            % Use peaks of fluorescence as control points.
            % Repeat a few times to walk towards far-off solutions.
            % FIXME: doesn't work so well aligning with actual molecules when
            % iterating multiple times. Why? Do we need multiple settings?
            a = [];
            for itr=1:3,
                a = alignSearch_cpt( donor_t, target_t, a );
            end
            newAlign(i) = a;
        
        elseif params.alignMethod==4,
            % Old, slow brute force method.
            newAlign(i) = alignSearch_weber( donor_t, target_t );
        end
        
        % Register acceptor side so that it is lined up with the donor.
        registered_t{i} = imtransform( target_t, newAlign(i).tform, ...
                           'bicubic', 'XData',[1 ncol], 'YData',[1 nrow]);
        
        total_t = total_t + registered_t{i};
        tform{i} = newAlign(i).tform;
        
        % Measure the "quality" of the alignment as the magnitude increase in
        % score compared to a "random" alignment.
        % FIXME: the threshold here is for total intensity, which may be much
        % brighter than the combination of any two channels. This could give 
        % low quality scores even when the alignment is good.
        quality(i) =  weberQuality(donor_t,registered_t{i},0.7*params.don_thresh);
    end
    
    % If the optimal alignment is not trivial, re-pick molecule locations and
    % derive alignment deviation score.
    if ~( all([newAlign.dx]==0) && all([newAlign.dy]==0) && all([newAlign.theta]==0) ),
        % Give a warning for poor quality alignment.
        if any( quality<1.1 & quality>0 ),
            warning('gettraces:lowConfidenceAlignment', ...
                    'Low confidence alignment. Parameters are out of range or data quality is poor.');
        end
        
        % Pick peaks from the aligned, total intensity image.
        [total_picks,rejected] = pickPeaks( total_t, params.don_thresh, ...
                                     params.nhoodSize, params.overlap_thresh );

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
        % The rmsd is then the distance between each channel and the donor.
        refinedPicks = getCentroids( image_t, picks, params.nhoodSize );
        residuals = refinedPicks-picks; 
        
        for i=1:nCh,
            dev = residuals(i:nCh:end,:)-residuals(indD:nCh:end,:);
            dev = dev(~rejected,:);  %remove overlapped peaks.
            newAlign(i).abs_dev = mean(  sqrt( dev(:,1).^2 + dev(:,2).^2 )  );
            newAlign(i).quality = quality(i);
        end
        align = newAlign;
    end
end

if ~params.quiet && nCh>1,
    disp('Weber contrast alignment quality score:');
    disp( quality(2:end) );
end

% Save output
stkData.rejectedTotalPicks = total_picks(rejected,:);
stkData.total_peaks = total_picks(~rejected,:);

good = repmat( ~rejected, [nCh,1] ); good=logical(good(:));
stkData.rejectedPicks = picks( ~good,: );
picks = picks( good,: );

stkData.fractionOverlapped = size(stkData.rejectedTotalPicks,1) / nPicked;

stkData.total_t = total_t;
stkData.alignStatus = align;

end %function getPeaks



function [picks,boolRejected,centroids] = pickPeaks( image_t, threshold, nhood, overlap_thresh )
% Localizes the peaks of fluorescence.
%   picks = locations (x,y) of all molecules selected.
%   rejectedPicks = locations of molecules that are too close to a neighbor
%                        and should be ignored in analysis.
%

[nrow,ncol] = size(image_t);

% Detect molecules as fluorescence maxima over local 3x3 regions,
% ignoring any that are below the detection threshold or near the edges.
% FIXME: may want to remove edges at the end, after overlap detection.
kernel = [0 1 0 ; 1 1 1; 0 1 0];

maxima = imregionalmax(image_t,kernel) & image_t>threshold;
[rows,cols] = find(maxima);
picks = [cols,rows];

ed = size(kernel,1);  %distance from the edge to avoid
edge = rows<=ed | rows>=nrow-ed | cols<=ed  | cols>=ncol-ed;
picks = picks(~edge,:);


% Find weighted center of each peak for more accurate PSF overlap detection.
% imdilate+regionprops works for this, but tends to merge nearby peaks.
% FIXME: consider using these centroid locations as picks. The only problem is
% that they msut be converted back to integers for indexing the image.
centroids = getCentroids( image_t, picks, nhood );
nMol = size(centroids,1);

% Detect maxima that are very close together and probably have overlapping
% point-spread functions.
if overlap_thresh==0 || nMol==0,
    boolRejected = false(1,nMol);  %no overlap rejection.
else
    [~,dist] = knnsearch( centroids, centroids, 'k',2 );
    boolRejected = dist(:,2)'<=overlap_thresh;
end



end %FUNCTION pickPeaks



function centroids = getCentroids( image_t, picks, nhood )
% Find weighted center of each peak for more accurate PSF overlap detection.
% imdilate+regionprops works for this, but tends to merge nearby peaks.
% FIXME: when operating on a flat (zero) signal, gives NaN values. This should
% instead just return the input, possibly with a warning.

if nargin<3,
    nhood=1;  %3x3 region.
end

nMol = size(picks,1);
centroids = zeros(nMol,2);

for i=1:nMol,
    % Extract a window region around the molecule.
    x_window = picks(i,1) + (-nhood:+nhood);
    y_window = picks(i,2) + (-nhood:+nhood);
    block = image_t( y_window, x_window );
    
    % Calculate an intensity-weighted average position of molecule w/i window.
    tot = sum(block(:));
    x = sum( x_window .* sum(block,1)/tot );
    y = sum( y_window .* sum(block,2)'/tot  );
    
    centroids(i,:) = [x y];
end

% NaN values may appear for regions that are entirely black (zero) due to
% division by zero. Replace these values with the input and give a warning.
badPicks = isnan(centroids);
if any( badPicks(:) ),
    centroids(badPicks) = picks(badPicks);
    warning('gettraces:getCentroids:NaN','NaN values found when searching for centroids (%d, %.0f%%). This can happen when a field is empty (zero). Using input molecule locations instead for these molecules.', ...
            sum(badPicks(:)), 100*sum(badPicks(:))/numel(badPicks) );
end

end %FUNCTION getCentroids




function align = alignSearch_cpt( ref_img, target_img, initAlign )
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
% A lower threshold is used to better detect dim channels (Cy7).
[~,reject,ref_peaks]    = pickPeaks( ref_img, 0.7*params.don_thresh, params.nhoodSize, params.overlap_thresh );
ref_peaks = ref_peaks(~reject,:);

[~,reject,target_peaks] = pickPeaks( target_img, 0.7*params.don_thresh, params.nhoodSize, params.overlap_thresh );
target_peaks = target_peaks(~reject,:);


% 2) For each peak in the reference, find the corresponding peak in the target.
% K-nearest neighbor search (Statistics toolbox), where K=1. Using K>1 could be
% used to verify there is no uncertainty (other points are far away).
% This will not work if the beads are too dense or if the distortion is too
% severe -- the "correct" choice has to be closer than nearby beads.
if nargin>=3 && ~isempty(initAlign),
    % Pre-align the target side if an initial guess is provided. This allows an
    % iterative search for the correct alignment if it is far off.
    temp = [target_peaks(:,1)-(ncol/2) target_peaks(:,2)-(nrow/2)];
    temp = tformfwd(initAlign.tform, temp);
    temp = [temp(:,1)+(ncol/2) temp(:,2)+(nrow/2)];
    
    [idx,dist] = knnsearch( temp, ref_peaks );
else
    [idx,dist] = knnsearch( target_peaks, ref_peaks );
end

target_peaks = target_peaks(idx,:);
sel = dist<params.nPixelsToSum;  %/2 may work better with dense, but close alignment.

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
assert( size(ref_peaks,1)>3, 'Not enough control points to create an alignment' );
tform = cp2tform( target_peaks, ref_peaks, 'nonreflective similarity' );

                      
% 4) Calculate distortion parameters from tform.
ss = tform.tdata.Tinv(2,1);
sc = tform.tdata.Tinv(1,1);
tx = tform.tdata.Tinv(3,1);
ty = tform.tdata.Tinv(3,2);
scale = sqrt(ss*ss + sc*sc);
theta = atan2(ss,sc)*180/pi;

align = struct( 'dx',tx,'dy',ty, 'theta',theta, 'sx',scale,'sy',scale, ...
                'abs_dev',0, 'tform',tform );

end %FUNCTION alignSearch_cpt




function [bestAlign,bestReg] = alignSearch_weber( donor_t, acceptor_t )
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
params.alignment.theta = -1:0.1:1;
params.alignment.dx    = -8:1:8;
params.alignment.dy    = params.alignment.dx;  %all dx and dy must be integers
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
    parfor_progress( ntheta,'Searching for the optimal alignment');
end

% Start the matlab thread pool if not already running, unless disabled.
if constants.enable_parfor,
    pool = gcp;
    M = pool.NumWorkers;
else
    M = 0;
end

% Reserve space for the best alignment for each possible rotation value.
% This is required for parfor to work correctly because all threads must be
% totally independent.
%scores = zeros(ntheta,ndx,ndy);
bestScores = zeros(ntheta,1);
bestAligns = cell(ntheta,1);
bestRegs = cell(ntheta,1);

parfor (t=1:ntheta,M)
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
            % FIXME: biased against moving too far because edges are zero?
            registered = zeros( size(acceptor_t) );

            registered( max(1,1-dy):min(nrow,nrow-dy), max(1,1-dx):min(ncol,ncol-dx) ) = ...
                rot_a( max(1,1+dy):min(nrow,nrow+dy), max(1,1+dx):min(ncol,ncol+dx) );

            % Calculate Weber contrast as a score. When the fields are not
            % aligned, the peak intensities go down and some peaks fall
            % below the threshold of background, lowering contrast.
            total = donor_t+registered;
            picks = total>params.don_thresh;    %pixels above background
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
        %waitbar( t/ntheta, h );
        parfor_progress();
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


if ~params.quiet && ntheta>1,
    parfor_progress(0);
end

end %FUNCTION alignSearch



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

N = 3; %number of repititions for averaging.
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


end %FUNCTION weberRandom




% --------------------- SAVE PICKED TRACES TO FILE --------------------- %

function getIntegrationWindows( stk_top, peaks )
% For each molecule location in "peaks", find the most intense pixels in
% its immediate neighborhood (defined by params.nPixelsToSum). These
% regions are used by integrateAndSave() to sum most of the intensity for
% each peak and generate fluorescence-time traces.
% To minimize the contribution of nearby molecules, the molecules closest to the
% peak center are added first and progressively out to the edge.
%
% FIXME: instead of just returning the N best pixels, return the entire sorted
% neighborhood and choose the highest N pixels later. This would give more
% flexibility and information?

hw = params.nhoodSize;  % distance from peak to consider (eg, 1=3x3 area)
squarewidth = 1+2*hw;   % width of neighborhood to examine.


% Get x,y coordinates of picked peaks
Npeaks = size(peaks,1);
x = peaks(:,1);
y = peaks(:,2);


% Define regions over which to integrate each peak --
% Done separately for each channel!
stkData.integrationEfficiency = zeros(Npeaks,squarewidth^2);
regions = zeros(params.nPixelsToSum,2,Npeaks);  %pixel#, dimension(x,y), peak#

imgReused = zeros( size(stk_top) );  %marks where pixels are re-used
idxs = zeros( params.nPixelsToSum, Npeaks );
stkData.fractionWinOverlap = zeros(Npeaks,1);

for m=1:Npeaks
    % Get a window of pixels around the intensity maximum (peak).
    nhood = stk_top( y(m)-hw:y(m)+hw, x(m)-hw:x(m)+hw );
    center = sort( nhood(:), 'descend' );
    
    % Find the most intense pixels, starting from the center and moving out.
    % This reduces the chance of getting intensity from nearby molecules.
    [A,B] = find( nhood>=center(params.nPixelsToSum), params.nPixelsToSum );

    % Estimate the fraction of intensity in each pixel,
    % relative to the total intensity in the full window region around the peak.
    % This is just an estimate and depends on the window size.
    % High molecule density can also distort this if there are overlapping PSFs.
    stkData.integrationEfficiency(m,:) = cumsum( center/sum(center) )';
    
    % Convert to coordinates in the full FOV image and save.
    regions(:,:,m) = [ A+y(m)-hw-1, B+x(m)-hw-1  ];
    
    % Note where pixels are being reused for later calculation.
    idxs(:,m) = sub2ind( size(stk_top), regions(:,1,m), regions(:,2,m) );
    imgReused(idxs(:,m)) = imgReused(idxs(:,m)) +1;
end
    
% If this entire neighborhood is empty, integrationEfficiency will be NaN.
% Give a warning, but leave the NaN.
badWindows = isnan(stkData.integrationEfficiency);
badWindows = any(badWindows'); %select if any entry per molecule is NaN.

if any(badWindows),
    warning('gettraces:getIntegrationWindows:NaN','NaN values found when getting integration windows (%d, %.0f%%). This can happen when a field is empty (zero).', ...
        sum(badWindows), 100*sum(badWindows)/Npeaks );
end
    
% For each peak, get the fraction of re-used pixels.
% FIXME: idxs could be used to construct regions to save time/memory.
for m=1:Npeaks,
    stkData.fractionWinOverlap(m) = sum(imgReused(idxs(:,m))-1) / params.nPixelsToSum;
end

stkData.regions = regions;


end %function getIntegrationWindows




function integrateAndSave( peaks, stk_fname )
% For each location in "peaks", sum a region that includes most of the
% intensity for that spot and repeat for each time point to get a
% fluorescence-time trace for each peak. Background subtraction, crosstalk
% correction, and calculation of derived signals (FRET traces) is all done
% here. Then the result is saved as a .rawtraces file with metadata.
%
% tic;

movie = stkData.movie;
nFrames = movie.nFrames;

% Start the progress bar before initial setup; indicate something is happening.
quiet = params.quiet;
if ~quiet,
    parfor_progress( nFrames/10,'Extracting traces from movie data');
end

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
if ~quiet && data.time(1)==1,
    disp('Time axis information is not present in movie!');
    a = inputdlg('Time resolution (ms):','No time axis in movie!');
    a = str2double(a);
    if ~isempty(a), %empty if user hit Cancel or entered an invalid number
        data.time = a*( 0:nFrames-1 );
    else
        disp('Using frames as the time axis');
    end
end

% Parallelize very large TIFF movies, where disk access is quick compared to
% image processing. 
if nTraces*nFrames/2000 > 1500 && constants.enable_parfor && isa(movie,'Movie_TIFF'),
    % Processing large TIFF movies is CPU limited. Use parfor to parallelize.
    pool = gcp;
    M = pool.NumWorkers;
else
    % For small datasets, do everything in the GUI thread (regular for loop).
    M = 0;
end

% Create a trace for each molecule across the entire movie.
% The estimated background image is also subtracted to help with molecules
% that do not photobleach during the movie.
traces = zeros(Npeaks,nFrames,'single');

idx = sub2ind( [movie.nY movie.nX], regions(:,1,:), regions(:,2,:) );
bg = single(stkData.background);
nPx = params.nPixelsToSum;

parfor (k=1:nFrames, M)
    % NOTE: 25% faster by converting to int16, with no change to sCMOS data.
    % But EMCCD have slight differences due to 15-bit overflows?
    frame = single(movie.readFrame(k)) - bg;
    
    if nPx>1,
        traces(:,k) = sum( frame(idx) );
    else
        traces(:,k) = frame(idx);
    end
    
    % Update waitbar. Using mod speeds up the loop, but isn't ideal because
    % indexes are executed somewhat randomly. Reasonably accurate despite this.
    if mod(k,10)==0 && ~quiet,
        parfor_progress();
    end
end
if ~quiet,
    parfor_progress('Correcting traces and calculating FRET...');
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

% If this is multi-color FRET data and the metadata doesn't specify how FRET
% should be calculated, ask the user for clarification.
% NOTE: this is a temporary field and will be replaced by something more general
% in a future version.
if data.isChannel('acceptor2'),
    result = questdlg('Can you assume there is no donor->acceptor2 FRET?', ...
                    '3-color FRET calculation','Yes','No','Cancel','No');
    if strcmp(result,'Cancel'),  return;  end
    data.fileMetadata(1).isTandem3 = double( strcmp(result,'Yes') );
end

% Subtract background and calculate FRET
data = correctTraces(data);


% ---- Metadata: save various metadata parameters from movie here.
parfor_progress('Saving traces...');

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

if ~quiet,
    parfor_progress(0);
end
% disp(toc);

end %function integrateAndSave




end %FUNCTION gettraces




