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
%                 that field is disregarded. OPTIONAL

%   Copyright 2007-2015 Cornell University All Rights Reserved.


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
[stkData.regions, stkData.integrationEfficiency, stkData.fractionWinOverlap] = ...
                                getIntegrationWindows( image_t, peaks, params );


%------ Extract traces using peak locations
if nargin>=3 && ischar(varargin{3}),
    outputFilename = varargin{3};
    integrateAndSave( peaks, outputFilename );
end
    
return;









%% --------------------- OPEN STK MOVIE (CALLBACK) --------------------- %

function [stkData] = OpenStk(filename)

% Load movie data from file.
movie = Movie.load(filename);
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
log_fid = fopen( fullfile(direct,'gettraces.log'), 'wt' );

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
align = struct('dx',{},'dy',{},'theta',{},'sx',{},'sy',{},'abs_dev',{},'quality',{});
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

        tform = fitgeotrans( d(~rejected,:), a(~rejected,:), 'NonreflectiveSimilarity' );
        ss = tform.T(2,1);
        sc = tform.T(1,1);
        scale = sqrt(ss*ss + sc*sc);
        theta = atan2(ss,sc)*180/pi;

        % Measure the "quality" of the alignment as the magnitude increase in
        % score compared to a "random" alignment.
        target_t = allFields(:,:,params.idxFields(i));
        quality(i) =  weberQuality(donor_t,target_t,0.7*params.don_thresh);

        align(i) = struct( 'dx',tform.T(3,1), 'dy',tform.T(3,2), 'theta',theta, ...
                'sx',scale, 'sy',scale, 'abs_dev',abs_dev, 'quality',quality(i) );
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
        % tform moves the acceptor field to be aligned with the donor.
        if params.alignMethod==2,
            % Nothing to search, just apply the alignment.
            newAlign(i) = params.alignment(i);
        
        elseif params.alignMethod==3,
            % Iterative closest point algorithm.
            newAlign(i) = icpalign( donor_t, target_t, params );
        
        elseif params.alignMethod==4,
            % Old, slow brute force method.
            error('weberalign method no longer available');
            %newAlign(i) = weberalign( donor_t, target_t, params );
        end
        
        % Register acceptor side so that it is lined up with the donor.
        % imref2d specifies the center of the image is the origin (0,0).
        R = imref2d( size(target_t), [-1 1]*nrow/2, [-1 1]*ncol/2 );
        registered_t{i} = imwarp( target_t,R, newAlign(i).tform,...
                                    'Interp','cubic', 'OutputView',R );
        
        total_t = total_t + registered_t{i};
        tform{i} = newAlign(i).tform;
        
        % Measure the "quality" of the alignment as the magnitude increase in
        % score compared to a "random" alignment.
        % FIXME: the threshold here is for total intensity, which may be much
        % brighter than the combination of any two channels. This could give 
        % low quality scores even when the alignment is good.
        quality(i) = weberQuality(donor_t,registered_t{i},0.7*params.don_thresh);
    end
    
    % If the optimal alignment is not trivial, re-pick molecule locations and
    % derive alignment deviation score.
    if ~( all([newAlign.dx]==0) && all([newAlign.dy]==0) && all([newAlign.theta]==0) ),
        % Give a warning for poor quality alignment.
        if any( quality<1.1 & quality>0 ),
            disp('Gettraces: Low confidence alignment. Parameters are out of range or data quality is poor.');
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




% --------------------- SAVE PICKED TRACES TO FILE --------------------- %
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
    wbh = parfor_progressbar(1.1*nFrames,'Extracting traces from movie data');
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

if isfield(params,'zeroMethod'),
    data.fileMetadata(1).zeroMethod = params.zeroMethod;
end

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
        wbh.iterate(10);
    end
end
if ~quiet,
    wbh.message = 'Correcting traces and calculating FRET...';
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
    
    % Spectral crosstalk correction.
    if numel(params.crosstalk)==1,
        data.acceptor = data.acceptor - params.crosstalk*data.donor;
    else
        % The order of operations matters here.
        nFluor = numel(data.idxFluor);
        for src=1:nFluor,
            for dst=1:nFluor,
                if src>=dst, continue; end  %only consider forward crosstalk
                ch1 = data.channelNames{src};
                ch2 = data.channelNames{dst};
                crosstalk = params.crosstalk(src,dst);
                data.(ch2) = data.(ch2) - crosstalk*data.(ch1);
            end
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

% Subtract background
data = correctTraces(data);

% Scale acceptor channel to correct for unequal brightness (gamma is not 1).
% Highly scaled (dim) channels can confuse the background subtraction method,
% so this must be done after background subtraction.
if isfield(params,'scaleAcceptor') && ~isempty(params.scaleAcceptor),
    data.fileMetadata(1).scaleAcceptor = params.scaleAcceptor;
    
    data.acceptor = data.acceptor*params.scaleAcceptor(1);
    for i=2:numel(params.scaleAcceptor),
        name = sprintf('acceptor%d',i);
        data.(name) = data.(name)*params.scaleAcceptor(i);
    end
end
data.recalculateFret();


% ---- Metadata: save various metadata parameters from movie here.
if ~quiet,
    wbh.iterate(0.1*nFrames);
    wbh.message = 'Saving traces...';
end

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
    close(wbh);
end
% disp(toc);

end %function integrateAndSave




end %FUNCTION gettraces




