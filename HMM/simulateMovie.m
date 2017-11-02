function output = simulateMovie(data, bgMovieFilename, outFilename, inputParams)
%% simulateMovie   Simulate wide-field smFRET movies from fluorescence data
% 
%    simulateMovie( traceData, bgMovieFilename, outputFilename, params )
%    Generates wide-field smFRET movies, with donor fluorescence intensity
%    on the left half of the field and acceptor fluorescence on the right half.
%    Experimental movies can be used as to approximate background noise.
%    Fluorescence intensities are spread over symmetric 2D Gaussian functions
%    and placed at random positions within the field-of-view.
%    
%    TRACES can be a Traces object or the path to a .traces object.
%    If any filenames are not given, the user will be prompted for them.
%
%    Recognized parameters:
%    - sigmaPSF: standard deviation of symmetric 2D Gaussian point-spread
%                   functions used to distribute fluorescence in movies
%    - density:  number of molecules placed per field.
%    - grid:     place fluorophores on a regular grid if true.
%    - alignX,alignY,alignTheta,alignScale: displacement of the acceptor 
%          relative to the donor (in pixels/degrees) to simulate misalingment.
%    - edgeBuffer: number of pixels around the edges to avoid placing PSFs.
%    - aduPhoton:  ADU units per photon conversion in the output movie.
%
% See also: simulate, batchKinetics, Movie, Tiff.

%   Copyright 2007-2017 Cornell University All Rights Reserved.

% TODO: choose parameter values like buffer size and the region over which
% Gaussians are distributed based on the size of the point-spread function.
% Also, re-implement random seed?


%% Parse input parameters
narginchk(0,4)

% Ask for input traces file if necessary.
if nargin<1 || isempty(data)
    data = loadTraces;
    if isempty(data), return; end  %user hit cancel
    
elseif ischar(data)
    data = loadTraces(data);
end

% Prompt user for background movie file if necessary.
if nargin<2 || isempty(bgMovieFilename),
    bgMovieFilename = getFile('*.tif;*.tiff;*.stk', 'Select background movie');
    if isempty(bgMovieFilename), return; end  %user hit cancel

elseif iscell(bgMovieFilename),
    bgMovieFilename = bgMovieFilename{1};
end
if ~exist(bgMovieFilename,'file')
    error('Invalid background movie filename');
end

% Prompt user for output movie file name.
if nargin<3 || isempty(outFilename)
    [f,p] = uiputfile('*.tif','Save movie');
    if isequal(p,0), return; end  %user hit cancel
    outFilename = fullfile(p,f);
end


%% Assign default values to parameters not specified.
% FIXME: should be in cascadeConstants.
% FIXME: these are specific to sCMOS.
params.sigmaPSF   = 0.8;
params.density    = 5500;
params.alignX     = 0;
params.alignY     = 0;
params.alignTheta = 0;
params.alignScale = 1;
params.aduPhoton  = 2;
params.edgeBuffer = 15;  %4+ceil(params.sigmaPSF*8);
params.grid       = false;

% Parameters actually specified (second argument) take precedence
if nargin>=4,
    params = mergestruct( params, inputParams );
end


constants = cascadeConstants;
if data.nTraces > 2000 && constants.enable_parfor,
    % Processing large TIFF movies is CPU limited. Use parfor to parallelize.
    pool = gcp;
    M = pool.NumWorkers;
else
    % For small datasets, do everything in the GUI thread (regular for loop).
    M = 0;
end


%% Simulate Wide-field Fluorescence Movies.

% --- 1. Load background movie as initial image to which peaks are added.
movie = Movie.load( bgMovieFilename );
stkX = movie.nX;
stkY = movie.nY;
stkNFrames = min(movie.nFrames,data.nFrames);

if data.nFrames>movie.nFrames,
    disp('Truncating traces to fit in the movie.');
elseif data.nTraces<params.density,
    warning('Not enough traces to fill the movie to desired density');
    params.density = data.nTraces;
end

assert( mod(stkX,2)==0, 'Movie widths must be even!' );



% --- 2. Select fluorophore centroid positions randomly.
border = params.edgeBuffer;

% If user requests a regular grid of positions, add some small
% variability (up to 1 pixel) to approximate random placement.
if params.grid,
    % Assign a minimal spacing between peaks (3 std's on each side).
    minSpacing = ceil(8*params.sigmaPSF);

    % Calculate a spacing that best approximates the requested density.
    % If the spacing is too small to avoid overlap, reassign.
    idealSpacing = sqrt(  ( (stkX/2-2*border)*(stkY-2*border) )/params.density  )-1;
    spacing = max(minSpacing, floor(idealSpacing));

    % Assign mean fluorophore positions.
    x = (1+border):spacing:(stkX/2-border);  nX = numel(x);
    y = (1+border):spacing:(stkY-border);    nY = numel(y);
    x = repmat( x, [nY 1] );
    y = repmat( y, [nX 1] )';

    % Position each fluorophore at a random position within the
    % central pixel to simulate random placement.
    donorPos = [y(:) x(:)];
    donorPos = donorPos + rand(size(donorPos))-0.5;

% Otherwise, place fluorophores at random positions. The edges of the
% field are avoided b/c they will not be processed by gettraces.
% NOTE: some PSFs will overlap, leading to artifacts when fluorescence
% traces are generated.
else
    donorPos = border + [ rand(params.density,1)*(stkY  -2*border) ...
                          rand(params.density,1)*(stkX/2-2*border) ];
end

% Construct a transformation to mimic field misalignment.
sc = params.alignScale * cosd(params.alignTheta);
ss = params.alignScale * sind(params.alignTheta);
T = [ sc -ss  0;
      ss  sc  0;
      params.alignX  params.alignY  1];
tform = affine2d(T);

% Tranform acceptor side positions accordingly.
% Assumes the field image is symmetric. The image is shifted so the
% center is 0 so that the rotation happens about the origin.
acceptorPos(:,[2,1]) = transformPointsForward( tform, ...
                                    donorPos(:,[2,1])-(stkY/2) ) + (stkY/2);

% Save peak locations for later comparison (Dy,Dx,Ay,Ax).
if nargout>0, output=[donorPos acceptorPos]; end
% figure;
% scatter( donorPos(:,2), donorPos(:,1), 'bo' ); hold on;
% scatter( acceptorPos(:,2), acceptorPos(:,1), 'ro' );
% title('Simulated molecule locations');

% Translate acceptor to neighboring field (camera) and combining locations and
% fluorescence into an interleaved list for each processing.
acceptorPos(:,2) = acceptorPos(:,2)+stkX/2;  %translate
peaks = [donorPos; acceptorPos];
nPeaks = size(peaks,1);
fluor = [data.donor(1:nPeaks/2,1:stkNFrames) ; data.acceptor(1:nPeaks/2,1:stkNFrames)];

% --- 3. Estimate Gaussian PDF for each peak.
wbh = waitbar(0,'Estimating point-spread functions...'); 

limit = ceil(4*params.sigmaPSF);  %half-width of window
bins = -limit:limit;  %Pixel edges within window for evaluating 2D Gaussians
nw = 2*limit+1;
psfs = zeros( nw,nw, nPeaks );

% Coordinates of pixel closest to actual molecule position.
origins = round(peaks);

% Peak locations relative to window centers.
cy = peaks(:,1)-origins(:,1);
cx = peaks(:,2)-origins(:,2);

% Crude but fast approximation to the area under the Gaussian in each pixel.
% [y,x] = meshgrid(-limit:limit,-limit:limit);
% for i=1:nPeaks,
%     % Estimate the Gaussian PDF over the pixel grid.
%     psfs(:,:,i) = exp( -0.5* ((x-cx(i)).^2+(y-cy(i)).^2)/(params.sigmaPSF^2) );
% end

% Use a Riemann integral to get the area under the 2D Gaussian PDF for each
% pixel. x,y define the partitions in the subpixel array for the sum.
delta = 0.1;  %sub-pixel integration step size
[y,x] = meshgrid( (0:delta:0.9)+delta/2, (0:delta:0.9)+delta/2 );  %subpixel bin centers
prefactor = -0.5/(params.sigmaPSF^2);

parfor (i=1:nPeaks,M)
    for xx=1:nw,
        dx2 = (bins(xx)+x-cx(i)).^2;   %#ok<PFBNS>
        for yy=1:nw,
            % 2D Gaussian evaluated at the partition points
            dy2 = (bins(yy)+y-cy(i)).^2;
            val = exp( prefactor*(dx2+dy2) );

            % Riemann sum to approximate the integral.
            psfs(yy,xx,i) = sum( val(:) );
        end
    end
end

% Normalization
psfs = psfs*(delta^2)/(2*pi* params.sigmaPSF^2);
% disp( mean(psfs,3) );

% --- 4. Distribute fluorescence intensities into each PSF
waitbar(0.15,wbh,'Distributing fluorescence information...');

% Create TIFF tag structure. FIXME: insure these match the movie.
tags.ImageLength = movie.nY;
tags.ImageWidth = movie.nX;
tags.Photometric = Tiff.Photometric.MinIsBlack;
tags.BitsPerSample = 16;
tags.SamplesPerPixel = 1;
tags.RowsPerStrip = movie.nY;  %one strip per image.
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tags.Software = [mfilename ' (' constants.software ')'];
%tags.ImageDescription = ...  %simulation settings, etc here?

% Generate DateTime tags to get exposure time, etc.
% FIXME: Tiff class doesn't seem to allow the ms part (.FFF).
% tags.ExporeTime = data.sampling/1000;  %Tiff class doesn't allow custom tags?
% times = now+data.time/24/60/60/1000; %fractions of a day
% time_str = datestr(times,'dd-mmm-yy HH:MM:SS.FFF');

% Create the TIFF file.
hTiff = Tiff(outFilename,'w8');  % w=tiff, w8=bigtiff

tic;

chunkSize = 1;  %number of frames to process at once.
pxrange = (-limit:limit);
selAll = data.total>50;

for k=1:chunkSize:stkNFrames,
    % Load the background movie data
    fidx = k:min(k+chunkSize-1,stkNFrames);  %frame number in whole movie.
    frames = movie.readFrames(fidx);
    % or gaussian noise background frame if no bg movie given.
    
    % Only add noise to regions where the donor is alive.
    % Helps speed up the poissrnd step, which is very slow.
    sel = selAll(:,k);
    
    % For each peak, add fluorescence onto the frames.
    % bsxfun multiplies the fluorescence (f) into each pixel of the PSF (psfs).
    % poissrnd simulates shot noise from these values (must be in photons!)
    fluorpsfs = bsxfun( @times, psfs, reshape(fluor(:,fidx),[1 1 nPeaks]) );
    fluorpsfs(sel) = poissrnd(fluorpsfs(sel));
    fluorpsfs = uint16( params.aduPhoton*fluorpsfs );
    
    for i=1:nPeaks,
        yl = origins(i,1)+pxrange;
        xl = origins(i,2)+pxrange;
        frames(yl,xl) = frames(yl,xl) + fluorpsfs(:,:,i);
    end

    % Write the finished frame to disk.
    for i=1:numel(fidx),
        %tags.DateTime = time_str(k,:);
        if k>1 || i>1,
            writeDirectory(hTiff);
        end
        setTag(hTiff, tags);
        write(hTiff, frames(:,:,i));
    end

    % Update waitbar
    if mod(k,10)==0
        waitbar(0.15+0.85*(k/stkNFrames),wbh);
    end
end

close(hTiff);
close(wbh);
fprintf('Saved simulated movie to %s.\n', bgMovieFilename);
disp(toc);

end


