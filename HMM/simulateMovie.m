function peakLocations = simulateMovie(data,bgMovieFilename,outFilename,params)
%% simulateMovie   Simulate wide-field smFRET movies from fluorescence data
% 
%    simulateMovie( traceData, bgMovieFilename, outputFilename, params )
%    Generates wide-field smFRET movies, with donor fluorescence intensity
%    on the left half of the field and acceptor fluorescence on the right half.
%    Experimental movies can be used as to approximate background noise.
%    Fluorescence intensities are spread over symmetric 2D Gaussian functions
%    and placed at random positions within the field-of-view.
%
%    simulateMovie( tracesFilename, ... )
%    Alternatively, a .traces file with fluorescence data can be given.
%
%    Recognized parameters:
%    - sigmaPSF: standard deviation of symmetric 2D Gaussian point-spread
%                functions used to distribute fluorescence in movies
%    - density: number of molecules placed per field.
%    - grid: place fluorophores on a regular grid if =1.
%    - randomSeed: seed for pseudorandom number generator.
%    - edAlignX,edAlignY,edAlignTheta: displacement of the acceptor
%                relative to donor (in pixels) to simulate misalingment.
%

% TODO: choose parameter values like buffer size and the region over which
% Gaussians are distributed based on the size of the point-spread function.


%% Parse input parameters

% No parameters given; nothing to do but exit.
assert( nargin>=3, 'Invalid parameter list');

% If a traces filename is given instead of the data load it.
if ischar(data)
    data = loadTraces(data);
end

if iscell(bgMovieFilename),
    bgMovieFilename = bgMovieFilename{1};
end
if ~( ischar(bgMovieFilename) && exist(bgMovieFilename,'file') )
    error('Invalid background movie filename');
end


%% Assign default values to parameters not specified.
% FIXME: if anything, these should be in cascadeConstants.
defaultParams.sigmaPSF = 0.9;
defaultParams.randomSeed = 0;
defaultParams.density = 5500;
defaultParams.alignX = 0;
defaultParams.alignY = 0;
defaultParams.alignTheta = 0;

% Parameters actually specified (second argument) take precedence
if nargin>=4,
    params = catstruct( defaultParams, params );
else
    params = defaultParams;
end

sigmaPSF = params.sigmaPSF;
if ~isfield(params,'edgeBuffer') || isempty(params.edgeBuffer),
    params.edgeBuffer = 4+ceil(sigmaPSF*4);
end



%% Simulate Wide-field Fluorescence Movies.
tic;

% --- 1. Load background movie as initial image to which peaks are added.
movie = Movie.load( bgMovieFilename );
stkX = movie.nX;
stkY = movie.nY;
stkNFrames = min(movie.nFrames,data.nFrames);

if data.nFrames>movie.nFrames,
    disp('Truncating traces to fit in the movie.');
end

assert( mod(stkX,2)==0, 'Movie widths must be even!' );



% --- 2. Select fluorophore centroid positions randomly.
border = 20;%ceil(8*sigmaPSF);%params.edgeBuffer;

% If user requests a regular grid of positions, add some small
% variability (up to 1 pixel) to approximate random placement.
if isfield(params,'grid') && ~isempty(params.grid) && params.grid,
    % Assign a minimal spacing between peaks (3 std's on each side).
    minSpacing = ceil(8*sigmaPSF);

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
% FIXME: add scale!
dx = params.alignX; dy = params.alignY; % translation factors
sx = 1; sy = 1; % scaling (magnification) factors
theta = params.alignTheta*pi/180; % degrees of rotation about center.
T = [ sx*cos(theta)     sin(theta)  0 ; ...
        -sin(theta)  sy*cos(theta)  0 ; ...
               dx             dy    1 ];
tform = maketform('affine',T);

% Tranform acceptor side positions accordingly.
% Assumes the field image is symmetric. The image is shifted so the
% center is 0 so that the rotation happens about the origin.
acceptorPos(:,[2,1]) = tformfwd( tform, donorPos(:,[2,1])-(stkY/2) ) + (stkY/2);


% Save peak locations for later comparison (Dy,Dx,Ay,Ax).
peakLocations = [donorPos acceptorPos];  %return parameter
% figure;
% scatter( donorPos(:,2), donorPos(:,1), 'bo' ); hold on;
% scatter( acceptorPos(:,2), acceptorPos(:,1), 'ro' );
% title('Simulated molecule locations');


% FIXME: here, combine donor/acceptor points and data so that they can
% be indexed as one big array to simplify the code below.
acceptorPos(:,2) = acceptorPos(:,2)+stkX/2;  %translate

peaks = zeros( size(donorPos,1)*2, 2 );
peaks(1:2:end,:) = donorPos;
peaks(2:2:end,:) = acceptorPos;
nPeaks = size(peaks,1);


% Crop data so we only keep the traces we care about. FIXME
fluor = zeros( nPeaks, stkNFrames );
fluor(1:2:end,:) = data.donor(1:nPeaks/2,1:stkNFrames);
fluor(2:2:end,:) = data.acceptor(1:nPeaks/2,1:stkNFrames);

% --- 3. Estimate Gaussian PDF for each peak.
% FIXME: create a grid with 10x the density and sum to better
% estimate the PSF. -limit:0.1:limit...
limit = ceil(4*sigmaPSF);  %half-width of window
nw = 2*limit+1;
psfs = zeros( nw,nw, nPeaks );

% Movie coordinates for each windows origin (0,0), which mol's center
origins = round(peaks);

% Peak locations relative to window centers.
cy = peaks(:,1)-origins(:,1);
cx = peaks(:,2)-origins(:,2);

[y,x] = meshgrid(-limit:limit,-limit:limit);

for i=1:nPeaks,
    % Estimate the Gaussian PDF over the pixel grid.
    psfs(:,:,i) = exp( -0.5* ((x-cx(i)).^2+(y-cy(i)).^2)/(sigmaPSF^2) );
end

% 2D Gaussian normalization factor
psfs = psfs/(2*pi* sigmaPSF^2);


% --- 4. Distribute fluorescence intensities into each PSF
% Create TIFF tag structure. FIXME: insure these match the movie.
constants = cascadeConstants();
tags.ImageLength = movie.nY;
tags.ImageWidth = movie.nX;
tags.Photometric = Tiff.Photometric.MinIsBlack;
tags.BitsPerSample = 16;
tags.SamplesPerPixel = 1;
tags.RowsPerStrip = movie.nY;  %one strip per image.
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tags.Software = ['simulateMovie (ver. movie' constants.version ')'];
%tags.ImageDescription = ...  %simulation settings, etc here?

% Generate DateTime tags to get exposure time, etc.
% FIXME: Tiff class doesn't seem to allow the ms part (.FFF).
% tags.ExporeTime = data.sampling/1000;  %Tiff class doesn't allow custom tags?
% times = now+data.time/24/60/60/1000; %fractions of a day
% time_str = datestr(times,'dd-mmm-yy HH:MM:SS.FFF');

% Create the TIFF file.
hTiff = Tiff(outFilename,'w8');  % w=tiff, w8=bigtiff

chunkSize = 100;  %number of frames to process at once.
wbh = waitbar(0,'Distributing fluorescence information...'); 

for k=1:chunkSize:stkNFrames,
    % Load the background movie data
    fidx = k:min(k+chunkSize-1,stkNFrames);
    frames = movie.readFrames(fidx);

    % For each peak, add fluorescence onto the frames
    f = reshape( fluor(:,fidx), [1 nPeaks chunkSize] );
    
    for i=1:nPeaks,
        yl = origins(i,1)+(-limit:limit);
        xl = origins(i,2)+(-limit:limit);
        frames(yl,xl,:) = frames(yl,xl,:) + uint16( ...
                                     bsxfun( @times, psfs(:,:,i), f(1,i,:) )  );
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
    waitbar(k/stkNFrames,wbh);
end

close(hTiff);
close(wbh);
disp(toc);


end
