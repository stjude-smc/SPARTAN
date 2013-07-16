function peakLocations = simulateMovie(varargin)
%% simulateMovie   Simulate wide-field smFRET movies from fluorescence data
% 
%    simulateMovie( donor, acceptor, bgMovieFilenames, params )
%    Generates wide-field smFRET movies, with donor fluorescence intensity
%    on the left half of the field and acceptor fluorescence on the right half.
%    Experimental movies can be used as to approximate background noise.
%    Fluorescence intensities are spread over symmetric 2D Gaussian functions
%    and placed at random positions within the field-of-view.
%
%    simulateMovie( tracesFilename, ... )
%    Use donor and acceptor fluorescence traces from the file specified in
%    tracesFilename. WARNING: fluorescence data should not have background
%    noise already!
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

nargs = numel(varargin);

% No parameters given; nothing to do but exit.
if nargin<3 || nargin>4,
    error('Invalid parameter list');
end

% First parameter is a filename; load traces data from file.
if ischar(varargin{1})
    tracesFilename = varargin{1};
    data = loadTraces( tracesFilename );
    
    bgMovieFilenames = varargin{2};
    params = varargin{3};
    
% First parameter is a matrix of numbers, specifying donor/acceptor traces.
else
    tracesFilename = [pwd filesep 'sim.traces'];
    data.donor    = varargin{1};
    data.acceptor = varargin{2};
    data.time = 1:size(data.donor,2);
    
    bgMovieFilenames = varargin{3};
    params = varargin{4};
end

% If only a single background movie filename is given, convert it to cell array.
if ischar(bgMovieFilenames),
    bgMovieFilenames = {bgMovieFilenames};
end


%% Assign default values to parameters not specified.
defaultParams.sigmaPSF = 0.8;
defaultParams.randomSeed = 0;
defaultParams.density = 300;
defaultParams.useAllMovies = 0;
defaultParams.alignX = 0;
defaultParams.alignY = 0;
defaultParams.alignTheta = 0;

% Parameters actually specified (second argument) take precedence
params = catstruct( defaultParams, params );

sigmaPSF = params.sigmaPSF;
if ~isfield(params,'edgeBuffer') || isempty(params.edgeBuffer),
    params.edgeBuffer = 4+ceil(sigmaPSF*4);
end


 
% Initialize random number generator. Insure values are not correlated
% if traces are generated with the same random seed.
rand( 'twister', params.randomSeed+1299 );



%% Alter fluorescence intensities to simulate spectral artifacts.

% Simulate the effect of donor to acceptor channel crosstalk
constants = cascadeConstants;
data.acceptor = data.acceptor + constants.crosstalk*data.donor;



%% Simulate Wide-field Fluorescence Movies.

[nTracesTotal,traceLen] = size(data.donor);
nTracesUsed = 0;
peakLocations = cell(0,1);

h = waitbar(0,'Distributing fluorescence information...'); tic;

for n=1:numel(bgMovieFilenames)

    % --- 1. Load background movie as initial image to which peaks are added.
    movie = Movie_STK( bgMovieFilenames{n} );
    assert( traceLen<=movie.nFrames, 'Background movie is too short!' );
    stk = movie.readFrames(1,traceLen);
    
    stkOutput = double( stk );  %use background movies
    [stkY,stkX,stkNFrames] = size(stkOutput);
    assert( mod(stkX,2)==0, 'Movie widths must be even!' );

    
    % --- 2. Select fluorophore centroid positions randomly.
    border = params.edgeBuffer;
    
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
    
    % Verify there are enough traces to fully populate the movie
    nPeaks = size(donorPos,1);
    if ~params.useAllMovies && nPeaks>(nTracesTotal-nTracesUsed)
        disp('Not enough traces to fill this movie');
        break;
    end
    
    
    % Construct a transformation to mimic field misalignment.
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
    peakLocations = [donorPos acceptorPos];
    figure;
    scatter( donorPos(:,2), donorPos(:,1), 'bo' ); hold on;
    scatter( acceptorPos(:,2), acceptorPos(:,1), 'ro' );
    title('Simulated molecule locations');
    

    % --- 3. Distribute fluorescence intensities into a point-spread function
    %        around each centroid.            
    
    % 2D Gaussian normalization factor
    A = 1/(2*pi* sigmaPSF^2);
    limit = ceil(4*sigmaPSF);
    
    % Iterate over peaks.
    for i=1:nPeaks,
        traceID = mod(i+nTracesUsed-1,nTracesTotal)+1;
        d = data.donor(traceID,:);
        a = data.acceptor(traceID,:);
        
        cy = donorPos(i,1);
        cx = donorPos(i,2);
        cya = acceptorPos(i,1);
        cxa = acceptorPos(i,2);
        
        % Define point-spread function
        for x=max(1,floor(cx)-limit):min(stkX/2,ceil(cx)+limit)
            for y=max(1,floor(cy)-limit):min(stkY,ceil(cy)+limit)
                % Donor fluorophore
                pixelDonor    = A.*exp(-0.5* ((x-cx)^2+(y-cy)^2)/(sigmaPSF^2)  ).*d;
                stkOutput(y,x,:) = stkOutput(y,x,:) + reshape(pixelDonor,[1 size(pixelDonor)]);
            end
        end
                
        for x=max(1,floor(cxa)-limit):min(stkX/2,ceil(cxa)+limit)
            for y=max(1,floor(cya)-limit):min(stkY,ceil(cya)+limit)
                % Acceptor fluorophore
                pixelAcceptor = A.*exp(-0.5* ((x-cxa)^2+(y-cya)^2)/(sigmaPSF^2)  ).*a;
                stkOutput(y,x+stkX/2,:) = stkOutput(y,x+stkX/2,:) + reshape(pixelAcceptor,[1 size(pixelAcceptor)]);
            end
        end

        if ~params.useAllMovies
            waitbar( (0.95*i+nTracesUsed)/nTracesTotal, h );
        else
            waitbar( n/numel(bgMovieFilenames), h );
        end
    end

    % --- 4. Save newly generated movie to file
    [p,name] = fileparts(tracesFilename);
    stkFilename = [p filesep strrep(name,'.traces','') num2str(n) '.tiff'];
    
    imwrite( uint16(stkOutput(:,:,1)), stkFilename, 'Compression','none'); %overwrite existing
    for k=2:stkNFrames,
        imwrite( uint16(stkOutput(:,:,k)), stkFilename, 'WriteMode','append',  'Compression','none');
    end
    
    % Also save the peak locations for later comparison
    save( strrep(stkFilename,'.tiff','.loc'), 'peakLocations', '-ASCII' );
    
    %
    if ~params.useAllMovies
        waitbar( (i+nTracesUsed)/nTracesTotal, h );
    else
        waitbar( n/numel(bgMovieFilenames), h );
    end
    nTracesUsed = nTracesUsed +nPeaks;
end

close(h); disp(toc);


end
