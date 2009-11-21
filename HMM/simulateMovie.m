function stkOutput = simulateMovie(varargin)
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

% TODO: provide options for simulating misalignment, choose parameter values
%  like buffer size and the region over which Gaussians are distributed
%  based on the size of the point-spread function.


%% Parse input parameters

nargs = numel(varargin);

% No parameters given; nothing to do but exit.
if nargin<3 || nargin>4,
    error('Invalid parameter list');
end

% First parameter is a filename; load traces data from file.
if ischar(varargin{1})
    tracesFilename = varargin{1};
    [donor,acceptor,ids,time] = loadTraces( tracesFilename );
    
    bgMovieFilenames = varargin{2};
    params = varargin{3};
    
% First parameter is a matrix of numbers, specifying donor/acceptor traces.
else
    tracesFilename = [pwd filesep 'sim.traces'];
    donor    = varargin{1};
    acceptor = varargin{2};
    time = 1:size(donor,2);
    
    bgMovieFilenames = varargin{3};
    params = varargin{4};
end

% If only a single background movie filename is given, convert it to cell array.
if ischar(bgMovieFilenames),
    bgMovieFilenames = {bgMovieFilenames};
end


%% Assign default values to parameters not specified.
if ~isfield(params,'sigmaPSF') || isempty(params.sigmaPSF)
    params.sigmaPSF = 0.8;
end
sigmaPSF = params.sigmaPSF;

if ~isfield(params,'randomSeed') || isempty(params.randomSeed),
    params.randomSeed = 0;
end

if ~isfield(params,'density') || isempty(params.density),
    params.density = 300;
end

if ~isfield(params,'edgeBuffer') || isempty(params.edgeBuffer),
    params.edgeBuffer = ceil(sigmaPSF*3);
end

 
% Initialize random number generator. Insure values are not correlated
% if traces are generated with the same random seed.
rand( 'twister', params.randomSeed+1299 );



%% Alter fluorescence intensities to simulate spectral artifacts.

% Simulate the effect of donor to acceptor channel crosstalk
constants = cascadeConstants;
acceptor = acceptor + constants.crosstalk*donor;



%% Simulate Wide-field Fluorescence Movies.

[nTracesTotal,traceLen] = size(donor);
nTracesUsed = 0;

h = waitbar(0,'Distributing fluorescence information...'); tic;

for n=1:numel(bgMovieFilenames)

    % --- 1. Load background movie as initial image to which peaks are added.
    stk = tiffread( bgMovieFilenames{n} );
    [stkY,stkX,stkNFrames] = size(stk);
    assert( mod(stkX,2)==0, 'Movie widths must be even!' );
    assert( traceLen<=stkNFrames, 'Background movie is too short!' );
    
    stkOutput = double( stk(:,:,1:traceLen) );  %use background movies

    
    % --- 2. Select fluorophore centroid positions randomly.
    border = params.edgeBuffer;
    
    % If user requests a regular grid of positions, add some small
    % variability (up to 1 pixel) to approximate random placement.
    if isfield(params,'grid') && params.grid,
        spacing = ceil(5*sigmaPSF);
        
        x = (1+border):spacing:(stkX/2-border);  nX = numel(x);
        y = (1+border):spacing:(stkY-border);    nY = numel(y);
        x = repmat( x, [nY 1] );
        y = repmat( y, [nX 1] )';
        
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
    if nPeaks>(nTracesTotal-nTracesUsed)
        disp('Not enough traces to fill this movie');
        break;
    end
    

    % --- 3. Distribute fluorescence intensities into a point-spread function
    %        around each centroid.
    A = 1/(2*pi* sigmaPSF^2); % 2D Gaussian normalization factor
    
    for i=1:nPeaks,
        traceID = i+nTracesUsed;
        d = donor(traceID,:);
        a = acceptor(traceID,:);
        
        cy = donorPos(i,1);
        cx = donorPos(i,2);
        limit = ceil(4*sigmaPSF);

        % Define point-spread function
        for x=max(1,floor(cx)-limit):min(stkX/2,ceil(cx)+limit)
            for y=max(1,floor(cy)-limit):min(stkY,ceil(cy)+limit)
                % Donor fluorophore
                pixelDonor    = A.*exp(-0.5* ((x-cx)^2+(y-cy)^2)/(sigmaPSF^2)  ).*d;
                stkOutput(y,x,:) = stkOutput(y,x,:) + reshape(pixelDonor,[1 size(pixelDonor)]);
                
                % Acceptor fluorophore
                % TODO: introduce optional misalignment error here.
                pixelAcceptor = A.*exp(-0.5* ((x-cx)^2+(y-cy)^2)/(sigmaPSF^2)  ).*a;
                stkOutput(y,x+stkX/2,:) = stkOutput(y,x+stkX/2,:) + reshape(pixelAcceptor,[1 size(pixelAcceptor)]);
            end
        end

        waitbar( (0.95*i+nTracesUsed)/nTracesTotal, h );
    end

    % --- 4. Save newly generated movie to file
    [p,name] = fileparts(tracesFilename);
    stkFilename = strrep(name,'.traces','');
    stkFilename = [stkFilename num2str(n) '.movie'];
    
    fid = fopen( stkFilename, 'w' );
    fwrite( fid, [stkX stkY traceLen], 'int16' );
    % fwrite( fid, time, 'float' );
    fwrite( fid, stkOutput, 'int16' );
    fclose(fid);
    
    waitbar( ((i+nTracesUsed)/nTracesTotal), h );
    nTracesUsed = nTracesUsed +nPeaks;
end

close(h); disp(toc);


end
