function [LL,A,mu,sigma,p0,ps] = mkMethod( ...
          observations, framerate, nStates, varargin )
% MKMETHOD  McKinney2006 method for naive discovery of FRET values
%
%    [LL,A,mu,sigma,p0] = mkMethod( OBSERVATIONS, framerate, nStates, OPTIONS )
%    Runs the Baum-Welch parameter optimization proceedure on each trace,
%    obtaining FRET values and transition probabilities, and idealizes
%    each trace using these parameters.
%
%    OBSERVATIONS is an MxN matrix of FRET trajectories
%    (M traces, N datapoints per trace) that are used as input data for
%    optimization.  nStates specifies the number of states to use when
%    optimizing.  For naieve analysis, choose a number larger than you expect.
%
%    Output parameters are cell arrays (Mx1) of output parameters observed
%    in each trace.
%
%    (OPTIONAL): added argument "idealize",filename.dwt will idealize the
%    data using the fitted values

%    Algorithm:
%    1. Initial value for mu: to 0 (1/N)*[0:N-1]
%    2. Initial value for A:  all rates 2 /sec
%    3. BWoptimize to find highest likelihood of data given model
%    4. idealize using parameter estimates from 4
%    5. save idealization

%    Useful things to do with the results:
%    1. Plot histogram of mean FRET values (mu{:})
%    2. Plot transition density plots, fit to 2D histograms
%       to discover number of FRET states
%    3. Cluster the resulting data based on parameter estimates

% Process optional arguments
if nargin>3,
    fields = {varargin{1:2:end}};
    values = {varargin{2:2:end}};
    
    options = cell2struct( values', fields );
end


%-------------- USER TUNABLE PARAMETERS ----------------

initialRate  = 2; %/sec
initialSigma = 0.09;

%----- 1. Initial parameter values

% Create initial A-matrix
p0_start = ones(1,nStates)/nStates;   % state initial probabilities

A_start = zeros(nStates);       % transition probability matrix
x = initialRate*(1/framerate);
A_start(:) = x;
A_start( find(eye(nStates)) ) = 1-((nStates-1)*x); %normalization

assert( all(sum(A_start,2)-1 <0.0001) ); %test normalization

% Create initial FRET distribution parameters
N = nStates;
% mu_start    = (1/N).*(0:N-1);  %ex, [0.0, 0.2, 0.4, 0.6, 0.8], N=5
% mu_start    = (1/(N-1)).*(0:N-1);
mu_start = [0.01 0.52 0.75];
sigma_start = repmat( initialSigma, 1, N );


%------ 2. Optimize kinetic parameter values
[nTraces,nFrames] = size(observations);

LL    = cell(nTraces,1);
A     = cell(nTraces,1);
mu    = cell(nTraces,1);
sigma = cell(nTraces,1);
p0    = cell(nTraces,1);
ps    = cell(nTraces,1);

h = waitbar(0, 'Optimizing parameters...');
for i=1:nTraces,
    trace = observations(i,:);
    
    [LL{i},A{i},mu{i},sigma{i},p0{i},ps{i}] = BWoptimize( ...
                trace, A_start, mu_start, sigma_start, p0_start );
    waitbar(i/nTraces,h);
end
close(h);


%------ 3. Idealize each trace using optimized parameter values

% Skip this step unless explicitly requested
if ~exist('options','var') || ~isfield(options,'idealize'),
    return;
end
dwtFilename = options.idealize;

%
idealization = cell(nTraces,1);
models       = cell(nTraces,1);

for i=1:nTraces,
    % Load trace, ignoring zero-FRET region at end
    trace = observations(i,:);
    NT = find(trace~=0,1,'last')+4;
    NT = min(NT,length(trace));
    trace = trace(1:NT);
    
    % Idealize using parameter values optimized for this trace
    models{i} = [mu{i}' sigma{i}'];
    idl = idealize( trace, models{i}, p0{i}, A{i} );
    idealization{i} = idl{1}; %first and only trace
end

% Add offsets to relate idealization back to raw data
offsets = (0:(nTraces-1))*nFrames;

% Save idealization to dwell-time file
sampling = 1000/framerate; %in ms
saveDWT(dwtFilename, idealization, offsets, models, sampling);




