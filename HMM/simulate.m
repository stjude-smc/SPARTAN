 function [dwt,fret,donor,acceptor] =  simulate( ...
                                dataSize, sampling, model, Q, varargin )
% SIMULATE   Simulate smFRET data
%
%    [IDL,FRET,DONOR,ACCEPTOR] = SIMULATE( SIZE, FRAMERATE, MODEL, Q )
%    Generates noiselss FRET trajectories (FRET) at 1ms time resolution
%    using the given FRET emission MODEL (col1=means, col2=stdevs) and
%    rates in the given rate matrix (Q(i,j)=k_i->j), averages to
%    FRAMERATE (/sec) to simulate the effect of missed events, and
%    adds Gaussian noise to produce the final FRET trajctory (DATA).
%
%    Fluorescence trajectories (DONOR, ACCEPTOR) are also simulated,
%    with specified mean total intensity, stdev across traces, 
%    and added fluorescence noise (only if specified).
%
%    OPTIONAL Parameter values:
%      'totalIntensity'     ->  for simulated fluorescence trajectories (10000)
%      'stdTotalIntensity'  ->  std of distribution (0)
%          *** if stdBackground is set, traces will have varying
%          signal-to-noise ratios! 
%      'stdPhoton'      ->  noise within fluorescence traces (multiplicitive)
%      'stdBackground'  ->  background noise (additive)
%
%      'randomSeed'     ->  for reproducable results
%      'kBleach'        ->  rate of acceptor photobleaching (0)

%  TODO:
%   - Simulate anisotropy noise
%   - Simulate effect of background photobleaching (either channel)
%   - Simulate effect of gamma~=1
%   - Simulate distribution of gamma across traces (stdGamma)
%   - Simulate effect of varying quantum yield in each state... 
%   - Simulate donor->accpetor (acceptor->donor) crosstalk. 
%

if nargin<1,
    simulate_gui;
    return;
end


% CONSTANTS

%donor bleaching is half as fast as apparent acceptor bleaching
kBleachDonorFactor = 0.5; 


%%
% PARSE REQUIRED PARAMETER VALUES
nTraces  = dataSize(1);
traceLen = dataSize(2);
framerate = 1/sampling;

mu    = model(:,1);
sigma = model(:,2);
nStates = length(mu);
 

% PARSE OPTIONAL PARAMETER VALUES: initial kinetic parameter values
optargs = struct( varargin{:} );

if isfield(optargs,'totalIntensity') && ~isempty(optargs.totalIntensity)
    totalIntensity = optargs.totalIntensity;
    assert( totalIntensity>=0, 'totalIntensity must be a positive number' );
else
    totalIntensity = 10000;
end

if isfield(optargs,'stdTotalIntensity') && ~isempty(optargs.stdTotalIntensity)
    stdTotalIntensity = optargs.stdTotalIntensity;
    assert( stdTotalIntensity>=0, 'stdBackground must be a positive number' );
else
    stdTotalIntensity = 0;
end


if isfield(optargs,'stdBackground') && ~isempty(optargs.stdBackground)
    stdBackground = optargs.stdBackground;
    assert( stdBackground>=0, 'stdBackground must be a positive number' );
else
    stdBackground = 0;
end

if isfield(optargs,'stdPhoton') && ~isempty(optargs.stdPhoton)
    stdPhoton = optargs.stdPhoton;
    assert( stdPhoton>=0, 'stdPhoton must be a positive number' );
else
    stdPhoton = 0;
end

if isfield(optargs,'gamma') && ~isempty(optargs.gamma)
    gamma = optargs.gamma;
    assert( gamma>0, 'gamma must be a positive number' );
else
    gamma = 1;
end

%random number generator seed value
if isfield(optargs,'randomSeed') && ~isempty(optargs.randomSeed)
    randomSeed = optargs.randomSeed;
else
    randomSeed  = 1272729;
end

% Simulated acceptor photobleaching rate (0 disables).
% Data will have an *apparent* acceptor bleaching rate of this value.
% Actual value used in simulation also depends on donor bleaching rate.
if isfield(optargs,'kBleach') && ~isempty(optargs.kBleach)
    kBleach = optargs.kBleach;
    assert(kBleach~=Inf,'kBleach must be finite');
else
    kBleach = 0;
end

kBleachDonor    = kBleachDonorFactor*kBleach;
kBleachAcceptor = kBleach - kBleachDonor;



%%

simFramerate = 1000; %1000/sec = 1ms frames
binFactor = simFramerate/framerate;

tic;


% Predict steady-state probabilities for use as initial probabilities.
if isfield(optargs,'startProb') && ~isempty(optargs.startProb)
    p0 = optargs.startProb;
    assert( abs(sum(p0)-1)<=1e-4, 'p0 is not normalized' );
    assert( numel(p0)==nStates, 'p0 is the wrong size' );
else
    % Calculate A matrix from rates using the matrix exponential method.
    Qorig = Q;
    Q( logical(eye(size(Q))) ) = 0;
    Q( logical(eye(size(Q))) ) = -sum( Q,2 );
    dt = 1/simFramerate;
    A = expm( Q.*dt );
    Q = Qorig;
    
    % If not specified, estimate the equilibrium probabilities.
    p0 = A^10000;
    p0 = p0(1,:);
    p0 = p0/sum(p0);
end

% FRET values are adjusted here so that they match the input values
% after the introduction of the gamma artifact later.
mu = mu ./ ( mu + gamma - gamma*mu );



%%

rand('twister',101+randomSeed);
savedSeed = rand(); % save a seed for later that is not dependant
                    % on how many timesrand() is called.

% Generate a set of uniform random numbers for choosing states
h = waitbar(0,'Simulating...');



%--- Draw lifetimes for each trace before photobleaching
rand('twister',savedSeed);

if kBleach > 0,
    pbTimes = -log(rand(1,nTraces))./kBleachAcceptor;
    pbTimes = ceil(pbTimes*simFramerate);
end



%--- Generate noiseless state trajectory at 1 ms
idl  = zeros( nTraces, traceLen*binFactor );  %state assignment at every time point
dwt  = cell(nTraces, 1);

for i=1:nTraces,
    
    % Choose the initial state
    curState = find( rand <= cumsum(p0), 1, 'first' );
    
    % Draw dwell times from exponential distribution.
    % This is important so that results change as little as possible
    % when a single parameter is modified.
    endTime = 1000*(traceLen*sampling); %in ms
    states = [];
    times  = [];
    
    randData = rand(floor(endTime),nStates);
    itr=1;
    
    while sum(times)<pbTimes(i) && sum(times)<endTime, %end when no more dwells are needed.
        
        assert( itr<size(randData,1), 'simulate: N. dwells exceeded RNG buffer size' );
        
        % Draw dwell time from exponential distribution
        tau = 1000./ Q( curState, : ); %in ms.
        choices = -tau .* log( randData(itr,:) );
        
        % Choose event that happens first
        [dwellTime,nextState] = min( choices );
        if sum(choices == dwellTime) > 1,
            %if more than one event are chosen, try again (never happens?)
            itr = itr+1;
            continue;
        end
        
        % Round to 1ms time resolution
        dwellTime = round(dwellTime);
        
        states(end+1) = curState;
        times(end+1) = dwellTime;
        
        curState = nextState;
        itr = itr+1;
    end
    
    % Truncate last dwell to fit into time window
    totalTime = sum(times);
    times(end) = times(end) - (totalTime-endTime);
    
    dwt{i} = [states' times'];
    
    % Sample the series at 1ms to produce "real" trace
    trace = [];
    nDwells = numel(states);

    for j=1:nDwells,
        dtime = round( times(j) );
        trace = [trace repmat( states(j), 1,dtime ) ];
    end
    
    % Truncate idealization from photobleaching times
    trace(pbTimes(i):end) = 1;  % set the state to blinking

    % Save results to output data.
    idl(i,:) = trace;
    
    waitbar(i/(nTraces*2),h);
end



%--- Generate and and time-average the noiseless FRET trajectories (slow)
noiseless_fret = zeros( nTraces, traceLen );

waitbar(0.5,h, 'Time-averaging...');

for i=1:nTraces,
    
    trace = mu( idl(i,:) );
    trace = binFretData( trace, binFactor );
    noiseless_fret(i,:) = trace;
    
    waitbar( (i+nTraces)/(2*nTraces), h);
end


%--- Simulate fluorescence traces, adding gaussian read noise

randn('state',111+randomSeed);

% Add read/background noise to fluorescence and recalculate FRET
if stdBackground~=0 || stdPhoton ~=0
    % Generate total intensity profile
    variation = stdTotalIntensity*randn( nTraces,1 );
    intensity = totalIntensity + variation;
    intensity = repmat( intensity, 1, traceLen );
    
    % Generate noiseless fluorescence traces (variable intensity)
    donor    = (1-noiseless_fret).*intensity;
    acceptor = intensity-donor;
    
    % Alter fluorescence traces to simulate non-equal collection efficiencies
    % for donor and acceptor fluorescence intensities. The factor gamma
    % is a correction for this effect.
    % Not that this alters total fluorescence intensities!
    donor = (1/gamma)*donor;
    
    % Simulate donor photobleaching
    if kBleach > 0,
        pbTimes = -log(rand(1,nTraces))./kBleachDonor;
        pbTimes = ceil(pbTimes*simFramerate);
        
        for i=1:nTraces,
            pbtime = (pbTimes(i)/binFactor);
            pbtime = round( min(pbtime,traceLen) );
            pbtime = max(1,pbtime);

            donor(i,pbtime:end) = 0;
            acceptor(i,pbtime:end) = 0;
        end
    end
    
    % Rescale fluorescence intensities so that the total intensity
    % (and SNR) are roughly the same as if there was not gamma correction.
    nNonZero = sum( donor(:)>0 );
    obsIntensity = ( sum(donor(:))+sum(acceptor(:)) )/nNonZero
    donor    = donor*    (totalIntensity/obsIntensity);
    acceptor = acceptor* (totalIntensity/obsIntensity);
    obsIntensity = ( sum(donor(:))+sum(acceptor(:)) )/nNonZero
        
    % Simulate the effect of photophysical noise. Here the variance
    % is propotional to the intensity of the signal.
    donor    = donor    .*  (1+ stdPhoton*randn(size(donor))    );
    acceptor = acceptor .*  (1+ stdPhoton*randn(size(acceptor)) );
    
    % Add background and read noise to fluorescence traces. Here the
    % variance is invariant of the signal.
    alive = donor~=0;
    donor    = donor    + stdBackground*randn(size(donor));
    acceptor = acceptor + stdBackground*randn(size(donor));
    
    % Recalculate FRET from fluorescence traces
    fret = acceptor./(donor+acceptor);
    fret(~alive) = 0; %set undefined values to zero
    assert( ~any(isnan( fret(:) )) );
    
    
% If fluorescence noise is not specified, add noise directly to fret
% and generate perfectly correlated donor/acceptor fluorescence traces.
% THIS DOES NOT WORK!
else
    error('Direct FRET simulation not supported');
    fret = noiseless_fret + sigma(end)*randn( size(noiseless_fret) );
    donor    = (1-fret).*totalIntensity;
    acceptor = totalIntensity-donor;
end



close(h);



disp(toc);


end %FUNCTION simulate




function output = binFretData( trace, factor )

assert( factor >= 1 );
assert( floor(factor)==factor, 'binFretData: only integer values accepted' );

traceLen = numel(trace);
nFrames2 = floor(traceLen/factor);

% Group each set of frames to be averaged as its own column vector.
% Then, each group is time averaged by taking the mean of each vector.
x = reshape(trace,[factor nFrames2]);
output = mean(x);


end %FUNCTION binFretData
