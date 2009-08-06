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
%          *** if stdFluorescence is set, traces will have varying
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


%%
% PARSE REQUIRED PARAMETER VALUES
nTraces  = dataSize(1);
traceLen = dataSize(2);
framerate = 1/sampling;

mu    = model(:,1);
sigma = model(:,2);
nStates = length(mu);

% Zero probabilities must be avoided, because they can lead
% to zero probability traces, which will cause an overflow/crash.
Q( Q<=0 ) = eps;
 

% PARSE OPTIONAL PARAMETER VALUES: initial kinetic parameter values
optargs = struct( varargin{:} );

if isfield(optargs,'totalIntensity') % total intensity of fluorescence
    totalIntensity = optargs.totalIntensity;
else
    totalIntensity = 10000;
end

if isfield(optargs,'stdTotalIntensity')
    stdTotalIntensity = optargs.stdTotalIntensity;
else
    stdTotalIntensity = 0;
end

if isfield(optargs,'stdPhoton')
    stdPhoton = optargs.stdPhoton;
else
    stdPhoton = 0;
end

if isfield(optargs,'stdBackground')
    stdBackground = optargs.stdBackground;
else
    stdBackground = 0;
end

%random number generator seed value
if isfield(optargs,'randomSeed') && optargs.randomSeed~=0
    randomSeed = optargs.randomSeed;
else
    randomSeed  = 1272729;
    %randomSeed = sum(100*clock);
end

% Simulated acceptor photobleaching rate (0 disables)
if isfield(optargs,'kBleach')
    kBleach = optargs.kBleach;
else
    kBleach = 0;
end



%%

simFramerate = 1000; %1ms
binFactor = simFramerate/framerate;
% roundTo = 25; %ms

tic;


% Construct transition probability matrix (A) for 1ms simulation
A = Q./simFramerate;
A( find(eye(nStates)) ) = 1-sum(A,2);


% Predict steady-state probabilities for use as initial probabilities.
if isfield(optargs,'startProb') && ~isempty(optargs.startProb)
    p0 = optargs.startProb;
else
    p0 = A^10000;
    p0 = p0(1,:);
    p0 = p0/sum(p0);
end


clear A;



%%

rand('twister',101+randomSeed);

s = rand();

% Generate a set of uniform random numbers for choosing states
h = waitbar(0,'Simulating...');

%--- Generate noiseless state trajectory at 1 ms
idl  = zeros( nTraces, traceLen*binFactor );  %state assignment at every time point
dwt  = cell(nTraces, 1);

for i=1:nTraces,
    
    % Choose the initial state
    curState = find( rand <= cumsum(p0), 1, 'first' );
    
    % Draw dwell times from exponential distribution
    % This is important so that results change as little as possible
    % when a single parameter is modified.
    endTime = 1000*(traceLen*sampling); %in ms
    states = [];
    times  = [];
    
    randData = rand(floor(endTime),nStates);
    itr=1;
    
    while sum(times)<endTime,
        
        assert( itr<floor(endTime), 'simulate: N. dwells exceeded RNG buffer size' );
        
        % Draw dwell time from exponential distribution
        tau = 1000./ Q( curState, : ); %in ms.
        choices = -tau .* log( randData(itr,:) );
        %choices(curState) = Inf;
        
        % Choose event that happens first
        [dwellTime,nextState] = min( choices );
        if sum(choices == dwellTime) > 1,
            %if more than one event are chosen, try again (never happens?)
            itr = itr+1;
            continue;
        end
        
        % Round to 1ms time resolution
        if exist('roundTo','var')
            dwellTime = round(dwellTime/roundTo)*roundTo;
            if dwellTime<25,  %FIXME: framerate dependant test!
                itr = itr+1;
                continue;
            end
        else
            dwellTime = round(dwellTime);
        end
        
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

    idl(i,:) = trace;
    
    waitbar(i/(nTraces*2),h);
end




%--- Chop each trace to simulate exponential photobleaching of acceptor dye
rand('twister',s);

if kBleach > 0,
    pbTimes = -log(rand(1,nTraces))./kBleach;
    pbTimes = ceil(pbTimes*simFramerate);

    for i=1:nTraces,
        idl(i,pbTimes(i):end) = 1;  % set the state to blinking
    end 
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
if stdPhoton~=0 || stdBackground~=0
    % Generate total intensity profile
    variation = stdTotalIntensity*randn( nTraces,1 );
    intensity = totalIntensity + variation;
    intensity = repmat( intensity, 1, traceLen );
    
    % Generate noiseless fluorescence traces (variable intensity)
    donor    = (1-noiseless_fret).*intensity;
    acceptor = intensity-donor;
    
    % Simulate donor photobleaching
    if kBleach > 0,
        pbTimes = -log(rand(1,nTraces))./(0.7*kBleach);
        pbTimes = ceil(pbTimes*simFramerate);
        
        for i=1:nTraces,
            pbtime = (pbTimes(i)/binFactor);
            pbtime = round( min(pbtime,traceLen) );
            pbtime = max(1,pbtime);

            donor(i,pbtime:end) = 0;
            acceptor(i,pbtime:end) = 0;
        end
    end
    
    % Simulate fluctuations in fluorescence intensity (shot noise,etc).
    donor    = donor    .* (1+ stdPhoton*randn(size(donor)) );
    acceptor = acceptor .* (1+ stdPhoton*randn(size(donor)) );
    
    % Simulate the effect of donor->acceptor channel crosstalk
    acceptor = acceptor + 0.075*donor;
    
    % Add background noise to fluorescence traces
    alive = donor~=0;
    donor    = donor    + stdBackground*randn(size(donor));
    acceptor = acceptor + stdBackground*randn(size(donor));
    
    % Generate final FRET trajectories
    fret = acceptor./(donor+acceptor);
    fret(~alive) = 0; %set undefined values to zero
    assert( ~any(isnan( fret(:) )) );
    %fret(fret>1) = 1;
    

% If fluorescence noise is not specified, add noise directly to fret
% and generate perfectly correlated donor/acceptor fluorescence traces.
% THIS DOES NOT WORK!
else
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
output = zeros( 1,nFrames2 );

for j=1:nFrames2

    s = factor*(j-1) +1;
    e = factor*j;

    output( j ) = mean( trace(s:e) );

end


end %FUNCTION binFretData
