 function [dwt,fret,donor,acceptor] =  simulate( ...
                                dataSize, sampling, model, varargin )
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
%  
%  TODO:
%   - Simulate anisotropy noise
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

mu    = model.mu;
Q = model.rates;
p0 = model.p0;
p0 = p0/sum(p0);
 

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
% if isfield(optargs,'randomSeed') && ~isempty(optargs.randomSeed)
%     randomSeed = optargs.randomSeed;
% else
%     randomSeed  = 1272729;
% end

% Simulated acceptor photobleaching rate (0 disables).
% Data will have a FRET lifetime of 1/kBleach, but the exact donor and
% acceptor lifetimes will be the values that give that apparent rate with
% the donor bleaching kBleachDonorFactor as fast as the acceptor.
if isfield(optargs,'kBleach') && ~isempty(optargs.kBleach)
    kBleach = optargs.kBleach;
    assert(kBleach~=Inf,'kBleach must be finite');
else
    kBleach = 0;
end

kBleachAcceptor = kBleach/(1+kBleachDonorFactor);
kBleachDonor    = kBleach - kBleachAcceptor;


%%

simFramerate = 1000; %1000/sec = 1ms frames
binFactor = round(simFramerate/framerate);

tic;



    
% FRET values are adjusted here so that they match the input values
% after the introduction of the gamma artifact later.
mu = mu ./ ( mu + gamma - gamma*mu );



%%

% Generate a set of uniform random numbers for choosing states
h = waitbar(0,'Simulating...');



%--- Draw lifetimes for each trace before photobleaching
if kBleach > 0,
    pbTimes = exprnd( 1/kBleachAcceptor, 1,nTraces );
    pbTimes = pbTimes*simFramerate;
end



% Start the matlab thread pool if not already running. perfor below will
% run the calculations of the available processors.
% If a random seed is used, the results will no longer be deterministic.
% Use for (not parfor) to get the predictable behavior.
% if isempty( gcp('nocreate') ),   parpool;   end


%--- Generate noiseless state trajectory at 1 ms
dwt  = cell(nTraces, 1);
noiseless_fret = zeros( nTraces, traceLen );
dt = 1000*sampling; %integration time (timestep) in ms.

% parfor i=1:nTraces,
for i=1:nTraces,   %use this instead to turn off multi-threading
    
    % Choose the initial state
    curState = find( rand <= cumsum(p0), 1, 'first' );
    
    
    endTime = 1000*(traceLen*sampling); %in ms
    states = [];
    times  = [];
    
    % Below is essentially the (direct) Gillipse algorithm. Dwell times are
    % drawn from exponential distributions and the first event to occur
    % takes the system to a new "current state", continuing until the end
    % of the trace (acceptor dye photobleaching).
    % Note that, as with the experiments, this process is biased toward
    % short dwells when the bleaching time is on the order of the dwell
    % times in any state (even if the last dwell isn't truncated).
    while sum(times)<pbTimes(i) && sum(times)<endTime, %end when no more dwells are needed.
        
        % Simulate the time until the next transition.
        l_tot = sum( Q(curState,:) ); %rate for any transition to occur (the sum of all rates)
        dwellTime = exprnd(1000./l_tot);        %in ms
        
        % Simulate the transition type with probabilities being the
        % fraction of transition of each type by rate.
        x = cumsum( Q(curState,:)./l_tot );
        nextState = find( rand<=x, 1, 'first' );
        
        states(end+1) = curState;
        times(end+1) = dwellTime;
        
        curState = nextState;
    end
    
    
    % Truncate last dwell to fit into time window
    totalTime = sum(times);
    times(end) = times(end) - ( totalTime-min(endTime-1,pbTimes(i)) );
    
    % Save the dwell-time information for .dwt file.
    % Times are rounded to 1ms time resolution. DWT files do not support
    % continuous time and this is usually enough to see most dwells.
    dwt{i} = [model.class(states) round(times)'];
    nDwells = numel(states);
    
    % Generate noiseless FRET traces with time averaging.
    e = 1+cumsum(times./dt);    % dwell end times (continuous frames)
    s = [1 e(1:end-1)+eps];     % dwell start times
    eh = floor(e);              % discrete frame in which dwell ends
    sh = floor(s);              % discrete frame in which dwell starts
    
    trace = zeros( 1,traceLen );
    
    for j=1:nDwells,
        dwellFretValue = mu( model.class(states(j)) );
        
        % 1. Isolated stretch (shortdwell within one frame).
        % Add the fret value to the current frame as a fraction of how much
        % time it takes up, so that the total "probability" (the first
        % number) sums to one for each frame at the end of the loop.        
        if sh(j)==eh(j),
            trace(sh(j)) = trace(sh(j)) + (e(j)-s(j))*dwellFretValue;
            continue;
        end
    
        % 2. If the dwell covers several frames completely (the most common case),
        % just fill them in (weight is unity).
        trace( (sh(j)+1):(eh(j)-1) ) = dwellFretValue;

        % 3. Dwell starts in the middle of a frame.
        trace(sh(j)) = trace(sh(j)) + ((sh(j)+1)-s(j))*dwellFretValue;

        % 4. Dwell ends in the middle of a frame.
        trace(eh(j)) = trace(eh(j)) + (e(j)-eh(j))*dwellFretValue;
    end
    
    noiseless_fret(i,:) = trace;
    
    if mod(i,25)==0,
        waitbar( 0.9*i/nTraces, h );  %crashes with parfor
    end
end


%--- Simulate fluorescence traces, adding gaussian read noise

waitbar( 0.9, h, 'Adding noise...' );

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
    
    % Simulate donor photobleaching.
    % NOTE: this is not correctly time averaged.
    if kBleach > 0,
        pbTimes = exprnd( 1/kBleachDonor, 1,nTraces );
        pbTimes = round(pbTimes*simFramerate/binFactor);
        pbTimes = min( max(1,pbTimes), traceLen );
        
        for i=1:nTraces,
            donor(i,pbTimes(i):end) = 0;
            acceptor(i,pbTimes(i):end) = 0;
        end
    end
    
    % Rescale fluorescence intensities so that the total intensity
    % (and SNR) are roughly the same as if there was not gamma correction.
    nNonZero = sum( donor(:)>0 );
    obsIntensity = ( sum(donor(:))+sum(acceptor(:)) )/nNonZero;
    donor    = donor*    (totalIntensity/obsIntensity);
    acceptor = acceptor* (totalIntensity/obsIntensity);
        
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
    %fret = noiseless_fret + sigma(end)*randn( size(noiseless_fret) );
    %donor    = (1-fret).*totalIntensity;
    %acceptor = totalIntensity-donor;
end



close(h);



disp(toc);


end %FUNCTION simulate


