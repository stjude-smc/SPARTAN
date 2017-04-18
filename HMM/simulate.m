 function [dwt,data] = simulate( dataSize, sampling, model, varargin )
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
%      'totalIntensity'     ->  for simulated fluorescence trajectories
%      'stdTotalIntensity'  ->  std of distribution (0)
%          *** if stdBackground is set, traces will have varying
%          signal-to-noise ratios! 
%      'stdPhoton'      ->  noise within fluorescence traces (multiplicitive)
%      'stdBackground'  ->  background noise (additive)
%
%      'kBleach'        ->  rate of acceptor photobleaching (0)
%  
%  See also: simulateMovie, batchKinetics, QubModel.

%   Copyright 2007-2017 Cornell University All Rights Reserved.

%  TODO:
%   - re-implement random seed, if possible.
%   - Simulate anisotropy noise
%   - Simulate distribution of gamma across traces (stdGamma)
%   - Simulate effect of varying quantum yield in each state... 
%   - Simulate donor->accpetor (acceptor->donor) crosstalk. 


% CONSTANTS
kBleachDonorFactor = 0.5; %donor relative to acceptor bleaching rate.


%%
% PARSE REQUIRED PARAMETER VALUES
nTraces  = dataSize(1);
traceLen = dataSize(2);
framerate = 1/sampling;

Q = model.rates;
p0 = model.p0;
p0 = p0/sum(p0);

% Default parameter values. FIXME: should be in cascadeConstants?
params = struct('totalIntensity',500, 'stdTotalIntensity',0, 'stdBackground',0, ...
                'stdPhoton',0, 'gamma',1, 'shotNoise',true, 'emNoise',false, ...
                'kBleach',0); 
params = mergestruct( params, struct(varargin{:}) );

% Simulated acceptor photobleaching rate (0 disables).
% Data will have a FRET lifetime of 1/kBleach, but the exact donor and
% acceptor lifetimes will be the values that give that apparent rate with
% the donor bleaching kBleachDonorFactor as fast as the acceptor.
kBleachAcceptor = params.kBleach/(1+kBleachDonorFactor);
kBleachDonor    = params.kBleach - kBleachAcceptor;


%%

simFramerate = 1000; %1000/sec = 1ms frames
binFactor = round(simFramerate/framerate);



%%

% Generate a set of uniform random numbers for choosing states
wbh = parfor_progressbar(1.25*nTraces,'Simulating state sequences...');


%--- Draw lifetimes for each trace before photobleaching
if params.kBleach > 0,
    pbTimes = exprnd( 1/kBleachAcceptor, 1,nTraces );
    pbTimes = pbTimes*simFramerate;
else
    pbTimes = (traceLen+1)*ones(1,nTraces);
end


% Start the matlab thread pool if not already running. perfor below will
% run the calculations of the available processors.
constants = cascadeConstants;
if nTraces*traceLen/1000 > 10 && constants.enable_parfor,
    % Processing large TIFF movies is CPU limited. Use parfor to parallelize.
    pool = gcp;
    M = pool.NumWorkers;
else
    % For small datasets, do everything in the GUI thread (regular for loop).
    M = 0;
end

simTraces = nargout<2;


%--- Generate noiseless state trajectory at 1 ms
dwt  = cell(nTraces, 1);
noiseless_fret = zeros( nTraces, traceLen );
dt = 1000*sampling; %integration time (timestep) in ms.

parfor (i=1:nTraces,M)
% for i=1:nTraces,
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
        dwellTime = (-1000./l_tot) .* log( rand(size(l_tot)) );  %in ms
        
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
    
    % Convert states to classes and merge dwells in the same class.
    % Times must be rounded to 1ms time resolution for the .dwt format.
    [classes,times] = mergedwells( model.class(states), times ); %#ok<PFBNS>
    dwt{i} = [classes round(times)];
    
    if mod(i,20)==0,  %fixme; not very useful with long traces.
        wbh.iterate(20); %#ok<PFBNS>
    end
    
    
    % Don't simulate FRET traces if not requested to save time.
    if simTraces, continue; end
    
    % Generate noiseless FRET traces with time averaging.
    e = 1+cumsum(times'./dt);    % dwell end times (continuous frames)
    s = [1 e(1:end-1)+eps];     % dwell start times
    eh = floor(e);              % discrete frame in which dwell ends
    sh = floor(s);              % discrete frame in which dwell starts
    
    trace = zeros( 1,traceLen );
    
    for j=1:numel(times),
        dwellFretValue = model.mu( classes(j) );
        
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
end



%--- Simulate fluorescence traces, adding gaussian read noise
% Don't simulate traces unless requested to save time.
if nargout>1,
    wbh.message = 'Simulating noise...';

    % Generate total intensity profile
    variation = params.stdTotalIntensity*randn( nTraces,1 );
    intensity = params.totalIntensity + variation;
    intensity = repmat( intensity, 1, traceLen );
    
    % Generate noiseless fluorescence traces (variable intensity)
    donor    = (1-noiseless_fret).*intensity;
    acceptor = intensity-donor;
    
    % Alter fluorescence traces to simulate non-equal collection efficiencies
    % for donor and acceptor fluorescence intensities. The factor gamma
    % is a correction for this effect.
    % Not that this alters total fluorescence intensities!
    donor = (1/params.gamma)*donor;
    
    % Simulate donor photobleaching.
    % NOTE: this is not correctly time averaged.
    if params.kBleach > 0,
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
    donor    = donor*    (params.totalIntensity/obsIntensity);
    acceptor = acceptor* (params.totalIntensity/obsIntensity);
        
    % Simulate the effect of photophysical noise. Here the variance
    % is propotional to the intensity of the signal.
    donor    = donor    .*  (1+ params.stdPhoton*randn(size(donor))    );
    acceptor = acceptor .*  (1+ params.stdPhoton*randn(size(acceptor)) );
        
    % Simulate Poisson shot "noise" or the variance of signal intensity simply
    % due to photon statistics. Note that this also discretizes the data.
    if params.shotNoise
        donor    = poissrnd(donor);
        acceptor = poissrnd(acceptor);
    end
    
    % Add background and read noise to fluorescence traces. Here the
    % variance is invariant of the signal.
    alive = donor~=0;
    donor    = donor    + params.stdBackground*randn(size(donor));
    acceptor = acceptor + params.stdBackground*randn(size(donor));
    
    % Recalculate FRET from fluorescence traces
    fret = acceptor./(donor+acceptor);
    fret(~alive) = 0; %set undefined values to zero
    %assert( ~any(isnan( fret(:) )) );
    fret( isnan(fret) ) = 0;
end

close(wbh);
drawnow;

% Construct output Traces object
data = TracesFret(nTraces, traceLen);
data.fret     = fret;
data.donor    = donor;
data.acceptor = acceptor;
data.time     = sampling*1000*(0:traceLen-1);
data.fileMetadata(1).wavelengths = [532 640];


end %FUNCTION simulate




function [newstates,newtimes] = mergedwells(states, times)
% Merge repeated dwells in the same state/class into a single longer dwell.

newstates = states(1);
newtimes  = times(1);

for i=2:numel(states),
    % Accumulate time if still in the same state.
    if states(i)==newstates(end),
        newtimes(end) = newtimes(end) + times(i);
        
    % Add the next dwell if there is a transition to a new state.
    else
        newstates(end+1) = states(i);
        newtimes(end+1)  = times(i);
    end
end

newstates = newstates';
newtimes  = newtimes';

end
