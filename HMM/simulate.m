 function [data,dwt] = simulate( nTraces, traceLen, sampling, model, varargin )
% SIMULATE   Stochastic simulation of fluorescence and FRET traces
%
%    DATA = SIMULATE( nTraces, nFrames, SAMPLING, MODEL )
%    Generates fluorescence and FRET traces (returned in the Traces object
%    DATA) from stochastic simulation of the dynamic process described by
%    the QubModel object MODEL. State sequences and exponentially 
%    distributed dwell times are simulated using the Gillespie algorithm.
%    These continuous-time dwells are used to generate FRET traces time
%    averaged to the time resolution specified in the SAMPLING argument
%    in seconds). nTraces/nFrames define the size of the simulated data.
%
%    [DATA,DWT] = SIMULATE(...) also returns the list of dwell states and
%    times for each trace (one cell array element per trace). Note that
%    states rather than classes are specified and times are not rounded,
%    unlike the typical dwell-time format (.dwt file). Use with caution.
%    
%    SIMULATE(..., PARAMS) specifies additional parameters as a struct:
%      totalIntensity    : mean donor+acceptor intensity per frame.
%      stdTotalIntensity : stdev of total intensity across molecules.
%      snr          : signal-to-background noise ratio.
%      shotNoise    : if true, simulate the Poisson-distributed 'noise'
%                        of with photon counting (true).
%      gamma        : scale donor intensity to simulate unequal apparent
%                         acceptor/donor brightness.
%      donLife      : average time before donor bleaching (seconds)
%      accLife      : average time before acceptor bleaching (seconds)
%      simExpBleach : if true, use exponentially-distributed bleaching
%                         times. If false, use fixed times for all traces.
%  
%  See also: gillespie, simphotons, simulateMovie, batchKinetics, QubModel.

%   Copyright 2007-2018 Cornell University All Rights Reserved.


%  TODO:
%   - Simulate distribution of gamma across traces (stdGamma)
%   - Simulate effect of varying quantum yield in each state... 
%   - Simulate donor->accpetor (acceptor->donor) crosstalk. 



%% Parse input parameters

narginchk(4,5);
nargoutchk(1,2);
saveDwt = nargout>=2;

% FIXME: should defaults be defined in cascadeConstants?
params = struct('totalIntensity',500, 'stdTotalIntensity',0, ...
                'snr',30, 'shotNoise',true, 'gamma',1, ...
                'accLife',Inf, 'donLife',Inf, 'simExpBleach',false); 
params = mergestruct( params, struct(varargin{:}) );

% Draw photobleaching times (in ms) from exponential distributions.
% endTime is the last point to simulate dwells.
if params.simExpBleach
    accBleachTimes = exprnd( 1000*params.accLife, 1,nTraces );
    donBleachTimes = exprnd( 1000*params.donLife, 1,nTraces );
else
    accBleachTimes = 1000*params.accLife*ones(1,nTraces);
    donBleachTimes = 1000*params.donLife*ones(1,nTraces);
end
endTimes = min( min(accBleachTimes,donBleachTimes), 1000*(traceLen*sampling)-1 );


% Paralellize across cores except simulations of small datasets.
if nTraces*traceLen/1000 > 10 && cascadeConstants('enable_parfor')
    pool = gcp;
    M = pool.NumWorkers;
else
    M = 0;
end



%% Simulate noiseless FRET trajectories

wbh = parfor_progressbar(1.25*nTraces,'Simulating state sequences...');

dwt = cell(nTraces, 1);
noiseless_fret = zeros( nTraces, traceLen );
dt = 1000*sampling; %integration time (timestep) in ms.

parfor (i=1:nTraces,M)
% for i=1:nTraces
    
    %--- Simulate state dwells using the (direct) Gillespie algorithm.;
    [states,dwelltimes] = gillespie( model, endTimes(i) );
    
    if saveDwt, dwt{i}=[to_col(double(states)) to_col(dwelltimes)]; end
    
    % Truncate last dwell to fit into time window
    dwelltimes(end) = dwelltimes(end) - (sum(dwelltimes)-endTimes(i));
    dwelltimes = to_col(dwelltimes);
    
    
    %--- Generate noiseless FRET traces with time averaging.
    e = 1+cumsum(dwelltimes'./dt);    % dwell end times (continuous frames)
    s = [1 e(1:end-1)+eps];     % dwell start times
    eh = floor(e);              % discrete frame in which dwell ends
    sh = floor(s);              % discrete frame in which dwell starts
    
    trace = zeros( 1,traceLen );
    fretValues = model.mu( model.class(states) );
    
    for j=1:numel(fretValues),
        dwellFretValue = fretValues(j);
        
        % 1. Isolated stretch (short dwell within one frame).
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
    
    if mod(i,10)==0, wbh.iterate(10); end
end



%% Simulate fluorescence traces, including experimental noise.
wbh.message = 'Simulating noise...';

% Simulate uneven excitation as Gaussian-distribution total intensities.
intensity = params.totalIntensity + params.stdTotalIntensity*randn( nTraces,1 );

% Generate noiseless fluorescence traces, scaled to simulate uneven
% apparent brightness of the two dyes (gamma).
donor    = bsxfun(@times, 1-noiseless_fret, intensity) /params.gamma;
acceptor = bsxfun(@times,   noiseless_fret, intensity);

% Simulate donor photobleaching by setting intensity to zero.
% FIXME: add time averaging by dividing by fraction of time in zero state
% in bleaching dwell. Should be easy to implement.
pbFrame = max(1, round(donBleachTimes/dt)+1 );
for i=1:nTraces
    if pbFrame(i)>traceLen, continue; end
    donor(i,pbFrame(i):end) = 0;
    acceptor(i,pbFrame(i):end) = 0;
end
dark = donor==0;

% Scale fluorescence intensities so that the observed mean total intensity
% (D+A) matches the value requested by the user. Necessary if gamma=/=1.
% obsIntensity = mean( donor(donor~=0) + acceptor(donor~=0) );
% donor    = donor*    (params.totalIntensity/obsIntensity);
% acceptor = acceptor* (params.totalIntensity/obsIntensity);

% Simulate Poisson-distributed photon counting statistics (shot noise).
if params.shotNoise
    donor    = poissrnd(donor);
    acceptor = poissrnd(acceptor);
end

% Add Gaussian-distributed background noise to fluorescence traces.
stdbg = params.totalIntensity/(sqrt(2)*params.snr);
donor    = donor    + stdbg*randn(size(donor));
acceptor = acceptor + stdbg*randn(size(donor));

% Calculate FRET efficiency traces from simulated fluorescence traces
fret = acceptor./(donor+acceptor);
fret( dark | isnan(fret) ) = 0;

% Construct output Traces object
data = TracesFret(nTraces, traceLen);
data.fret     = fret;
data.donor    = donor;
data.acceptor = acceptor;
data.time     = sampling*1000*(0:traceLen-1);
data.fileMetadata(1).wavelengths = [532 640];


close(wbh);
drawnow;


end %FUNCTION simulate



