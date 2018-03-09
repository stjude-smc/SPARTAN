 function [data,dwt] = simphotons( dataSize, sampling, model, varargin )
% SIMPHOTONS   Simulate fluorophore photophysics
%
%    DATA = SIMPHOTONS(SIZE, SAMPLING, MODEL) performs a stochastic simulation
%    of MODEL (QubModel object) using the Gillespie algorithm, where each dwell
%    in states of class 2 (or higher) adds a photon to the frame in the output
%    trace at the time it was sampled. This function is designed for simulating
%    photophysical processes by simulating photon emission series directly.
%    SAMPLING is the output time resolution in seconds.
%
%      DATA.donor contains raw fluorescence counts.
%      DATA.fret has the raw fluorescence counts normalized to 1.
%      DATA.acceptor also includes background noise.
%
%    OPTIONAL Parameter values:
%      'snr'  ->  signal to background noise ratio of output traces
%  
%  See also: simulate, batchKinetics, QubModel.

%   Copyright 2007-2018 Cornell University All Rights Reserved.



% PARSE REQUIRED PARAMETER VALUES
nTraces  = dataSize(1);
traceLen = dataSize(2);

Q = model.rates;
p0 = to_row(model.p0);
p0 = p0/sum(p0);

% Default parameter values. FIXME: should be in cascadeConstants?
params = struct('stdBackground',0); 
params = mergestruct( params, struct(varargin{:}) );



%% Simulate noiseless fluorescence traces
% FIXME: for some reason, parfor is SLOWER here.
tic;


% Initialize variables for Gillespie algorithm
wbh = parfor_progressbar(nTraces,'Simulating photon series...');
dwt  = cell(nTraces, 1);
noiseless_traces = zeros( nTraces, traceLen );
endTime = 1000*(traceLen*sampling); %in ms
nStates = numel(p0);


% Pre-calculate state time constants for Gillespie direct method
Qtau    = cell(1,nStates);  %mean dwell time for each state
Qcumsum = cell(1,nStates);  %cumsum of probability of each possible exit from a state
for s=1:nStates
    Qtau{s}   = -1000 / sum( Q(s,:) );
    Qcumsum{s} = cumsum(  Q(s,:) ./ sum(Q(s,:))  );
end


for i=1:nTraces
    
    traceStates = [];
    traceTimes  = [];
    cumTime = 0;
    
    % Sample initial state from initial probabilities distribution (p0)
    curState = find( rand <= cumsum(p0), 1, 'first' );
    
    %--- Simulate state dwells using the (direct) Gillespie algorithm.
    while cumTime<endTime,
        
        % Randomly sample the time until the next transition.
        dwellTime = Qtau{curState} .* log(rand);  %in ms
        
        % Randomly sample final state with probabilities calculated as the
        % fraction of all possible rate constants exiting current state.
        nextState = find( rand<=Qcumsum{curState}, 1 );  %'first' is default
        
        traceStates(end+1) = curState;
        traceTimes(end+1) = dwellTime;

        cumTime = cumTime + dwellTime;
        curState = nextState;
        
    end %while not enough dwells to fill trace
    
    
    %---- Sum photon arrivals into time bins to get fluorescence traces.
    % FIXME: photons should only be emitted from a specific transition, not any
    % dwell in the singlet excited state.
    
    % 1. Get explicit arrival times of each photon in fractional frames
    traceClasses = model.class(traceStates);
    eventTimes = cumsum(traceTimes);
    eventTimes = eventTimes(traceClasses>1) ./ (1000*sampling);

    % 2. Assign each photon to the appropriate time bin.
    frameIdx = floor(eventTimes+1);

    % 3. Sum events for each frame to get total fluorescence intensity.
    if ~isempty(frameIdx)
        noiseless_traces(i,:) = histc( frameIdx, 1:traceLen );
    end


    % Save dwell-time series.
    % FIXME: .dwt files only support CLASS series, but these are mostly
    % useless for this type of simulation...
    %if nargout>1
    %    dwt{i} = [traceClasses' traceTimes'];
    %end

    
    if mod(i,5)==0,  %fixme; not very useful with long traces.
        wbh.iterate(5);
    end
    
end %for each trace



%% Simulate noise

% Add background and read noise to fluorescence traces
mti = mean( noiseless_traces(:,1) );
stdbg = mti/params.snr;  %/(sqrt(2)??
noisy_traces = noiseless_traces + stdbg*randn(size(noiseless_traces));


% Construct output Traces object.
% FIXME: consider setting the 'fret' channel normalized to the expected
% intensity if there are no non-radiative pathways (from the 1->2 rate constant)
% to give a visual of the amount of quenching.
data = TracesFret(nTraces, traceLen);
data.donor    = noisy_traces;
data.acceptor = noiseless_traces;
data.fret = 0.8 * noiseless_traces ./ mti;  %for display purposes only
data.time     = sampling*1000*(0:traceLen-1);
data.fileMetadata(1).wavelengths = [532 640];


disp(toc);
close(wbh);
drawnow;



end %FUNCTION simulate




