function [data,dwt] = simphotons( dataSize, sampling, model, varargin )
% SIMPHOTONS   Simulate fluorophore photophysics
%
%   DATA = SIMPHOTONS(SIZE, SAMPLING, MODEL) simulates a photon emission
%   photophysical process using the Gillespie algorithm.
%
%     SIZE is the number of traces and frames to simulate, respectively.
%     SAMPLING is the time binning resolution (in milliseconds).
%     MODEL is a QubModel object describing the photophysical behavior.
%       Transitions of states 2->1 contribute to photon counting.
%     DATA is a Traces object, with fluorescence counts in 'donor' and a
%     normalized trace in 'fret'.
%
%   ... = SIMPHOTONS(..., PARAMS) specifies parameters as a struct:
%     'snr'       -> signal to background noise ratio of output traces
%     'detection' -> percent of emitted photons that are recorded
%  
%   See also: gillespie, simulate, batchKinetics, QubModel.

%   Copyright 2018 Cornell University All Rights Reserved.



% PARSE REQUIRED PARAMETER VALUES
nTraces  = dataSize(1);
traceLen = dataSize(2);

% Default parameter values. FIXME: should be in cascadeConstants?
params = struct('stdBackground',0, 'detection',22); 
params = mergestruct( params, struct(varargin{:}) );

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



%% Simulate noiseless fluorescence traces
tic;

wbh = parfor_progressbar(nTraces,'Simulating photon series...');
noiseless_traces = zeros( nTraces, traceLen );
% dwt = cell(nTraces, 1);


parfor (i=1:nTraces,M)
% for i=1:nTraces
    
    %--- Simulate state dwells using the (direct) Gillespie algorithm.
    [traceStates,traceTimes] = gillespie( model, 1000*(traceLen*sampling) );
    
    
    %---- Sum photon arrivals into time bins to get fluorescence traces.
    
    % 1. Find photon emission events (transition from state 2 to 1)
    events = traceStates(1:end-1)==2 & traceStates(2:end)==1;
    
    % 2. Get arrival times of each photon in fractional frames.
    eventTimes = cumsum(traceTimes);
    eventTimes = eventTimes(events) ./ (1000*sampling);
    
    % 2. Remove a fraction of photons to simulate detection efficiency
    nEvents = numel(eventTimes);
    idxKeep = randi( nEvents, 1, floor(nEvents*params.detection/100) );
    eventTimes = eventTimes(idxKeep);

    % 3. Assign each photon to the appropriate time bin.
    frameIdx = floor(eventTimes+1);

    % 4. Sum events for each frame to get total fluorescence intensity.
    if ~isempty(frameIdx)
        noiseless_traces(i,:) = histc( frameIdx, 1:traceLen );
    end


    % Save dwell-time series.
    % FIXME: .dwt files only support CLASS series, but these are mostly
    % useless for this type of simulation...
    %if nargout>1
    %    dwt{i} = [model.class(traceStates)' traceTimes'];
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
data.donor = noisy_traces;
data.fret  = 0.8 * noiseless_traces ./ mti;  %for display purposes only
data.time  = sampling*1000*(0:traceLen-1);
data.fileMetadata(1).wavelengths = [532 640];


disp(toc);
close(wbh);
drawnow;



end %FUNCTION simulate




