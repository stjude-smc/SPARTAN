function [meanTPS,stdTPS] = transitionsPerSecond( dwtFilenames )
% transitionsPerSecond  Calculate average transition rates
%
%   TPS = transitionsPerSecond( FILES ) calculates the average number of
%   transitions per second observed in each of the specified dwell-time FILES.
%   Transitions to and from the dark state (blinking) are ignored.
%
%   See also: percentTime, makeplots.

%   Copyright 2007-2015 Cornell University All Rights Reserved.



%% USER TUNABLE PARAMTERS

% Only consider non-zero FRET states. zero-FRET state assumed to be state 1.
REMOVEZERO = true;


%% Get input filenames from user if not specified.
if nargin<1,
    % Request filenames from user
    dwtFilenames = getFiles('*.dwt','Choose an idealization file:');
    
else
    % If only a single file is specified, turn it into a cell array
    if ischar(dwtFilenames),
        dwtFilenames = {dwtFilenames};
    end
end

nFiles = numel(dwtFilenames);


%% Calculate transitions/sec for each datafile using bootstrap sampling.
bootfun = @(N,time) sum(N)/sum(time);  %transitions per second

meanTPS = zeros(0,0);
stdTPS = zeros(0,0);

for i=1:nFiles,    
    % Load DWT data
    [dwells,sampling] = loadDWT( dwtFilenames{i} );
    nTraces = length(dwells);

    nEvents = zeros( nTraces, 1 );
    totalTime = zeros( nTraces, 1 );

    for trace=1:nTraces,
        states = dwells{trace}(:,1);
        times  = dwells{trace}(:,2);

        if REMOVEZERO,
            % Remove dwells in lowest FRET state (assuming it is the dark state)
            times  = times(states>1);
            states = states(states>1);

            % Combine dwells that are now in the same state by converting into
            % an idealization and then back to a dwell-time sequence.
            if ~isempty(times),
                idl = dwtToIdl( [states times] );
                newDwt = RLEncode(idl);
                times  = newDwt(:,2);
            end
        end

        % Calculate number of events and total time in this trace.
        nEvents(trace)   = numel(times)-1;
        totalTime(trace) = sum(times)*sampling/1000; %in seconds.

    end %for each trace

    % Calculate bootstrap samples to estimate standard error.
    meanTPS(i,:) = bootfun(nEvents,totalTime);
    stdTPS(i,:) = std(  bootstrp(1000, bootfun, nEvents,totalTime)  );
    
end %for each file



%% Display the result.
if nFiles<2, return; end

figure;
nStates = size(meanTPS,2);
errorbar( repmat(1:nFiles,nStates,1)', meanTPS, stdTPS );

xlabel('File number');
ylabel('Transitions per second');
xlim([0.5 nFiles+0.5]);
set(gca,'XTick',1:nFiles);



end %FUNCTION transitionsPerSecond.





