function [meanTPS,stdTPS] = transitionsPerSecond( dwtFilenames )
% transitionsPerSecond  Calculate average transition rates
%
%   TPS = transitionsPerSecond( FILES )
%   Calculates the average number of transitions per second observed in
%   each of the specified dwell-time FILES (see loadDWT.m).

%   Calculates the percentage time spent in each non-zero FRET state
%   using a QuB idealization data file (DWT).  States are listed in columns
%   from low to high FRET.  Each dataset is a row.
%   
% NOTE: unlike percentTime, transitions to the zero-FRET state are not removed.

% CAUTION: assumes first model state (1) is a zero-FRET state and removes it!


rand('twister',sum(100*clock));


% USER TUNABLE PARAMTERS
bootstrapN = 10000;
removeBlinks = 1;


% Get input filenames from user if not specified.
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

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end



% Calculate transitions/sec for each datafile using bootstrap sampling.
meanTPS = zeros(0,0);
stdTPS = zeros(0,0);

for i=1:nFiles,
    [meanTPS(i,:),stdTPS(i,:)] = TPS( dwtFilenames{i}, bootstrapN, removeBlinks );
end %for each file

end %FUNCTION transitionsPerSecond.





function [meanTPS,stdTPS] = TPS( dwtfilename, bootstrapN, removeBlinks )

% Load DWT data
[dwells,sampling] = loadDWT( dwtfilename );
nTraces = length(dwells);

% Calculate number of events & total dwell time for each trace.
traceNevents = zeros( nTraces, 1 );
traceTotalTime = zeros( nTraces, 1 );

for i=1:nTraces %for each trace
    
    states = dwells{i}(:,1);
    times  = dwells{i}(:,2);
    
    if removeBlinks,
        % Remove dwells in lowest FRET state (assuming it is the dark state)
        times  = times( states>1 );
        states = states( states>1 );
        
        % Combine dwells that are now in the same state by converting into
        % an idealization and then back to a dwell-time sequence.
        idl = dwtToIdl( {[states times]}, sum(times), 0 );
        newDwt = RLEncode(idl);
        states = newDwt(:,1);
        times  = newDwt(:,2);
    end
    
    % Calculate data needed to calculate TPS.
    traceNevents(i)   = numel(times)-1;
    traceTotalTime(i) = sum(times)*sampling/1000; %in seconds.
    
end %for each trace


% Calculate percent time in each state for each subset
bootstrapTPS = zeros(bootstrapN,1);

for s=1:bootstrapN, %for each subset

    % Generate random bootstrap datasets (drawn with replacement)
    bootstrapSet = floor(rand(nTraces,1)*nTraces)+1;
    
    % Calculate TPS for bootstrap dataset
    bootstrapTPS(s) = sum( traceNevents(bootstrapSet)     )  /  ...
                      sum( traceTotalTime(bootstrapSet)   );
    
end %for each subset


% Calculate mean and std PT across subsets
meanTPS = mean( bootstrapTPS );
stdTPS  = std( bootstrapTPS );

end %FUNCTION TPS


