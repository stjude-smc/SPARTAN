function [meanPT,stdPT] = percentTime( filenames )
% PERCENTTIME  Stable state probabilities
%
%   Calculates the percentage time spent in each non-zero FRET state
%   using a QuB idealization data file (DWT).  States are listed in columns
%   from low to high FRET.  Each dataset is a row.
%   
% CAUTION: assumes first model state (1) is a zero-FRET state and removes it!

rand('twister',sum(100*clock));

bootstrapN = 10000;


if nargin<1,
    % Request filenames from user
    filenames = getFiles('*.dwt','Choose an idealization file:');
    
else
    % If only a single file is specified, turn it into a cell array
    if ischar(filenames),
        filenames = {filenames};
    end
end

nFiles = numel(filenames);

if nFiles<1,
    disp('No files specified, exiting.');
    return;
end


% Load dwell lifetime info from idealization data
meanPT = zeros(0,0);
stdPT = zeros(0,0);

for i=1:nFiles,
    [meanPT(i,:),stdPT(i,:)] = calcPT( filenames{i}, bootstrapN );
end %for each file



end % FUNCTION percentTime





function [meanPT,stdPT] = calcPT( dwtfilename, bootstrapN )

% Load DWT data
[dwells,sampling,offsets,fretModel] = loadDWT( dwtfilename );
nStates = numel(fretModel/4);
nTraces = length(dwells);

% Calculate percent time of each trace seperately
tracePT = zeros( nTraces,nStates );

for i=1:nTraces %for each trace
    
    states = dwells{i}(:,1);
    times  = dwells{i}(:,2);

    % Add times to total
    for j=1:nStates,
        tracePT(i,j) = tracePT(i,j) + sum( times(states==j) );
    end

end %for each trace


% Allocate space for bootstrapped results
bootstrapPT = zeros(bootstrapN,nStates-1);

% Calculate percent time in each state for each subset
for s=1:bootstrapN, %for each subset
    
    % Generate random bootstrap datasets (drawn with replacement)
    bootstrapSet = floor(rand(nTraces,1)*nTraces)+1;
    
    % Sum total times across all traces in bootstrap dataset
    totals = sum( tracePT(bootstrapSet,:) );

    % Remove 0-FRET state
    totals = totals(2:end);
    
    % Calculate and save percentages
    bootstrapPT(s,:) = 100*totals/sum(totals);
    
end %for each subset



% Calculate mean and std PT across subsets
meanPT = mean( bootstrapPT, 1 );
stdPT  = std( bootstrapPT, 0, 1 );


end % FUNCTION calcPT





