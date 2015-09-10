function [meanPT,stdPT] = percentTime( filenames, truncateLength )
% PERCENTTIME  Stable state probabilities
%
%   [MEAN,STD] = percentTime( FILESNAMES, LEN )
%
%   Calculates the percentage time spent in each non-zero FRET state
%   using a QuB idealization data file (DWT). FILENAMES is a cell array of
%   .dwt file names of datasets to process. For each file, bootstrapping is
%   used to estimate the statistical error in the measurement by finding
%   the variance in percent time of many random subsets of the data. MEAN
%   is the average (over traces) time spent in each state, with states
%   listed in columns from low to high FRET and one file per row. STD is
%   the standard deviation across these bootstrap sets.
%
%   LEN (optional) is the length of traces to use, which can be used to
%   truncate all traces to a specified length. This is useful to ensure
%   that the percent-time measurements here match with plots made with
%   makeplots (statehist, TD, and contour plots), which all use a specific
%   window size (if that option is enabled). It is also useful when the
%   data are biased because by differing rates of photobleaching in each
%   state, where the system collects in the state with the slowest
%   photobleaching rate.
%   
% CAUTION: assumes first model state (1) is a zero-FRET state and removes it!

%   Copyright 2007-2015 Cornell University All Rights Reserved.


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

% If no truncation length is set, don't truncate at all
if nargin<2,
    truncateLength = [];
end


% Load dwell lifetime info from idealization data
meanPT = zeros(0,0);
stdPT = zeros(0,0);

for i=1:nFiles,
    [meanPT(i,:),stdPT(i,:)] = calcPT( filenames{i}, bootstrapN, truncateLength );
end %for each file


% Make titles for each file
titles = strrep(filenames,'_',' ');
for i=1:nFiles,
    [~,titles{i}] = fileparts( titles{i} );
end

% Plot the results
figure;
nStates = size(meanPT,2);
errorbar( repmat(1:nFiles,nStates,1)', meanPT, stdPT );
legend(titles);
xlabel('File number');
ylabel('Fraction occupancy');

end % FUNCTION percentTime





function [meanPT,stdPT] = calcPT( dwtfilename, bootstrapN, truncateLength )

% Load DWT data
[dwells,sampling,offsets,fretModel] = loadDWT( dwtfilename );
nStates = numel(fretModel)/2;
nTraces = length(dwells);

% Calculate percent time of each trace seperately
tracePT = zeros( nTraces,nStates );

for i=1:nTraces %for each trace
    
    % Convert trace to idealization for easier handling.
    len = sum(dwells{i}(:,2));
    idl = dwtToIdl( dwells(i), len, 0 );
    
    % Truncate it to the specific length
    if ~isempty(truncateLength) && truncateLength<len,
        idl = idl( :, 1:truncateLength );
    end

    % Add times to total
    for j=1:nStates,
        tracePT(i,j) = tracePT(i,j) + sum( idl==j );
    end

end %for each trace


% Allocate space for bootstrapped results
bootstrapPT = zeros(bootstrapN,nStates-1);

% Calculate percent time in each state for each subset
for s=1:bootstrapN, %for each subset
    
    % Generate random bootstrap datasets (drawn with replacement)
    if s==1,
        bootstrapSet = 1:nTraces;
    else
        bootstrapSet = floor(rand(nTraces,1)*nTraces)+1;
    end
    
    % Sum total times across all traces in bootstrap dataset
    totals = sum( tracePT(bootstrapSet,:) );

    % Remove 0-FRET state
    totals = totals(2:end);
    
    % Calculate and save percentages
    bootstrapPT(s,:) = 100*totals/sum(totals);
    
end %for each subset



% Calculate mean and std PT across subsets.
% Mean is just the mean of all data, not across bootstrap sets.
meanPT = bootstrapPT(1,:);
stdPT  = std( bootstrapPT, 0, 1 );


end % FUNCTION calcPT





