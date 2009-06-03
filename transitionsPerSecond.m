function [meanTPS,stdTPS] = transitionsPerSecond( dwtFilenames )
% PERCENTTIME  Stable state probabilities
%
%   Calculates the percentage time spent in each non-zero FRET state
%   using a QuB idealization data file (DWT).  States are listed in columns
%   from low to high FRET.  Each dataset is a row.
%   
% CAUTION: assumes first model state (1) is a zero-FRET state and removes it!

rand('twister',sum(100*clock));

if nargout>1
    bootstrapN = 1000;
else
    bootstrapN = 1;
end


if nargin<1,
    
    dwtFilenames = cell(1,0);
    
    % Request filenames from user
    while 1,
        [file path]=uigetfile('*.dwt','Choose an idealization file:');
        if file==0, break; end  %user hit "cancel"

        dwtFilenames{end+1} = [path filesep file];
        disp( dwtFilenames{end} );
    end
else
    if ischar(dwtFilenames),
        filenames = {filenames};
    end
end

nFiles = numel(dwtFilenames);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end


meanTPS = zeros(0,0);
stdTPS = zeros(0,0);

% Load dwell lifetime info from idealization data
for i=1:nFiles,
    
    [meanTPS(i,:),stdTPS(i,:)] = TPS( dwtFilenames{i}, bootstrapN );
    
    % Make sure all idealizations have the same number of states
%     assert( size(output,1)==0 || numel(totals) == size(output,2), ...
%       'Model mismatch' );
    
end %for each file




function [meanTPS,stdTPS] = TPS( dwtfilename, bootstrapN )

% Load DWT data
[dwells,sampling,offsets,fretModel] = loadDWT( dwtfilename );
nTraces = length(dwells);
nStates = numel(fretModel)/2;

% Allocate space for bootstrapped results
bootstrapTPS = zeros(bootstrapN,1);

% Calculate percent time in each state for each subset
for s=1:bootstrapN, %for each subset

    % Generate random bootstrap datasets (drawn with replacement)
    if s==1,
        bootstrapSet = 1:nTraces;
    else
        bootstrapSet = floor(rand(nTraces,1)*nTraces)+1;
    end
    
    %
    bootstrapTPS(s) = totalTransitions( dwells(bootstrapSet), sampling );
    
end %for each subset


% Calculate mean and std PT across subsets
meanTPS = mean( bootstrapTPS );
% meanTPS = bootstrapTPS(1);
stdTPS  = std( bootstrapTPS );






function tps = totalTransitions( dwells, sampling )
totalTime = 0;
totalTransitions = 0;

for i=1:numel(dwells)

    states = dwells{i}(:,1);
    times  = double(dwells{i}(:,2)) * sampling/1000;
    nDwells = numel(states);
    
    for j=2:nDwells,
        if states(j-1)>1 && states(j)>1
            totalTransitions = totalTransitions + 1;
        end
    end
    
    totalTime = totalTime + sum( times(states>1) );
end

tps = totalTransitions/totalTime;




