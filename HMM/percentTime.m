function [meanPT,stdPT] = percentTime( filenames, truncateLength )
% PERCENTTIME  Stable state probabilities
%
%   [MEAN,STD] = percentTime( FILESNAMES ) calculates the percentage time spent 
%   in each non-zero FRET state (MEAN) and standard error of the measurement
%   (STD) estimated with bootstrap samples from the dwell-time information in
%   each file in the cell array FILENAMES. For both outputs, states are listed
%   across columns and files across rows.
%
%   [...] = percentTime( FILESNAMES, LEN ) truncates the traces to LEN frames.
%   This is useful to match the window size of makeplots, particularly to avoid
%   potential bias cuased by differing photobleaching rates in each state.
%   
%   NOTE: the zero-FRET state (1) is ignored.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


%% OPTIONS

% Only consider non-zero FRET states. zero-FRET state assumed to be state 1.
REMOVEZERO = true;


%% Get filename from user if not specified.
if nargin<1,
    filenames = getFiles('*.dwt','Choose an idealization file:');    
else
    % If only a single file is specified, turn it into a cell array
    if ischar(filenames),
        filenames = {filenames};
    end
end

nFiles = numel(filenames);


%% Load dwell-times and calculate percent time in each state for each file.
bootfun = @(times) 100*sum(times)/sum(times(:));

meanPT = zeros(0);
stdPT = zeros(0);

for i=1:nFiles,    
    % Load dwell-time information and convert to state assignment matrix.
    [dwt,~,~,model] = loadDWT(filenames{i});
    nStates = size(model,1);
    
    idl = dwtToIdl(dwt);
    [nTraces,len] = size(idl);

    % Truncate the idealization if necessary
    if nargin>=2,
        idl = idl( :, 1:min(len,truncateLength) );
    end

    % Calculate percent time of each trace seperately
    tracePT = zeros(nTraces, nStates);

    for state=1:nStates,
        tracePT(:,state) = sum(idl==state,2);
    end

    % Remove zero state from consideration:
    if REMOVEZERO,
        tracePT = tracePT(:,2:end);
    end

    % Calculate bootstrap samples to estimate standard error.
    meanPT(i,:) = bootfun(tracePT);
    stdPT(i,:)  = std(  bootstrp(1000, bootfun, tracePT)  );
end


%% Plot the results
if nFiles<2, return; end

figure;
nStates = size(meanPT,2);
errorbar( repmat(1:nFiles,nStates,1)', meanPT, stdPT );

% Construct titles with the state number and FRET values.
states = (1:nStates)+REMOVEZERO;
fret = model(states,1);

titles = cell(nStates,1);
for i=1:nStates,
    titles{i} = sprintf('State %d (%.2f)\t',states(i),fret(i));
end
legend(titles);

xlabel('File number');
ylabel('Fraction occupancy');
xlim([0.5 nFiles+0.5]);
set(gca,'XTick',1:nFiles);



end % FUNCTION percentTime




