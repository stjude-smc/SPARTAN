function [dwells,sampling,model] = loadDwelltimes( dwtfilename, dropLastDwell )
% Returns an cell array (Nx1, N states) which contains
% a list of the dwell times (in ms) for each state.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


if nargin<2,
    dropLastDwell=0;
end

% Load dwell times
[idl,sampling,~,model] = loadDWT( dwtfilename );
nStates = numel(model)/2;
nTraces = numel(idl);

% Add dwells to collection
dwells = cell(nStates,1);

for i=1:nTraces,
    d = idl{i};
    
    % Skip last dwell if parameter set
    if dropLastDwell,
        d = d(1:end-1,:);
    end
    
    states = d(:,1);
    times  = double( d(:,2) );
    
    % Add all dwell times in state j to list (dwells{j})
    for j=1:nStates
        dwells{j} = [dwells{j} ; times(states==j)];
    end
end

% Convert dwell times from frames to msec
for j=1:nStates
    dwells{j} = dwells{j}.*sampling;
end
