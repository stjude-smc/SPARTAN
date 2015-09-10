function Z = keq( varargin )
% KEQ  Calcualtes equilibrium constant between each pair of states
%
%   Z = KEQ( dwtFilename ) loads the idealization in <dwtFilename> and
%   returns the equilibrium contants between each pair of states. 
%   Z(2,3) is [TotalTime instate 2]/[TotalTime instate 3].
%

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% ask user for data file
if nargin<1,
    f = getFiles('*.dwt');
    dwtFilename = f{1};
    
elseif ischar( varargin{1} )
    dwtFilename = varargin{1};
        
end


% Load idealization
[dwells,sampling,model] = loadDwelltimes( dwtFilename, 0 );
nStates = length(dwells);

totalTimes = zeros( nStates, 1 );

% Sum dwell times
for i=1:nStates,    
    totalTimes(i) = sum( dwells{i} );
end

% Calculate Keq
Z = zeros( nStates, nStates );

for i=1:nStates,    
    for j=1:nStates,
        Z(i,j) = totalTimes(i)/totalTimes(j);
    end
end
