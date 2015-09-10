function [path,offsets] = tIdealize( traces, level, options )
% Assigns to one of two states based on threshold level,
% ignoring 1-frame dwells on either side.
% Last dwell is ignored because it includes dark time.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


[nTraces,len] = size(traces);

path = cell(0,1);

% Second threshold level to determine photobleaching point.
% if nargin<3,
    bleachingLevel = 0.2;
% end

% Apply an initial filtering to the data = 2-point running average
% Neccessary to avoid noise transitions
% if nargin<3,
    windowSize = 3;
% end

traces = runningAverage( traces, windowSize );
% traces = medianfilter( traces, windowSize );

% Threshold the data
% 0 = dark, 1 = mid, 2 = high
idl = zeros( size(traces) );

idl( traces>level ) = 3;
idl( traces<=level ) = 2;


% Run-length encode to produce idealization
nOutput = 1;
skipped = [];

for i=1:nTraces,
    
    traceLen = find( traces(i,:)>bleachingLevel, 1, 'last' )-1;
    
    % Run-length encode
    if isempty(traceLen) || traceLen<2,
        skipped = [skipped i];
    else
        path{nOutput} = RLEncode( idl(i,1:traceLen) );
        nOutput = nOutput+1;
    end
end

% Add offsets to relate idealization back to raw data
offsets = (0:(nTraces-1))*len;
offsets = offsets( ~ismember(1:nTraces,skipped) );

end %function tIdealize



function output = runningAverage( input, windowSize )
% Filter runs down the columns, but <input> is in the rows...

% Pad input with the first value.
% By default, effectively padded w/ zeros, leading to a decreased level...
paddedInput = [repmat(input(:,1), 1,windowSize) input];

weights = ones(1,windowSize)/windowSize;
output = filter( weights, 1, paddedInput' )';

% Remove padded region
output = output( :, 1+windowSize:end );

end %function runningAverage

