function result = rleFilter( x, tau )
% RLE_FILTER  Run-length-encoding filter
%
%   VALUES = RLE_FILTER( SIGNAL, TAU )
%   Filters logical array SIGNAL, leaving only runs more than TAU in length

%   Copyright 2007-2015 Cornell University All Rights Reserved.

% Parts adopted from:
% http://home.online.no/~pjacklam/matlab/doc/mtt/index.html

if isempty(x), result=x; return; end

% Create a run-length-encoded version of signal
i = [ find(x(1:end-1) ~= x(2:end)) numel(x) ]; %run end positions
len = diff([ 0 i ]);
val = x(i);

% Remove points without required number of neighbors
val( len<=tau ) = 0;

% Re-expand the RLE signal
%i = cumsum(len);
j = zeros(1, i(end));
j(i(1:end-1)+1) = 1;
j(1) = 1;
result = val(cumsum(j));

end

