function [alive,thresholds] = thresholdTotal( total, thresholds )
%thresholdTotal   Idealize total fluorescence intensity using thresholds
%
%   [IDL,THRESH] = thresholdTotal(TOTAL) idealizes the total fluorescence 
%   intensity traces in the rows of the matrix TOTAL using an automatically-
%   calculated thresholds above baseline noise to detect donor dark states.
%   IDL is false where dark (blinking/bleached), true when it is fluorescent.
%   The tresholds used are returned in the vector THRESH. 
%
%   IDL = thresholdTotal(..., THRESH) uses the supplied thresholds in the 
%   vector THRESH instead of calculating them automatically.
%
%   See also: skmTotal, calcLifetime, TracesFret.recalculateFret.

%   Copyright 2015 Cornell University All Rights Reserved.

narginchk(1,2);
nargoutchk(1,2);

const = cascadeConstants;

% Process input arguments
[nTraces,nFrames] = size(total);
if nargin<2,
    thresholds = zeros(nTraces,1);
end
assert( isnumeric(thresholds) && numel(thresholds)==nTraces, 'Invalid threshold values');

% Determine where the donor bleached to restrict search for blinking events.
lt = max(1, calcLifetime(total) );

alive = true( size(total) );  %true where donor is alive.

for i=1:nTraces,
    % Set FRET to zero after donor photobleaching.
    alive( i, lt(i):end ) = false;

    % Set FRET to zero in areas where the donor is dark (blinking).
    % FIXME: Should manual thresholds be saved in traceMetadata?
    s = lt(i)+5;
    range = s:min(s+const.NBK, nFrames);
    if numel(range)>=10,
        if nargin<2, 
            thresholds(i) = const.blink_nstd*std1(total(i,range));
        end
        darkRange = total( i, 1:lt(i) ) <= thresholds(i);
        alive(i,darkRange) = false;
    end
end

end %function thresholdTotal
