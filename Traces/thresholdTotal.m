function [alive,thresholds] = thresholdTotal( total, thresholds )
%skmTotal   Idealizes total fluorescence intensity using SKM donor blinking.
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

constants = cascadeConstants;
nBkg = constants.NBK;
nStd = constants.blink_nstd;

[nTraces,nFrames] = size(total);
if nargin<2,
    thresholds = zeros(nTraces,1);
end

alive = true( size(total) );  %true where alive.

lt = max(1, calcLifetime(total) );  %photobleaching event time

for i=1:nTraces,
    % Set FRET to zero after donor photobleaching.
    alive( i, lt(i):end ) = false;

    % Set FRET to zero in areas where the donor is dark (blinking).
    % FIXME: Should manual thresholds be saved in traceMetadata?
    s = lt(i)+5;
    range = s:min(s+nBkg, nFrames);
    if numel(range)>=10,
        if nargin<2, 
            thresholds(i) = nStd*std(total(i,range));
        end
        darkRange = total( i, 1:lt(i) ) <= thresholds(i);
        alive(i,darkRange) = false;
    end
end

end %function thresholdTotal
