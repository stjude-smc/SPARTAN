function data = bgsub( data, constants )
% BGSUB  Automatically subtract background labels after donor bleaching
%
%   BGSUB( DATA ) subtracts background fluorscence from all channels, 
%   using frames just after donor photobleaching.  Does not recalc. FRET!

%   Copyright 2007-2015 Cornell University All Rights Reserved.

if nargin<2,
    constants = cascadeConstants;
end

[nTraces,len] = size(data.donor);


% Determine the type of experiment so we know how to process the data, and
% insure that the fields make sense.
assert( data.isChannel('donor'), 'Unregonized trace data (no donor channel)' );

if data.isChannel('acceptor2') && ~data.isChannel('acceptor'),
    error('Found acceptor2 but not acceptor1');
end

isFret       = data.isChannel('acceptor');
isThreeColor = data.isChannel('acceptor2');

if data.isChannel('donor2'),
    warning([mfilename ':multiDonor'],'This function is not designed for multiple donors');
end

if data.isChannel('factor'),
    disp('Warning: handling of factor signals is in an early stage and may be incomplete.');
end


% Calculate donor lifetime
lt = calcLifetime(data.total,constants.TAU,constants.NSTD);

if data.isChannel('acceptorDirect')
    ltA = calcLifetime(data.acceptorDirect,constants.TAU,constants.NSTD);
    bgA = zeros(nTraces,1);
else
    ltA = [];
end


% Subtract fluorescence intensity so the baseline after photobleaching in
% zero. For traces that do not photobleach, no correction is made, but the
% baseline will be close because an estimated background image is
% subtracted from each frame in gettraces.
% FIXME: consider making this a method in Traces class
bg = zeros(nTraces,3);

for m=1:nTraces,
    
    % Subtract acceptorDirect signal
    if ~isempty(ltA)
        s = ltA(m)+5;  %ignore the frames around the photobleaching event
        range = s:min(s+constants.NBK,len);
        nRange = numel(range);
        if nRange >= 10
            bgA(m,1) = sum( data.acceptorDirect(m,range) )/nRange;
        end
    end

    s = lt(m)+5;  %ignore the frames around the photobleaching event
    range = s:min(s+constants.NBK,len);
    nRange = numel(range);
    
    % Make background correction
    if nRange >= 10
        bg(m,1) = sum( data.donor(m,range) )/nRange;

        if isFret,
            bg(m,2) = sum( data.acceptor(m,range) )/nRange;
        end

        if isThreeColor,
            bg(m,3) = sum( data.acceptor2(m,range) )/nRange;
        end
    end
end

data.donor = bsxfun(@minus, data.donor, bg(:,1));

if isFret,
    data.acceptor = bsxfun(@minus, data.acceptor, bg(:,2));
end

if ~isempty(ltA)
    data.acceptorDirect = bsxfun(@minus, data.acceptorDirect, bgA(:,1));
end

if isThreeColor,
    data.acceptor2 = bsxfun(@minus, data.acceptor2, bg(:,3));
end


