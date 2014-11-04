function data = correctTraces( data, constants )
% CORRECTTRACES  Makes simple adjustments to traces
%
%   [D,A,F] = CORRECTTRACES( DATA )
%   Subtracts background fluorscence from both channels, using 100
%   frames after donor photobleaching.  Also calculates FRET
%   efficiency, with E=0 where donor is blinking or photobleached.
%   

% TODO: gamma correction (whole pipeline), background drift correction?
% Consider splitting this into several functions for different types of data,
% doing background subtraction for all channels, but calculating FRET only where
% it makes sense.

if nargin<2,
    constants = cascadeConstants;
end


% Determine the type of experiment so we know how to process the data, and
% insure that the fields make sense.
assert( data.isChannel('donor'), 'Unregonized trace data (no donor channel)' );

if data.isChannel('acceptor2') && ~data.isChannel('acceptor'),
    error('Found acceptor2 but not acceptor1');
end

isFret       = data.isChannel('acceptor');
isThreeColor = data.isChannel('acceptor2');

if data.isChannel('donor2'),
    warning('correctTraces:multiDonor','This function is not designed for multiple donors');
end

if data.isChannel('factor'),
    disp('Warning: handling of factor signals is in an early stage and may be incomplete.');
end



% Calculate donor lifetime
lt = calcLifetime(data.total,constants.TAU,constants.NSTD);



% Subtract fluorescence intensity so the baseline after photobleaching in
% zero. For traces that do not photobleach, no correction is made, but the
% baseline will be close because an estimated background image is
% subtracted from each frame in gettraces.
[nTraces,len] = size(data.donor);

for m=1:nTraces,

    s = lt(m)+5;  %ignore the frames around the photobleaching event
    range = s:min(s+constants.NBK,len);
    
    if numel(range)<10,
        continue; %not enough baseline to calculate background. skip trace.
    end

    % Make background correction
    data.donor(m,:) = data.donor(m,:)- mean( data.donor(m,range) );
    
    if isFret,
        data.acceptor(m,:) = data.acceptor(m,:) - mean( data.acceptor(m,range) );
    end
    
    if isThreeColor,
        data.acceptor2(m,:) = data.acceptor2(m,:) - mean( data.acceptor2(m,range) );
    end
end



% Calculate FRET efficiencies. For three-color, these could be the fraction
% of all energy (emitted directly or indirectly from the donor) that is
% emitted by one specific acceptor. To get true FRET efficiencies, we will need
% to integrate data from alternating excitations (ALEX).
total = data.total;

if isFret,
    data.fret  = data.acceptor  ./ total;
else
    return;  %nothing more to do.
end

if isThreeColor,
    data.fret2 = data.acceptor2 ./ total;
end



% Sets FRET value to 0 when intensity is below a calculated threshold
% based on the number of standard devations above background
for m=1:nTraces,
    data.fret(m, lt(m):end) = 0;
    
    if isThreeColor,
        data.fret2(m, lt(m):end) = 0;
    end

    s = lt(m)+5;
    range = s:min(s+constants.NBK,len);
    if numel(range)<10, continue; end 

    % Set FRET to 0 when Cy3 is blinking or photobleached:
    % ie, when FRET is below a calculated threshold = 4*std(background)
    darkRange = total(m,1:lt(m)) <= constants.blink_nstd*std(total(m,range));
    data.fret(m,darkRange)  = 0;
    
    if isThreeColor,
        data.fret2(m,darkRange) = 0;
    end
end





