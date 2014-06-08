function data = correctTraces( data, constants )
% CORRECTTRACES  Makes simple adjustments to traces
%
%   [D,A,F] = CORRECTTRACES( DATA )
%   Subtracts background fluorscence from both channels, using 100
%   frames after donor photobleaching.  Also calculates FRET
%   efficiency, with E=0 where donor is blinking or photobleached.
%   

% TODO: gamma correction (whole pipeline), background drift correction?

if nargin<2,
    constants = cascadeConstants;
end

% If the acceptor is a cell array, this is a three-color FRET experiment
% and we have two acceptor signals.
assert( isfield(data,'donor') && isfield(data,'acceptor') && ...
        ~isempty(data.donor) && ~isempty(data.acceptor) );
isThreeColor = isfield(data,'acceptor2');

if isfield(data,'donor2'),
    warning('correctTraces:multiDonor','This function is not designed for multiple donors');
end


% Calculate donor lifetime
lt = calcLifetime(data.total,constants.TAU,constants.NSTD);



% Subtract fluorescence intensity so the baseline after photobleaching in
% zero. For traces that do not photobleach, no correction is made, but the
% baseline will be close because an estimated background image is
% subtracted from each frame in gettraces.
for m=1:data.nTraces,

    s = lt(m)+5;  %ignore the frames around the photobleaching event
    range = s:min(s+constants.NBK,len);
    
    if numel(range)<10,
        continue; %not enough baseline to calculate background. skip trace.
    end

    % Make background correction
    data.donor(m,:)    = data.donor(m,:)    - mean( data.donor(m,range) );
    data.acceptor(m,:) = data.acceptor(m,:) - mean( data.acceptor(m,range) );
    
    if isThreeColor,
        data.acceptor2(m,:) = data.acceptor2(m,:) - mean( data.acceptor2(m,range) );
    end
end

% Calculate FRET efficiencies. For three-color, these could be the fraction
% of all energy (emitted directly or indirectly from the donor) that is
% emitted by one specific acceptor. In that case, FRET1+FRET2=1 always.
% Here I calculate FRET as the fraction of energy transferred to one
% acceptor as if the other acceptor were not there -- we only consider the
% not transferred to the other acceptor. In this case, FRET1+FRET2 can be
% more than 1. This has the advantage that the FRET signals can be
% semi-independent.
total = data.total;

if isThreeColor,
    %fret  = acceptor  ./ (donor+acceptor);
    %fret2 = acceptor2 ./ (donor+acceptor2);
    data.fret  = data.acceptor  ./ total;
    data.fret2 = data.acceptor2 ./ total;
else
    data.fret  = data.acceptor ./ total;
end



% Sets FRET value to 0 when intensity is below a calculated threshold
% based on the number of standard devations above background
for m=1:Ntraces,
    data.fret(m, lt(m):end) = 0;
    data.fret2(m, lt(m):end) = 0;

    s = lt(m)+5;
    range = s:min(s+constants.NBK,len);
    if numel(range)<10, continue; end 

    % Set FRET to 0 when Cy3 is blinking or photobleached:
    % ie, when FRET is below a calculated threshold = 4*std(background)
    darkRange = total(m,1:lt(m)) <= constants.blink_nstd*std(total(m,range));
    data.fret(m,darkRange)  = 0;
    data.fret2(m,darkRange) = 0;
end



