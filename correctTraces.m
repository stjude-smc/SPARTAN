function [donor,acceptor,fret] = correctTraces( ...
                                      donor,acceptor, constants, indexes )
% CORRECTTRACES  Makes simple adjustments to traces
%
%   [D,A,F] = CORRECTTRACES( DONOR, ACCEPTOR, CONST, INDEXES )
%   Subtracts background fluorscence from both channels, using 100
%   frames after donor photobleaching.  Also calculates FRET
%   efficiency, with E=0 where donor is blinking or photobleached.
%
%   INDEXES specifies the indexes of traces to load. Useful if only a small
%   number of traces are needed from a large file.
%   

% TODO: gamma correction (whole pipeline), background drift correction?

if nargin<2,
    error('CorrectTraces: not all required parameters given!')
end
if nargin<3,
    constants = cascadeConstants;
end

% If the acceptor is a cell array, this is a three-color FRET experiment
% and we have two acceptor signals.
isThreeColor = 0;

if iscell(acceptor),
    assert( numel(acceptor)==2 );
    a = acceptor{1};
    acceptor2 = acceptor{2};
    acceptor = a;
    isThreeColor=1;
end

% Extract just the traces requested (default=all)
if nargin>=4 && ~isempty(indexes),
    donor = donor(indexes,:);
    acceptor = acceptor(indexes,:);
    if isThreeColor,
        acceptor2 = acceptor2(indexes);
    end
end

% Calculate donor lifetime
if isThreeColor,
    total = donor+acceptor+acceptor2;
else
    total = donor+acceptor;
end

lt = calcLifetime(total,constants.TAU,constants.NSTD);


[Ntraces,len] = size(donor);

% Subtract fluorescence intensity so the baseline after photobleaching in
% zero. For traces that do not photobleach, no correction is made, but the
% baseline will be close because an estimated background image is
% subtracted from each frame in gettraces.
for m=1:Ntraces,

    s = lt(m)+5;  %ignore the frames around the photobleaching event
    range = s:min(s+constants.NBK,len);
    
    if numel(range)<10,
        continue; %not enough baseline to calculate background. skip trace.
    end

    % Make background correction    
    donor(m,:) = donor(m,:)       - mean( donor(m,range) );
    acceptor(m,:) = acceptor(m,:) - mean( acceptor(m,range) );
    
    if isThreeColor,
        acceptor2(m,:) = acceptor2(m,:) - mean( acceptor2(m,range) );
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
if isThreeColor,
    total = donor+acceptor+acceptor2;
    %fret  = acceptor  ./ (donor+acceptor);
    %fret2 = acceptor2 ./ (donor+acceptor2);
    fret  = acceptor  ./ total;
    fret2 = acceptor2 ./ total;
else
    total = donor+acceptor;
    fret  = acceptor ./ total;
    fret2 = zeros(size(donor));
end



% Sets FRET value to 0 when intensity is below a calculated threshold
% based on the number of standard devations above background
for m=1:Ntraces,

    fret(m, lt(m):end) = 0;
    fret2(m, lt(m):end) = 0;

    s = lt(m)+5;
    range = s:min(s+constants.NBK,len);
    if numel(range)<10, continue; end 

    % Set FRET to 0 when Cy3 is blinking or photobleached:
    % ie, when FRET is below a calculated threshold = 4*std(background) 
    stdbg = std( total(m,range) );
    fret(m, total(m,1:lt(m))<=constants.blink_nstd*stdbg)=0;
    fret2(m, total(m,1:lt(m))<=constants.blink_nstd*stdbg)=0;
end


% For three-color FRET, acceptors are returned as they are sent in: as a
% cell array.
if isThreeColor,
    acceptor = {acceptor,acceptor2};
    fret = {fret,fret2};
end


