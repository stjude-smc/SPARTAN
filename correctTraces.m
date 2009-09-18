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

% Extract just the traces requested (default=all)
if nargin>=4 && ~isempty(indexes),
    donor = donor(indexes,:);
    acceptor = acceptor(indexes,:);
end

% Calculate donor lifetime
total = donor+acceptor;
lt = calcLifetime(total,constants.TAU,constants.NSTD);


[Ntraces,len] = size(donor);

% make mean background fluorescence zero.
% no correction for traces with < 10 points of background
for m=1:Ntraces,

    s = lt(m)+5;
    range = s:min(s+constants.NBK,len);
    if numel(range)<10, continue; end

    % Make background correction
    md = mean( donor(m,range) );
    ma = mean( acceptor(m,range) );
    
    donor(m,:) = donor(m,:) - md;
    acceptor(m,:) = acceptor(m,:) - ma;
end

total = donor+acceptor;
fret  = acceptor ./ total;


% Sets FRET value to 0 when FRET is below a calculated threshold
% based on the number of standard devations above background
for m=1:Ntraces,

    fret(m, lt(m):end) = 0;

    s = lt(m)+5;
    range = s:min(s+constants.NBK,len);
    if numel(range)<10, continue; end 

    % Set FRET to 0 when Cy3 is blinking or photobleached:
    % ie, when FRET is below a calculated threshold = 4*std(background) 
    stdbg = std( total(m,range) );
    fret(m, total(m,1:lt(m))<=constants.blink_nstd*stdbg)=0;
end


