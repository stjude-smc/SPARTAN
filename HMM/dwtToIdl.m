function [idl,dwtIDs] = dwtToIdl( dwt, offsets, traceLen, nTraces )
% dwtToIdl  Convert dwell-time list to state assignment traces (idealization)
%
%   IDL = dwtToIdl( DWT, OFFSETS, LEN, NTRACES ) creates the state assignment
%   matrix IDL from the cell array of state-dwelltime pairs DWT.
%   LEN is the trace length in frames. OFFSETS are the zero-based offsets
%   indicating the start (in frames) of each idealized segment.
%   States are numbered 1..N, with zero indicating no idealization information.
%   If IDL will ever be used with FRET data, it is strongly recommended to use
%   this version. The versions below may not result in matching sized matrices.
% 
%   IDL = dwtToIdl( DWT, OFFSETS, LEN ) guesses the number of traces from the
%   offsets, but may be short if the last trace(s) are not idealized.
% 
%   IDL = dwtToIdl( DWT, OFFSETS ) uses the stride between traces as the trace
%   length.
%
%   IDL = dwtToIdl( DWT ) creates an idealization with a minimum size to
%   contain the dwell-time information, with no unidealized traces.
%   Most useful when IDL is not used in conjunction with coganate FRET data.
% 
%   [IDL,IDS] = dwtToIdl(...) also returns the DWT segement number associated
%   with each FRET trace (row).
% 
%   See also idlToDwt, loadDWT, saveDWT.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


if ~iscell(dwt),
    dwt = {dwt};
end

% If only the dwell-times are given, create a state assignment matrix that is
% the minimum size to contain the idealization, assuming all traces are
% idealized. Often used when using the idealization alone without FRET data.
if nargin==1,
    empty = cellfun(@isempty,dwt);
    traceLen = max(  cellfun( @(x) sum(x(:,2)), dwt(~empty) )  );
    offsets = (0:1:numel(dwt)-1)*traceLen;
    nTraces = numel(dwt);

% Guess the trace length from the interval between idealized traces.
elseif nargin==2
    traceLen = median( diff(offsets) );
    nTraces = numel(dwt);

% If the number of traces is not specified, create an idealization that is
% long enough to cover all idealized traces. 
% WARNING: If there are some traces at the end that are not idealized,
% idl will not be the same size as fret!
elseif nargin==3,
    nTraces = round(offsets(end)/traceLen)+1;
    assert( numel(dwt)<=nTraces );
end


idl = zeros(nTraces,traceLen);
dwtIDs = zeros( 1, nTraces );


for dwtID=1:numel(dwt)
    traceID = floor(offsets(dwtID)/traceLen)+1;
    dwtIDs(traceID)=dwtID;
    trace = dwt{dwtID};
    
    if isempty(trace), continue; end
    
    states = floor( trace(:,1) );
    times  = floor( trace(:,2) );

    % For all dwells in this trace, get the start and end times.
    ends = cumsum(times);
    starts = [1; ends(1:end-1)+1];
    
    if ends(end)>traceLen || offsets(dwtID)+traceLen > nTraces*traceLen
        error('Idealization is longer than trace length');
    end
    
    % Add dwells to idealization output
    for j=1:numel(states),
        idl(traceID, starts(j):ends(j)) = states(j);
    end
end




