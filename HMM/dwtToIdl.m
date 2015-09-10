function [idl,dwtIDs] = dwtToIdl( dwt, traceLen, offsets, nTraces )
% dwtToIdl  Convert dwell-time list to state assignment traces (idealization)
%
%    IDL = dwtToIdl( DWT, LEN, OFFSETS, NTRACES ) creates an idealization (state
%    assignment trace, IDL) from a cell array of state-dwelltime pairs (DWT).
%    LEN is the trace length (in frames). OFFSETS are the zero-based offsets
%    indicating the start (in frames) of each idealized segment.
%    States are numbered 1..N, with zero indicating no idealization information.
% 
%    IDL = dwtToIdl( DWT, LEN, OFFSETS ) will create an idealization that
%    includes all idealized traces, 
% 
%    IDL = dwtToIdl( DWT, LEN ) returns an idealization assuming that each
%    segment (element in DWT) is at the beginning of sequential traces
% 
%    [IDL,IDS] = dwtToIdl( ... ) returns a list of the DWT segement number
%    associated with each trace.
% 
%    See also idlToDwt, loadDWT, saveDWT.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% If no offsets are provided, we assume all traces have an idealization,
% or that any empty idealizations have an empty element in dwt.
if nargin<3,
    offsets = traceLen*( 0:1:numel(dwt)-1 );
end

if nargin<4,
    % If the number of traces is not specified, create an idealization that is
    % long enough to cover all idealized traces. 
    % WARNING: If there are some traces at the end that are not idealized,
    %idl will not be the same size as fret!
    nTraces = ceil(offsets(end)/traceLen)+1;
    assert( numel(dwt)<=nTraces );
end

idl = zeros(nTraces,traceLen);
dwtIDs = zeros( 1, nTraces );


for dwtID=1:numel(dwt)
    traceID = floor(offsets(dwtID)/traceLen)+1;
    dwtIDs(traceID)=dwtID;
    
    states = floor( dwt{dwtID}(:,1) );
    times  = floor( dwt{dwtID}(:,2) );

    % For all dwells in this trace, get the start and end times.
    ends = cumsum(times);
    starts = [1; ends(1:end-1)+1];
    
    if ends(end)>traceLen,
        error('Idealization is longer than trace length');
    end
    
    % Add dwells to idealization output
    for j=1:numel(states),
        idl(traceID, starts(j):ends(j)) = states(j);
    end
end




