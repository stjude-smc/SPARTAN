function [idl,dwtIDs] = dwtToIdl( dwt, traceLen, offsets )
% dwtToIdl  Converts a list of dwell times to state assignment (idealization)
%
%    IDL = dwtToIdl( DWT, LEN, OFFSETS )
% 
%    Expands the cell array list of state-dwelltime pairs (DWT) into an
%    idealization (state number assignment at each point in time. LEN is
%    the trace length (in frames). OFFSETS are the zero-based offsets
%    indicating the start (in frames) of each trace. See loadDWT and
%    saveDWT for more information. States are numbered 1..N, with zero
%    indicating areas without idealization information in DWT.
%
%    dwtIDs is an array of indexes into dwt that correspond to each trace.
%    This assumes at most one dwt per trace.
%
%    See also idlToDwt, loadDWT, saveDWT.
%
% NOTE: if the last trace(s) are not idealized, idl will have a different
% (slightly smaller) size than the trace data because the last trace(s) are
% not represented!
%


nTraces = ceil(offsets(end)/traceLen)+1;

idl = zeros(nTraces,traceLen);
dwtIDs = zeros( 1, nTraces );


for dwtID=1:numel(dwt)
    traceID = floor(offsets(dwtID)/traceLen)+1;
    dwtIDs(traceID)=dwtID;
    
    states = dwt{dwtID}(:,1);
    times  = dwt{dwtID}(:,2);

    % For all dwells in this trace, get the start and end times.
    ends = cumsum(times);
    starts = [1; ends(1:end-1)+1];
    
    for j=1:numel(states),
        idl(traceID, starts(j):ends(j)) = states(j);
    end
end




