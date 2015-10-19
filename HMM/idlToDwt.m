function [dwt,offsets] = idlToDwt( idl )
% idlToDwt  Converts state assignment (idealization) to a list of dwell times
%
%    [DWT,OFFSETS] = idlToDwt( IDL )
% 
%    Creates a cell array of dwell-time vector pairs. See loadDWT/saveDWT
%    for more details. IDL is a set of row vectors for each idealized
%    trace. Data that are not idealized are zero to indicate "no
%    information", the first state is 1, etc. All times are in frames.
%
%    See also dwtToIdl, loadDWT, saveDWT.
%
% NOTE: this is intended for state # / dwell-time data, not FRET
% idealizations. But it will still work as long as 0s are placed
% wherever the data aren't actually idealized.
%
% NOTE: we assume that there is at most one idealization per trace!

%   Copyright 2007-2015 Cornell University All Rights Reserved.


assert( nargin==1 && isnumeric(idl) );

[~,traceLen] = size(idl);


% Find un-idealized portions (zero state) as segment separators.
% Note that this includes the last one, which marks the end of the file!
segments = find( any(idl,2) )'; %select traces that have an idealization
offsets = (segments-1)*traceLen;

nSeg = numel(segments);
dwt = cell(1,nSeg);

for i=1:nSeg,
    seg = segments(i);
    
    % "Compress" the idealization into dwell-time pairs by converting 
    % runs of the same state index (or FRET value) into state/time pairs.
    d = RLEncode( idl(seg,:) );
    states = d(:,1);
    times  = d(:,2); %in frames
    
    % Save the dwt, removing zero-state dwells, which are actually regions
    % that were not idealized.
    dwt{i} = [ states(states>0) times(states>0) ];
    
end %for each trace







end %idlToDwt





