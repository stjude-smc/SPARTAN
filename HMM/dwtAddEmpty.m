function [dwtOut,offsetsOut] = dwtAddEmpty( dwtIn, offsetsIn, nFrames, nTraces )
% dwtAddEmpty  re-insert empty elements into dwell-time cell array.
%
%   

narginchk(4,4);
nargoutchk(0,2);

dwtOut = cell( nTraces, 1 );
idxTrace = (offsetsIn/nFrames)+1;

for i=1:numel(idxTrace)
    dwtOut{ idxTrace(i) } = dwtIn{i};
end
offsetsOut = ((1:nTraces)-1)*nFrames;

end