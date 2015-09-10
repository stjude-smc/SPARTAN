function varargout = squashIndexes( starts, ends )
% SQUASHINDEXES  Remove empty spans in a list of ranges
%
%   [starts,ends] = squashIndexes( STARTS, ENDS )
%   STARTS and ENDS are vectors of indexes (integers), where each pair of
%   START/END elements defines a range of space (just like START:END). This
%   function removes gaps in these ranges such that they represent a sequential
%   span in space.
%
%   [list] = squashIndexes( LIST )
%   does the same thing, where LIST is an Nx2 array where the first column are
%   the start values and the last column are the end points.
%
%   Example:
%       squashIndexes( [1 10; 11 20; 55 60] )
%   returns at sequence of start/end indexes without gaps:
%       [1 10; 11 20; 21 26]

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Check input arguments
if nargin==1,
    assert( size(starts,2)==2 );
    ends = starts(:,2);
    starts = starts(:,1);
end

assert( exist('ends','var') && numel(starts)==numel(ends) && numel(starts)==max(size(starts)), ...
        'Malformed input arguments' );
    
%
for i=1:numel(starts),
    if i==1,
        difference = starts(1)-1;
    else
        difference = (starts(i)-ends(i-1)-1);
    end
    starts(i:end) = starts(i:end)-difference;
    ends(i:end)   = ends(i:end)  -difference;
end

% Return result in the same form as input
if nargin==1 && nargout<=1,
    varargout{1} = [starts ends];
else
    varargout = {starts,ends};
end

end %FUNCTION
