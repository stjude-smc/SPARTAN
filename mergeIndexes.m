function list = mergeIndexes( starts, ends )
% MERGEINDEXES  Creates a list of indexes from ranges
%
%   values = mergeIndexes( STARTS, ENDS )
%   STARTS and ENDS are vectors of indexes (integers), where each pair of
%   START/END elements defines a range of space (just like START:END). This
%   function operates like the colon operator (eg, 1:10), creating a list of
%   numbers starting with each START value and ending at the END value. This
%   is repeated for each entry in the vectors, creating a linear sequence
%   of values.
%
%   values = squashIndexes( LIST )
%   does the same thing, where LIST is an Nx2 array where the first column are
%   the start values and the last column are the end points.
%
%   Example:
%       mergeIndexes( [1 5; 10 12; 50 53] )
%   returns
%       [ 1  2  3  4  5    10  11  12    50  51  52  53 ]
%


% Check input arguments
if nargin==1,
    assert( size(starts,2)==2 );
    ends = starts(:,2);
    starts = starts(:,1);
end

assert( exist('ends','var') && numel(starts)==numel(ends) && numel(starts)==max(size(starts)), ...
        'Malformed input arguments' );

% For each start/end pair, add the range of values to the list.
list = [];

for i=1:numel(starts),
    list = [list colon( starts(i), ends(i) )];
end

% Remove duplicates and sort the results.
% list = sort( unique(list) );

end %FUNCTION
