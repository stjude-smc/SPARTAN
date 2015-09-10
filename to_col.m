function vector = to_col(vector)
% Convert a vector to a column vector.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


assert( size(vector,1)==1 | size(vector,2)==1, 'Not a vector' );
vector = reshape( vector, [numel(vector) 1] );

end

