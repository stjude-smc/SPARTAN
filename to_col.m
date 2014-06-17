function vector = to_col(vector)
% Convert a vector to a column vector.

assert( size(vector,1)==1 | size(vector,2)==1, 'Not a vector' );
vector = reshape( vector, [numel(vector) 1] );

end

