
function vector = to_row(vector)
% Convert a vector to a column vector.

assert( size(vector,1)==1 | size(vector,2)==1, 'Not a vector' );
vector = reshape( vector, [1 numel(vector)] );

end
