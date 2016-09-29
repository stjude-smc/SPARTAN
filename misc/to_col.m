function varargout = to_col(varargin)
% Convert a vector to a column vector.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% assert( size(vector,1)<=1 | size(vector,2)<=1, 'Not a vector' );

for i=1:nargin,
    varargout{i} = reshape( varargin{i}, [numel(varargin{i}) 1] );  %#ok<AGROW>
end

end

