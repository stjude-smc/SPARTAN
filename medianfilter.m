function output = medianfilter(input,n)
%MEDFILT1D  One dimensional median filter.
%
% Inputs:
%   x     - row vector
%   n     - order of the filter

% Use fast compiled version if available
if exist('medianfilterx','file') == 3
    output = medianfilterx(input')';
    return;
end
    
nx = size(input,2);

% if rem(n,2)~=1    % n even
%     m = n/2;
% else
    m = (n-1)/2;
% end
assert(rem(n,2)==1, 'n must be odd!!');

output = zeros( size(input) );


% Create a matrix of indeces into the data as such:
% 1 2 3 4 5 ... nx
% 2 3 4 5 6 ... nx+1
% . . . . . . . .
% N N+1 . . ... nx+N-1
% 
% On the columns, these are the indexes of the windows
ind = repmat( 1:nx, n,1 ) + repmat( (0:n-1)', 1,nx ); 


for row=1:size(input,1)
    
    % Pad array with zeros so edges are defined
    X = [zeros(m,1); input(row,:)'; zeros(m,1)];

    % Do the median filtering
%     output(row,:) = median( X(ind), 1 );

    vals = sort( X(ind), 1 );
    output(row,:) = vals( floor(n/2)+1, 1:nx );
    %xx = X(ind);
%     output(row,:) = fastmedian( X, 9 );
end

% end function MEDFILT1D


