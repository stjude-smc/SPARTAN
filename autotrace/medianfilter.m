function output = medianfilter(input,n)
%MEDIANFILTER  One dimensional median filter.
%
%    OUTPUT = MEDIANFILTER( INPUT, TAU )
%    Median filters the INPUT data over rows with the window size TAU.
%    This type of filter is good at smoothing noise while maintaining sharp
%    edges when there are large changes in the values. This is primarily
%    used in calcLifetime/traceStat for detecting photobleaching events.
%
%    NOTE: this is fairly slow, so I made a compiled C (MEX) version with
%    the same name. If binary folder is first in the path (it should be),
%    the faster version will take precedence.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Verify input arguments.
assert( nargin==2, 'Invalid input arguments' );
assert( rem(n,2)==1, 'Window size must be odd!!' );

try
    if any(isnan(input(:))),
        % FIXME: add support for NaN in mex code?
        %warning('NaN values not supported');
        input(isnan(input)) = 0;
    end
    output = medianfilterx(double(input),n);
    return;
catch
    %disp('Falling back on the matlab version');
end

nx = size(input,2);  %trace length
m = (n-1)/2;         %half the window size (padding size)

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
    % Pad array with data so edges are defined
    X = [repmat(input(row,1),m,1); input(row,:)'; repmat(input(row,end),m,1)];

    % For each window, sort the values and pick the middle (median) one.
    vals = sort( X(ind), 1 );
    output(row,:) = vals( floor(n/2)+1, 1:nx );
end

% end function MEDIANFILTER


