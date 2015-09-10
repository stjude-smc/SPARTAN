function output = RLEncode( x )
% RLEncode   run-length encoding
% 
%   OUTPUT = RLEncode( INPUT )
%   produces a run-length encoded version of the 1xN vector INPUT. OUTPUT
%   is a Nx2 array, where the first column are the values of each run and
%   the second column are the lengths of each run.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


x = reshape( x, [1 numel(x)] );

i = [ find(x(1:end-1) ~= x(2:end)) length(x) ]; %run end positions
len = diff([ 0 i ]);
val = x(i);

output = [val' len'];

end %function RLEncode
