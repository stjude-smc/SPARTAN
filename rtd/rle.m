function [starts,ends] = rle( x, v )
% Returns a run-length encoded version fo the input signal,
% which is assumed to be a vector.
% 

% Create a run-length-encoded version of signal
ends = find(x(1:end-1) ~= x(2:end)); % length(x) ]; %run end positions
len = diff([ 0 ; ends ]);

starts = [1 ; cumsum(len)+1];  % start positions of all runs
starts = starts( x(starts)==v );  %select only positive runs

ends = [ends ; length(x)]; %add terminating step
ends = ends( x(ends)==v );  %select only positive runs

end





