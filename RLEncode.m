function output = RLEncode( x )

x = reshape( x, [1 numel(x)] );

i = [ find(x(1:end-1) ~= x(2:end)) length(x) ]; %run end positions
len = diff([ 0 i ]);
val = x(i);

output = [val' len'];

end %function RLEncode
