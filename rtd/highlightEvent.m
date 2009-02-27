function highlightEvent( starts, ends, limits, color )

N = length(starts);

for i=1:N,
    
    X = [starts(i) starts(i) ends(i) ends(i)];
    Y = [limits(1) limits(end) limits(end) limits(1)];
    
    patch( X, Y, color, 'EdgeColor',color );
    
end
