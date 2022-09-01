function j =gumble_sample(p)
[~,j]    = max(p-log( -log( rand( size(p) ) ) ));
