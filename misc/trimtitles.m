function titles = trimtitles(titles)
% Remove common elements at the beginning and end of sets of strings
%
%   STR = trimtitles(STR) removes any text at the beginning and end of all
%   elements in the cell array STR that are identical in all elements.
%   Useful for creating figure titles from long filenames.
%
%   See also: makeplots, frethistComparison, avgFretTime

%   Copyright 2016 Cornell University All Rights Reserved.

% FIXME: What about titles that end up empty because they are all common?

narginchk(1,1);
nargoutchk(0,1);

% Remove underscores that can be interpreted as subscripts.
titles = strrep(titles,'_',' ');

% Nothing more to do if there is only one title.
ntitles = numel(titles);
if ntitles<2, return; end

% Remove common characters from the beginning
titlecat = char(titles);
mask = sum(  titlecat - repmat(titlecat(1,:),[ntitles 1])  );
first = find(mask>0,  1, 'first');
titles = cellstr( titlecat(:,first:end) );

% Remove from the end by right justifying the text
titlecat = strjust( char(titles),'right' );
mask = sum(  titlecat - repmat(titlecat(1,:),[ntitles 1])  );
last = find(mask~=0, 1, 'last');
titles = cellstr( strjust(titlecat(:,1:last),'left') );

% Insert a placeholer for any empty elements
e = cellfun(@isempty, titles);
[titles{e}] = deal('-');

end



