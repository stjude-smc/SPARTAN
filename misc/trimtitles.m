function titles = trimtitles(files)
%trimtitles Generate short figure titles for a cell array of file names
%
%   STR = trimtitles(FILES) strips the path and extension, and any text at the
%   beginning and end in common, from each elements in the cell array FILES
%   to generate short figure titles for figure legends.
%
%   See also: makeplots, frethistComparison, avgFretTime, dwellhist

%   Copyright 2016 Cornell University All Rights Reserved.

% FIXME: only chop at word bounderies.


% Verify input arguments
narginchk(1,1);
nargoutchk(0,1);

if ischar(files),
    files={files};
elseif ~iscell(files) || ~all(cellfun(@ischar,files))
    error('Invalid input. Must be a cell array of strings');
end
nFiles = numel(files);


% Extract file names and remove underscores (interpreted as subscripts).
[~,titles] = cellfun(@fileparts, files, 'UniformOutput',false);
titles = strrep(titles,'_',' ');

if nFiles<2, return; end  %nothing more to do.


% Remove common characters from the beginning
titlecat = char(titles);
mask = sum(  abs( titlecat - repmat(titlecat(1,:),[nFiles 1]) )  );
first = find(mask~=0,  1, 'first');
if ~isempty(first),
    titles = cellstr( titlecat(:,first:end) );
else
    % All file names are identical; use numbers instead.
    titles = cellfun(@num2str,num2cell(1:nFiles),'UniformOutput',false);
    return;
end


% Remove from the end by right justifying the text
titlecat = strjust( char(titles),'right' );
mask = sum(  abs( titlecat - repmat(titlecat(1,:),[nFiles 1]) )  );
last = find(mask~=0, 1, 'last');
titles = cellstr( strjust(titlecat(:,1:last),'left') );


% Insert placeholers for any empty elements.
e = cellfun(@isempty,titles);
[titles{e}] = deal('-');


end



