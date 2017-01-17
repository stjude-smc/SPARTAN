function titles = trimtitles(files)
%trimtitles Generate short figure titles for a cell array of file names
%
%   STR = trimtitles(FILES) strips the path and extension, and any text at the
%   beginning and end in common, from each elements in the cell array FILES
%   to generate short figure titles for figure legends.
%
%   See also: makeplots, frethistComparison, avgFretTime, dwellhist

%   Copyright 2016 Cornell University All Rights Reserved.


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
titles = strtrim( strrep(titles,'_',' ') );

if nFiles<2, return; end  %nothing more to do.


% Split each title into words, delimited by any symbol.
temp = cell(numel(titles), 0);
for i=1:nFiles,
    line = titles{i};
    
    idx = [1 find(~arrayfun(@isalnum,line)) numel(line)+1];
    for j=1:numel(idx)-1,
        temp{i,j} = line(idx(j):idx(j+1)-1);
    end
end

% Prevent errors with unique.
[temp{cellfun(@isempty,temp)}] = deal('');

% Remove identical words at the beginning.
while ~isempty(temp) && numel( unique(temp(:,1)) )==1
    temp(:,1) = [];
end

% Use numbers if all the files are identical.
if isempty(temp),
    titles = cellfun(@num2str, num2cell(1:4), 'UniformOutput',false);
    return;
end

% Remove identical words at end by right aligning cell array
e = sum( cellfun(@isempty,temp), 2);
for i=1:size(temp,1),
    temp(i,:) = circshift(temp(i,:), e(i), 2);
end
while numel( unique(temp(:,end)) )==1
    temp(:,end) = [];
end

% Recombine words to full strings for each file.
titles = cell(size(titles));
for i=1:numel(titles),
    titles{i} = strtrim( [temp{i,:}] );
end

% Insert placeholers for any empty elements.
e = cellfun(@isempty,titles);
[titles{e}] = deal('-');


end



function out = isalnum(str)
% Check if all characters are alpha-numeric or underscore.
out = all( ismember(upper(str),'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789') );
end




