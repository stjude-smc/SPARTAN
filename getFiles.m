function files = getFiles( filter, prompt )
% GETFILES   User prompt for filenames
%
%   FILES = GETFILES( FILTER, PROMPT )
%   Asks the user to select files, where the listed files are filtered
%   according to FILTER (for example, '*.txt'). The text in PROMPT is shown
%   in the title bar to the user. A file-open dialog will appear repeatedly
%   to allow the user to select multiple files, until the user pressed
%   "Cancel". The filenames (FILES) selected are returned in the cell array.
%

persistent filterIndex;

files = {};

if nargin<2,
    prompt = 'Select a file, hit cancel when finished:';
end
if nargin<1 || isempty(filter),
    filter = {'*.txt;*.traces','All Traces Files (*.txt,*.traces)'; ...
              '*.traces','Binary Traces Files (*.traces)'; ...
              '*.txt','Text Files (*.txt)'; ...
              '*.dwt;*.traces','Dwell-Time Files (*.dwt)'; ...
              '*.*','All Files (*.*)'};
end

while 1,
    % re-order filter list to make the selection go to the top.
    % (for convenience when selection many files).
    % If the user enters a custom filter, we cannot get it so do nothing
    % (this happens when filterIndex>number of defined filters).
    if ~isempty(filterIndex) && filterIndex<=size(filter,1) && filterIndex>0,
        ind = 1:size(filter,1);
        filter = [ filter(filterIndex,:); filter(ind~=filterIndex,:) ];
    end
    
    [f,p,filterIndex] = uigetfile(filter, prompt, 'MultiSelect','on');
    
    if iscell(f),
        files = [files strcat(p,f)];
    else
        if f==0, break; end  %user hit cancel
        files{end+1} = [p f];
    end
    
    % re-order filter list to make the selection go to the top.
    % (for convenience when selection many files).
    % If the user enters a custom filter, we cannot get it so do nothing
    % (this happens when filterIndex>number of defined filters).
    if filterIndex<=size(filter,1),
        ind = 1:size(filter,1);
        filter = [ filter(filterIndex,:); filter(ind~=filterIndex,:) ];
    end
end
