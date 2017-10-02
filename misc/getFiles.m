function files = getFiles( filter, prompt )
% GETFILES   User prompt for filenames
%
%   FILES = GETFILES( FILTER, PROMPT )
%   Asks the user to select files, where the listed files are filtered
%   according to FILTER (for example, '*.txt'). The text in PROMPT is shown
%   in the title bar to the user. A file-open dialog will appear repeatedly
%   to allow the user to select multiple files, until the user pressed
%   "Cancel". The filenames (FILES) selected are returned in the cell array.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


persistent filterIndex;
persistent last;  %file list created on last call to getFiles.


% Return the previously requested file list if running 'getFiles -last'.
if nargin>=1 && ischar(filter) && strcmpi(filter,'-last')
    files=last;
    return;
end


%%
files = {};

if nargin<2,
    prompt = 'Select files, hit cancel when finished:';
end
if nargin<1 || isempty(filter),
    filter = {'*.traces','Binary Traces Files (*.traces)'; ...
              '*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.txt','Text Files (*.txt)'; ...
              '*.dwt','Dwell-Time Files (*.dwt)'; ...
              '*.qmf','QuB Model Files (*.qmf)'; ...
              '*.*','All Files (*.*)'};
end

while 1,
    % Get file name from user.
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
    if filterIndex<=size(filter,1) && filterIndex>0,
        ind = 1:size(filter,1);
        filter = [ filter(filterIndex,:); filter(ind~=filterIndex,:) ];
    end
end

last = files;