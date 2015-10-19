function files = getFileGroups( filter, prompt )
% GETFILES   User prompt for filenames in groups
%
%   FILES = GETFILES( FILTER, PROMPT )
%   Asks the user to select files, where the listed files are filtered
%   according to FILTER (for example, '*.traces'). The text in PROMPT is 
%   shown in the title bar to the user. A file-open dialog will appear
%   repeatedly to allow the user to select multiple files, until the user
%   pressed "Cancel".
%   
%   FILES is a cell array with each element containing all the files
%   selected each time the user clicked "Open" in the dialog. This allows
%   groups of files within each directory to be explicitly grouped in the
%   output. This differs from getFiles, where all selected files are
%   returned as one large group.
%
%   See also: getFiles, getFile.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


persistent filterIndex;

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
        files{end+1} = strcat(p,f);
    else
        if f==0, break; end  %user hit cancel
        files{end+1} = {[p f]};
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
