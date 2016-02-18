function filename = getFile( filter, prompt )
% GETFILE   User prompt for a filename
%
%   FILES = GETFILE( FILTER, PROMPT )
%   Asks the user to select a file, where the listed files are filtered
%   according to FILTER (for example, '*.txt'). The text in PROMPT is shown
%   in the title bar to the user. A file-open dialog will appear repeatedly
%   to allow the user to select multiple files, until the user pressed
%   "Cancel". The filenames (FILES) selected are returned in the cell array.

%   Copyright 2015 Cornell University All Rights Reserved.


persistent filterIndex;

filename = [];

if nargin<2,
    prompt = 'Select a file, hit cancel when finished:';
end
if nargin<1 || isempty(filter),
    filter = {'*.traces','Binary Traces Files (*.traces)'; ...
              '*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.txt','Text Files (*.txt)'; ...
              '*.dwt','Dwell-Time Files (*.dwt)'; ...
              '*.qmf','QuB Model Files (*.qmf)'; ...
              '*.*','All Files (*.*)'};
end

% Get file name from user.
[f,p,filterIndex] = uigetfile(filter, prompt, 'MultiSelect','off');

if f==0, return; end  %user hit cancel
filename = fullfile(p,f);

% re-order filter list to make the selection go to the top.
% (for convenience when selection many files).
% If the user enters a custom filter, we cannot get it so do nothing
% (this happens when filterIndex>number of defined filters).
if filterIndex<=size(filter,1) && filterIndex>0,
    ind = 1:size(filter,1);
    filter = [ filter(filterIndex,:); filter(ind~=filterIndex,:) ];
end
