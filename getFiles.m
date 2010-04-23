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

files = {};

if nargin<2,
    prompt = 'Select a file, hit cancel when finished:';
end
if nargin<1,
    filter = {'*.txt','Text Files (*.txt)'; ...
              '*.txt;*.traces','All Traces Files (*.txt,*.traces)'; ...
              '*.traces','Binary Traces Files (*.traces)'; ...
              '*.dwt;*.traces','Dwell-Time Files (*.dwt)'; ...
              '*.*','All Files (*.*)'};
end

while 1,
    [f,p,filterIndex] = uigetfile(filter, prompt, 'MultiSelect','on');
    
    if iscell(f),
        files = [files strcat(p,f)];
    else
        if f==0, break; end  %user hit cancel
        files{end+1} = [p f];
    end
end
