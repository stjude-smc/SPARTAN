function files = getFiles( filter, prompt )
% GETFILES   User prompt for getting filenames
%


files = {};

if nargin<2,
    prompt = 'Select a file:';
end
if nargin<1,
    filter = '*.txt';
end

while 1,
    [f,p] = uigetfile(filter, prompt, 'MultiSelect','on');
    
    if iscell(f),
        files = [files strcat(p,f)];
    else
        if f==0, break; end  %user hit cancel
        files{end+1} = [p f];
    end
end
