function resizeTraces( finalLen, files, mode )
% resizeTraces   Change length of a .traces files
%
%   resizeTraces( TRACE_LEN )
%   prompts the user for files and resizes all traces to TRACE_LEN.
%
%   resizeTraces( TRACE_LEN, FILES )
%   resizes all traces to TRACE_LEN in the specified files.
%
%   If TRACE_LEN is empty (or the user hits cancel when asked), all traces will
%   be truncated to the same length. To extend instead:
%      resizeTraces( ..., 'extend' );
%
%   If traceLen > actual, the last data value is repeated at the end of every
%   trace.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Parse input arguments, asking the user for any that are missing.
if nargin<1,
    f = inputdlg('Final trace length (optional):');
    finalLen = str2double(f);
end

if nargin<2,
    filter = {'*.traces','Binary Traces Files (*.traces)'; ...
              '*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.*','All Files (*.*)' };
    files = getFiles(filter);
end
if ischar(files),
    files = {files};
end
if numel(files)==0, return; end

% If no final trace length given, choose a size so all traces match.
% By default truncate, or extend if the option is given.
if isempty(finalLen),
    if nargin<3 || strcmp(mode,'truncate')
        finalLen = min( sizeTraces(files,2) );
    elseif strcmp(mode,'extend')
        finalLen = max( sizeTraces(files,2) );
    else
        error('Invalid resize mode');
    end
end


% For each file in the user-selected directory
for i=1:numel(files),    
    % Read traces file
    data = loadTraces( files{i} );
    originalLen = data.nFrames;
    
    % Modify, if they do not match the target trace length, and save.
    if finalLen == originalLen,
        disp( ['Skipping ' files{i}] );
        continue;  %nothing to do
    elseif finalLen > originalLen,
        saveTraces(files{i}, data.extend(finalLen) );
    else
        saveTraces(files{i}, data.truncate(finalLen) );
    end
    
    fprintf('Resized %.0f to %.0f: %s\n',originalLen,finalLen,files{i});
end




