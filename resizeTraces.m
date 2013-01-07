function resizeTraces( traceLen, files )
% resizeTraces   Change length of a .traces files
%
%   resizeTraces( TRACE_LEN )
%   loads all .traces files in the directory specified by the user,
%   resizes them to the TRACE_LEN, saves them back to their original
%   filenames. This ensures that all files will be the same length.
%
%   If traceLen < actual, the traces will be cropped.
%   If traceLen > actual, the last data value will be appended
% 

if nargin<1,
    f = inputdlg('How many frames to keep?');
    if isempty(f), return; end
    traceLen = str2double(f);
end

if nargin<2,
    filter = {'*.traces;*.rawtraces','Binary Traces Files (*.traces)'; ...
              '*.*','All Files (*.*)'};
    files = getFiles(filter);
end

constants = cascadeConstants();


% For each file in the user-selected directory
for i=1:numel(files),
    
    % ---- Read traces file
    filename = files{i};
    data = loadTraces( filename );
    [nTraces,actualLen] = size( data.donor );
    
    % ---- Reconstruct time axis (assuming linear, contiguous acquisition)
    dt = data.time(2)-data.time(1);
    data.time = cumsum( [data.time(1) repmat(dt,1,traceLen-1)] );
    
    % ---- Undo crosstalk correction (otherwise it will be done twice).
    % This only applies to old format files!
    if isfield(data,'traceMetaData') && ~isempty(data.traceMetadata),
        data.acceptor = data.acceptor + constants.crosstalk*data.donor;
    end
    
    % ---- Modify traces, if they do not match the target trace length.
    if traceLen == actualLen,
        %no changes needed
        continue;
    
    elseif traceLen > actualLen,
        % Expand traces using last value
        delta = traceLen-actualLen;
        data.donor    = [data.donor    repmat( data.donor(:,end), 1, delta )];
        data.acceptor = [data.acceptor repmat( data.acceptor(:,end), 1, delta )];
        data.fret     = [data.fret     zeros( nTraces, delta )];
    
    else
        % Truncate traces
        data.donor    = data.donor(:,1:traceLen);
        data.acceptor = data.acceptor(:,1:traceLen);
        data.fret     = data.fret(:,1:traceLen);
    end
    
    disp( sprintf('Resizing %.0f to %.0f: %s',actualLen,traceLen,filename) );
    
    % ---- Save the results.
    saveTraces( filename,'traces',data );
    
    
end

disp('Finished.');



