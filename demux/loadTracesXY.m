function tracesXY = loadTracesXY(traces_files)
    % Check input arguments
    if nargin < 1
        filter = {'*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
                  '*.traces','Binary Traces Files (*.traces)';};
        traces_files = getFiles(filter);
    end

    % Ensure traces_files is a cell array
    if ischar(traces_files) || isstring(traces_files)
        traces_files = {traces_files};
    end
    
    % loadTraces from each file
    X = [];
    Y = [];
    nX = [];
    nY = [];
    tracesXY = struct('X', [], 'Y', [], 'nX', [], 'nY', []);

    for i = 1:numel(traces_files)
        data = loadTraces(traces_files{i});

        X = [X, [data.traceMetadata.donor_x]];
        Y = [Y, [data.traceMetadata.donor_y]];
    
        % Check if image dimensions are specified
        if isfield(data.fileMetadata, 'nX') && isfield(data.fileMetadata, 'nY')
            nX = [nX, data.fileMetadata.nX];
            nY = [nY, data.fileMetadata.nY];
        else
            error('File metadata must specify nX and nY for image dimensions.');
        end
    
        tracesXY = struct('X', X, 'Y', Y, ...
                          'nX', max(nX), 'nY', max(nY));
    end
end
