function tracesXY = loadTracesXY(traces_files)
    % Camera pixel size
    cam_px_size   = 6.5; % µm
    magnification = 60;
    binning       = 2;
    px_size = cam_px_size*binning/magnification;

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
    x = [];
    y = [];
    nx = [];
    ny = [];
    for i = 1:numel(traces_files)
        data = loadTraces(traces_files{i});

        x = [x, [data.traceMetadata.donor_x] * px_size];
        y = [y, [data.traceMetadata.donor_y] * px_size];
    
        % Check if image dimensions are specified
        if isfield(data.fileMetadata, 'nX') && isfield(data.fileMetadata, 'nY')
            nx = [nx, data.fileMetadata.nX * px_size];
            ny = [ny, data.fileMetadata.nY * px_size];
        else
            error('File metadata must specify nX and nY for image dimensions.');
        end
    
        tracesXY = struct('x', x, 'y', y, 'nx', max(nx), 'ny', max(ny));
    end

end
