function selectPrintedSpots(files, PLT)
% Select molecules with positions that are within printed spots of the
% field of view in multiplexed surface patterning experiments.
% This implementation uses edge detection and circle fitting for spot selection.
% 
% Parameters:
%   files: A string or cell array of input file names.
%   PLT: Plotting behavior
%       - false: No plotting
%       - [] (default): Plot and save using input file basename + '_demux.png'
%       - string: Save plot to the specified file name
%       - cell array: Save plots to corresponding file names for each input file.

% Parameters
suffix = {'A', 'B', 'C', 'D'}; % Add to split file names

% Default value for PLT
if nargin < 2
    PLT = [];
end

% Check input arguments
if nargin < 1
    files = getFiles('.rawtraces');
elseif ischar(files)
    files = {files};
elseif ~iscell(files)
    error('Input must be a string or cell array of strings');
end

% Validate PLT if it is a cell array
if iscell(PLT) && numel(PLT) ~= numel(files)
    error('PLT must have the same length as files if provided as a cell array');
end

for i = 1:numel(files)
    % Load trace data and extract metadata
    data = loadTraces(files{i});
    [p, f, e] = fileparts(files{i});
    x = [data.traceMetadata.donor_x];
    y = [data.traceMetadata.donor_y];

    % Check if image dimensions are specified
    if isfield(data.fileMetadata, 'nX') && isfield(data.fileMetadata, 'nY')
        nX = data.fileMetadata.nX;
        nY = data.fileMetadata.nY;
    else
        error('File metadata must specify nX and nY for image dimensions.');
    end

    % Determine whether to plot and save, based on PLT
    if isequal(PLT, false)
        ax = []; % No plotting
        plt_name = ''; % No output file for plotting
    elseif isempty(PLT)
        figure;
        ax = gca(); % Create a figure and get the axes
        plt_name = fullfile(p, [f, '_demux.png']); % Default plot file name
    elseif ischar(PLT)
        figure;
        ax = gca(); % Create a figure and get the axes
        plt_name = PLT; % Use the specified file name
    elseif iscell(PLT)
        figure;
        ax = gca(); % Create a figure and get the axes
        plt_name = PLT{i}; % Use the file name corresponding to the current file
    else
        error('Invalid value for PLT');
    end

    % Perform edge detection and circle fitting
    edge_coords = find_spot_edges(data, ax);
    circles = fit_circles(data, edge_coords, ax);

    % Save the plot if axes were created
    if ~isempty(ax)
        margin = 80;
        xlim([1 - margin, nY + margin]); % Horizontal limits
        ylim([1 - margin, nX + margin]); % Vertical limits
        title([f e], 'Interpreter', 'none');
        saveas(gcf(), plt_name); % Save the plot as an image
        close(gcf());
    end

    % Process each circle
    for j = 1:size(circles, 1)
        % Extract circle parameters
        cx = circles(j, 1);
        cy = circles(j, 2);
        r = circles(j, 3);
    
        % Compute distances from all points to the circle center
        distances = sqrt((x - cx).^2 + (y - cy).^2);
    
        % Create a boolean mask for points within the circle
        in_circle = distances <= r;
    
        % Extract the subset of traces for this circle
        subset = data.getSubset(in_circle);
    
        % Skip saving if the subset is empty
        if r == 0
            fprintf('Warning: no spot detected in quadrant %s of %s\n', suffix{j}, [f e]);
            continue;
        end
    
        % Save the subset to the corresponding output file
        outname = fullfile(p, [f, '_', suffix{j}, e]);
        saveTraces(outname, subset);
    end
end
end
