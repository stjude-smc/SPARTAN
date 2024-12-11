function demux_spots_internal(input, output, PLT, keep_plot)
% Select molecules with positions that are within printed spots of the
% field of view in multiplexed surface patterning experiments.
% This implementation uses edge detection and circle fitting for spot selection.
% 
% Parameters:
%   input: input file name
%   output: cell array with four output file names
%   PLT: Plotting behavior
%       - false: No plotting
%       - [] (default): Plot and save using input file basename + '_demux.png'
%       - string: Save plot to the specified file name
%   keep_plot: whether to keep or close the plot after plotting

    suffix = {'A', 'B', 'C', 'D'};

    % Load trace data and extract metadata
    data = loadTraces(input);
    [p, f, e] = fileparts(input);
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
        if ~keep_plot
            close(gcf());
        end
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
        saveTraces(output{j}, subset);
    end

end
