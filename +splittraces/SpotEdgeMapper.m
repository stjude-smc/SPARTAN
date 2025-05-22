classdef SpotEdgeMapper < handle
    properties
        traces_xy = struct('X', [], 'Y', [], 'nX', [], 'nY', []);
        traces_files = [];  % Open traces files
        thumbnail;      % Downscaled and processed thumbnail
        edges;          % Canny edge detection result
        edge_coordinates;   % Detected edge coordinates
        auto_update = false;
        ax_in;
        ax_out;
        px_size;
        mask;
    end

    properties (Dependent)
        % Edge detection settings
        Downscale1;     % Downscale step 1
        Dilation;       % Dilation radius
        Downscale2;     % Downscale step 2
        LowThreshold;   % Lower threshold for Canny
        HighThreshold;  % Higher threshold for Canny
        Sigma;          % Sigma for Canny
        % Output settings
        CreateSubdir;   % Whether a new directory should be created for every input file
        CombineDatasets;  % Whether output traces should be combined into a single file
    end

    properties (Access = private)
        % Default values for edge detection settings, if GUI didn't provide them to constructor
        params = struct( ...
            'Downscale1', 4, ...
            'Dilation', 6, ...
            'Downscale2', 2, ...
            'LowThreshold', 0.05, ...
            'HighThreshold', 0.6, ...
            'Sigma', 7, ...
            'CreateSubdir', true, ...
            'PlotResult', true, ...
            'SavePlot', true, ...
            'CombineDatasets', true ...
        );
        FOVrectHandle;
        spot_layout;
    end

    methods (Access = public)
        % Constructor
        function self = SpotEdgeMapper(ax_in, ax_out, px_size, spot_layout)
            self.ax_in = ax_in;
            self.ax_out = ax_out;
            self.px_size = px_size;
            self.spot_layout = spot_layout;

            % Configure axes
            self.configure_axes(self.ax_out);
            self.create_FOV_rectangle(self.ax_out);
        end

        % loading trace XY coordinates
        function load_traces_XY(self, traces_files)
            % Ensure traces_files is a cell array
            if ischar(traces_files) || isstring(traces_files)
                traces_files = {traces_files};
            end

            self.traces_files = convert_to_char_cell(traces_files);
            
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
            
                self.traces_xy = struct('X', X, 'Y', Y, ...
                                  'nX', max(nX), 'nY', max(nY));
            end

            self.compute_edges(self.traces_xy);
        end

        % Enable/disable automatic update
        function enable_auto_update(self, active)
            self.auto_update = active;
            self.compute_edges(self.traces_xy);
        end

        function labeledImage = createLabeledImage(self)
            Spots = self.spot_layout.Spots;
            nX = self.traces_xy.nX;
            nY = self.traces_xy.nY;

            % Initialize the labeled image
            labeledImage = zeros(nY, nX);

            % Create a circular mask for each spot
            [X, Y] = meshgrid(1:nX, 1:nY); % Generate grid for image coordinates

            ids = Spots.keys;
            for k = 1:numel(ids)
                spot = Spots(ids{k});
                xCenter = spot.position(1) / self.px_size;
                yCenter = spot.position(2) / self.px_size;
                spotID = spot.id;
                spotSize = spot.size / self.px_size;

                % Create a circular mask for the current spot
                mask = (X - xCenter).^2 + (Y - yCenter).^2 <= (spotSize / 2)^2;

                % Fill the labeled image with the spot ID
                labeledImage(mask) = spotID;
            end

            % Downscale the labeled image to match the thumbnail and detected edges
            labeledImage = imresize(labeledImage, 1/(self.params.Downscale1*self.params.Downscale2), 'nearest');

        end

        % Compute edges
        function compute_edges(self, traces_xy)
            if nargin < 2
                traces_xy = self.traces_xy;
            end

            if self.auto_update & ~isempty(traces_xy.X)
                % Extract metadata
                x = traces_xy.X;
                y = traces_xy.Y;
                nX = traces_xy.nX;
                nY = traces_xy.nY;

                % Initialize binary image
                location_img = zeros(nY, nX);
                for i = 1:length(x)
                    location_img(y(i), x(i)) = 1;
                end

                % Mask and resize
                self.thumbnail = imresize(location_img, 1 / self.params.Downscale1, 'bilinear');


                % Dilation
                self.thumbnail = imdilate(self.thumbnail, strel('disk', round(self.params.Dilation / self.params.Downscale1)));

                % Downscale again
                self.thumbnail = imresize(self.thumbnail, 1 / self.params.Downscale2, 'bilinear');
                self.thumbnail = self.thumbnail ./ max(self.thumbnail(:));
                
                self.display_thumbnail(self.ax_in);

                % Canny edge detection
                self.edges = edge(self.thumbnail, 'Canny', [self.params.LowThreshold, self.params.HighThreshold], self.params.Sigma);

                % Labeled image acts as a mask. Everything outside circles
                % is set to 0, edges within circle get a circle number assigned
                labeled_edges = self.edges .* self.createLabeledImage();

                % Extract edge coordinates and save as cell array (one element per spot)
                self.edge_coordinates = extract_edge_coordinates(labeled_edges);

                % Rescale edge coordinates back into original image size
                scale_factor = self.params.Downscale1 * self.params.Downscale2;

                keys = self.edge_coordinates.keys();
                for i = 1:numel(keys)
                    key = keys{i};
                    self.edge_coordinates(key) = self.edge_coordinates(key) * scale_factor;
                end

                self.display_edges(self.edge_coordinates, self.ax_out);

                self.fit_circles(self.ax_out);
            end
        end

        % Generic getter
        function value = get(self, propName)
            value = self.params.(propName);
        end

        % Generic setter
        function set(self, propName, value)
            self.params.(propName) = value;
            self.compute_edges(self.traces_xy); % Recompute edges whenever a parameter changes
        end

        % Save settings
        function save_settings(self)
            try
                [file, path] = uiputfile('edge_detection_settings.json', 'Save edge detection settings');

                if ischar(file)
                    filepath = fullfile(path, file);
                    jsonStr = jsonencode(self.params);
                    fid = fopen(filepath, 'w');
                    fwrite(fid, jsonStr, 'char');
                    fclose(fid);
                end
            catch ME
                disp(['Error saving edge detection setting: ', ME.message]);
            end
        end

        % Load settings
        function load_settings(self, app)
            [file, path] = uigetfile('*.json', 'Load edge detection settings');

            if ischar(file)
                try
                    filepath = fullfile(path, file);
                    jsonStr = fileread(filepath);
                    self.params = jsondecode(jsonStr);

                    % update GUI inputs
                    app.ddDownscale1.Value = [num2str(self.params.Downscale1), 'x'];
                    app.sDilation.Value = self.params.Dilation;
                    app.ddDownscale2.Value = [num2str(self.params.Downscale2), 'x'];
                    app.sLowThreshold.Value = self.params.LowThreshold;
                    app.sHighThreshold.Value = self.params.HighThreshold;
                    app.sSigma.Value = self.params.Sigma;

                    self.compute_edges(self.traces_xy);

                catch ME
                    disp(['Error loading edge detection settings: ', ME.message]);
                end
            end
        end

        % Display thumbnail in axes
        function display_thumbnail(self, ax)
            if nargin < 2
                figure();
                ax = gca();
            end

            if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
                imshow(self.thumbnail, 'Parent', ax);
            end
        end

        % Display edges in axes
        function display_edges(self, edge_coordinates, ax)
            if nargin < 3
                figure();
                ax = gca();
            end

            % Check if the axes object is valid
            if isempty(ax) || ~isgraphics(ax, 'axes')
                error('Invalid axes object provided.');
            end

            cla(ax);
            self.configure_axes(self.ax_out);
            self.create_FOV_rectangle(self.ax_out);

            % Generate a colormap for different object IDs
            keys = edge_coordinates.keys();  % Get all keys
            num_objects = numel(keys);
            cmap = lines(num_objects);       % Use the 'lines' colormap for distinct colors

            hold(ax, 'on');  % Enable holding for multiple scatter plots

            % Loop through each object and plot its coordinates
            for i = 1:num_objects
                object_id = keys{i};                       % Get current object ID
                coordinates = edge_coordinates(object_id); % Get coordinates for current object

                % Scatter plot for the current object's edge coordinates
                scatter(ax, coordinates(:, 1) * self.px_size, ...
                        coordinates(:, 2) * self.px_size, ...
                        3, cmap(i, :), 'filled');
            end

            hold(ax, 'off');  % Release hold
        end


        function circles = fit_circles(self, ax, px_size)
            % Fit circles to edges of spots using RANSAC and optionally plot them.
            % Output:
            %   circles - Dictionary with keys as object IDs and values as [cx, cy, r]

            if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
                hold(ax, 'on');
                axis(ax, 'equal');
            end

            if nargin < 3  % px_size not provided
                px_size = self.px_size;
            end

            % Initialize the output dictionary
            keys = self.edge_coordinates.keys();  % Get all keys
            num_objects = numel(keys);
            circles = containers.Map('KeyType', 'int32', 'ValueType', 'any');

            % Loop through each object
            Spots = self.spot_layout.Spots;
            for i = 1:num_objects
                % Extract edge coordinates for the current object
                obj_id = keys{i};
                edges = self.edge_coordinates(obj_id);

                % Default values for center and radius
                margin = self.params.Downscale1 * self.params.Downscale2 / 2;
                cx = Spots(obj_id).position(1) / self.px_size;
                cy = Spots(obj_id).position(2) / self.px_size;
                r = 0.5 * Spots(obj_id).size / self.px_size;

                % Refine center and radius if have enough data points
                if size(edges, 1) >= 10
                    % Fit the circle using the RANSAC algorithm
                    [cx, cy, r] = fit_circle(edges, margin, cx, cy, r);
                end

                % Store the circle in the dictionary
                circles(obj_id) = [cx, cy, r];

                % Optional plotting if ax is provided
                if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
                    % Convert coordinates to µm
                    cx_plot = cx * px_size;
                    cy_plot = cy * px_size;
                    r_plot = r * px_size;

                    % Plot the circular patch
                    rectangle(ax, 'Position', [cx_plot - r_plot, cy_plot - r_plot, 2*r_plot, 2*r_plot], ...
                              'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 0.5);

                    % Plot the center
                    plot(ax, cx_plot, cy_plot, 'g+', 'MarkerSize', 10, 'LineWidth', 2, ...
                         'DisplayName', ['Object ', num2str(obj_id)]);
                end
            end
        end

        function split_traces(self)
            [split_files, spot_id_list, created_dirs] = self.split_dataset_traces();

            % If combine checkbox is checked, combine datasets per each spot
            if self.params.CombineDatasets
                self.combine_datasets(split_files, spot_id_list, created_dirs);
            end
        end

        function [split_files, spot_id_list, created_dirs] = split_dataset_traces(self)
            % Get all spot IDs from the layout (order preserved)
            spot_id_list = cell2mat(self.spot_layout.Spots.keys); % e.g., [2 4 5 9 ...]
            n_input = numel(self.traces_files);
            n_spots = numel(spot_id_list);

            % Preallocate output cell array: rows=input files, cols=spot IDs
            split_files = cell(n_input, n_spots);
            created_dirs = {};

            for i = 1:n_input
                fn = self.traces_files{i};
                [parent_path, name, extension] = fileparts(fn);

                % Optionally create subdirectory for this input file
                if self.get('CreateSubdir')
                    out_dir = fullfile(parent_path, name);
                    if ~exist(out_dir, 'dir')
                        mkdir(out_dir);
                    end
                    created_dirs{end+1} = out_dir;
                else
                    out_dir = parent_path;
                end

                % Plotting setup
                doPlot = self.get('PlotResult') || self.get('SavePlot');
                if doPlot
                    fig = figure();
                    ax = gca();
                    self.configure_axes(ax);
                    hold(ax, 'on');
                else
                    ax = [];
                end

                % Load traces file
                data = loadTraces(fn);

                % Extract coordinates for all points
                x = [data.traceMetadata.donor_x];
                y = [data.traceMetadata.donor_y];
                traces_xy = struct('X', x, 'Y', y, 'nX', data.fileMetadata.nX, 'nY', data.fileMetadata.nY);

                % Add FOV rectangle
                if ~isempty(ax)
                    self.create_FOV_rectangle(ax, 1);
                end

                % Compute edges and fit circles for this file
                self.compute_edges(traces_xy);
                circles = self.fit_circles(ax, 1);

                if doPlot
                    scatter(ax, x, y, 5, [0.4, 0.4, 0.4], 'filled');
                    xlabel(ax, 'x, px');
                    ylabel(ax, 'y, px');
                end

                % For each spot in the layout, try to split and save
                cmap = lines(n_spots);
                for j = 1:n_spots
                    spot_id = spot_id_list(j);
                    if ~isKey(self.spot_layout.Spots, spot_id)
                        split_files{i, j} = '';
                        continue;
                    end

                    % Only process if a circle was fitted for this spot in this file
                    if ~isKey(circles, spot_id)
                        split_files{i, j} = '';
                        continue;
                    end

                    circle = circles(spot_id);
                    cx = circle(1);
                    cy = circle(2);
                    r = circle(3) + 5; % Add margin

                    distances = sqrt((x - cx).^2 + (y - cy).^2);
                    in_spot = distances <= r;

                    if ~any(in_spot)
                        split_files{i, j} = '';
                        continue;
                    end

                    subset = data.getSubset(in_spot);

                    if self.get('CreateSubdir')
                        out_fn = fullfile(out_dir, sprintf('spot_%02d%s', spot_id, extension));
                    else
                        out_fn = fullfile(out_dir, sprintf('%s_%02d%s', name, spot_id, extension));
                    end

                    saveTraces(out_fn, subset);

                    split_files{i, j} = out_fn;

                    if doPlot
                        scatter(ax, [subset.traceMetadata.donor_x], [subset.traceMetadata.donor_y], 5, cmap(j,:), 'filled');
                    end
                end

                if doPlot
                    title(ax, name, 'Interpreter', 'none');
                    legend(ax, 'off');
                    hold(ax, 'off');

                    if self.get('SavePlot')
                        if self.get('CreateSubdir')
                            plot_fn = fullfile(out_dir, 'split_traces.png');
                        else
                            plot_fn = fullfile(parent_path, [name, '.png']);
                        end
                        saveas(fig, plot_fn);
                        disp(['Plot saved as: ', plot_fn]);
                        if ~self.get('PlotResult')
                            close(fig);
                        end
                    else
                        if ~self.get('PlotResult')
                            close(fig);
                        end
                    end
                end
            end
        end

        function combine_datasets(self, split_files, spot_id_list, created_dirs)
            % Prompt user for combined output prefix
            defaultPrefix = 'combined';
            % Flatten, remove empties, and find common directory
            all_files = split_files(:);
            all_files = all_files(~cellfun(@isempty, all_files));
            if isempty(all_files)
                disp('No files to combine.');
                return;
            end
            p = commonDir(all_files);

            % there is just one input file
            if (numel(self.traces_files) == 1 && numel(created_dirs) > 0)
                p = fileparts(p); % move one step up
            end


            fileFilter = {
                '*.rawtraces', 'Raw Traces files (*.rawtraces)';
                '*.traces',    'Traces files (*.traces)';
                '*.*',         'All Files (*.*)'
            };

            [filename, pathname] = uiputfile(fileFilter, ...
                'Select prefix to save combined traces', ...
                fullfile(p, defaultPrefix));
            if isequal(filename,0) || isequal(pathname,0)
                disp('User canceled file selection.');
                return;
            end
            prefix_fullpath = fullfile(pathname, filename);


            keep_files = questdlg('Do you want to keep intermediate files?', ...
                'Save intermediates?', ...
                'Keep', 'Delete', 'Keep');

            n_spots = numel(spot_id_list);
            n_input = size(split_files, 1);

            % Use the extension of the first non-empty file for each spot
            for j = 1:n_spots
                spot_id = spot_id_list(j);
                files_to_combine = split_files(:, j);
                files_to_combine = files_to_combine(~cellfun(@isempty, files_to_combine));
                if isempty(files_to_combine)
                    fprintf('No traces to combine for spot %d\n', spot_id);
                    continue;
                end

                % Determine extension to use
                [~, ~, ext] = fileparts(files_to_combine{1});
                [folder, base, ~] = fileparts(prefix_fullpath);
                outputPath = fullfile(folder, sprintf('%s_%02d%s', base, spot_id, ext));

                % Combine the files
                combineDatasets(files_to_combine, outputPath);

                disp(['Combined spot ' num2str(spot_id) ' -> ' outputPath]);
            end

            if strcmp(keep_files, 'Delete')
                % Flatten the cell array and remove empties
                files_to_delete = split_files(:);
                files_to_delete = files_to_delete(~cellfun(@isempty, files_to_delete));

                % Delete intermediate files
                for k = 1:numel(files_to_delete)
                    if exist(files_to_delete{k}, 'file')
                        delete(files_to_delete{k});
                    end
                end

                % Delete empty folders that contained intermediate files
                unique_dirs = unique(created_dirs); % Remove duplicates
                for i = 1:numel(unique_dirs)
                    d = unique_dirs{i};
                    if isfolder(d)
                        files = dir(d);
                        files = files(~ismember({files.name},{'.','..'}));
                        if isempty(files)
                            rmdir(d);
                            fprintf('Deleted empty folder: %s\n', d);
                        end
                    end
                end
            end

        end
    end

    methods (Access = private)
        % Configure axes
        function configure_axes(self, ax)
            box(ax, 'on');
            xlabel(ax, 'x, µm');
            ylabel(ax, 'y, µm');
            axis(ax, 'equal');
            set(ax, 'Color', 'k'); % Black background color
            set(ax, 'YDir', 'reverse');
        end

        % Create Field of View (FOV) rectangle
        function create_FOV_rectangle(self, ax, px)
            if nargin < 3
                px = self.px_size;
            end

            if isempty(self.traces_xy.X)
                nX = 250;
                nY = 250;
            else
                nX = self.traces_xy.nX * px;
                nY = self.traces_xy.nY * px;
            end

            rectangle(ax, ...
                'Position', [0, 0, nX, nY], ...
                'FaceColor', 'w', 'EdgeColor', 'none', ...
                'HitTest', 'off');
            axis(self.ax_out, 'tight');
        end
    end
end



% Helper functions

function edge_coordinates = extract_edge_coordinates(labeled_edges)
    % Find unique object IDs
    object_ids = unique(labeled_edges);
    object_ids(object_ids == 0) = [];  % Remove background (ID = 0)
    
    % Initialize a dictionary to hold coordinates for each object ID
    edge_coordinates = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    
    % Loop through each object ID
    for i = 1:numel(object_ids)
        object_id = object_ids(i);
        
        % Find coordinates of pixels belonging to the current object ID
        [Y, X] = find(labeled_edges == object_id);

        % Combine X and Y coordinates into a 2D array and store in the dictionary
        edge_coordinates(object_id) = [X, Y];  % Columns: [X, Y]
    end
end




function [cx, cy, r] = fit_circle(edges, margin, cx0, cy0, r0)
% Fit a circle to a set of edge points using RANSAC for robustness.
% Input:
%   edges - Nx2 array of edge coordinates [x, y]
%   cx0, cy0, r0 - starting values for center and radius, in px
% Output:
%   [cx, cy, r] - Fitted circle parameters (center and radius)

    % Restrictions on circle radius
    min_r = 0.75*r0;
    max_r = 1.25*r0;

    % Parameters for RANSAC
    max_iterations = 1000; % Maximum number of RANSAC iterations
    inlier_threshold = 4;  % Distance threshold to count a point as an inlier (in µm)
    min_inliers = 10;      % Minimum number of inliers for a valid model

    % Initialize best circle parameters and inlier count
    best_cx = 0;
    best_cy = 0;
    best_r = 0;
    max_inliers = 0;

    % Extract x and y coordinates
    x = edges(:, 1);
    y = edges(:, 2);

    % Set the random number generarator seed for reproducibility
    rng(7845);

    % RANSAC loop
    for i = 1:max_iterations
        % Randomly select 3 points (minimum required to define a circle)
        idx = randperm(size(edges, 1), 3);
        pts = edges(idx, :);

        % Fit a circle to the 3 points
        [cx_tmp, cy_tmp, r_tmp] = circle_from_three_points(pts);

        % Skip if the radius is invalid
        if isnan(r_tmp) || r_tmp <= 0
            continue;
        end

        % Compute distances of all points to the fitted circle
        distances = abs(sqrt((x - cx_tmp).^2 + (y - cy_tmp).^2) - r_tmp);

        % Count inliers within the threshold
        inliers = distances <= inlier_threshold;
        num_inliers = sum(inliers);

        % Update the best circle if the current one has more inliers...
        if num_inliers > max_inliers && num_inliers >= min_inliers
            % ... and has radius within the reasonable range
            if min_r < r_tmp && r_tmp < max_r
                best_cx = cx_tmp;
                best_cy = cy_tmp;
                best_r = r_tmp;
                max_inliers = num_inliers;
            end
        end
    end


    % Check how far away the new circle moved
    circle_shift = sqrt((cx0 - best_cx).^2 + (cy0 - best_cy).^2);

    % Check if a valid circle was found
    if max_inliers >= min_inliers && circle_shift < 0.8*r0
        cx = best_cx;
        cy = best_cy;
        r = best_r + margin;
    else
        % Return default values if no valid circle was found
        fprintf('Warning: RANSAC failed to find a valid circle.\n');
        cx = cx0;
        cy = cy0;
        r = r0;
    end
end

function [cx, cy, r] = circle_from_three_points(pts)
% Compute the circle passing through three points
% Input:
%   pts - 3x2 array of [x, y] coordinates
% Output:
%   [cx, cy, r] - Circle center and radius

    % Extract points
    x1 = pts(1, 1); y1 = pts(1, 2);
    x2 = pts(2, 1); y2 = pts(2, 2);
    x3 = pts(3, 1); y3 = pts(3, 2);

    % Compute the perpendicular bisectors of (x1, y1)-(x2, y2) and (x2, y2)-(x3, y3)
    A = [x2 - x1, y2 - y1; x3 - x2, y3 - y2];
    b = 0.5 * [(x2^2 - x1^2 + y2^2 - y1^2); (x3^2 - x2^2 + y3^2 - y2^2)];

    % Solve for the circle center
    if abs(det(A)) < 1e-10 % Check for degeneracy
        cx = NaN; cy = NaN; r = NaN;
        return;
    end
    center = A \ b;

    % Extract center coordinates
    cx = center(1);
    cy = center(2);

    % Compute the radius
    r = sqrt((x1 - cx)^2 + (y1 - cy)^2);
end



function char_cell = convert_to_char_cell(input)
    if ischar(input) || isstring(input)
        char_cell = cellstr(input); % Convert string or char to cell array of char vectors
    elseif iscell(input)
        char_cell = cellfun(@char, input, 'UniformOutput', false); % Convert each string in cell to char vector
    else
        error('Input must be a string, char, or cell array of strings.');
    end
end

function dirName = commonDir( files )
    % Find the common directory containing all give FILES.

    % Take out only path names from files
    nFiles = numel(files);
    for i=1:nFiles,
        files{i} = fileparts(files{i});
    end

    % Concatinate all pathnames into a single string matrix.
    files = char(files);

    % Find any differences
    diffs = zeros(1,size(files,2));

    for i=1:nFiles
        diffs = diffs | files(i,:)~=files(1,:);
    end

    % No differences - we got the folder!
    if ~any(diffs)
        dirName = files(1, :);
        return;
    end

    lastDiff = find(diffs);
    if isempty(lastDiff), lastDiff = size(files,2); end

    % Find directory name by going back to the last path seperating character.
    dirName = fileparts( files(1,1:lastDiff) );

end
