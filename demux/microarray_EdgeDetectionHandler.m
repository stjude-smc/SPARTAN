classdef microarray_EdgeDetectionHandler < handle
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
        Spots = struct('patch', {}, 'position', {}, 'text', {}, 'id', {});
        mask;
    end

    properties (Dependent)
        Downscale1;     % Downscale step 1
        Dilation;       % Dilation radius
        Downscale2;     % Downscale step 2
        LowThreshold;   % Lower threshold for Canny
        HighThreshold;  % Higher threshold for Canny
        Sigma;          % Sigma for Canny
    end

    properties (Access = private)
        params = struct( ...
            'Downscale1', 4, ...
            'Dilation', 8, ...
            'Downscale2', 4, ...
            'LowThreshold', 0.05, ...
            'HighThreshold', 0.3, ...
            'Sigma', 5 ...
        );
    end

    methods
        % Constructor
        function self = microarray_EdgeDetectionHandler(ax_in, ax_out, px_size)
            self.ax_in = ax_in;
            self.ax_out = ax_out;
            self.px_size = px_size;
        end

        % loading trace XY coordinates
        function load_traces_XY(self, traces_files)
            % Ensure traces_files is a cell array
            if ischar(traces_files) || isstring(traces_files)
                traces_files = {traces_files};
            end

            self.traces_files = traces_files;
            
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

            self.compute_edges();
        end

        % Enable/disable automatic update
        function enable_auto_update(self, active)
            self.auto_update = active;
            self.compute_edges();
        end

        function set_spots(self, Spots)
            self.Spots = struct('id', {Spots.id}, 'position', {Spots.position}, 'size', {Spots.size});
        end

        function labeledImage = createLabeledImage(self)
            Spots = self.Spots;
            nX = self.traces_xy.nX;
            nY = self.traces_xy.nY;

            % Initialize the labeled image
            labeledImage = zeros(nY, nX);

            % Create a circular mask for each spot
            [X, Y] = meshgrid(1:nX, 1:nY); % Generate grid for image coordinates

            for i = 1:length(Spots)
                % Get spot center and ID
                xCenter = Spots(i).position(1) / self.px_size;
                yCenter = Spots(i).position(2) / self.px_size;
                spotID = Spots(i).id;
                spotSize = Spots(i).size / self.px_size;

                % Create a circular mask for the current spot
                mask = (X - xCenter).^2 + (Y - yCenter).^2 <= (spotSize / 2)^2;

                % Fill the labeled image with the spot ID
                labeledImage(mask) = spotID;
            end

            % Downscale the labeled image to match the thumbnail and detected edges
            labeledImage = imresize(labeledImage, 1/(self.params.Downscale1*self.params.Downscale2), 'nearest')

        end

        % Compute edges
        function compute_edges(self)
            if self.auto_update
                % Validate input traces
                if isempty(self.traces_xy.X)
                    error('Invalid input: `traces` is empty');
                end

                % Extract metadata
                x = self.traces_xy.X;
                y = self.traces_xy.Y;
                nX = self.traces_xy.nX;
                nY = self.traces_xy.nY;

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

                % Extract edge coordinates
                self.edge_coordinates = extract_edge_coordinates(labeled_edges);

                % Rescale edge coordinates back into original image size
                scale_factor = self.params.Downscale1 * self.params.Downscale2 * self.px_size;
                self.edge_coordinates = cellfun(@(coords) coords * scale_factor, self.edge_coordinates, 'UniformOutput', false);

                show_edge_coordinates(self.edge_coordinates, self.ax_out);

                %figure();
                fit_circles(self.edge_coordinates, self.ax_out);
            end
        end

        % Generic getter
        function value = get(self, propName)
            value = self.params.(propName);
        end

        % Generic setter
        function set(self, propName, value)
            self.params.(propName) = value;
            self.compute_edges(); % Recompute edges whenever a parameter changes
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

                    self.compute_edges();

                catch ME
                    disp(['Error loading edge detection settings: ', ME.message]);
                end
            end
        end

        % Display thumbnail in axes
        function display_thumbnail(self, ax)
            if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
                imshow(self.thumbnail, 'Parent', ax);
            end
        end

        % Display edges in axes
        function display_edges(self, ax)
            if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
                imshow(self.edges, 'Parent', ax);
            end
        end

        % Run demultiplexing
        function demux(self)
            % TODO
            edge_coords = find_spot_edges(data, ax);
            circles = fit_circles(data, edge_coords, ax);
        end

    end
end



% Helper functions

function edge_coordinates = extract_edge_coordinates(labeled_edges)
    % Find unique object IDs
    object_ids = unique(labeled_edges);
    object_ids(object_ids == 0) = [];  % Remove background (ID = 0)
    
    % Initialize a cell array to hold coordinates for each object ID
    edge_coordinates = cell(numel(object_ids), 1);
    
    % Loop through each object ID
    for i = 1:numel(object_ids)
        object_id = object_ids(i);
        
        % Find coordinates of pixels belonging to the current object ID
        [Y, X] = find(labeled_edges == object_id);

        % Combine X and Y coordinates into a 2D array
        edge_coordinates{i} = [X, Y];  % Columns: [X, Y]
    end
end

function show_edge_coordinates(edge_coordinates, ax)
    if nargin < 2
        figure();
        ax = gca();
    end

    % Check if the axes object is valid
    if isempty(ax) || ~isgraphics(ax, 'axes')
        error('Invalid axes object provided.');
    end
    
    % Generate a colormap for different object IDs
    num_objects = numel(edge_coordinates);
    cmap = lines(num_objects);  % Use the 'lines' colormap for distinct colors
    
    % Clear the axes before plotting
    cla(ax);
    hold(ax, 'on');  % Enable holding for multiple scatter plots
    
    % Loop through each object and plot its coordinates
    for object_id = 1:num_objects
        coordinates = edge_coordinates{object_id};  % Get coordinates for current object
        
        % Scatter plot for the current object's edge coordinates
        scatter(ax, coordinates(:, 1), coordinates(:, 2), 3, cmap(object_id, :), 'filled');
    end
    
    % Add labels and title
    xlabel(ax, 'x / µm');
    ylabel(ax, 'y / µm');
    title(ax, 'Detected edges');
    legend(ax, 'off');
    
    hold(ax, 'off');  % Release hold
end


function circles = fit_circles(edge_coords, ax)
    % Fit circles to edges of arbitrary objects using RANSAC and optionally plot them.
    % Input:
    %   traces       - Struct containing image metadata (e.g., dimensions)
    %   edge_coords  - Cell array where each cell is Nx2 array of edge pixel coordinates
    %   ax           - Optional axes object for plotting
    % Output:
    %   circles      - Mx3 array of [cx, cy, r] for the M circles (one for each object)

    if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
        hold(ax, 'on');
        axis(ax, 'equal');
    end

    % Initialize the output array
    num_objects = numel(edge_coords); % Number of objects
    circles = zeros(num_objects, 3); % To store [cx, cy, r] for each circle

    % Loop through each object
    for obj_id = 1:num_objects
        % Extract edge coordinates for the current object
        edges = edge_coords{obj_id};

        % Check if edges are empty
        if isempty(edges)
            % Return a zero-sized circle
            circles(obj_id, :) = [0, 0, 0];
            continue;
        end

        % Fit the circle using the RANSAC algorithm
        [cx, cy, r] = fit_circle(edges);
        circles(obj_id, :) = [cx, cy, r];

        % Optional plotting if ax is provided
        if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
            % Plot the circular patch
            rectangle(ax, 'Position', [cx - r, cy - r, 2*r, 2*r], ...
                      'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 0.5);

            % Plot the center
            plot(ax, cx, cy, 'g+', 'MarkerSize', 10, 'LineWidth', 2, ...
                 'DisplayName', ['Object ', num2str(obj_id)]);
        end
    end

    % Add legend if plotting
    if nargin > 2 && ~isempty(ax) && isgraphics(ax, 'axes')
        legend(ax, 'show');
        hold(ax, 'off');
    end
end



function [cx, cy, r] = fit_circle(edges)
% Fit a circle to a set of edge points using RANSAC for robustness.
% Input:
%   edges - Nx2 array of edge coordinates [x, y]
% Output:
%   [cx, cy, r] - Fitted circle parameters (center and radius)

    % Restrictions on circle radius %FIXME
    min_r = 20;  % µm
    max_r = 120; % µm

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

    % Check if a valid circle was found
    if max_inliers >= min_inliers
        cx = best_cx;
        cy = best_cy;
        r = best_r;
    else
        % Return default values if no valid circle was found
        fprintf('Warning: RANSAC failed to find a valid circle.\n');
        cx = 0;
        cy = 0;
        r = 0;
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
