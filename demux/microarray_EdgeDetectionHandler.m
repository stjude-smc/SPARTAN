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

        function fit_circles(self, data, edge_coords)
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
    end
end