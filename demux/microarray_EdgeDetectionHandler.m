classdef microarray_EdgeDetectionHandler < handle
    properties
        traces;         % Input traces
        thumbnail;      % Downscaled and processed thumbnail
        edges;          % Canny edge detection result
        pixel_coords;   % Detected edge coordinates
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
        function self = EdgeDetectionHandler(traces)
            if nargin > 0
                self.traces = traces;
            end
        end

        % Compute edges
        function compute_edges(self)
            % Validate input traces
            if ~isfield(self.traces, 'traceMetadata') || ~isfield(self.traces, 'fileMetadata')
                error('Invalid input: `traces` must contain `traceMetadata` and `fileMetadata` fields.');
            end

            % Extract metadata
            x = [self.traces.traceMetadata.donor_x];
            y = [self.traces.traceMetadata.donor_y];
            nX = self.traces.fileMetadata.nX;
            nY = self.traces.fileMetadata.nY;

            % Initialize binary image
            location_img = zeros(nY, nX);
            for i = 1:length(x)
                location_img(y(i), x(i)) = 1;
            end

            % Mask and resize
            msk = make_spot_mask(location_img, 160, 150); % Placeholder for mask generation
            self.thumbnail = imresize(location_img .* msk, 1 / self.downscale_step1, 'bilinear');

            % Dilation
            self.thumbnail = imdilate(self.thumbnail, strel('disk', round(self.dilation_radius / self.downscale_step1)));

            % Downscale again
            self.thumbnail = imresize(self.thumbnail, 1 / self.downscale_step2, 'bilinear');
            self.thumbnail = self.thumbnail ./ max(self.thumbnail(:));

            % Canny edge detection
            self.edges = edge(self.thumbnail, 'Canny', [self.low_threshold, self.high_threshold], self.sigma);

            % Extract edge coordinates
            [row_coords, col_coords] = find(self.edges);
            self.pixel_coords = [row_coords, col_coords] * self.downscale_step1 * self.downscale_step2;
        end

        % Generic getter
        function value = get(self, propName)
            value = self.params.(propName);
        end

        % Generic setter
        function set(self, propName, value)
            disp("Setter called")
            self.params.(propName) = value;
            'Calling `compute_edges()`'
            % FIXME self.compute_edges(); % Recompute edges whenever a parameter changes
        end

        % Save settings
        function save_settings(self)
            try
                % FIXME edge_detect_settings = self.prepare_edge_detect_settings();
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
                title(ax, 'Thumbnail');
            end
        end

        % Display edges in axes
        function display_edges(self, ax)
            if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
                imshow(self.edges, 'Parent', ax);
                title(ax, 'Edges');
            end
        end
    end
end