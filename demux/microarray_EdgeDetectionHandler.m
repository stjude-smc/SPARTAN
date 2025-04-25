classdef microarray_EdgeDetectionHandler < handle
    properties
        traces;         % Input traces
        thumbnail;      % Downscaled and processed thumbnail
        edges;          % Canny edge detection result
        pixel_coords;   % Detected edge coordinates

        % Edge detection parameters
        downscale_step1 = 4; % Downscale step 1
        dilation_radius = 2; % Dilation radius
        downscale_step2 = 4; % Downscale step 2
        low_threshold = 0.05; % Lower threshold for Canny
        high_threshold = 0.3; % Higher threshold for Canny
        sigma = 5; % Sigma for Canny
    end

    methods (Access = public)
        % Constructor
        function self = microarray_EdgeDetectionHandler(traces)
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
            %FIXME msk = make_spot_mask(location_img, 160, 150); % Placeholder for mask generation
            %self.thumbnail = imresize(location_img .* msk, 1 / self.downscale_step1, 'bilinear');
            self.thumbnail = imresize(location_img, 1 / self.downscale_step1, 'bilinear');

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

        % Save settings
        function save_settings(self, filepath)
            settings = struct('downscale_step1', self.downscale_step1, ...
                              'dilation_radius', self.dilation_radius, ...
                              'downscale_step2', self.downscale_step2, ...
                              'low_threshold', self.low_threshold, ...
                              'high_threshold', self.high_threshold, ...
                              'sigma', self.sigma);
            jsonStr = jsonencode(settings);
            fid = fopen(filepath, 'w');
            fwrite(fid, jsonStr, 'char');
            fclose(fid);
        end

        % Load settings
        function load_settings(self, filepath)
            jsonStr = fileread(filepath);
            settings = jsondecode(jsonStr);
            self.downscale_step1 = settings.downscale_step1;
            self.dilation_radius = settings.dilation_radius;
            self.downscale_step2 = settings.downscale_step2;
            self.low_threshold = settings.low_threshold;
            self.high_threshold = settings.high_threshold;
            self.sigma = settings.sigma;

            % Recompute edges with new settings
            self.compute_edges();
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