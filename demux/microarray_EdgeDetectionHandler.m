classdef microarray_EdgeDetectionHandler < handle
    properties
        traces_xy = struct('X', [], 'Y', [], 'nX', [], 'nY', []);
        traces_files = [];  % Open traces files
        thumbnail;      % Downscaled and processed thumbnail
        edges;          % Canny edge detection result
        pixel_coords;   % Detected edge coordinates
        auto_update = false;
        ax_in;
        ax_out;
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
        function self = microarray_EdgeDetectionHandler(ax_in, ax_out)
            self.ax_in = ax_in;
            self.ax_out = ax_out;
        end

        % Enable/disable automatic update
        function enable_auto_update(self, active)
            self.auto_update = active;
            self.compute_edges();
        end

        function set_traces(self, traces_xy, traces_files)
            self.traces_xy = traces_xy;
            self.traces_files = traces_files;
            self.compute_edges();
        end

        function set_spots(self, spots)

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
                %msk = make_spot_mask(location_img, 160, 150); % Placeholder for mask generation
                %self.thumbnail = imresize(location_img .* msk, 1 / self.Downscale1, 'bilinear');
                self.thumbnail = imresize(location_img, 1 / self.params.Downscale1, 'bilinear');


                % Dilation
                self.thumbnail = imdilate(self.thumbnail, strel('disk', round(self.params.Dilation / self.params.Downscale1)));


                % Downscale again
                self.thumbnail = imresize(self.thumbnail, 1 / self.params.Downscale2, 'bilinear');
                self.thumbnail = self.thumbnail ./ max(self.thumbnail(:));
                
                self.display_thumbnail(self.ax_in);

                % Canny edge detection
                self.edges = edge(self.thumbnail, 'Canny', [self.params.LowThreshold, self.params.HighThreshold], self.params.Sigma);

                self.display_edges(self.ax_out);

                % Extract edge coordinates
                [row_coords, col_coords] = find(self.edges);
                self.pixel_coords = [row_coords, col_coords] * self.params.Downscale1 * self.params.Downscale2;
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
    end
end