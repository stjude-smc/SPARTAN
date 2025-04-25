classdef microarrayDesigner_ArrayLayout < handle
    properties (Access = public)
        traces_X = [];
        traces_Y = [];
        nX = 0;
        nY = 0;
        px_size;
        ax              matlab.ui.control.UIAxes;
        spot_size = 0;
    end

    properties (Access = private)
        Spots = struct('patch', {}, 'position', {}, 'text', {}, 'id', {}); % Circle data
        AvailableIDs = []; % Pool of reusable IDs
        NextID = 1; % ID for the next circle

        FOVrectHandle;
        scatterHandle;
    end

    methods (Access = public)
        % Constructor
        function self = microarrayDesigner_ArrayLayout(ax, spot_size, px_size)
            self.px_size = px_size;
            self.spot_size = self.validate_spot_size(spot_size);
            self.ax = ax;

            % Configure axes
            self.configure_axes();
            self.create_FOV_rectangle();
        end

        % Draw traces and manage Z-order
        function draw_traces(self, traces_xy)
            if isempty(traces_xy.X)
                return
            end

            hold(self.ax, 'on'); % Preserve existing objects
            self.update_FOV_rectangle(traces_xy);
            self.update_scatter_plot(traces_xy);

            % Arrange objects in the correct Z-order
            self.adjust_z_order();
        end

        % Add a new spot (with optional ID)
        function add_spot(self, x, y, id)
            if nargin < 4  % If id is not provided
                id = [];  % Default to an empty value
            end

            id = self.determine_spot_id(id);  % Determine the ID based on availability
            [xCircle, yCircle] = self.compute_circle_coordinates(x, y);

            % Create the circle patch
            patchHandle = self.create_patch(xCircle, yCircle, id);

            % Add text label for the circle
            textHandle = self.create_text(x, y, id);

            % Store circle data
            self.store_spot_data(patchHandle, textHandle, [x, y], id);
        end

        % Update spot size
        function update_spot_size(self, new_size)
            new_size = self.validate_spot_size(new_size);
            self.spot_size = new_size;

            for i = 1:length(self.Spots)
                pt = self.Spots(i).position;
                [xCircle, yCircle] = self.compute_circle_coordinates(pt(1), pt(2));
                set(self.Spots(i).patch, 'XData', xCircle, 'YData', yCircle);
            end
        end

        % Load layout
        function spot_size = load_layout(self)
            spot_size = 'none';
            [file, path] = uigetfile('*.json', 'Load microarray layout');

            if ischar(file)
                try
                    layoutData = self.read_layout_file(fullfile(path, file));
                    spot_size = self.process_loaded_layout(layoutData);
                catch ME
                    disp(['Error loading layout: ', ME.message]);
                end
            end
        end

        % Save layout
        function save_layout(self)
            try
                layoutData = self.prepare_layout_data();
                [file, path] = uiputfile('microarray_layout.json', 'Save microarray layout');

                if ischar(file)
                    self.write_layout_file(fullfile(path, file), layoutData);
                end
            catch ME
                disp(['Error saving layout: ', ME.message]);
            end
        end
    end

    methods (Access = private)
        % Configure axes
        function configure_axes(self)
            box(self.ax, 'on');
            xlabel(self.ax, 'x, µm');
            ylabel(self.ax, 'y, µm');
            axis(self.ax, 'equal');
            set(self.ax, 'Color', 'k'); % Black background color
        end

        % Create Field of View (FOV) rectangle
        function create_FOV_rectangle(self)
            self.FOVrectHandle = rectangle(self.ax, ...
                'Position', [0, 0, 250, 250], ...
                'FaceColor', 'w', 'EdgeColor', 'none', ...
                'HitTest', 'off');
            axis(self.ax, 'tight');
        end

        % Update FOV rectangle
        function update_FOV_rectangle(self, traces_xy)
            px = self.px_size;

            if isempty(self.FOVrectHandle) || ~isvalid(self.FOVrectHandle)
                self.FOVrectHandle = rectangle(self.ax, 'Position', [0, 0, traces_xy.nX * px, traces_xy.nY * px], ...
                                               'FaceColor', 'w', 'EdgeColor', 'none', ...
                                               'HitTest', 'off');
            else
                set(self.FOVrectHandle, 'Position', [0, 0, traces_xy.nX * px, traces_xy.nY * px]);
            end
        end

        % Update scatter plot
        function update_scatter_plot(self, traces_xy)
            px = self.px_size;

            if isempty(self.scatterHandle) || ~isvalid(self.scatterHandle)
                self.scatterHandle = scatter(self.ax, traces_xy.X * px, traces_xy.Y * px, 3, 'k', ...
                                             'filled', ...
                                             'HitTest', 'off');
            else
                set(self.scatterHandle, 'XData', traces_xy.X * px, 'YData', traces_xy.Y * px);
            end
        end

        % Adjust Z-order
        function adjust_z_order(self)
            uistack(self.scatterHandle, 'bottom');
            uistack(self.FOVrectHandle, 'bottom');

            for i = 1:length(self.Spots)
                if isvalid(self.Spots(i).patch)
                    uistack(self.Spots(i).patch, 'top');
                end
                if isvalid(self.Spots(i).text)
                    uistack(self.Spots(i).text, 'top');
                end
            end
        end



        % Validate spot size
        function spot_size = validate_spot_size(self, spot_size)
            if spot_size <= 0
                error('Spot size must be positive and non-zero.');
            end
        end

        % Compute circle coordinates
        function [xCircle, yCircle] = compute_circle_coordinates(self, x, y)
            theta = linspace(0, 2 * pi, 100);
            xCircle = x + self.spot_size * cos(theta) / 2;
            yCircle = y + self.spot_size * sin(theta) / 2;
        end

        % Determine spot ID
        function id = determine_spot_id(self, id)
            if nargin < 2 || isempty(id)
                if ~isempty(self.AvailableIDs)
                    id = self.AvailableIDs(1);
                    self.AvailableIDs(1) = [];
                else
                    id = self.NextID;
                    self.NextID = self.NextID + 1;
                end
            end
        end

        % Create patch
        function patchHandle = create_patch(self, xCircle, yCircle, id)
            patchHandle = patch('XData', xCircle, 'YData', yCircle, ...
                                'FaceColor', [0.6 0.7 1.0], 'EdgeColor', [0 0 0.8], ...
                                'FaceAlpha', 0.4, ...
                                'Parent', self.ax, ...
                                'ButtonDownFcn', @(src, event) self.handle_spot_click(src, event));
        end

        % Create text
        function textHandle = create_text(self, x, y, id)
            textHandle = text(self.ax, x, y, num2str(id), ...
                              'HorizontalAlignment', 'center', ...
                              'VerticalAlignment', 'middle', ...
                              'Color', 'k', 'FontSize', 25, ...
                              'FontWeight', 'bold', ...
                              'ButtonDownFcn', @(src, event) self.handle_spot_click(src, event));
        end

        % Store spot data
        function store_spot_data(self, patchHandle, textHandle, position, id)
            self.Spots(end+1).patch = patchHandle;
            self.Spots(end).position = position;
            self.Spots(end).text = textHandle;
            self.Spots(end).id = id;
        end

        % Read layout file
        function layoutData = read_layout_file(self, filepath)
            jsonStr = fileread(filepath);
            layoutData = jsondecode(jsonStr);
        end

        % Process loaded layout
        function spot_size = process_loaded_layout(self, layoutData)
            spot_size = layoutData.spot_size;

            % Delete existing spots
            for i = 1:length(self.Spots)
                delete(self.Spots(i).patch);
                delete(self.Spots(i).text);
            end
            self.Spots = struct('patch', {}, 'position', {}, 'text', {}, 'id', {}); % Circle data

            % Reconstruct the spots
            for i = 1:length(layoutData.Spots)
                x = layoutData.Spots(i).position(1);
                y = layoutData.Spots(i).position(2);
                id = layoutData.Spots(i).id;
                self.add_spot(x, y, id);
            end

            self.update_spot_size(spot_size);
        end

        % Prepare layout data for saving
        function layoutData = prepare_layout_data(self)
            roundedPositions = arrayfun(@(spot) round(spot.position), self.Spots, 'UniformOutput', false);
            layoutData.Spots = struct('id', {self.Spots.id}, 'position', roundedPositions);
            layoutData.spot_size = self.spot_size;
        end

        % Write layout file
        function write_layout_file(self, filepath, layoutData)
            jsonStr = jsonencode(layoutData);
            fid = fopen(filepath, 'w');
            fwrite(fid, jsonStr, 'char');
            fclose(fid);
        end

        % Handle click on a spot
        function handle_spot_click(self, src, event)
            parentFigure = ancestor(self.ax, 'figure');
            modifiers = get(parentFigure, 'CurrentModifier');
            isCtrl = ismember('control', modifiers);

            if isCtrl
                self.delete_spot(src);
            else
                self.move_spot(src);
            end
        end

        % Delete spot
        function delete_spot(self, src)
            for i = length(self.Spots):-1:1
                if self.Spots(i).patch == src || self.Spots(i).text == src
                    delete(self.Spots(i).patch);
                    delete(self.Spots(i).text);

                    self.AvailableIDs(end+1) = self.Spots(i).id;
                    self.AvailableIDs = sort(self.AvailableIDs);

                    self.Spots(i) = [];
                    break;
                end
            end
        end

        % Move spot
        function move_spot(self, src)
            pt = get(self.ax, 'CurrentPoint');
            x = pt(1, 1);
            y = pt(1, 2);

            for i = length(self.Spots):-1:1
                if self.Spots(i).patch == src || self.Spots(i).text == src
                    [xCircle, yCircle] = self.compute_circle_coordinates(x, y);
                    set(self.Spots(i).patch, 'XData', xCircle, 'YData', yCircle);

                    self.Spots(i).position = [x, y];
                    set(self.Spots(i).text, 'Position', [x, y]);
                    break;
                end
            end
        end
    end
end