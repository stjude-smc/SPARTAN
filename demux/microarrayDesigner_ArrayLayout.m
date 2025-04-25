classdef microarrayDesigner_ArrayLayout < handle
    properties (Access = public)
        traces_X = [];
        traces_Y = [];
        nX = 0;
        nY = 0;
        px_size;
        ax              matlab.ui.control.UIAxes;
        spot_size = 0;

        Spots = struct('patch', {}, 'position', {}, 'text', {}, 'id', {}); % Circle data
        AvailableIDs = []; % Pool of reusable IDs
        NextID = 1; % ID for the next circle
    end

    methods (Access = public)
        % constructor
        function self = microarrayDesigner_ArrayLayout(ax, spot_size)
            cam_px_size   = 6.5; % µm
            magnification = 60;
            cam_binning   = 2;
            self.px_size = cam_px_size*cam_binning/magnification;
            self.ax = ax;
            self.spot_size = spot_size;
        end

        function draw_traces(self)
        end

        function add_spot(self, x, y)
            spotFillColor   = [0.6 0.7 1.0];
            spotEdgeColor   = [0 0 0.8];
            spotAlpha       = 0.4;

            % Determine the ID for the new circle
            if ~isempty(self.AvailableIDs)
                id = self.AvailableIDs(1); % Reuse the first available ID
                self.AvailableIDs(1) = []; % Remove it from the pool
            else
                id = self.NextID; % Use the next sequential ID
                self.NextID = self.NextID + 1; % Increment for future use
            end
    
    
            % Create the circle patch
            theta = linspace(0, 2*pi, 100);
            xCircle = x + self.spot_size * cos(theta);
            yCircle = y + self.spot_size * sin(theta);
            patchHandle = patch('XData', xCircle, 'YData', yCircle, ...
                                'FaceColor', spotFillColor, 'EdgeColor', spotEdgeColor, ...
                                'FaceAlpha', spotAlpha, ...
                                'Parent', self.ax, ...
                                'ButtonDownFcn', @(src, event) self.handle_spot_click(src, event));
    
            % Add text label for the circle
            textHandle = text(self.ax, x, y, num2str(id), ...
                              'HorizontalAlignment', 'center', ...
                              'VerticalAlignment', 'middle', ...
                              'Color', 'k', 'FontSize', 25, ...
                              'FontWeight', 'bold', ...
                              'ButtonDownFcn', @(src, event) self.handle_spot_click(src, event));
    
            % Store circle data
            self.Spots(end+1).patch = patchHandle;
            self.Spots(end).position = [x, y];
            self.Spots(end).text = textHandle;
            self.Spots(end).id = id;
        end

    end

    methods (Access = private)

        function handle_spot_click(self, src, event)
            % Get the parent UIFigure of the axes
            parentFigure = ancestor(self.ax, 'figure'); % Retrieves the parent UIFigure
            modifiers = get(parentFigure, 'CurrentModifier'); % Query modifier keys
            isCtrl = ismember('control', modifiers);

            if isCtrl  % Delete the spot
                for i = length(self.Spots):-1:1
                    if self.Spots(i).patch == src || self.Spots(i).text == src
                        % Delete the patch and its associated data
                        delete(self.Spots(i).patch);
                        delete(self.Spots(i).text);

                         % Add the circle's ID to the pool of reusable IDs
                        self.AvailableIDs(end+1) = self.Spots(i).id;
                        self.AvailableIDs = sort(self.AvailableIDs);

                        self.Spots(i) = []; % Remove from list
                        break;
                    end
                end
            else  % Move the spot

                % Get mouse click position in axes coordinates
                pt = get(self.ax, 'CurrentPoint');
                x = pt(1, 1); % X-coordinate of the mouse click
                y = pt(1, 2); % Y-coordinate of the mouse click

                for i = length(self.Spots):-1:1
                    if self.Spots(i).patch == src || self.Spots(i).text == src
                        % Update the patch's position
                        theta = linspace(0, 2*pi, 100);
                        xCircle = x + self.spot_size * cos(theta);
                        yCircle = y + self.spot_size * sin(theta);
                        set(self.Spots(i).patch, 'XData', xCircle, 'YData', yCircle);

                        % Update the stored position
                        self.Spots(i).position = [x, y];

                        % Move the text label
                        set(self.Spots(i).text, 'Position', [x, y]);

                        break;
                    end
                end

            end
        end

    end
end