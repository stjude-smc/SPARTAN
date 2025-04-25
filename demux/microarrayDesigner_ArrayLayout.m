classdef microarrayDesigner_ArrayLayout < handle
    properties (Access = public)
        traces_X = [];
        traces_Y = [];
        nX = 0;
        nY = 0;
        px_size;
        ax;

        Spots = struct('patch', {}, 'position', {}, 'text', {}, 'id', {}); % Circle data
        AvailableIDs = []; % Pool of reusable IDs
        NextID = 1; % ID for the next circle
    end

    methods (Access = public)
        % constructor
        function obj = microarrayDesigner_ArrayLayout(axes)
            cam_px_size   = 6.5; % µm
            magnification = 60;
            cam_binning   = 2;
            obj.px_size = cam_px_size*cam_binning/magnification;
            obj.ax = axes;
        end

        function draw_traces(self)
        end

        function add_spot(self, app, x, y)
            'Adding a circle '
            spotFillColor   = [0.6 0.7 1.0];
            spotEdgeColor   = [0 0 0.8];
            spotAlpha       = 0.4;

            ss = app.sSpotSize.Value;
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
            xCircle = x + ss * cos(theta);
            yCircle = y + ss * sin(theta);
            patchHandle = patch('XData', xCircle, 'YData', yCircle, ...
                                'FaceColor', spotFillColor, 'EdgeColor', spotEdgeColor, ...
                                'FaceAlpha', spotAlpha, ...
                                'Parent', app.axMicroarray);
    
            % Add text label for the circle
            textHandle = text(app.axMicroarray, x, y, num2str(id), ...
                              'HorizontalAlignment', 'center', ...
                              'VerticalAlignment', 'middle', ...
                              'Color', 'k', 'FontSize', 25, ...
                              'FontWeight', 'bold');
    
            % Store circle data
            self.Spots(end+1).patch = patchHandle;
            self.Spots(end).position = [x, y];
            self.Spots(end).text = textHandle;
            self.Spots(end).id = id;
        end

        function remove_spot(self)
        end

        function drag_spot(self)
        end
    end
end