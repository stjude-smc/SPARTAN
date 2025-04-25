classdef microarrayDesigner_ArrayLayout
    properties (Access = public)
        traces_X = [];
        traces_Y = [];
        nX = 0;
        nY = 0;
        px_size;
        ax;
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

        function draw_traces()
        end

        function add_spot()
        end

        function remove_spot()
        end

        function drag_spot()
        end
    end
end