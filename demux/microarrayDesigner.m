classdef microarrayDesigner < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        % GUI layout elements (public)
        axCannyIn                    matlab.ui.control.UIAxes
        axCannyOut                   matlab.ui.control.UIAxes
        axMicroarray                 matlab.ui.control.UIAxes
        btnDemuxTraces               matlab.ui.control.Button
        btnLoadLayout                matlab.ui.control.Button
        btnLoadSettings              matlab.ui.control.Button
        btnLoadTraces                matlab.ui.control.Button
        btnSaveLayout                matlab.ui.control.Button
        btnSaveSettings              matlab.ui.control.Button
        chkAutoUpdate                matlab.ui.control.CheckBox
        chkCreateSubdir              matlab.ui.control.CheckBox
        ddDownscale1                 matlab.ui.control.DropDown
        ddDownscale2                 matlab.ui.control.DropDown
        sDilation                    matlab.ui.control.Spinner
        sHighThreshold               matlab.ui.control.Spinner
        sLowThreshold                matlab.ui.control.Spinner
        sSigma                       matlab.ui.control.Spinner
        sSpotSize                    matlab.ui.control.Spinner

        % Public data
        px_size;
        SpotSize = 50; % Size of the circles

        % Coordinates of loaded traces
        traces_xy = struct('X', [], 'Y', [], 'nX', [], 'nY', []);
        traces_files = [];  % Open traces files


        % Microarray layout manager - consolidates microarray-related methods
        array_layout                 microarrayDesigner_ArrayLayout
        edge_detect                  microarray_EdgeDetectionHandler
    end

    properties (Access = private)
        % Auto-reflow width
        onePanelWidth = 700;

        % GUI layout elements
        UIFigure                     matlab.ui.Figure
        EdgedetectionPanel           matlab.ui.container.Panel
        GridLayout                   matlab.ui.container.GridLayout
        GridLayout2                  matlab.ui.container.GridLayout
        GridLayout3                  matlab.ui.container.GridLayout
        GridLayout4                  matlab.ui.container.GridLayout
        GridLayout5                  matlab.ui.container.GridLayout
        GridLayout6                  matlab.ui.container.GridLayout
        GridLayout6_2                matlab.ui.container.GridLayout
        GridLayout7                  matlab.ui.container.GridLayout
        GridLayout8                  matlab.ui.container.GridLayout
        lblAxMicroarrayLayout        matlab.ui.control.Label
        lblDilationpxSpinner         matlab.ui.control.Label
        lblDownscalestep1DropDown    matlab.ui.control.Label
        lblDownscalestep2DropDown    matlab.ui.control.Label
        lblHighthreshold             matlab.ui.control.Label
        lblLowthresholdSpinner       matlab.ui.control.Label
        lblSigmaSpinner              matlab.ui.control.Label
        lblSpotsizemSpinner          matlab.ui.control.Label
        LeftPanel                    matlab.ui.container.Panel
        OutputPanel                  matlab.ui.container.Panel
        RightPanel                   matlab.ui.container.Panel

        % Internal data
        % Camera pixel size
        cam_px_size   = 6.5; % µm
        magnification = 60;
        cam_binning   = 2;
    end
    

    % Callbacks that handle component events
    methods (Access = private)

    %%
    % Callbacks
    %

        % Button pushed function: btnLoadTraces
        function btnLoadTracesButtonPushed(app)
            filter = {'*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.traces','Binary Traces Files (*.traces)';};
            app.traces_files = getFiles(filter);
            app.traces_xy = loadTracesXY(app.traces_files)
            app.array_layout.draw_traces(app.traces_xy);
        end

        function load_layout_pushed(app)
            ss = app.array_layout.load_layout();
            if isnumeric(ss)
                app.sSpotSize.Value = ss;
            end
        end

        % Button down function: axMicroarray
        function axMicroarrayButtonDown(app)
            % Get current point
            pt = get(app.axMicroarray, 'CurrentPoint');
            x = pt(1, 1);
            y = pt(1, 2);
    
            % Check for modifier keys
            modifiers = get(app.UIFigure, 'CurrentModifier');    
            if ismember('shift', modifiers)
                app.array_layout.add_spot(x, y);
            end
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {515, 515};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {422, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

    %%
    % Component initialization
    %

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 909 515];
            app.UIFigure.Name = 'Microarray Designer';
            app.UIFigure.SizeChangedFcn = @(src, event) app.updateAppLayout();

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {422, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.LeftPanel);
            app.GridLayout2.ColumnWidth = {'1x'};
            app.GridLayout2.RowHeight = {'1x', '10x', '2x'};

            % Create axMicroarray
            app.axMicroarray = uiaxes(app.GridLayout2);
            app.axMicroarray.Layout.Row = 2;
            app.axMicroarray.Layout.Column = 1;
            app.axMicroarray.ButtonDownFcn = @(src, event) app.axMicroarrayButtonDown();

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GridLayout2);
            app.GridLayout3.ColumnWidth = {'1x', '1x', '1x'};
            app.GridLayout3.Layout.Row = 3;
            app.GridLayout3.Layout.Column = 1;

            % Create btnLoadTraces
            app.btnLoadTraces = uibutton(app.GridLayout3, 'push');
            app.btnLoadTraces.ButtonPushedFcn = @(src, event) app.btnLoadTracesButtonPushed();
            app.btnLoadTraces.Layout.Row = 1;
            app.btnLoadTraces.Layout.Column = 1;
            app.btnLoadTraces.Text = 'Load traces';

            % Create SpotsizemSpinnerLabel
            app.lblSpotsizemSpinner = uilabel(app.GridLayout3);
            app.lblSpotsizemSpinner.HorizontalAlignment = 'right';
            app.lblSpotsizemSpinner.Layout.Row = 1;
            app.lblSpotsizemSpinner.Layout.Column = 2;
            app.lblSpotsizemSpinner.Text = 'Spot size, µm';

            % Create sSpotSize
            app.sSpotSize = uispinner(app.GridLayout3);
            app.sSpotSize.Step = 5;
            app.sSpotSize.Limits = [10 250];
            app.sSpotSize.RoundFractionalValues = 'on';
            app.sSpotSize.HorizontalAlignment = 'left';
            app.sSpotSize.Layout.Row = 1;
            app.sSpotSize.Layout.Column = 3;
            app.sSpotSize.Value = 100;
            app.sSpotSize.ValueChangedFcn = @(src, event) app.array_layout.update_spot_size(src.Value);

            % Create bLoadLayout
            app.btnLoadLayout = uibutton(app.GridLayout3, 'push');
            app.btnLoadLayout.Layout.Row = 2;
            app.btnLoadLayout.Layout.Column = 2;
            app.btnLoadLayout.Text = 'Load layout';
            app.btnLoadLayout.ButtonPushedFcn = @(src, event) app.load_layout_pushed()

            % Create btnSaveLayout
            app.btnSaveLayout = uibutton(app.GridLayout3, 'push');
            app.btnSaveLayout.Layout.Row = 2;
            app.btnSaveLayout.Layout.Column = 3;
            app.btnSaveLayout.Text = 'Save layout';
            app.btnSaveLayout.ButtonPushedFcn = @(src, event) app.array_layout.save_layout()

            % Create lblAxMicroarrayLayout
            app.lblAxMicroarrayLayout = uilabel(app.GridLayout2);
            app.lblAxMicroarrayLayout.HorizontalAlignment = 'center';
            app.lblAxMicroarrayLayout.Layout.Row = 1;
            app.lblAxMicroarrayLayout.Layout.Column = 1;
            app.lblAxMicroarrayLayout.Text = {'Microarray layout'; 'Add spots with Shift+click, remove with Ctrl+click'};

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.RightPanel);
            app.GridLayout4.ColumnWidth = {'1x'};
            app.GridLayout4.RowHeight = {'6x', '1x'};

            % Create EdgedetectionPanel
            app.EdgedetectionPanel = uipanel(app.GridLayout4);
            app.EdgedetectionPanel.BorderType = 'none';
            app.EdgedetectionPanel.Title = 'Edge detection';
            app.EdgedetectionPanel.Layout.Row = 1;
            app.EdgedetectionPanel.Layout.Column = 1;

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.EdgedetectionPanel);
            app.GridLayout5.ColumnWidth = {'2x', '3x'};
            app.GridLayout5.RowHeight = {'7x', '7x', '2x'};

            % Create axCannyIn
            app.axCannyIn = uiaxes(app.GridLayout5);
            title(app.axCannyIn, 'Canny edges input')
            xlabel(app.axCannyIn, 'X')
            ylabel(app.axCannyIn, 'Y')
            zlabel(app.axCannyIn, 'Z')
            app.axCannyIn.Layout.Row = 1;
            app.axCannyIn.Layout.Column = 1;

            % Create axCannyOut
            app.axCannyOut = uiaxes(app.GridLayout5);
            title(app.axCannyOut, 'Canny edges output')
            xlabel(app.axCannyOut, 'X')
            ylabel(app.axCannyOut, 'Y')
            zlabel(app.axCannyOut, 'Z')
            app.axCannyOut.Layout.Row = 2;
            app.axCannyOut.Layout.Column = 1;

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.GridLayout5);
            app.GridLayout6.ColumnWidth = {'2x', '1x'};
            app.GridLayout6.RowHeight = {'1x', '1x', '1x'};
            app.GridLayout6.Layout.Row = 1;
            app.GridLayout6.Layout.Column = 2;

            % Create lblDownscalestep1DropDown
            app.lblDownscalestep1DropDown = uilabel(app.GridLayout6);
            app.lblDownscalestep1DropDown.HorizontalAlignment = 'right';
            app.lblDownscalestep1DropDown.Layout.Row = 1;
            app.lblDownscalestep1DropDown.Layout.Column = 1;
            app.lblDownscalestep1DropDown.Text = 'Downscale step 1';

            % Create ddDownscale1
            app.ddDownscale1 = uidropdown(app.GridLayout6);
            app.ddDownscale1.Items = {'1x', '2x', '4x', '8x'};
            app.ddDownscale1.Layout.Row = 1;
            app.ddDownscale1.Layout.Column = 2;
            app.ddDownscale1.Value = '4x';
            app.ddDownscale1.ValueChangedFcn = @(src, event) app.edge_detect.set('Downscale1', str2double(event.Value(1)));

            % Create lblDilationpxSpinner
            app.lblDilationpxSpinner = uilabel(app.GridLayout6);
            app.lblDilationpxSpinner.HorizontalAlignment = 'right';
            app.lblDilationpxSpinner.Layout.Row = 2;
            app.lblDilationpxSpinner.Layout.Column = 1;
            app.lblDilationpxSpinner.Text = 'Dilation, px';

            % Create sDilation
            app.sDilation = uispinner(app.GridLayout6);
            app.sDilation.Limits = [1 16];
            app.sDilation.HorizontalAlignment = 'left';
            app.sDilation.Layout.Row = 2;
            app.sDilation.Layout.Column = 2;
            app.sDilation.Value = 2;
            app.sDilation.ValueChangedFcn = @(src, event) app.edge_detect.set('Dilation', event.Value);

            % Create lblDownscalestep2DropDown
            app.lblDownscalestep2DropDown = uilabel(app.GridLayout6);
            app.lblDownscalestep2DropDown.HorizontalAlignment = 'right';
            app.lblDownscalestep2DropDown.Layout.Row = 3;
            app.lblDownscalestep2DropDown.Layout.Column = 1;
            app.lblDownscalestep2DropDown.Text = 'Downscale step 2';

            % Create ddDownscale2
            app.ddDownscale2 = uidropdown(app.GridLayout6);
            app.ddDownscale2.Items = {'1x', '2x', '4x', '8x'};
            app.ddDownscale2.Layout.Row = 3;
            app.ddDownscale2.Layout.Column = 2;
            app.ddDownscale2.Value = '4x';
            app.ddDownscale2.ValueChangedFcn = @(src, event) app.edge_detect.set('Downscale2', str2double(event.Value(1)));

            % Create GridLayout6_2
            app.GridLayout6_2 = uigridlayout(app.GridLayout5);
            app.GridLayout6_2.ColumnWidth = {'2x', '1x'};
            app.GridLayout6_2.RowHeight = {'1x', '1x', '1x'};
            app.GridLayout6_2.Layout.Row = 2;
            app.GridLayout6_2.Layout.Column = 2;

            % Create lblLowthresholdSpinner
            app.lblLowthresholdSpinner = uilabel(app.GridLayout6_2);
            app.lblLowthresholdSpinner.HorizontalAlignment = 'right';
            app.lblLowthresholdSpinner.Layout.Row = 1;
            app.lblLowthresholdSpinner.Layout.Column = 1;
            app.lblLowthresholdSpinner.Text = 'Low threshold';

            % Create sLowThreshold
            app.sLowThreshold = uispinner(app.GridLayout6_2);
            app.sLowThreshold.Step = 0.05;
            app.sLowThreshold.Limits = [0 95];
            app.sLowThreshold.HorizontalAlignment = 'left';
            app.sLowThreshold.Layout.Row = 1;
            app.sLowThreshold.Layout.Column = 2;
            app.sLowThreshold.Value = 0.05;
            app.sLowThreshold.ValueChangedFcn = @(src, event) app.edge_detect.set('LowThreshold', event.Value);

            % Create lblHighthreshold
            app.lblHighthreshold = uilabel(app.GridLayout6_2);
            app.lblHighthreshold.HorizontalAlignment = 'right';
            app.lblHighthreshold.Layout.Row = 2;
            app.lblHighthreshold.Layout.Column = 1;
            app.lblHighthreshold.Text = 'High threshold';

            % Create sHighThreshold
            app.sHighThreshold = uispinner(app.GridLayout6_2);
            app.sHighThreshold.Step = 0.05;
            app.sHighThreshold.Limits = [0 95];
            app.sHighThreshold.HorizontalAlignment = 'left';
            app.sHighThreshold.Layout.Row = 2;
            app.sHighThreshold.Layout.Column = 2;
            app.sHighThreshold.Value = 0.3;
            app.sHighThreshold.ValueChangedFcn = @(src, event) app.edge_detect.set('HighThreshold', event.Value);

            % Create lblSigmaSpinner
            app.lblSigmaSpinner = uilabel(app.GridLayout6_2);
            app.lblSigmaSpinner.HorizontalAlignment = 'right';
            app.lblSigmaSpinner.Layout.Row = 3;
            app.lblSigmaSpinner.Layout.Column = 1;
            app.lblSigmaSpinner.Text = 'Sigma';

            % Create sSigma
            app.sSigma = uispinner(app.GridLayout6_2);
            app.sSigma.Step = 0.2;
            app.sSigma.Limits = [1 10];
            app.sSigma.HorizontalAlignment = 'left';
            app.sSigma.Layout.Row = 3;
            app.sSigma.Layout.Column = 2;
            app.sSigma.Value = 5;
            app.sSigma.ValueChangedFcn = @(src, event) app.edge_detect.set('Sigma', event.Value);

            % Create GridLayout7
            app.chkAutoUpdate = uicheckbox(app.GridLayout5);
            app.chkAutoUpdate.Text = 'Auto update edges';
            app.chkAutoUpdate.Value = false;
            app.chkAutoUpdate.Layout.Row = 3;
            app.chkAutoUpdate.Layout.Column = 1;
            app.chkAutoUpdate.ValueChangedFcn = @(src, event) app.edge_detect.enable_auto_update(event.Value);

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.GridLayout5);
            app.GridLayout7.RowHeight = {'1x'};
            app.GridLayout7.Layout.Row = 3;
            app.GridLayout7.Layout.Column = 2;

            % Create btnSaveSettings
            app.btnSaveSettings = uibutton(app.GridLayout7, 'push');
            app.btnSaveSettings.Layout.Row = 1;
            app.btnSaveSettings.Layout.Column = 1;
            app.btnSaveSettings.Text = 'Save settings';
            app.btnSaveSettings.ButtonPushedFcn = @(src, event) app.edge_detect.save_settings();

            % Create btnLoadSettings
            app.btnLoadSettings = uibutton(app.GridLayout7, 'push');
            app.btnLoadSettings.Layout.Row = 1;
            app.btnLoadSettings.Layout.Column = 2;
            app.btnLoadSettings.Text = 'Load settings';
            app.btnLoadSettings.ButtonPushedFcn = @(src, event) app.edge_detect.load_settings(app);

            % Create OutputPanel
            app.OutputPanel = uipanel(app.GridLayout4);
            app.OutputPanel.BorderType = 'none';
            app.OutputPanel.Title = 'Output';
            app.OutputPanel.Layout.Row = 2;
            app.OutputPanel.Layout.Column = 1;

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.OutputPanel);
            app.GridLayout8.ColumnWidth = {'2x', '1x'};
            app.GridLayout8.RowHeight = {'1x'};

            % Create chkCreateSubdir
            app.chkCreateSubdir = uicheckbox(app.GridLayout8);
            app.chkCreateSubdir.Text = 'Create a subdirectory';
            app.chkCreateSubdir.Layout.Row = 1;
            app.chkCreateSubdir.Layout.Column = 1;
            app.chkCreateSubdir.Value = true;

            % Create btnDemuxTraces
            app.btnDemuxTraces = uibutton(app.GridLayout8, 'push');
            app.btnDemuxTraces.Layout.Row = 1;
            app.btnDemuxTraces.Layout.Column = 2;
            app.btnDemuxTraces.Text = 'Demux traces';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = microarrayDesigner()
            app.px_size = compute_pixel_size(6.5, 60, 2);

            % Create UIFigure and components
            createComponents(app);

            % Register the app with App Designer
            registerApp(app, app.UIFigure);

            app.array_layout = microarrayDesigner_ArrayLayout(app.axMicroarray, app.sSpotSize.Value, app.px_size);
            app.edge_detect = microarray_EdgeDetectionHandler();

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end


% Compute pixel size
function px_size = compute_pixel_size(cam_px_size, magnification, cam_binning)
    px_size = cam_px_size * cam_binning / magnification;
end
