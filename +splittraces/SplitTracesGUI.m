classdef SplitTracesGUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        % GUI layout elements (public)
        axCannyIn                    matlab.ui.control.UIAxes
        axCannyOut                   matlab.ui.control.UIAxes
        axMicroarray                 matlab.ui.control.UIAxes
        btnSplitTraces               matlab.ui.control.Button
        btnLoadLayout                matlab.ui.control.Button
        btnLoadSettings              matlab.ui.control.Button
        btnLoadTraces                matlab.ui.control.Button
        btnSaveLayout                matlab.ui.control.Button
        btnSaveSettings              matlab.ui.control.Button
        chkAutoUpdate                matlab.ui.control.CheckBox
        chkCreateSubdir              matlab.ui.control.CheckBox
        chkPlotResult                matlab.ui.control.CheckBox
        chkSavePlot                  matlab.ui.control.CheckBox
        chkCombineTraces             matlab.ui.control.CheckBox
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
    end

    properties (Access = private)
        % Helper classes
        spot_layout splittraces.SpotLayout;
        edge_mapper splittraces.SpotEdgeMapper;

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
            app.edge_mapper.load_traces_XY(getFiles(filter));
            app.spot_layout.draw_traces(app.edge_mapper.traces_xy);
        end

        function load_layout_pushed(app)
            ss = app.spot_layout.load_layout();
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
                app.spot_layout.add_spot(x, y);
            end
        end

        function autoUpdateChanged(app, value)
            app.btnSplitTraces.Enable = matlab.lang.OnOffSwitchState(value);
            app.edge_mapper.enable_auto_update(value);
        end


    %%
    % Component initialization
    %

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'on';
            app.UIFigure.Position = [100 100 1000 600];
            app.UIFigure.Name = 'Split Microarray Traces';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'2x', '3x'};
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
            app.sSpotSize.ValueChangedFcn = @(src, event) app.spot_layout.update_spot_size(src.Value);

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
            app.btnSaveLayout.ButtonPushedFcn = @(src, event) app.spot_layout.save_layout()

            % Create lblAxMicroarrayLayout
            app.lblAxMicroarrayLayout = uilabel(app.GridLayout2);
            app.lblAxMicroarrayLayout.HorizontalAlignment = 'center';
            app.lblAxMicroarrayLayout.Layout.Row = 1;
            app.lblAxMicroarrayLayout.Layout.Column = 1;
            app.lblAxMicroarrayLayout.Text = {'Microarray layout'; 'Add spots with Shift+click, click to move, delete with Ctrl+click'};

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.RightPanel);
            app.GridLayout4.ColumnWidth = {'1x'};
            app.GridLayout4.RowHeight = {'1x', 60};

            % Create EdgedetectionPanel
            app.EdgedetectionPanel = uipanel(app.GridLayout4);
            app.EdgedetectionPanel.BorderType = 'none';
            app.EdgedetectionPanel.Title = 'Edge detection';
            app.EdgedetectionPanel.Layout.Row = 1;
            app.EdgedetectionPanel.Layout.Column = 1;

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.EdgedetectionPanel);
            app.GridLayout5.ColumnWidth = {'1x', '1x'};
            %app.GridLayout5.RowHeight = {'14x', '5x', '2x'};
            app.GridLayout5.RowHeight = {'1x', 120, 50};

            % Create axCannyIn
            app.axCannyIn = uiaxes(app.GridLayout5);
            title(app.axCannyIn, 'Canny edges input');
            axis(app.axCannyIn, 'equal');
            axis(app.axCannyIn, 'off');
            app.axCannyIn.Layout.Row = 1;
            app.axCannyIn.Layout.Column = 1;

            % Create axCannyOut
            app.axCannyOut = uiaxes(app.GridLayout5);
            title(app.axCannyOut, 'Detected edges');
            axis(app.axCannyOut, 'equal');
            app.axCannyOut.Layout.Row = 1;
            app.axCannyOut.Layout.Column = 2;

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.GridLayout5);
            app.GridLayout6.ColumnWidth = {'2x', '1x'};
            app.GridLayout6.RowHeight = {'1x', '1x', '1x'};
            app.GridLayout6.Layout.Row = 2;
            app.GridLayout6.Layout.Column = 1;

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
            app.ddDownscale1.ValueChangedFcn = @(src, event) app.edge_mapper.set('Downscale1', str2double(event.Value(1)));

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
            app.sDilation.Value = 6;
            app.sDilation.ValueChangedFcn = @(src, event) app.edge_mapper.set('Dilation', event.Value);

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
            app.ddDownscale2.Value = '2x';
            app.ddDownscale2.ValueChangedFcn = @(src, event) app.edge_mapper.set('Downscale2', str2double(event.Value(1)));

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
            app.sLowThreshold.ValueChangedFcn = @(src, event) app.edge_mapper.set('LowThreshold', event.Value);

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
            app.sHighThreshold.Value = 0.6;
            app.sHighThreshold.ValueChangedFcn = @(src, event) app.edge_mapper.set('HighThreshold', event.Value);

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
            app.sSigma.Value = 7;
            app.sSigma.ValueChangedFcn = @(src, event) app.edge_mapper.set('Sigma', event.Value);

            % Create GridLayout7
            app.chkAutoUpdate = uicheckbox(app.GridLayout5);
            app.chkAutoUpdate.Text = 'Auto update edges';
            app.chkAutoUpdate.Value = true;
            app.chkAutoUpdate.Layout.Row = 3;
            app.chkAutoUpdate.Layout.Column = 1;
            app.chkAutoUpdate.ValueChangedFcn = @(src, event) app.autoUpdateChanged(event.Value)

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
            app.btnSaveSettings.ButtonPushedFcn = @(src, event) app.edge_mapper.save_settings();

            % Create btnLoadSettings
            app.btnLoadSettings = uibutton(app.GridLayout7, 'push');
            app.btnLoadSettings.Layout.Row = 1;
            app.btnLoadSettings.Layout.Column = 2;
            app.btnLoadSettings.Text = 'Load settings';
            app.btnLoadSettings.ButtonPushedFcn = @(src, event) app.edge_mapper.load_settings(app);

            % Create OutputPanel
            app.OutputPanel = uipanel(app.GridLayout4);
            app.OutputPanel.BorderType = 'none';
            app.OutputPanel.Title = 'Output';
            app.OutputPanel.Layout.Row = 2;
            app.OutputPanel.Layout.Column = 1;

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.OutputPanel);
            app.GridLayout8.ColumnWidth = {'13x', '8x', '8x', '11x', '8x'};
            app.GridLayout8.RowHeight = {'1x'};

            % Create chkCreateSubdir
            app.chkCreateSubdir = uicheckbox(app.GridLayout8);
            app.chkCreateSubdir.Text = 'Create subdirectory';
            app.chkCreateSubdir.Layout.Row = 1;
            app.chkCreateSubdir.Layout.Column = 1;
            app.chkCreateSubdir.Value = true;
            app.chkCreateSubdir.ValueChangedFcn = @(src, event) app.edge_mapper.set('CreateSubdir', event.Value);

            % Create chkPlotResult
            app.chkPlotResult = uicheckbox(app.GridLayout8);
            app.chkPlotResult.Text = 'Show plots';
            app.chkPlotResult.Layout.Row = 1;
            app.chkPlotResult.Layout.Column = 2;
            app.chkPlotResult.Value = true;
            app.chkPlotResult.ValueChangedFcn = @(src, event) app.edge_mapper.set('PlotResult', event.Value);

            % Create chkSavePlot
            app.chkSavePlot = uicheckbox(app.GridLayout8);
            app.chkSavePlot.Text = 'Save plots';
            app.chkSavePlot.Layout.Row = 1;
            app.chkSavePlot.Layout.Column = 3;
            app.chkSavePlot.Value = true;
            app.chkSavePlot.ValueChangedFcn = @(src, event) app.edge_mapper.set('SavePlot', event.Value);

            % Create chkCombineTraces
            app.chkCombineTraces = uicheckbox(app.GridLayout8);
            app.chkCombineTraces.Text = 'Combine datasets';
            app.chkCombineTraces.Layout.Row = 1;
            app.chkCombineTraces.Layout.Column = 4;
            app.chkCombineTraces.Value = true;
            app.chkCombineTraces.ValueChangedFcn = @(src, event) app.edge_mapper.set('CombineDatasets', event.Value);

            % Create btnSplitTraces
            app.btnSplitTraces = uibutton(app.GridLayout8, 'push');
            app.btnSplitTraces.Layout.Row = 1;
            app.btnSplitTraces.Layout.Column = 5;
            app.btnSplitTraces.Text = 'Split traces';
            app.btnSplitTraces.ButtonPushedFcn = @(src, event) app.edge_mapper.split_traces();

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)
        % Called when microarray layout has been changed, can trigger additional functions
        function spots_changed(app)
            app.edge_mapper.compute_edges();
        end

        % Construct app
        function app = SplitTracesGUI()
            app.px_size = compute_pixel_size(6.5, 60, 2);

            % Create UIFigure and components
            createComponents(app);

            % Register the app with App Designer
            registerApp(app, app.UIFigure);

            % Create microarrya layout handler
            app.spot_layout = splittraces.SpotLayout(app.axMicroarray, app, app.sSpotSize.Value, app.px_size);

            % Create edge detection handler and propagate settings from GUI
            app.edge_mapper = splittraces.SpotEdgeMapper(app.axCannyIn, app.axCannyOut, app.px_size, app.spot_layout);
            
            % Propagate values in GUI to the backend
            app.autoUpdateChanged(app.chkAutoUpdate.Value);  % "Auto update" checkbox
            app.edge_mapper.set('Downscale1', str2double(app.ddDownscale1.Value(1)));
            app.edge_mapper.set('Dilation', app.sDilation.Value);
            app.edge_mapper.set('Downscale2', str2double(app.ddDownscale2.Value(1)));
            app.edge_mapper.set('LowThreshold', app.sLowThreshold.Value);
            app.edge_mapper.set('HighThreshold', app.sHighThreshold.Value);
            app.edge_mapper.set('Sigma', app.sSigma.Value);
            app.edge_mapper.set('CreateSubdir', app.chkCreateSubdir.Value);
            app.edge_mapper.set('PlotResult', app.chkPlotResult.Value);
            app.edge_mapper.set('SavePlot', app.chkSavePlot.Value);
            app.edge_mapper.set('CombineDatasets', app.chkCombineTraces.Value);

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

