function layoutData = microarrayLayoutDesigner()
    % Persistent global variables
    spotIDs = [];
    traces_files = [];
    traces_xy = struct('x', [], 'y', [], 'nx', [], 'ny', []);
    %%

    % Spot colors
    inactiveFillColor = [0.9 0.9 0.9];
    inactiveEdgeColor = [0.6 0.6 0.6];
    activeFillColor   = [0.6 0.7 1.0];
    activeEdgeColor   = [0 0 0.8];
    AlphaActive       = 0.4;
    AlphaInactive     = 0.1;

    dlg = uifigure('Name', 'Microarray Layout Designer', 'Position', [100 100 800 500]);
 
    %% 
    % UI input fields
    wl = 70;  % width of label field
    wn = 60;   % width of numeric field
    g  = 10;   % gap
    c1 = 10 + g;
    c2 = c1 + wl + g;
    c3 = c2 + wn + g;
    c4 = c3 + wl + g;
    c5 = c4 + wn + g;

    uilabel(dlg, 'Position', [c1 420 wl 22], 'Text', 'Rows:');
    rowField = uispinner(dlg, 'Position', [c2 420 wn 22], 'Limits', [1 6], 'Value', 2);
    
    uilabel(dlg, 'Position', [c3 420 wl 22], 'Text', 'Columns:');
    colField = uispinner(dlg, 'Position', [c4 420 wn 22], 'Limits', [1 6], 'Value', 2);
    
    uilabel(dlg, 'Position', [c1 390 wl 22], 'Text', 'V pitch (µm):');
    vSpaceField = uispinner(dlg, 'Position', [c2 390 wn 22], 'Limits', [50 250], 'Step', 5, 'Value', 150);
    
    uilabel(dlg, 'Position', [c3 390 wl 22], 'Text', 'H pitch (µm):');
    hSpaceField = uispinner(dlg, 'Position', [c4 390 wn 22], 'Limits', [50 250], 'Step', 5, 'Value', 150);
    
    uilabel(dlg, 'Position', [c1 360 c4-c1 22], 'Text', 'Spot size (µm):');
    spotField = uispinner(dlg, 'Position', [c4 360 wn 22], 'Limits', [20 200], 'Step', 5, 'Value', 75);
    
    uilabel(dlg, 'Position', [c1 320 wl 22], 'Text', 'dX (µm):');
    shiftXField = uispinner(dlg, 'Position', [c2 320 wn 22], 'Limits', [-500 500], 'Step', 2, 'Value', 30);
    
    uilabel(dlg, 'Position', [c3 320 wl 22], 'Text', 'dY (µm):');
    shiftYField = uispinner(dlg, 'Position', [c4 320 wn 22], 'Limits', [-500 500], 'Step', 2, 'Value', 30);
    
    uilabel(dlg, 'Position', [c1 290 c4-c1 22], 'Text', 'Grid rotation (deg):');
    rotationField = uispinner(dlg, 'Position', [c4 290 wn 22], 'Limits', [-10 10], 'Step', 0.25, 'Value', 0);
    %%

    % Save/Load/Trace Buttons
    bw = c3-c1-g;  % button width
    uibutton(dlg, 'Text', 'Save layout', 'Position', [c1 120 bw 30], 'ButtonPushedFcn', @(btn, event) saveLayout());
    uibutton(dlg, 'Text', 'Load layout', 'Position', [c3 120 bw 30], 'ButtonPushedFcn', @(btn, event) loadLayout());
    uibutton(dlg, 'Text', 'Load traces', 'Position', [c1 80 bw 30], 'ButtonPushedFcn', @(btn, event) load_traces_callback());
    clear_btn = uibutton(dlg, 'Text', 'Clear traces', 'Position', [c3 80 bw 30], 'ButtonPushedFcn', @(btn, event) clear_traces_callback());

    process_btn = uibutton(dlg, 'Text', 'Demux traces', 'Position', [c3, 40, bw, 30], 'ButtonPushedFcn', @(btn, event) demux_traces());
    

    %%
    ax = uiaxes(dlg, 'Position', [c5 10 800-c5-2*g 470]);
    ax.Toolbar.Visible = 'off';

    fields = [rowField, colField, spotField, vSpaceField, hSpaceField, rotationField, shiftXField, shiftYField];
    for f = fields
        f.ValueChangedFcn = @(~, ~) updatePreview();
    end

    updatePreview();

    function updatePreview()
        rows = rowField.Value;
        cols = colField.Value;
        spot = spotField.Value;
        vSpacing = vSpaceField.Value;
        hSpacing = hSpaceField.Value;
        theta = deg2rad(rotationField.Value);
        dx = shiftXField.Value;
        dy = shiftYField.Value;

        if numel(traces_xy.x) == 0
            process_btn.Enable = 'off';
            clear_btn.Enable = 'off';
        else
            process_btn.Enable = 'on';
            clear_btn.Enable = 'on';
        end

        if isempty(spotIDs) || size(spotIDs, 1) ~= rows || size(spotIDs, 2) ~= cols
            spotIDs = zeros(rows, cols);
        end

        % Compute center of grid (unrotated)
        cx = (cols - 1) * hSpacing / 2;
        cy = (rows - 1) * vSpacing / 2;

        cla(ax);
        hold(ax, 'on');
        title(ax, 'Microarray layout');
        subtitle(ax, 'Assign spots with click or Ctrl+click');
        axis(ax, 'equal');
        box(ax, 'on');
        xlabel(ax, 'x, µm');
        ylabel(ax, 'y, µm');

        if numel(traces_xy.x) > 0
            set(ax, 'Color', 'k');
            rectangle('Position', [0, 0, traces_xy.nx, traces_xy.ny], 'FaceColor', 'w', 'EdgeColor', 'none', 'Parent', ax);
            scatter(ax, traces_xy.x, traces_xy.y, 3, 'k', 'filled');
        else
            set(ax, 'Color', 'w');
        end

        % Draw rotation center marker (crosshair)
        plot(ax, cx + dx, cy + dy, '+', 'Color', 'w', 'MarkerSize', 12, 'LineWidth', 2.5);
        plot(ax, cx + dx, cy + dy, '+', 'Color', [0.5 0 0.5], 'MarkerSize', 12, 'LineWidth', 0.5);

        for r = 1:rows
            for c = 1:cols
                x0 = (c - 1) * hSpacing - cx;
                y0 = (r - 1) * vSpacing - cy;
                x = cos(theta) * x0 - sin(theta) * y0 + cx + dx;
                y = sin(theta) * x0 + cos(theta) * y0 + cy + dy;
                theta_circ = linspace(0, 2 * pi, 50);
                xc = (spot / 2) * cos(theta_circ) + x;
                yc = (spot / 2) * sin(theta_circ) + y;

                if spotIDs(r, c) > 0
                    faceCol = activeFillColor;
                    edgeCol = activeEdgeColor;
                    alphaVal = AlphaActive;
                else
                    faceCol = inactiveFillColor;
                    edgeCol = inactiveEdgeColor;
                    alphaVal = AlphaInactive;
                end

                fill(ax, xc, yc, faceCol, 'EdgeColor', edgeCol, 'LineWidth', 1, 'FaceAlpha', alphaVal, 'EdgeAlpha', alphaVal + 0.2, 'ButtonDownFcn', @(src, event) onClick(r, c));

                if spotIDs(r, c) > 0
                    text(ax, x, y, num2str(spotIDs(r, c)), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 28, ...
                        'Color', 'k', ...
                        'BackgroundColor', 'w', ...
                        'Margin', 2, ...
                        'ButtonDownFcn', @(src, event) onClick(r, c));

                end
            end
        end
        
        dx_rot = abs(sin(theta) * (spot + hSpacing * (cols - 1)/2));
        dy_rot = abs(sin(theta) * (spot + vSpacing * (rows - 1)/2));
        xmin = min(-5, -2*spot/3 + dx - dx_rot);

        nx = traces_xy.nx;
        if isempty(nx), nx = 0; end
        ny = traces_xy.ny;
        if isempty(ny), ny = 0; end

        xmax = max(5 + nx, hSpacing * (cols - 1) + 2*spot/3 + dx + dx_rot);
        ymin = min(-5, -2*spot/3 + dy - dy_rot);
        ymax = max(5 + ny, vSpacing * (rows - 1) + 2*spot/3 + dy + dy_rot);
        ax.XLim = [xmin, xmax];
        ax.YLim = [ymin, ymax];
        hold(ax, 'off');
    end

    function onClick(r, c, ~)
        clickType = dlg.SelectionType;
        modifiers = dlg.CurrentModifier;
        if strcmp(clickType, 'normal') || strcmp(clickType, 'open')
            if ismember('control', modifiers) && spotIDs(r, c) > 0
                spotIDs(r, c) = spotIDs(r, c) - 1;
            else
                spotIDs(r, c) = spotIDs(r, c) + 1;
            end
        end
        updatePreview();
    end

    function saveLayout()
        layoutData.spotIDs = spotIDs;
        layoutData.rows = rowField.Value;
        layoutData.cols = colField.Value;
        layoutData.spotSize = spotField.Value;
        layoutData.verticalSpacing = vSpaceField.Value;
        layoutData.horizontalSpacing = hSpaceField.Value;
        layoutData.rotation = rotationField.Value;
        layoutData.shiftX = shiftXField.Value;
        layoutData.shiftY = shiftYField.Value;

        [file, path] = uiputfile('microarray_layout.json', 'Save microarray layout');
        if ischar(file)
            jsonStr = jsonencode(layoutData);
            fid = fopen(fullfile(path, file), 'w');
            fwrite(fid, jsonStr, 'char');
            fclose(fid);
        end
    end

    function loadLayout()
        [file, path] = uigetfile('*.json', 'Load microarray layout');
        if ischar(file)
            jsonStr = fileread(fullfile(path, file));
            layoutData = jsondecode(jsonStr);

            rowField.Value = layoutData.rows;
            colField.Value = layoutData.cols;
            spotField.Value = layoutData.spotSize;
            vSpaceField.Value = layoutData.verticalSpacing;
            hSpaceField.Value = layoutData.horizontalSpacing;
            rotationField.Value = layoutData.rotation;
            shiftXField.Value = layoutData.shiftX;
            shiftYField.Value = layoutData.shiftY;
            spotIDs = layoutData.spotIDs;
            updatePreview();
        end
    end

    function load_traces_callback()
        filter = {'*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
          '*.traces','Binary Traces Files (*.traces)';};
        traces_files = getFiles(filter);
        traces_xy = loadTracesXY(traces_files);
        updatePreview();
    end

    function clear_traces_callback()
        traces_xy = struct('x', [], 'y', [], 'nx', [], 'ny', []);
        traces_files = [];
        updatePreview();
    end

    function demux_traces()
        % This should call demux_spots with specific parameters
    end
end
