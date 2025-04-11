function interactiveCircles()
    spotSize = 50;
    % Create a figure and axes
    figure('WindowButtonDownFcn', @mouseClickCallback);
    ax = axes;
    imshow(rand(500, 500), 'Parent', ax); % Display a random image
    hold on;

    % Store circle data
    circles = struct('patch', {}, 'position', {}, 'text', {}, 'id', {});
    availableIDs = []; % Pool of reusable IDs
    nextID = 1; % ID for the next circle, if no reusable IDs exist
    
    % Callback for mouse clicks
    function mouseClickCallback(~, ~)
        % Get current point
        pt = get(ax, 'CurrentPoint');
        x = pt(1, 1);
        y = pt(1, 2);
        
        % Check for modifier keys
        modifiers = get(gcf, 'CurrentModifier');
        isShift = ismember('shift', modifiers);
        isCtrl = ismember('control', modifiers);
        
        if isShift
            % Add a new circle
            addCircle(x, y);
        elseif isCtrl
            % Delete circle if clicked near
            deleteCircle(x, y);
        else
            % Drag existing circle
            dragCircle(x, y);
        end
    end

    % Add a new circle
    function addCircle(x, y)
        % Determine the ID for the new circle
        if ~isempty(availableIDs)
            id = sort(availableIDs(1)); % Reuse the first available ID
            availableIDs(1) = []; % Remove it from the pool
            availableIDs = sort(availableIDs);
        else
            id = nextID; % Use the next sequential ID
            nextID = nextID + 1; % Increment for future use
        end
        
        % Create the circle patch
        theta = linspace(0, 2*pi, 100);
        xCircle = x + spotSize * cos(theta);
        yCircle = y + spotSize * sin(theta);
        patchHandle = patch('XData', xCircle, 'YData', yCircle, ...
                            'FaceColor', 'w', 'EdgeColor', 'none', ...
                            'Parent', ax, 'FaceAlpha', 0.6);
        
        % Add text label for the circle
        textHandle = text(x, y, num2str(id), ...
                          'HorizontalAlignment', 'center', ...
                          'VerticalAlignment', 'middle', ...
                          'Color', 'k', 'FontSize', 25, ...
                          'FontWeight', 'bold');
        
        % Store circle data
        circles(end+1).patch = patchHandle;
        circles(end).position = [x, y];
        circles(end).text = textHandle;
        circles(end).id = id;
    end

    % Delete a circle
    function deleteCircle(x, y)
        for i = length(circles):-1:1
            % Check distance from circle center
            pos = circles(i).position;
            dist = sqrt((x - pos(1))^2 + (y - pos(2))^2);
            if dist < spotSize % Within radius
                % Delete patch and text
                delete(circles(i).patch);
                delete(circles(i).text);
                
                % Add the circle's ID to the pool of reusable IDs
                availableIDs(end+1) = circles(i).id;
                availableIDs = sort(availableIDs);
                
                % Remove the circle from the list
                circles(i) = [];
                return;
            end
        end
    end

    % Drag a circle
    function dragCircle(x, y)
        for i = length(circles):-1:1
            % Check distance from circle center
            pos = circles(i).position;
            dist = sqrt((x - pos(1))^2 + (y - pos(2))^2);
            if dist < spotSize % Within radius
                % Update position
                circles(i).position = [x, y];
                
                % Update patch coordinates
                theta = linspace(0, 2*pi, 100);
                xCircle = x + spotSize * cos(theta);
                yCircle = y + spotSize * sin(theta);
                set(circles(i).patch, 'XData', xCircle, 'YData', yCircle);
                
                % Update text position
                set(circles(i).text, 'Position', [x, y]);
                return;
            end
        end
    end
end