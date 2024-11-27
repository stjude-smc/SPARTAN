function circles = fit_circles(traces, edge_coords, ax)
    % Find the circle centers from edge coordinates and optionally plot them.
    % Input:
    %   traces       - Struct containing image metadata (e.g., dimensions)
    %   edge_coords  - Nx2 array of edge pixel coordinates
    %   ax           - Optional axes object for plotting
    % Output:
    %   circle_centers - 4x2 array of [cx, cy] for the circles in each quadrant

    % Get image center
    nX = traces.fileMetadata.nX;
    nY = traces.fileMetadata.nY;
    cx_img = nX / 2;
    cy_img = nY / 2;

    % Initialize quadrants
    tl = []; % Top-left quadrant
    tr = []; % Top-right quadrant
    bl = []; % Bottom-left quadrant
    br = []; % Bottom-right quadrant

    % Split edge_coords into quadrants
    for i = 1:size(edge_coords, 1)
        x = edge_coords(i, 1);
        y = edge_coords(i, 2);
        if x <= cx_img && y <= cy_img
            tl = [tl; x, y]; % Top-left
        elseif x <= cx_img && y > cy_img
            tr = [tr; x, y]; % Top-right
        elseif x > cx_img && y <= cy_img
            bl = [bl; x, y]; % Bottom-left
        else
            br = [br; x, y]; % Bottom-right
        end
    end

    % Fit a circle for each quadrant
    circles = zeros(4, 3); % To store [cx, cy, r] for each circle
    quadrants = {tl, tr, bl, br};
    for q = 1:4
        edges = quadrants{q};
        
        % Check if edges are empty
        if isempty(edges)
            % Return a zero-sized circle
            circles(q, :) = [0, 0, 0];
            continue;
        end

        % Fit the circle
        [cx, cy, r] = fit_circle(edges);
        circles(q, :) = [cx, cy, r];


        % Optional plotting if ax is provided
        if nargin > 2 && ~isempty(ax) && isgraphics(ax, 'axes')
            % Calculate the radius based on the average distance from center to edge points
            x = edges(:, 1);
            y = edges(:, 2);
            %r = mean(sqrt((x - cx).^2 + (y - cy).^2));

            % Plot the circular patch
            rectangle(ax, 'Position', [cy - r, cx - r, 2*r, 2*r], ...
                      'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 0.5);

            % Plot the center
            plot(ax, cy, cx, 'g+', 'MarkerSize', 10, 'LineWidth', 2, ...
                 'DisplayName', ['Center ', num2str(q)]);
        end
    end

end

function [cx, cy, r] = fit_circle(edges)
% Fit a circle to a set of edge points using least-squares
% Input:
%   edges - Nx2 array of edge coordinates [x, y]
% Output:
%   [cx, cy, r] - Fitted circle parameters (center and radius)

    % Extra padding around the fitted circle
    padding = 10; % 10px ~ 2um

    % Fit a circle to a set of edge points using least-squares
    % Input:
    %   edges - Nx2 array of edge coordinates [x, y]
    % Output:
    %   [cx, cy] - Fitted circle center

    % Extract x and y coordinates
    x = edges(:, 1);
    y = edges(:, 2);

    % Set up the linear system for least-squares
    A = [2*x, 2*y, ones(size(x))];
    b = x.^2 + y.^2;

    % Solve the normal equations A * [cx; cy; c] = b
    params = A \ b;

    % Extract circle center
    cx = params(1);
    cy = params(2);
    r = sqrt(params(3) + cx^2 + cy^2) + padding;
end
