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
    cx_img = floor(nX / 2);
    cy_img = floor(nY / 2);

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
            bl = [bl; x, y]; % Bottom-left
        elseif x > cx_img && y <= cy_img
            tr = [tr; x, y]; % Top-right
        else % x > cx_img && y > cy_img
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
% Fit a circle to a set of edge points using RANSAC for robustness.
% Input:
%   edges - Nx2 array of edge coordinates [x, y]
% Output:
%   [cx, cy, r] - Fitted circle parameters (center and radius)

    % Restrictions on circle radius
    min_r = 205; % px ~90um
    max_r = 280; % px ~120um

    % Parameters for RANSAC
    max_iterations = 1000; % Maximum number of RANSAC iterations
    inlier_threshold = 12;  % Distance threshold to count a point as an inlier (in pixels)
    min_inliers = 10;      % Minimum number of inliers for a valid model

    % Initialize best circle parameters and inlier count
    best_cx = 0;
    best_cy = 0;
    best_r = 0;
    max_inliers = 0;

    % Extract x and y coordinates
    x = edges(:, 1);
    y = edges(:, 2);

    % Set the random number generarator seed for reproducibility
    rng(7845);

    % RANSAC loop
    for i = 1:max_iterations
        % Randomly select 3 points (minimum required to define a circle)
        idx = randperm(size(edges, 1), 3);
        pts = edges(idx, :);

        % Fit a circle to the 3 points
        [cx_tmp, cy_tmp, r_tmp] = circle_from_three_points(pts);

        % Skip if the radius is invalid
        if isnan(r_tmp) || r_tmp <= 0
            continue;
        end

        % Compute distances of all points to the fitted circle
        distances = abs(sqrt((x - cx_tmp).^2 + (y - cy_tmp).^2) - r_tmp);

        % Count inliers within the threshold
        inliers = distances <= inlier_threshold;
        num_inliers = sum(inliers);

        % Update the best circle if the current one has more inliers...
        if num_inliers > max_inliers && num_inliers >= min_inliers
            % ... and has radius within the reasonable range
            if min_r < r_tmp && r_tmp < max_r
                best_cx = cx_tmp;
                best_cy = cy_tmp;
                best_r = r_tmp;
                max_inliers = num_inliers;
            end
        end
    end

    % Check if a valid circle was found
    if max_inliers >= min_inliers
        cx = best_cx;
        cy = best_cy;
        r = best_r;
    else
        % Return default values if no valid circle was found
        fprintf('Warning: RANSAC failed to find a valid circle.\n');
        cx = 0;
        cy = 0;
        r = 0;
    end
end

function [cx, cy, r] = circle_from_three_points(pts)
% Compute the circle passing through three points
% Input:
%   pts - 3x2 array of [x, y] coordinates
% Output:
%   [cx, cy, r] - Circle center and radius

    % Extract points
    x1 = pts(1, 1); y1 = pts(1, 2);
    x2 = pts(2, 1); y2 = pts(2, 2);
    x3 = pts(3, 1); y3 = pts(3, 2);

    % Compute the perpendicular bisectors of (x1, y1)-(x2, y2) and (x2, y2)-(x3, y3)
    A = [x2 - x1, y2 - y1; x3 - x2, y3 - y2];
    b = 0.5 * [(x2^2 - x1^2 + y2^2 - y1^2); (x3^2 - x2^2 + y3^2 - y2^2)];

    % Solve for the circle center
    if abs(det(A)) < 1e-10 % Check for degeneracy
        cx = NaN; cy = NaN; r = NaN;
        return;
    end
    center = A \ b;

    % Extract center coordinates
    cx = center(1);
    cy = center(2);

    % Compute the radius
    r = sqrt((x1 - cx)^2 + (y1 - cy)^2);
end
