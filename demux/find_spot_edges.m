function pixel_coords = find_spot_edges(traces, ax)
    % Scale factors
    s1 = 2;
    s2 = 4;

    if ~isfield(traces, 'traceMetadata') || ~isfield(traces, 'fileMetadata')
        error('Invalid input: `traces` must contain `traceMetadata` and `fileMetadata` fields.');
    end

    x = [traces.traceMetadata.donor_x];
    y = [traces.traceMetadata.donor_y];

    % Get dimensions from file metadata
    nX = traces.fileMetadata.nX;
    nY = traces.fileMetadata.nY;

    % Initialize a binary image with zeros
    location_img = zeros(nX, nY);

    % Set ones at the specified x and y locations
    for i = 1:length(x)
        location_img(x(i), y(i)) = 1;
    end

    % Mask and resize
    msk = make_spot_mask(location_img, 160, 150);
    thumbnail = imresize(location_img .* msk, 1/s1, 'bilinear');

    % Dilation to increase signal in circles
    thumbnail = imdilate(thumbnail, strel('disk', round(4))); 

    % Downscale once again
    thumbnail = imresize(thumbnail, 1/s2, 'bilinear');
    thumbnail = thumbnail ./ max(thumbnail(:));

    % Apply Canny edge detection
    edges = edge(thumbnail, 'Canny', [0.05, 0.3], 6);

    % Mask canny edges
    msk = make_spot_mask(edges, 160/s1/s2, 140/s1/s2) - make_spot_mask(edges, 170/s1/s2, 70/s1/s2);
    edges = edges .* msk;

    % Assuming `edges` is your binary image
    [row_coords, col_coords] = find(edges);

    % Combine the row and column coordinates into a single array (optional)
    pixel_coords = [row_coords, col_coords] * s1 * s2;

    
    % Optional plotting if ax is provided
    if nargin > 1 && ~isempty(ax) && isgraphics(ax, 'axes')
        hold(ax, 'on'); % Ensure we can plot multiple elements

        % Plot detected edges
        imshow(imdilate(location_img, strel('disk', round(4))));
        plot(ax, pixel_coords(:, 2), pixel_coords(:, 1), 'y.');
        axis(ax, 'equal');
    end
end
