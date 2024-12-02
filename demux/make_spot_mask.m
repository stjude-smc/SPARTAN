function mask_img = make_spot_mask(img, spot_distance, spot_dia)
    % make_spot_mask - Creates a binary mask with four circular spots.
    %
    % Syntax:
    % mask_img = make_spot_mask(img, spot_distance, spot_dia)
    %
    % Inputs:
    % img           - Input image to determine the mask size
    % spot_distance - Distance between the spots (in micrometers)
    % spot_dia      - Diameter of the spots (in micrometers)
    %
    % Output:
    % mask_img      - Binary mask with four circular spots
    
    % Get the size of the input image
    [nX, nY] = size(img);
    
    % Define the pixel size in micrometers
    pixel_size_um = 13 / 60; % 13 um pixel, 60x magnification
    
    % Convert diameter and distance to pixels
    circle_radius_px = spot_dia / pixel_size_um / 2;
    offset_px = spot_distance / pixel_size_um / 2;
    
    % Image center
    cx = nX / 2;
    cy = nY / 2;
    
    % Compute center positions for the circles (in pixels)
    circle_centers = [
        cx - offset_px, cy - offset_px;  % Bottom-left corner
        cx - offset_px, cy + offset_px;  % Bottom-right corner
        cx + offset_px, cy - offset_px;  % Top-left corner
        cx + offset_px, cy + offset_px;  % Top-right corner
    ];
    
    % Initialize the mask image
    mask_img = zeros(nX, nY);
    
    % Draw the circles manually
    [xx, yy] = meshgrid(1:nX, 1:nY); % Create a grid for the image
    for i = 1:size(circle_centers, 1)
        % Calculate the distance of each pixel to the circle's center
        distances = sqrt((xx - circle_centers(i, 1)).^2 + (yy - circle_centers(i, 2)).^2);
        
        % Set pixels within the circle's radius to 1
        mask_img(distances <= circle_radius_px) = 1;
    end
end
