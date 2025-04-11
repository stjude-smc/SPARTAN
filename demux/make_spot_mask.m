function mask_img = make_spot_mask(nX, nY, cX, cY, Dia)
    % make_spot_mask - Creates a binary mask with one circular spot.
    
    % Initialize the mask image
    mask_img = zeros(nX, nY);
    
    % Draw the circles manually
    [YY, XX] = meshgrid(1:nY, 1:nX); % Create a grid for the image

    % Calculate the distance of each pixel to the circle's center
    distances = sqrt((XX - cX).^2 + (YY - cY).^2);
        
    % Set pixels within the circle's radius to 1
    mask_img(distances <= Dia/2) = 1;
end
