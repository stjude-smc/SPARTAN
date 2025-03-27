function thumbnail = prep_traces_img(traces_xy)
    % Initialize a binary image with zeros
    location_img = zeros(int32(traces_xy.ny), int32(traces_xy.nx));

    % Set ones at the specified x and y locations
    for i = 1:length(traces_xy.x)
        location_img(traces_xy.y(i), traces_xy.x(i)) = 1;
    end

    % Resize
    thumbnail = imresize(location_img .* msk, 1/s1, 'bilinear');

    % Dilation to increase signal in circles
    thumbnail = imdilate(thumbnail, strel('disk', round(8/s1))); 

    % Downscale once again
    thumbnail = imresize(thumbnail, 1/s2, 'bilinear');
    thumbnail = thumbnail ./ max(thumbnail(:));
end