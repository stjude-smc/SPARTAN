function thumbnail = prep_traces_img(traces_xy)
    s1 = 4;
    s2 = 4;
    
    % Initialize a binary image with zeros
    location_img = zeros(int32(traces_xy.nY), int32(traces_xy.nX));

    % Set ones at the specified x and y locations
    for i = 1:length(traces_xy.X)
        location_img(traces_xy.Y(i), traces_xy.X(i)) = 1;
    end

    % Resize
    thumbnail = imresize(location_img, 1/s1, 'bilinear');
    %thumbnail = imresize(location_img .* msk, 1/s1, 'bilinear');

    % Dilation to increase signal in circles
    thumbnail = imdilate(thumbnail, strel('disk', round(8/s1))); 

    % Downscale once again
    thumbnail = imresize(thumbnail, 1/s2, 'bilinear');
    thumbnail = thumbnail ./ max(thumbnail(:));
end