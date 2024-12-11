function demux_optimize_canny_thresholds(traces_file, s1, s2)
    % Function to interactively tweak Canny edge thresholds and sigma
    arguments
        % Scale factors
        traces_file (1, :) char = getFile({'*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
                                           '*.traces','Binary Traces Files (*.traces)' ; ...
                                           '*.*','All Files (*.*)'});
        s1 (1, 1) double = 4;
        s2 (1, 1) double = 4;
    end

    % Initial values for Canny edges
    low_thresh = 0.05;
    high_thresh = 0.3;
    sigma = 5; % Initial sigma value
    fraction = 1;

    all_traces = loadTraces(traces_file);

    if ~isfield(all_traces, 'traceMetadata') || ~isfield(all_traces, 'fileMetadata')
        error('Invalid input: `traces` must contain `traceMetadata` and `fileMetadata` fields.');
    end

    thumbnail = prep_thumbnail(fraction);




    % Create a figure with sliders
    figure('Name', 'Interactive Canny Edge Detection');
    
    % Display the original thumbnail
    subplot(1, 2, 1);
    thumb_ax = imshow(thumbnail, []);
    title(sprintf('Thumbnail, fraction = %.2f', fraction));

    % Calculate initial Canny edges
    edges = edge(thumbnail, 'Canny', [low_thresh, high_thresh], sigma);

    % Mask for Canny edges
    msk = make_spot_mask(edges, 160/s1/s2, 140/s1/s2) - make_spot_mask(edges, 170/s1/s2, 70/s1/s2);

    % Display the initial edges
    subplot(1, 2, 2);

    edge_ax = imshow(edges - edges * 0.7 .* (~msk) + 0.3*msk);


    title(sprintf('Canny Edges [%.2f, %.2f], Sigma: %.2f', low_thresh, high_thresh, sigma));
    
    % Add sliders for data fraction
    uicontrol('Style', 'text', 'Position', [40, 50, 100, 20], 'String', 'Data fraction');
    frac_slider = uicontrol('Style', 'slider', 'Position', [40, 30, 100, 20], ...
        'Min', 0.1, 'Max', 1, 'Value', fraction, ...
        'Callback', @(src, event) update_thumbnail());

    % Add sliders for low threshold
    uicontrol('Style', 'text', 'Position', [160, 50, 100, 20], 'String', 'Low Threshold');
    low_slider = uicontrol('Style', 'slider', 'Position', [160, 30, 100, 20], ...
        'Min', 0, 'Max', 1, 'Value', low_thresh, ...
        'Callback', @(src, event) update_edges());

    % Add sliders for high threshold
    uicontrol('Style', 'text', 'Position', [280, 50, 100, 20], 'String', 'High Threshold');
    high_slider = uicontrol('Style', 'slider', 'Position', [280, 30, 100, 20], ...
        'Min', 0, 'Max', 1, 'Value', high_thresh, ...
        'Callback', @(src, event) update_edges());
    
    % Add sliders for sigma
    uicontrol('Style', 'text', 'Position', [400, 50, 100, 20], 'String', 'Sigma');
    sigma_slider = uicontrol('Style', 'slider', 'Position', [400, 30, 100, 20], ...
        'Min', 1, 'Max', 20, 'Value', sigma, ...
        'Callback', @(src, event) update_edges());
    
    % Nested function to update edges
    function update_edges()
        % Get the current slider values
        low_thresh = get(low_slider, 'Value');
        high_thresh = get(high_slider, 'Value');
        sigma = get(sigma_slider, 'Value');
        
        % Apply Canny edge detection
        edges = edge(thumbnail, 'Canny', [low_thresh, high_thresh], sigma);
        
        % Recompute the mask based on the new thumbnail size
        msk = make_spot_mask(edges, 160/s1/s2, 140/s1/s2) - make_spot_mask(edges, 170/s1/s2, 70/s1/s2);
    
        % Update the displayed edge image
        set(edge_ax, 'CData', edges - edges * 0.7 .* (~msk) + 0.3*msk);
        title(edge_ax.Parent, sprintf('Canny Edges [%.2f, %.2f], Sigma: %.2f', low_thresh, high_thresh, sigma));
    end


    function update_thumbnail()
        fraction = get(frac_slider, 'Value');
        thumbnail = prep_thumbnail(fraction);
        set(thumb_ax, 'CData', thumbnail);
        title(thumb_ax.Parent, sprintf('Thumbnail, fraction = %.2f', fraction));
    
        % Ensure mask and edges are updated to match the new thumbnail
        update_edges();
    end


    function thumbnail = prep_thumbnail(fraction)
        % Select a fraction of all traces
        mask = false(1, all_traces.nTraces);
        indices = randperm(all_traces.nTraces, floor(fraction*all_traces.nTraces));
        mask(indices) = true;
        traces = all_traces.getSubset(mask);
        % END DEBUG code

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
        thumbnail = imdilate(thumbnail, strel('disk', round(8/s1))); 
    
        % Downscale once again
        thumbnail = imresize(thumbnail, 1/s2, 'bilinear');
        thumbnail = thumbnail ./ max(thumbnail(:));
    end


end


