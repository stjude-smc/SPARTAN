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
    
    % Fraction slider
    uicontrol('Style', 'text', 'Position', [20, 50, 100, 20], 'String', 'Data fraction');
    frac_slider = uicontrol('Style', 'slider', 'Position', [20, 30, 100, 20], ...
        'Min', 0.1, 'Max', 1, 'Value', fraction, ...
        'Callback', @(src, event) update_thumbnail());

    % Low threshold slider
    uicontrol('Style', 'text', 'Position', [120, 50, 100, 20], 'String', 'Low Threshold');
    low_slider = uicontrol('Style', 'slider', 'Position', [120, 30, 100, 20], ...
        'Min', 0, 'Max', 1, 'Value', low_thresh, ...
        'Callback', @(src, event) update_edges());

    % High threshold slider
    uicontrol('Style', 'text', 'Position', [220, 50, 100, 20], 'String', 'High Threshold');
    high_slider = uicontrol('Style', 'slider', 'Position', [220, 30, 100, 20], ...
        'Min', 0, 'Max', 1, 'Value', high_thresh, ...
        'Callback', @(src, event) update_edges());

    % Sigma slider
    uicontrol('Style', 'text', 'Position', [320, 50, 100, 20], 'String', 'Sigma');
    sigma_slider = uicontrol('Style', 'slider', 'Position', [320, 30, 100, 20], ...
        'Min', 1, 'Max', 20, 'Value', sigma, ...
        'Callback', @(src, event) update_edges());

    % s1 slider
    uicontrol('Style', 'text', 'Position', [420, 50, 100, 20], 'String', 's1');
    s1_slider = uicontrol('Style', 'slider', 'Position', [420, 30, 100, 20], ...
        'Min', 2, 'Max', 4, 'Value', s1, 'SliderStep', [1, 1], ...
        'Callback', @(src, event) update_s1());    

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
        location_img = zeros(nY, nX);
    
        % Set ones at the specified x and y locations
        for i = 1:length(x)
            location_img(y(i), x(i)) = 1;
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

    function update_s1()
        % Get slider value and snap to 2 or 4
        s1 = round(get(s1_slider, 'Value'));
        
        % Ensure the slider snaps only to valid values
        if s1 < 3
            s1 = 2;
        else
            s1 = 4;
        end
        
        % Update the slider position (in case of snapping)
        set(s1_slider, 'Value', s1);
    
        % Recreate thumbnail and axes
        thumbnail = prep_thumbnail(fraction);
    
        % Recreate the axes
        clf; % Clear the figure
        subplot(1, 2, 1);
        thumb_ax = imshow(thumbnail, []);
        title('Original Thumbnail');
    
        % Calculate initial Canny edges
        edges = edge(thumbnail, 'Canny', [low_thresh, high_thresh], sigma);
    
        % Recompute the mask for the new thumbnail size
        msk = make_spot_mask(edges, 160/s1/s2, 140/s1/s2) - make_spot_mask(edges, 170/s1/s2, 70/s1/s2);
    
        % Display the edges
        subplot(1, 2, 2);
        edge_ax = imshow(edges - edges * 0.7 .* (~msk) + 0.3 * msk);
        title(sprintf('Canny Edges [%.2f, %.2f], Sigma: %.2f', low_thresh, high_thresh, sigma));
    
        % Redraw sliders (recreating the figure deletes them)
        draw_sliders();
    end
    
    function draw_sliders()
        % Fraction slider
        uicontrol('Style', 'text', 'Position', [20, 50, 100, 20], 'String', 'Data fraction');
        frac_slider = uicontrol('Style', 'slider', 'Position', [20, 30, 100, 20], ...
            'Min', 0.1, 'Max', 1, 'Value', fraction, ...
            'Callback', @(src, event) update_thumbnail());
    
        % Low threshold slider
        uicontrol('Style', 'text', 'Position', [120, 50, 100, 20], 'String', 'Low Threshold');
        low_slider = uicontrol('Style', 'slider', 'Position', [120, 30, 100, 20], ...
            'Min', 0, 'Max', 1, 'Value', low_thresh, ...
            'Callback', @(src, event) update_edges());
    
        % High threshold slider
        uicontrol('Style', 'text', 'Position', [220, 50, 100, 20], 'String', 'High Threshold');
        high_slider = uicontrol('Style', 'slider', 'Position', [220, 30, 100, 20], ...
            'Min', 0, 'Max', 1, 'Value', high_thresh, ...
            'Callback', @(src, event) update_edges());
    
        % Sigma slider
        uicontrol('Style', 'text', 'Position', [320, 50, 100, 20], 'String', 'Sigma');
        sigma_slider = uicontrol('Style', 'slider', 'Position', [320, 30, 100, 20], ...
            'Min', 1, 'Max', 20, 'Value', sigma, ...
            'Callback', @(src, event) update_edges());
    
        % s1 slider
        uicontrol('Style', 'text', 'Position', [420, 50, 100, 20], 'String', 's1');
        s1_slider = uicontrol('Style', 'slider', 'Position', [420, 30, 100, 20], ...
            'Min', 2, 'Max', 4, 'Value', s1, 'SliderStep', [1, 1], ...
            'Callback', @(src, event) update_s1());    
    end

end


