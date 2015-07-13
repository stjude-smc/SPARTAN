function centroids = getCentroids( image_t, picks, nhood )
% Find weighted center of each peak for more accurate PSF overlap detection.
% imdilate+regionprops works for this, but tends to merge nearby peaks.
% FIXME: when operating on a flat (zero) signal, gives NaN values. This should
% instead just return the input, possibly with a warning.

if nargin<3,
    nhood=1;  %3x3 region.
end

nMol = size(picks,1);
centroids = zeros(nMol,2);

for i=1:nMol,
    % Extract a window region around the molecule.
    x_window = picks(i,1) + (-nhood:+nhood);
    y_window = picks(i,2) + (-nhood:+nhood);
    block = image_t( y_window, x_window );
    
    % Calculate an intensity-weighted average position of molecule w/i window.
    tot = sum(block(:));
    x = sum( x_window .* sum(block,1)/tot );
    y = sum( y_window .* sum(block,2)'/tot  );
    
    centroids(i,:) = [x y];
end

% NaN values may appear for regions that are entirely black (zero) due to
% division by zero. Replace these values with the input and give a warning.
badPicks = isnan(centroids);
if any( badPicks(:) ),
    centroids(badPicks) = picks(badPicks);
    warning('gettraces:getCentroids:NaN','NaN values found when searching for centroids (%d, %.0f%%). This can happen when a field is empty (zero). Using input molecule locations instead for these molecules.', ...
            sum(badPicks(:)), 100*sum(badPicks(:))/numel(badPicks) );
end

end %FUNCTION getCentroids
