function centroids = getCentroids( image_t, picks, nhood )
% Find weighted center of each peak for more accurate PSF overlap detection.
% imdilate+regionprops works for this, but tends to merge nearby peaks.
% FIXME: when operating on a flat (zero) signal, gives NaN values. This should
% instead just return the input, possibly with a warning.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


if nargin<3,
    nhood=1;  %3x3 region.
end

nMol = size(picks,1);
centroids = zeros(nMol,2);
win = -nhood:nhood;
blocks = zeros(numel(win)^2,nMol);

% Extract a window region around the molecule.
for i=1:nMol,
    block = image_t( picks(i,2)+win, picks(i,1)+win );
    blocks(:,i) = block(:);
end

% Substract the minimum of each block so the floor is near zero
blocks = bsxfun(@minus, blocks, min(blocks) );  %blocks-min(blocks)
t = sum(blocks);

% Calculate center of mass for each peak within its local window
[ii,jj]=ndgrid(win); ii=ii(:); jj=jj(:);
centroids(:,1) = sum( bsxfun(@times, jj, blocks) )./t;   %sum( (jj.*blocks) )./t;
centroids(:,2) = sum( bsxfun(@times, ii, blocks) )./t;   %sum( (ii.*blocks) )./t;

% Add above deviations to peak locations to get centroids
centroids = centroids + picks;

% NaN values may appear for regions that are entirely black (zero) due to
% division by zero. Replace these values with the input and give a warning.
badPicks = isnan(centroids);
if any( badPicks(:) ),
    centroids(badPicks) = picks(badPicks);
    warning('gettraces:getCentroids:NaN','NaN values found when searching for centroids (%d, %.0f%%). This can happen when a field is empty (zero). Using input molecule locations instead for these molecules.', ...
            sum(badPicks(:)), 100*sum(badPicks(:))/numel(badPicks) );
end

end %FUNCTION getCentroids
