function [picks,boolRejected,centroids] = pickPeaks( image_t, threshold, nhood, overlap_thresh )
% Localizes the peaks of fluorescence.
%   picks = locations (x,y) of all molecules selected.
%   rejectedPicks = locations of molecules that are too close to a neighbor
%                        and should be ignored in analysis.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


[nrow,ncol] = size(image_t);

% Detect molecules as fluorescence maxima over local 3x3 regions,
% ignoring any that are below the detection threshold or near the edges.
% FIXME: may want to remove edges at the end, after overlap detection.
kernel = [0 1 0 ; 1 1 1; 0 1 0];

maxima = imregionalmax(image_t,kernel) & image_t>threshold;
[rows,cols] = find(maxima);
picks = [cols,rows];

ed = size(kernel,1);  %distance from the edge to avoid
edge = rows<=ed | rows>=nrow-ed | cols<=ed  | cols>=ncol-ed;
picks = picks(~edge,:);


% Find weighted center of each peak for more accurate PSF overlap detection.
% imdilate+regionprops works for this, but tends to merge nearby peaks.
% FIXME: consider using these centroid locations as picks. The only problem is
% that they msut be converted back to integers for indexing the image.
centroids = getCentroids( image_t, picks, nhood );
nMol = size(centroids,1);

% Detect maxima that are very close together and probably have overlapping
% point-spread functions.
if overlap_thresh==0 || nMol<2
    boolRejected = false(1,nMol);  %no overlap rejection.
else
    [~,dist] = knnsearch( centroids, centroids, 'k',2 );
    boolRejected = dist(:,2)'<=overlap_thresh;
end



end %FUNCTION pickPeaks




