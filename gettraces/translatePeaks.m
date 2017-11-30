function [picks,rejects] = translatePeaks( picks, fieldSize, tform )
% translatePeaks  apply transformations to peak locations
%
%    PEAKS = translatePeaks(PEAKS, SZ, TFORM) applies the transformation TFORM
%    (from maketform) to the matrix PEAKS (x coordinates in first column, y in
%    second), centered around an image of size SZ (nY,nX). Any transformed peak
%    locations that fall outside the image are removed automatically.
%
%    [PEAKS,REJ] = translatePeaks(PEAKS, SZ, TFORM) also returns a logical 
%    array marking transformed peaks that were removed because they fell
%    outside the image boundaries.
%
%    PEAKS = translatePeaks(PEAKS, SZ, QUAD) directly translates the peak
%    locations to neighboring fields in a side-by-side montaged image.
%    Possible values: 1=UL (no translation), 2=UR, 3=BL, 4=BR.

%   Copyright 2007-2017 Cornell University All Rights Reserved.


% Check input arguments
narginchk(2,3);
nargoutchk(1,2);
assert( size(picks,2)==2 );

% Define the size of the field-of-view.
nrow = fieldSize(1);
ncol = fieldSize(2);
rejects = false( size(picks,1),1 );

if nargin<3 || isempty(tform)
    % Nothing to do if no transform provided
    return;

elseif isnumeric(tform) || ischar(tform)
    % Translate into target field
    switch tform
        case {1,'L','UL'}  %nothing to do.
        case {2,'R','UR'}
                picks(:,1) = picks(:,1) + ncol;
        case {3,'B','BL'}
            picks(:,2) = picks(:,2) + nrow; % LL y
        case {4,'BR'}
            picks(:,1) = picks(:,1) + ncol; % LR x
            picks(:,2) = picks(:,2) + nrow; % LR y
        otherwise
            error('Invalid quadrant');
    end
    
else
    % Apply transformation to the peak locations.
    % We have to first center the peak locations at (0,0) so the rotation is
    % about the center of the image, then put it back afterward.
    % Peak (maxima) locations must be integers, so they are rounded.
    picks = transformPointsInverse(  tform,  [picks(:,1)-(ncol/2) picks(:,2)-(nrow/2)]  );
    picks = round( [picks(:,1)+(ncol/2) picks(:,2)+(nrow/2) ] );
    
    % Mark any peaks that now fall outside the field limits.
    rejects = picks(:,1)<3      | picks(:,2)<3       | ...
              picks(:,1)>ncol-2 | picks(:,2)>nrow-2;
end


end  %function splitPeaks
