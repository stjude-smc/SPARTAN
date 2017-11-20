function align = icpalign( ref_img, target_img, params )
% icpalign Align images using the iterative closest point algorithm
%
%   ALIGN = icpalign( REFERRENCE, TARGET )
%   detects molecules as maxima above a threshold spots in each image
%   separately. For each point in the reference (fixed points), find the
%   closest point in the target ("moving" points) and use these pairs as
%   control points to estimate a transformation matrix. Target points are
%   moved using the transformation toward (hopefully) the reference points.
%   The process is iterated several times until convergence, allowing the
%   algorithm to "walk" toward the correction solution.
%
%   The algorithm is robust to uncertainty in coordinates, missing points
%   (donor-only or high-FRET), and spurious points, but it is best suited
%   to beads or patterned arrays that have high SNR and identical points
%   on all channels.
%
% See also: WEBERALIGN

%   Copyright 2007-2015 Cornell University All Rights Reserved.


[nrow,ncol] = size(ref_img);
assert( all(size(ref_img)==size(target_img)) );

% 1) Detect beads as brights spots above background.
% A lower threshold is used to better detect dim channels (Cy7).
% NOTE: All coordinates are X,Y, but indexing in images is Y,X
% FIXME: use per-channel thresholds to better handle intensity mismatches (Cy7).
[~,reject,fixed] = pickPeaks( ref_img, 0.5*params.don_thresh, params.nhoodSize, params.overlap_thresh );
fixed = fixed(~reject,:);

[~,reject,moving] = pickPeaks( target_img, 0.5*params.don_thresh, params.nhoodSize, params.overlap_thresh );
moving = moving(~reject,:);

% Translate coordinates so the center of the image is (0,0), so that
% rotation and scale will be about the center of the image.
% NOTE: a correct imref2d must be used with tform as a result.
fixed  = [ fixed(:,1)-(ncol/2)  fixed(:,2)-(nrow/2)  ];
moving = [ moving(:,1)-(ncol/2) moving(:,2)-(nrow/2) ];


% 2) Iterative closest point algorithm to estimate alignment transformation
% Parameters
maxitr = 200;   %maximum number of iterations to converge
gradtol = 1e-6; %minimum change in error before we assume convergence.

% Initialization
registered = moving;
err = zeros(maxitr,1);

for itr=1:maxitr,
    % For each "fixed" spot, find the nearest neighbor.
    % registered(idx) or moving(idx) gives list of spots corresponding to fixed.
    [idx,dist] = knnsearch( registered, fixed );
    
    % Initial error estimate before optimization:
    if itr==1,
        delta = (fixed-moving(idx,:)).^2;
        err0 = sum(  sqrt( delta(:,1)+delta(:,2) )  )/size(delta,1);
    end
    
    % Use point pairs to estimate an optimal transformation.
    % Points with a very large deviation, which may have no corresponding point
    % on the other side at all, are removed.
    sel = dist<params.nPixelsToSum;
    m = moving(idx,:);
    tform = fitgeotrans( m(sel,:), fixed(sel,:), 'nonreflectivesimilarity' );

    % Forward transform (about center of image) to see how we did.
    registered = transformPointsForward( tform, moving );

    % Estimated error: residual distance between nearest neighbors.
    delta = (fixed-registered(idx,:)).^2;
    err(itr) = sum(  sqrt( delta(:,1)+delta(:,2) )  )/size(delta,1);
    
    % Termination: test for convergance
    if itr>1 && abs(err(itr)-err(itr-1)) < gradtol,
        break;
    end
end
err = [err0; err(1:itr)];

% Display minimization of error and convergence of fit.
% figure;
% plot( 0:itr, err );
% title('ICP algorithm iterations');
% xlabel('Iteration number');
% ylabel('Average residual error (pixels)');

if itr==maxitr && maxitr>1,
    warning('ICP:ExceededMaxIterations','Did not converge within maximum number of iterations');
end

if err(itr+1)-min(err) > 1,
    warning('ICP:Diverged','Fitting diverged; may not have detected peaks on all channels properly');
end

% Display final alignment.
% figure;
% scatter( fixed(:,1),fixed(:,2),'bo' );
% hold on;
% scatter( moving(:,1),moving(:,2),'ro' );
% scatter( registered(:,1),registered(:,2),'gx' );
% hold off;
% legend( {'Fixed','Moving','Registered'} );


% Document the degree of residual error. FIXME
if ~params.quiet,
    fprintf('ICP converged after %d iterations (residual error %.2f px)\n',itr,err(itr+1));
    
    u = unique(idx);
    n = histc(idx,u); %count how many times each index is used
    fprintf('ICP: %.0f (%0.0f%%) were used; of those %.0f (%.0f%%) are ambiguous\n\n', ...
             numel(idx), 100*numel(idx)/numel(dist), sum(n>1), 100*sum(n>1)/numel(n) );
end

% 4) Calculate distortion parameters from tform.
% These parameters describe the distortion that created the misalignment,
% which is the inverse of the tform to align the acceptor (moving) to donor (fixed).
tinv = invert(tform);
ss = tinv.T(2,1);
sc = tinv.T(1,1);
tx = tinv.T(3,1);
ty = tinv.T(3,2);
scale = sqrt(ss*ss + sc*sc);
theta = atan2(ss,sc)*180/pi;

align = struct( 'dx',tx,'dy',ty, 'theta',theta, 'sx',scale,'sy',scale, ...
                'abs_dev',0, 'tform',tform );


end %FUNCTION alignSearch_cpt

