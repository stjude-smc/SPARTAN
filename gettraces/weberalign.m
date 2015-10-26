function [bestAlign,bestReg] = weberalign( donor, acceptor, params )
% weberalign Exhaustive alignment search, optimizing image Weber contrast.
%
%   ALIGN = weberalign( REFERENCE, TARGET )
%   performes a full search for a transformation (rotation+translation) of
%   the acceptor side relative to the donor. For each combination of dx,
%   dy, and dtheta (across a range in each), transform the acceptor side
%   image, sum the donor and acceptor images, and calculate a contrast
%   score. The one with the brightest, sharpest peaks wins.
%
% See also: icpalign, weberQuality.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


[nrow,ncol] = size(donor);

% Search parameters (dx and dy must be whole numbers):
dx_range = -8:1:8;
dy_range = -8:1:8;
theta_range = -1:0.1:1;

ntheta = numel(theta_range);
ndx = numel(dx_range);
ndy = numel(dy_range);

if ~params.quiet && ntheta>1,
    wbh = parfor_progress( ntheta,'Searching for the optimal alignment');
end

% Start the matlab thread pool if not already running, unless disabled.
if constants.enable_parfor,
    pool = gcp;
    M = pool.NumWorkers;
else
    M = 0;
end

% Reserve space for the best alignment for each possible rotation value.
% This is required for parfor to work correctly because all threads must be
% totally independent.
%scores = zeros(ntheta,ndx,ndy);
bestScores = zeros(ntheta,1);
bestAligns = cell(ntheta,1);
bestRegs = cell(ntheta,1);

parfor (t=1:ntheta, M)
    theta = theta_range(t);
    scoreTemp = zeros(ndx,ndy); %contrast scores for all translations.
    
    % Rotate the image, removing excess around the edges ('crop').
    % Rotation is done first for speed since the rotated image can be
    % reused for the translations. imrotate blurs out the image, which
    % lowers the contrast ratio, so there is a slight bias toward no
    % theta=0. Rotating both fields by some large amount outside the search
    % range can remove the bias, but in practice this isn't necessary.
    % There may also be a minor bias toward low FRET using this since the
    % acceptor side will have weaker intensity.
    rot_a = imrotate( acceptor, theta, 'bicubic', 'crop' );
    
    for i=1:ndx,
        dx = dx_range(i); %#ok<PFBNS>
            
        for j=1:ndy,
            dy = dy_range(j); %#ok<PFBNS>
            
            % Translate the image, also removing the excess.
            % FIXME: biased against moving too far because edges are zero?
            registered = zeros( size(acceptor) );

            registered( max(1,1-dy):min(nrow,nrow-dy), max(1,1-dx):min(ncol,ncol-dx) ) = ...
                rot_a( max(1,1+dy):min(nrow,nrow+dy), max(1,1+dx):min(ncol,ncol+dx) );

            % Calculate Weber contrast as a score. When the fields are not
            % aligned, the peak intensities go down and some peaks fall
            % below the threshold of background, lowering contrast.
            total = donor+registered;
            picks = total>params.don_thresh;    %pixels above background
            Ib = mean( total(~picks) );  %background intensity
            S = ( mean(total(picks)) -Ib ) / Ib;  %Weber contrast score
            
            % If this is the best alignment so far, save it.
            if all( S>scoreTemp(:) ),
                bestAligns{t} = struct('dx',dx,'dy',dy,'theta',theta,'sx',1,'sy',1,'abs_dev',0);
                bestRegs{t} = registered;
                bestScores(t) = S;
            end
            
            scoreTemp(i,j) = S;  %fixme for direct indexing.
        end
    end
    
    %scores(t,:,:) = scoreTemp;
    
    if ~params.quiet && ntheta>1,
        parfor_progress(wbh);
    end
end


% Over all possible rotations, find the best one.
[~,bestIdx] = max(bestScores);
bestAlign = bestAligns{bestIdx};
bestReg   = bestRegs{bestIdx};
        

% Create a transformation matrix from the alignment parameters.
scale=1; %FIXME: should be optimized also.
sc = scale*cosd(bestAlign.theta);
ss = scale*sind(bestAlign.theta);
T = [ sc           -ss            0;
      ss            sc            0;
      bestAlign.dx  bestAlign.dy  1];
bestAlign.tform = affine2d(T);


if ~params.quiet && ntheta>1,
    close(wbh);
end

end %FUNCTION alignSearch
