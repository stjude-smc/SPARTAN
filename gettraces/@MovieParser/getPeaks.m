function stkData = getPeaks(stkData, params)
% Localizes the peaks of molecules from a summed image of the two channels
% (in FRET experiments). The selection is made on the total fluorescence
% intensity image (summing all channels into a single image) to minimize
% bias (e.g., molecules with a dim donor and bright acceptor would be
% missed if we selected just on the donor side). The error in alignment is
% estimated by looking in the immediate neighborhood of each spot for the
% actual intensity maximum, which could be different on each side if the
% two fields are not aligned or if there are optical distortions.
%
% (If params.alignTranslate and alignRotate are set:) If the fields are not
% closely aligned, alignSearch() will try many possible alignments over a
% range of values, returning the best one. This process is pretty slow and
% the rotation adds noise to peak locations (picks), so rotation is turned
% off by default.
%
%    picks   = peak center positions for each channel listed in order:
%            For dual-color (FRET), it is Cy3,Cy5,Cy3,Cy5,...
%            For quad-color, UL,LL,LR,UR is the order.
%    total_t = total intensity image that was used for selecting peaks.
%    align   = alignment information as struct(dx,dy,theta,mag,abs_dev).
%            If software alignment is used, these numbers correspond
%            to the alignment applied.
%    total_picks = peak locations in the total intensity "channel" used for
%            picking in the first place.
%    fractionOverlapped = fraction of traces removed because they were too
%            close to a neighboring peak.
%    rejectedPicks = locations of peaks that will not be considered because
%            they are too close to a neighbor.
%    rejectedTotalPicks = rejected peaks in total intensity image.

% A note on notation: i, indD, indA, etc are indexes into the list of channels
% as they will appear in the output data (donor,acceptor). params.idxFields and
% quadrants identify the physical position of each channel on the camera chip.
% When looking into the image, use idxFields.


%% Process input arguments

% If the threshold for detecting intensity peaks is not given, calculate it
% automatically from the std of background regions at the end of the movie.
if ~params.don_thresh
    if ~isfield(params,'thresh_std')
        thresh_std = cascadeConstants('gettracesThresholdStd');
    else
        thresh_std = params.thresh_std;
    end

    % Calculate threshold from variance in background intensity at end of movie.
    params.don_thresh = thresh_std*stkData.stdbg;
%     params.don_thresh = thresh_std*mean(stkData.stdbg(params.idxFields));  %improved version
end


% Used by all methods
nCh = numel(params.chNames);
quality = zeros(nCh,1);
align = struct('dx',{},'dy',{},'theta',{},'sx',{},'sy',{},'abs_dev',{},'quality',{});
indD = find( strcmp(params.chNames,'donor') ); %donor channel to align to.
fields = stkData.stk_top(params.idxFields);



%% Single color or no multicolor without software alignment

if params.alignMethod==1 || isscalar(params.geometry)
    % Sum intensity from all channels without alignment
    total = sum( cat(3,fields{:}), 3 );
    
    % Find peaks of intensity in unregistered total intensity image.
    [total_picks,rejected] = pickPeaks( total, params.don_thresh, ...
                                     params.nhoodSize, params.overlap_thresh );
    picks = repmat( total_picks, [1 1 nCh] );
end 



%% Estimate local misalignment with software alignment is disabled.

if params.alignMethod==1 && ~isscalar(params.geometry)  % && numel(picks)>0,
    
    % Measure apparent deviation from perfect alignment
    refinedPicks = zeros( size(picks) );
    for i=1:nCh
        refinedPicks(:,:,i) = getCentroids( fields{i}, picks(:,:,i), params.nhoodSize );
    end
    
    for i=1:nCh
        if i==indD, continue; end %don't try to align donor to itself.
        
        % Estimate of misalignment error from the deviation of the position of
        % the highest intensity pixel in each field vs total intensity image.
        dev = refinedPicks(~rejected,:,i) - total_picks(~rejected,:);
        abs_dev = mean(  sqrt( dev(:,1).^2 + dev(:,2).^2 )  );

        % Use the picked peak locations to create a simple transformation
        % (including translation, rotation, and scaling) from donor to
        % each of the other fields.
        tform = fitgeotrans( refinedPicks(~rejected,:,indD), refinedPicks(~rejected,:,i), ...
                                                'NonreflectiveSimilarity' );
        ss = tform.T(2,1);
        sc = tform.T(1,1);
        scale = sqrt(ss*ss + sc*sc);
        theta = atan2(ss,sc)*180/pi;

        % Contrast score relative to a random alignment (scrambled inputs)
        quality(i) = weberQuality(fields{indD}, fields{i}, 0.7*params.don_thresh);

        align(i) = struct( 'dx',tform.T(3,1), 'dy',tform.T(3,2), 'theta',theta, ...
                'sx',scale, 'sy',scale, 'abs_dev',abs_dev, 'quality',quality(i) );
    end
end



%% Apply software alignment, if requested.
% alignMethod: 1=disable, 2=auto (ICP), 3=load from file, 4=memorize.

if ~isscalar(params.geometry) && params.alignMethod>1 % && numel(picks)>0,
    
    newAlign = struct('dx',{},'dy',{},'theta',{},'sx',{},'sy',{},'abs_dev',{},'tform',{});
    
    donor = fields{indD};  %reference field for alignment
    total = donor;         %total intensity of registered channel images
    [nrow,ncol] = size(donor);
    
    for i=1:nCh,
        if i==indD, continue; end %don't try to align donor to itself.
        target = fields{i};
        
        % Load memorized alignment
        if params.alignMethod==2 || params.alignMethod==4,
            newAlign(i) = params.alignment(i);
        
        % Run iterative closest points algorithm.
        elseif params.alignMethod==3,
            newAlign(i) = icpalign( donor, target, params );
        end
        
        % Register acceptor side so that it is lined up with the donor.
        % imref2d specifies the center of the image is the origin (0,0).
        R = imref2d( size(target), [-1 1]*ncol/2, [-1 1]*nrow/2 );
        registered = imwarp( target, R, newAlign(i).tform,...
                                    'Interp','cubic', 'OutputView',R );
        total = total + registered;
        
        % Contrast score relative to a random alignment (scrambled inputs).
        % FIXME: the threshold here is for total intensity, which may be much
        % brighter than the combination of any two channels. This could give 
        % low quality scores even when the alignment is good.
        quality(i) = weberQuality(donor, registered, 0.7*params.don_thresh);
    end

    % Pick peaks from the aligned total intensity image.
    [total_picks,rejected] = pickPeaks( total, params.don_thresh, ...
                                 params.nhoodSize, params.overlap_thresh );

    % Project molecule locations using the transformation above
    nPicked = size(total_picks,1);
    picks  = zeros( nPicked, 2, nCh );
    remove = false( nPicked,1 );

    for i=1:nCh,
        [picks(:,:,i),r] = translatePeaks( total_picks,  size(total), newAlign(i).tform );
        remove = remove | r;
    end

    % Remove peaks that fall outside of image after the transformation
    total_picks = total_picks(~remove,:,:);
    rejected    = rejected(~remove);
    picks       = picks(~remove,:,:);

    % Estimate of residual misalignment from the deviation of the position  of
    % the highest intensity pixel in each field vs registered total intensity image.
    for i=1:nCh
        goodPicks = picks(~rejected,:,i);  %non-overlapping PSFs
        refinedPicks = getCentroids( fields{i}, goodPicks, params.nhoodSize );
        
        dev = refinedPicks-goodPicks;
        newAlign(i).abs_dev = mean(  sqrt( dev(:,1).^2 + dev(:,2).^2 )  );
        newAlign(i).quality = quality(i);
    end
    align = newAlign;
end



%% Save output
stkData.rejectedTotalPicks = total_picks( rejected,: );
stkData.total_peaks        = total_picks( ~rejected,: );
stkData.rejectedPicks = picks( rejected,:,: );
stkData.peaks         = picks( ~rejected,:,: );

stkData.fractionOverlapped = sum(rejected)/numel(rejected);
stkData.total_t = total;
stkData.alignStatus = align;

% Reset any stale data from later steps
[stkData.regionIdx, stkData.integrationEfficiency, stkData.fractionWinOverlap, ...
 stkData.bgMask, stkData.psfWidth] = deal([]);
    

end %function getPeaks





