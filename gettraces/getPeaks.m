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
% 

constants = cascadeConstants;
image_t = stkData.stk_top - stkData.background;


% If the threshold for detecting intensity peaks is not given, calculate it
% automatically from the std of background regions at the end of the movie.
% FIXME: This must be recalcualted after software alignment, if applicable.
if ~params.don_thresh
    if ~isfield(params,'thresh_std')
        thresh_std = constants.gettracesThresholdStd;
    else
        thresh_std = params.thresh_std;
    end

    % Calculate threshold from variance in background intensity at end of movie.
    % FIXME: this does not work well when background levels change during the
    % movie, for example with injection of Cy5-labeled tRNA.
    endBG = sort( stkData.endBackground(:) );
    endBG_lowerHalf = endBG( 1:floor(numel(endBG)*0.75) );
    params.don_thresh = thresh_std*std( endBG_lowerHalf );
else
    params.don_thresh = params.don_thresh-mean2(stkData.background);
end


% A note on notation: i, indD, indA, etc are indexes into the list of channels
% as they will appear in the output data (donor,acceptor). params.idxFields and
% quadrants identify the physical position of each channel on the camera chip.
% When looking into the image, use idxFields.
quadrants = params.idxFields;
channelNames = params.chNames;
nCh = numel(channelNames);  %# of channels TO USE.


% Define each channel's dimensions and sum fields together.
fields = subfield(image_t, params.geometry);
allFields = cat(3, fields{:});

% Sum fields to get a total intensity image for finding molecules.
% For now, we assume everything is aligned. FIXME: if an alignment file is
% loaded, there's no need to check first.
total_t = sum( allFields(:,:,quadrants), 3 );
[nrow,ncol] = size(total_t); %from now on, this is the size of subfields.



%---- 1. Pick molecules as peaks of intensity from summed (D+A) image)
[total_picks,rejected] = pickPeaks( total_t, params.don_thresh, ...
                                     params.nhoodSize, params.overlap_thresh );
nPicked = size(total_picks,1);

assert( all(total_picks(:))>0, 'bad peak locations' );


% Predict the locations of peaks in all channels.
if params.geometry==1,
    picks = total_picks;

elseif params.geometry>1,
    % For each channel, predict peak locations as straight translations 
    % from the peak locations in the total intensity image.
    picks = zeros(nCh*nPicked, 2);
    
    for i=1:nCh,  
        picks(i:nCh:end,:) = translatePeaks( total_picks, size(total_t), quadrants(i) );
    end
end

%%%%% Estimation of misalignment (for dual or multi-color)
align = struct('dx',{},'dy',{},'theta',{},'sx',{},'sy',{},'abs_dev',{},'quality',{});
indD = find( strcmp(channelNames,'donor') ); %donor channel to align to.
quality = zeros(nCh,1);

if params.geometry>1 && numel(picks)>0,
    % Refine peak locations and how much they deviate. This helps determine
    % if realignment is needed. FIXME: this could be improved (?) by using
    % centroid locations for translatePeaks above?
    refinedPicks = getCentroids( image_t, picks, params.nhoodSize );
    r_mod = [ mod(refinedPicks(:,1),ncol) mod(refinedPicks(:,2),nrow) ];
    residuals = refinedPicks-picks;
    
    % For each channel, find a crude alignment using control points. This
    % helps determine if software alignment is needed.
    d = r_mod(indD:nCh:end,:);  %donor (reference) points
    donor_t = allFields( :,:, params.idxFields(indD) ); %target field to align to
    
    for i=1:nCh,
        if i==indD, continue; end %don't try to align donor to itself.
        
        % Calculate mean deviations from donor->acceptor fields.
        dev = residuals(i:nCh:end,:)-residuals(indD:nCh:end,:);
        abs_dev = mean(  sqrt( dev(:,1).^2 + dev(:,2).^2 )  );

        % Use the picked peak locations to create a simple transformation
        % (including translation, rotation, and scaling) from donor to
        % each of the other fields.
        a = r_mod(i:nCh:end,:);

        tform = fitgeotrans( d(~rejected,:), a(~rejected,:), 'NonreflectiveSimilarity' );
        ss = tform.T(2,1);
        sc = tform.T(1,1);
        scale = sqrt(ss*ss + sc*sc);
        theta = atan2(ss,sc)*180/pi;

        % Measure the "quality" of the alignment as the magnitude increase in
        % score compared to a "random" alignment.
        target_t = allFields(:,:,params.idxFields(i));
        quality(i) =  weberQuality(donor_t,target_t,0.7*params.don_thresh);

        align(i) = struct( 'dx',tform.T(3,1), 'dy',tform.T(3,2), 'theta',theta, ...
                'sx',scale, 'sy',scale, 'abs_dev',abs_dev, 'quality',quality(i) );
    end
end


%%%%% Optional software alignment algorithm (for dual-color only!)
% If the alignment is close, no need to adjust.
% Just give a warning unless asked to do software alignment in settings.
if params.geometry>1 && params.alignMethod>1 && numel(picks)>0,
    % FIXME (?): the user may expect the alignment to be applied even if
    % the deviation is small when a specific alignment is loaded!
    
    % 
    %---- 2. Transform acceptor side so that the two channels align properly.
    
    % Try out all possible alignments within a range and find the one with
    % the best donor-acceptor intensity overlap. The quality score is the
    % mean aligned peak magnitude vs random alignment. If the score is low,
    % reject it and just say "we don't know".
    registered_t = cell(nCh,1);
    newAlign = struct('dx',{},'dy',{},'theta',{},'sx',{},'sy',{},'abs_dev',{},'tform',{});
    tform = cell(nCh,1);
    
    donor_t = allFields( :,:, params.idxFields(indD) ); %target field to align to
    total_t = donor_t;
    
    for i=1:nCh,
        if i==indD, continue; end %don't try to align donor to itself.
        target_t = allFields(:,:,params.idxFields(i));
        
        % Search for an optimal alignment of the selected field vs donor.
        % tform moves the acceptor field to be aligned with the donor.
        if params.alignMethod==2 || params.alignMethod==4,
            % Nothing to search, just apply the alignment.
            newAlign(i) = params.alignment(i);
        
        elseif params.alignMethod==3,
            % Iterative closest point algorithm.
            newAlign(i) = icpalign( donor_t, target_t, params );
        end
        
        % Register acceptor side so that it is lined up with the donor.
        % imref2d specifies the center of the image is the origin (0,0).
        R = imref2d( size(target_t), [-1 1]*nrow/2, [-1 1]*ncol/2 );
        registered_t{i} = imwarp( target_t,R, newAlign(i).tform,...
                                    'Interp','cubic', 'OutputView',R );
        
        total_t = total_t + registered_t{i};
        tform{i} = newAlign(i).tform;
        
        % Measure the "quality" of the alignment as the magnitude increase in
        % score compared to a "random" alignment.
        % FIXME: the threshold here is for total intensity, which may be much
        % brighter than the combination of any two channels. This could give 
        % low quality scores even when the alignment is good.
        quality(i) = weberQuality(donor_t,registered_t{i},0.7*params.don_thresh);
    end
    
    % If the optimal alignment is not trivial, re-pick molecule locations and
    % derive alignment deviation score.
    if ~( all([newAlign.dx]==0) && all([newAlign.dy]==0) && all([newAlign.theta]==0) ),
        % Give a warning for poor quality alignment.
        if any( quality<1.1 & quality>0 ),
            disp('Gettraces: Low confidence alignment. Parameters are out of range or data quality is poor.');
        end
        
        % Pick peaks from the aligned, total intensity image.
        [total_picks,rejected] = pickPeaks( total_t, params.don_thresh, ...
                                     params.nhoodSize, params.overlap_thresh );

        % For each channel, predict peak locations as straight translations 
        % from the peak locations in the total intensity image.
        nPicked = size(total_picks,1);
        picks = zeros( nCh*nPicked, 2 );
        remove = zeros( nPicked,1 );

        for i=1:nCh,
            [picks(i:nCh:end,:),r] = translatePeaks( total_picks, ...
                                 size(total_t), quadrants(i), tform{i} );
            remove = remove | r;
        end
        
        % Remove peaks that were moved outside the field boundries in at
        % least one of the fields due to the software alignment.
        total_picks  = total_picks(~remove,:);
        rejected     = rejected(~remove);
        nPicked = size(total_picks,1);

        % Remove molecules if the peaks in any field are out of range.
        z = repmat( remove', [nCh,1] ); z = z(:);
        picks = picks(~z,:);


        %---- 4. Re-estimate coordinates of acceptor-side peaks to verify alignment.
        % Then normalize so that the "expected" location of each peak is (0,0).
        % The rmsd is then the distance between each channel and the donor.
        refinedPicks = getCentroids( image_t, picks, params.nhoodSize );
        residuals = refinedPicks-picks; 
        
        for i=1:nCh,
            dev = residuals(i:nCh:end,:)-residuals(indD:nCh:end,:);
            dev = dev(~rejected,:);  %remove overlapped peaks.
            newAlign(i).abs_dev = mean(  sqrt( dev(:,1).^2 + dev(:,2).^2 )  );
            newAlign(i).quality = quality(i);
        end
        align = newAlign;
    end
end

% Save output
stkData.rejectedTotalPicks = total_picks(rejected,:);
stkData.total_peaks = total_picks(~rejected,:);

good = repmat( ~rejected, [nCh,1] ); good=logical(good(:));
stkData.rejectedPicks = picks( ~good,: );
picks = picks( good,: );

stkData.fractionOverlapped = size(stkData.rejectedTotalPicks,1) / nPicked;

stkData.total_t = total_t;
stkData.alignStatus = align;
stkData.peaks = picks;

end %function getPeaks


