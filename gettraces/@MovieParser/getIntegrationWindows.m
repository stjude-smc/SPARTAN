function stkData = getIntegrationWindows(stkData, params)
% For each molecule location in "peaks", find the most intense pixels in
% its immediate neighborhood (defined by params.nPixelsToSum). These
% regions are used by integrateAndSave() to sum most of the intensity for
% each peak and generate fluorescence-time traces.
% To minimize the contribution of nearby molecules, the molecules closest to the
% peak center are added first and progressively out to the edge.

%   Copyright 2007-2017 Cornell University All Rights Reserved.


nPx = params.nPixelsToSum;
hw  = params.nhoodSize;
Npeaks = size(stkData.peaks,1);

% Define regions over which to integrate each peak
nCh = size(stkData.peaks,3);
[stkData.regionIdx,stkData.bgMask] = deal( cell(size(nCh,1)) );
intEff = 0;

for i=1:nCh
    field = stkData.stk_top{ params.idxFields(i) };
    [idxs,eff] = findRegions(field, stkData.peaks(:,:,i), nPx, hw);
    stkData.regionIdx{i} = idxs;
    intEff = intEff + eff;

    % Give a warning for any empty neighborhoods (eff is NaN)
    if any(isnan(eff(:))),
        warning('gettraces:getIntegrationWindows:NaN','Empty integration windows found. Zero-value field?');
    end

    % For each peak, get the fraction of re-used pixels.
    nUsed = diff( find([true;diff(sort(idxs(:)))~=0;true]) );
    stkData.fractionWinOverlap(i) = sum(nUsed>1) /Npeaks /nPx;
    
    % Create a mask of background areas (no PSF intensity)
    bgMask = true( size(field) );
    bgMask(idxs) = false;
    idxRej = findRegions(field, stkData.rejectedPicks(:,:,i), nPx, hw);
    bgMask(idxRej) = false;
    
    stkData.bgMask{i} = imerode(bgMask, ones(3));  % Remove PSF tails
end

% Calculate the fraction of total fluorescence intensity captured by summing
% the nPx most intense pixels ("integration efficiency").
intEff = bsxfun(@minus, intEff, min(intEff));   %avoid negative values
intEff = cumsum(intEff);
intEff = bsxfun(@rdivide, intEff, max(intEff));  %normalize
stkData.integrationEfficiency = 100*nanmean( intEff(nPx,:) );

% Estimate PSF size as number of pixels to get ~70% of total intensity.
stkData.psfWidth = mean( sum(intEff<0.7)+1 );

end %function getIntegrationWindows



function [idxs,eff] = findRegions(stk_top, peaks, nPx, hw)
% Get the action integration regions.

Npeaks = size(peaks,1);
squarewidth = 1+2*hw;   % width of neighborhood to examine.
eff  = zeros(squarewidth^2, Npeaks);  %integration efficiency
idxs = zeros(nPx, Npeaks);
win = 1:squarewidth;
peaks = peaks-hw-1;  %peak location relative to window edges

for m=1:Npeaks
    x = peaks(m,1);
    y = peaks(m,2);
    
    % Get a window of pixels around the intensity maximum (peak).
    nhood = stk_top( y+win, x+win );
    sortpx = sort( nhood(:), 'descend' );
    
    % Estimate fraction of PSF in window vs neighborhood (fraction collected).
    if nargout>1, eff(:,m)=sortpx; end
    
    % Find the nPx most intense pixels.
    [A,B] = find( nhood>=sortpx(nPx), nPx );
    
    % Convert to coordinates in the full FOV image and save linear indices.
    %idxs(:,m) = sub2ind( size(stk_top), A+y-hw-1, B+x-hw-1 );
    idxs(:,m) = (A+y) + ((B+x)-1).*size(stk_top,1);
end

end %FUNCTION findRegions

