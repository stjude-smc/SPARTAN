function stkData = getIntegrationWindows(stkData, params)
% For each molecule location in "peaks", find the most intense pixels in
% its immediate neighborhood (defined by params.nPixelsToSum). These
% regions are used by integrateAndSave() to sum most of the intensity for
% each peak and generate fluorescence-time traces.
% To minimize the contribution of nearby molecules, the molecules closest to the
% peak center are added first and progressively out to the edge.
%
% FIXME: instead of just returning the N best pixels, return the entire sorted
% neighborhood and choose the highest N pixels later. This would give more
% flexibility and information?

%   Copyright 2007-2015 Cornell University All Rights Reserved.


stk_top = stkData.stk_top-stkData.background;
Npeaks = size(stkData.peaks,1);

% Define regions over which to integrate each peak
fractionWinOverlap = zeros(Npeaks,1);
[regions,integrationEfficiency,idxs,imgReused] = findRegions(stk_top, ...
                        stkData.peaks, params.nPixelsToSum, params.nhoodSize);
    
% If this entire neighborhood is empty, integrationEfficiency will be NaN.
% Give a warning, but leave the NaN.
badWindows = isnan(integrationEfficiency);
badWindows = any(badWindows'); %select if any entry per molecule is NaN.

if any(badWindows),
    warning('gettraces:getIntegrationWindows:NaN','NaN values found when getting integration windows (%d, %.0f%%). This can happen when a field is empty (zero).', ...
        sum(badWindows), 100*sum(badWindows)/Npeaks );
end
    
% For each peak, get the fraction of re-used pixels.
% FIXME: idxs could be used to construct regions to save time/memory.
for m=1:Npeaks,
    fractionWinOverlap(m) = sum(imgReused(idxs(:,m))-1) / params.nPixelsToSum;
end

stkData.regions = regions;
stkData.integrationEfficiency = integrationEfficiency;
stkData.fractionWinOverlap = fractionWinOverlap;


% Also get are "used" by overlapping spots
[~,~,~,mask] = findRegions(stk_top, stkData.rejectedPicks, ...
                         params.nPixelsToSum, params.nhoodSize);
stkData.bgMask = imgReused+mask==0;


end %function getIntegrationWindows



function [regions,int,idxs,mask] = findRegions(stk_top, peaks,nPx, hw)
% Get the action integration regions.

squarewidth = 1+2*hw;   % width of neighborhood to examine.
Npeaks  = size(peaks,1);

regions = zeros(nPx, 2, Npeaks);  %pixel#, dimension(x,y), peak#
int     = zeros(Npeaks,squarewidth^2);
idxs    = zeros(nPx, Npeaks);
mask    = zeros( size(stk_top) );  %marks where pixels are re-used

x = peaks(:,1);
y = peaks(:,2);

for m=1:Npeaks
    % Get a window of pixels around the intensity maximum (peak).
    nhood = stk_top( y(m)-hw:y(m)+hw, x(m)-hw:x(m)+hw );
    center = sort( nhood(:), 'descend' );
    
    % Find the most intense pixels, starting from the center and moving out.
    % This reduces the chance of getting intensity from nearby molecules.
    [A,B] = find( nhood>=center(nPx), nPx );

    % Estimate the fraction of intensity in each pixel,
    % relative to the total intensity in the full window region around the peak.
    % This is just an estimate and depends on the window size.
    % High molecule density can also distort this if there are overlapping PSFs.
    int(m,:) = cumsum( center/sum(center) )';
    
    % Convert to coordinates in the full FOV image and save.
    regions(:,:,m) = [ A+y(m)-hw-1, B+x(m)-hw-1  ];
    
    % Note where pixels are being reused for later calculation.
    idxs(:,m) = sub2ind( size(stk_top), regions(:,1,m), regions(:,2,m) );
    mask(idxs(:,m)) = mask(idxs(:,m)) +1;
end

end

