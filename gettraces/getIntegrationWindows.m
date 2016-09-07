function stkData = getIntegrationWindows(stkData, params)
% For each molecule location in "peaks", find the most intense pixels in
% its immediate neighborhood (defined by params.nPixelsToSum). These
% regions are used by integrateAndSave() to sum most of the intensity for
% each peak and generate fluorescence-time traces.
% To minimize the contribution of nearby molecules, the molecules closest to the
% peak center are added first and progressively out to the edge.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


stk_top = stkData.stk_top-stkData.background;
Npeaks = size(stkData.peaks,1);

% Define regions over which to integrate each peak
[idxs,eff] = findRegions(stk_top, stkData.peaks, params.nPixelsToSum, ...
                                                 params.nhoodSize);
stkData.regionIdx = idxs;
stkData.integrationEfficiency = eff;

% Give a warning for any empty neighborhoods (eff is NaN)
if any(isnan(eff(:))),
    warning('gettraces:getIntegrationWindows:NaN','Empty integration windows found. Zero-value field?');
end
    
% For each peak, get the fraction of re-used pixels.
nUsed = diff( find([true;diff(sort(idxs(:)))~=0;true]) );
stkData.fractionWinOverlap = sum(nUsed>1) /Npeaks /params.nPixelsToSum;

% Create a mask of background areas (no PSF intensity)
stkData.bgMask = true( size(stk_top) );
stkData.bgMask(idxs) = false;
idxs = findRegions(stk_top, stkData.rejectedPicks, params.nPixelsToSum, ...
                                                   params.nhoodSize);
stkData.bgMask(idxs) = false;

end %function getIntegrationWindows



function [idxs,eff] = findRegions(stk_top, peaks, nPx, hw)
% Get the action integration regions.

Npeaks = size(peaks,1);
squarewidth = 1+2*hw;   % width of neighborhood to examine.
eff  = zeros(Npeaks,squarewidth^2);  %integration efficiency
idxs = zeros(nPx, Npeaks);

for m=1:Npeaks
    x = peaks(m,1);
    y = peaks(m,2);
    
    % Get a window of pixels around the intensity maximum (peak).
    nhood = stk_top( y+(-hw:hw), x+(-hw:hw) );
    center = sort( nhood(:), 'descend' );
    
    % Estimate fraction of PSF in window vs neighborhood (fraction collected).
    if nargout>1,
        eff(m,:) = cumsum( center/sum(center) )';
    end
    
    % Find the nPx most intense pixels.
    [A,B] = find( nhood>=center(nPx), nPx );
    
    % Convert to coordinates in the full FOV image and save linear indices.
    %idxs(:,m) = sub2ind( size(stk_top), A+y-hw-1, B+x-hw-1 );
    idxs(:,m) = (A+y-hw-1) + ((B+x-hw-1)-1).*size(stk_top,1);
end

end %FUNCTION findRegions

