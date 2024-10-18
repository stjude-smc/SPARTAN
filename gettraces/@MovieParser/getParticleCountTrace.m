function stkData = getParticleCountTrace(stkData)
% Count number of distinguishable spots in each frame, creating a
% trejectory over time. This uses the same algorithm as getPeaks().
% Arbitrarily uses the first channel in the list (usually donor).


% If the threshold for detecting intensity peaks is not given, calculate it
% automatically from the std of background regions at the end of the movie.
if stkData.params.autoThresh
    stkData.params.don_thresh = stkData.params.thresh_std*stkData.stdbg;
%     stkData.params.don_thresh = stkData.params.thresh_std*mean(stkData.stdbg(idxFields));  %improved version
end
params = stkData.params;
stkData.chExtractor.verify();

% Mask regions outside the ROI (if any)
if ~isempty(stkData.roi)
    ROI = images.roi.Polygon('Position',stkData.roi);
    mask = ROI.createMask( stkData.chExtractor.nY, stkData.chExtractor.nX );
else
    mask = ones(stkData.chExtractor.nY, stkData.chExtractor.nX);
end

% if stkData.nChannels>1
%     errordlg('This function only works with single-channel data!');
%     return;
% end

% Count number of particles in every frame
bg = stkData.chExtractor.background{1};
output = zeros(stkData.nFrames,2);
output(:,1) = stkData.chExtractor.timeAxis/1000;  %convert to seconds
wbh = waitbar(0,'Counting particles...');

for i=1:stkData.nFrames    
    frame = stkData.chExtractor.read(i);
    field = (double(frame{1})-bg) .* mask;
    total_picks = pickPeaks( field, params.don_thresh, params.nhoodSize, params.overlap_thresh );
    output(i,2) = size(total_picks,1);
    waitbar(i/stkData.nFrames, wbh);
end
close(wbh);

% Save the counts to file with same name as source and .txt extension
[p,f] = fileparts( stkData.chExtractor.movie.filename );
save( fullfile(p, [f '_' mfilename '.txt']), 'output', '-ASCII' );

% Create a figure to show the result
figure;
plot( output(:,1), output(:,2) );
ylabel('Counts');
xlabel('Time (s)');
title(f);

end %function getParticleCountTrace





