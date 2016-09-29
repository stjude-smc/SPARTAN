function [mean_crosstalk] = crosstalkcorrect(files,mean_crosstalk)
% crosstalkcorrect subtract spectral crosstalk from fluorescence signals.
%
%  Fluorescence emission from one fluorophore can be detected on channels of
%  other fluorophores as an elevated baseline when the first fluorophore.
%  This elevated baseline is only present when the first fluorophore is
%  fluorescent and not after it photobleaches. For FRET traces, the crosstalk
%  value can be calculated as the residual acceptor intensity, if not zero,
%  after the acceptor photobleaches.
%
%  A single, average crosstalk value per file is used to make the correction.
%
%       CROSSTALK = crosstalkcorrect;   %will ask for files to process
%       CROSSTALK = crosstalkcorrect( FILENAME );
% 
%   To instead use a known crosstalk value:
%
%       crosstalkcorrect( FILENAME, CROSSTALK );
%
%   If more than one file is given, file names are chosen automatically as
%   "xxx_ccorr.traces" without an opportunity to change the name.
%   
%   NOTE: this may not work correctly with >15% crosstalk. Consider using an
%   initial correction with an approximate value to get it close and use the
%   script again to make a fine-tuned correction.
%
% See also: scaleacceptor, gammacorrect.

%   Copyright 2014-2016 Cornell University All Rights Reserved.


% If no files given, obtain a list from the user.
if nargin<1 || isempty(files),
    files = getFiles;
end

if ~iscell(files),
    files = {files};
end

nFiles = numel(files);


% If a single crosstalk value is given, apply it to all traces.
if nargin>=2 && isscalar(mean_crosstalk) && nFiles>1,
    % If a single gamma value is given, apply it to all traces.
    mean_crosstalk = repmat(mean_crosstalk,[nFiles,1]);
elseif nargin>1,
    assert( numel(mean_crosstalk)==nFiles );
else
    % If not value is given, it will be calculated. Allocate space only.
    mean_crosstalk = zeros( nFiles,1 );
end


% For each file, estimate the mean donor->acceptor fluorescence crosstalk, then
% subtract the crosstalk, recalculate FRET, and save the result.
% The remaining crosstalk should be close to zero.
for i=1:nFiles
    % Load the file and estimate crosstalk
    data = loadTraces( files{i} );
    
    if isChannel(data,'acceptor2'),
        disp('Warning: this function may not work correctly with multi-color FRET');
    end
    
    % Estimate the crosstalk if not given by the user.
    if nargin<2 && ~isempty(mean_crosstalk),
        mean_crosstalk(i) = calc_crosstalk(data);
    end
    
    % Add values to existing crosstalk matrix
    idxD  = strcmpi(data.channelNames,'donor');
    idxA1 = strcmpi(data.channelNames,'acceptor');
    crosstalk = cat(3, data.traceMetadata.crosstalk);
    crosstalk(idxD,idxA1,:) = crosstalk(idxD,idxA1,:) + mean_crosstalk(i);
    
    % Make the crosstalk correction and recalculate FRET
    data = correctTraces(data, crosstalk);
    data.recalculateFret();
     
    % Save resulting data
    [p,f,e] = fileparts( files{i} );
    outFilename = fullfile(p, [f '_ccorr' e]);
    
    if nFiles==1,
        [f,p] = uiputfile(outFilename,'Save corrected file');
        if ~ischar(f), return; end
        outFilename = fullfile(p,f);
    end
    
    saveTraces( outFilename, data );

end %for each file



end %function crosstalkcorrect





%%
function [mean_crosstalk_file] = calc_crosstalk(data)
% Find the fraction of acceptor fluorescence resulting from crosstalk from the
% donor channel due to incomplete spectral separation. The output is the mean 
% crosstalk estimate from all traces in the given file.
%


constants = cascadeConstants;
nTraces = size(data.donor,1);

% Calculate crosstalk values as the fraction of intensity on the acceptor
% channel relative to the donor, after acceptor photobleaching. If there is no
% crosstalk, the acceptor intensity should be zero (and crosstalk is also zero).
crosstalk = NaN( nTraces,1 ); %NaN corresponds to no data.

for i = 1:nTraces
    fret = data.fret(i,:);
    
    % Find the point at which the donor photobleaches.
    donor_lt =  find( fret~=0, 1,'last' );
    if isempty(donor_lt),  continue;  end
    
    % Find regions where the acceptor is clearly dark. We want to avoid the
    % blinking regions because they may not be dark partially quenched.
    fretRange = rleFilter( fret>=0.2, constants.rle_min );  % ignore 1-frame "events"
    acc_lt = find(fretRange, 1, 'last');
    
    if isempty(acc_lt),
        % If there is no acceptor, we can use the whole trace.
        range = 1:donor_lt-1;
    else
        range = acc_lt+1:donor_lt-1;
    end
    
    % Avoid donor blinking, which confuses the calculation.
    if any( fret(range)==0 ),  continue;  end
    
    % Estimate crosstalk as the fraction of intensity in the acceptor channel
    % after photobleaching (may be negative). Ignore short traces.
    if numel(range) >= 8
        crosstalk(i) = mean(data.acceptor(i,range)) / mean(data.donor(i,range));
    end
end

% Remove NaN values (from traces where crosstalk could not be estimated).
crosstalk = crosstalk( ~isnan(crosstalk) );

% Remove crosstalk values that are far out of the valid range.
% Without this, the std() may not be useful (outliers increase std).
crosstalk = crosstalk( crosstalk<1 & crosstalk>-1 );

% Ignore any crosstalk estimates that are two standard deviations from the mean
% and return the mean crosstalk estimate.
% Consider taking the std of the middle 90% or something.
% median_crosstalk = median(crosstalk);
% std_crosstalk = std(crosstalk);
% crosstalk = crosstalk(crosstalk < (median_crosstalk + 2*std_crosstalk) & crosstalk > (median_crosstalk - 2*std_crosstalk));


% Show what fraction of traces could be used to calculate gamma.
assert( ~isempty(crosstalk), 'No useful traces found for calculating gamma' );

percentUsed = 100*numel(crosstalk)/nTraces;
if percentUsed<5,
    warning('Only a few traces (%d, %.0f%%) could be used for correction! May not be accurate.', ...
            numel(crosstalk),percentUsed );
else
    fprintf('\n%.0f%% of traces were used to calculate crosstalk.\n',percentUsed);
end

% Return an average value for apparent crosstalk value.
mean_crosstalk_file = median(crosstalk);


end %function calc_crosstalk






