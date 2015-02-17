function mean_gamma = gammacorrect(files,mean_gamma)
% gammacorrect scale acceptor intensity by estimated gamma values.
%
%   Corrects for unequal sensitivity and/or quantum yield of donor and acceptor
%   fluorophores in two-color FRET traces. A gamma value (ratio of donor to
%   acceptor apparent brightness) is calculated from traces where the acceptor
%   photobleaches first as the ratio of the magnitude of change in donor vs
%   acceptor fluorophores -- these ideally should be identical. See TJ Ha, 2004.
%   The output is saved to a new file specified by the user.
%
%       GAMMA = gammacorrect;   %will ask for files to process
%       GAMMA = gammacorrect( FILENAME );
% 
%   To use a known gamma value, rather than calculate from the data:
%
%       gammacorrect( FILENAME, GAMMA );
%
%   If more than one file is given, file names are chosen automatically as
%   "xxx_gammacorrect.traces" without an opportunity to change the name.
% 
% See also: scaleacceptor, crosstalkcorrect.
% 


% If no files given, obtain a list from the user.
if nargin<1 || isempty(files),
    filter = {'*.*traces','Traces files (*.traces,*.rawtraces)'; ...
              '*.txt','Old format traces files (*.txt)'; ...
              '*.*','All files (*.*)'};
    files = getFiles(filter);
end

if ~iscell(files),
    files = {files};
end

nFiles = numel(files);



if nargin>1 && isscalar(mean_gamma) && nFiles>1,
    % If a single gamma value is given, apply it to all traces.
    mean_gamma = repmat(mean_crosstalk,[nFiles,1]);
elseif nargin>1,
    assert( numel(mean_gamma)==nFiles );
else
    % If not value is given, it will be calculated. Allocate space only.
    mean_gamma = ones( nFiles,1 );
end


% For each file, estimate the value of gamma and adjust the donor intensity so
% that gamma=1. Then save the result to _gammacorrect.traces files.
for i=1:nFiles

    % Load fluorescence data
    data = loadTraces( files{i} );
    
    if isChannel(data,'acceptor2'),
        warning('This function may not work correctly with multi-color FRET');
    end
  
    % Estimate the mean gamma value across all traces in the file
    if nargin<2 && ~isempty(mean_gamma),
        mean_gamma(i) = calc_gamma(data);
    end
    
    % Scale donor intensity by gamma estimate (so final gamma=1)
    data.acceptor = data.acceptor*mean_gamma(i);
    data.recalculateFret();
     
    % Save resulting data
    [p,f,e] = fileparts( files{i} );
    outFilename = fullfile(p, [f '_gcorr' e]);
    
    if nFiles==1,
        [f,p] = uiputfile(outFilename,'Save corrected file');
        if ~ischar(f), return; end
        outFilename = fullfile(p,f);
    end
    
    saveTraces( outFilename, data );

end %for each file


end %function gammacorrect





%%
function mean_gamma_file = calc_gamma(data)
%

constants = cascadeConstants;

nTraces = size(data.donor,1);
gamma = NaN( nTraces,1 );  %NaN corresponds to no data.

% Loop through traces and find where the donor and acceptor photobleach. If
% the donor photobleaches within 10 frames of the acceptor, do not
% calculate gamma for those traces. Otherwise, calculate gamma as the
% change in donor intensity between the mean over 10 frames before acc photobleaching 
% (skipping the frame immediately before photobleaching because it is time-averaged) 
% and the mean over 10 frames after acc photobleaching (skipping the frame
% immediately after photobleaching because it is time-averaged) divided by the
% change in acceptor intensity from the mean over 10 frames before acc
% photobleaching (skipping the frame immediately before photobleaching) and zero.
for i = 1:nTraces
    fret = data.fret(i,:);
    
    % Find the point at which the donor photobleaches.
    donor_lt =  find( fret~=0, 1,'last' );
    if isempty(donor_lt),  continue;  end
    
    % Determine the points at which the donor and acceptor photobleach.
    % To reduce spurious "events", ignore runs of <5 frames above background.
    fretRange = fret >= 0.2;  %constants.min_fret;
    fretRange = rleFilter( fretRange, constants.rle_min );
    acc_lt = find( fretRange, 1, 'last');
    
    % Use only when the region is long enough for calculation.
    if ~isempty(acc_lt) && acc_lt > 21 && (donor_lt-acc_lt) > 21,
        Id = data.donor(i,:);
        Ia = data.acceptor(i,:);
    
        % Frames before (pre) and after (post) acceptor bleaching.
        pre_range  = acc_lt-21:acc_lt-1;  %before acceptor bleaching
        post_range = acc_lt+1:donor_lt-1;  %after  acceptor bleaching
        
        % Ignore some corner cases that give strange results, including donor
        % or acceptor blinking during the window.
        if any( fret(pre_range(1):post_range(end))==0 ) || ...
           any( fret(post_range)>0.2 ),
            continue;
        end
        
        % Estimate gamma
        delta_acc   = mean( Ia(pre_range)  )  -  mean( Ia(post_range) );
        delta_donor = mean( Id(post_range) )  -  mean( Id(pre_range)  );
        gamma(i) = delta_donor / delta_acc;
    end

end %for each trace


% Remove NaN values (from traces where crosstalk could not be estimated).
gamma = gamma( ~isnan(gamma) );

% Remove crosstalk values that are far out of the valid range.
% Without this, the std() may not be useful (outliers increase std).
gamma = gamma( gamma>0 & gamma<10 );

% Ignore any gamma estimates that are two standard deviations from the mean
% and return the mean gamma estimate.
% Consider taking the std of the middle 90% or something.
% median_gamma = median(gamma);
% std_gamma = std(gamma);
% gamma = gamma(gamma < (median_gamma + 2*std_gamma) & gamma > (median_gamma - 2*std_gamma));

% Show what fraction of traces could be used to calculate gamma.
assert( ~isempty(gamma), 'No useful traces found for calculating gamma' );

percentUsed = 100*numel(gamma)/nTraces;
if percentUsed<5,
    warning('Only a few traces (%d, %.0f%%) could be used for correction! May not be accurate.', ...
            numel(gamma),percentUsed );
else
    fprintf('\n%.0f%% of traces were used to calculate gamma.\n',percentUsed);
end

% Return an average value for apparent gamma.
mean_gamma_file = median(gamma);


end %function calc_gamma









