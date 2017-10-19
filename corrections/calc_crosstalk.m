function mean_crosstalk_file = calc_crosstalk(input)
% CALC_CROSSTALK  estimate fraction of donor to acceptor spectral overlap.
%
%   CT = calc_crosstalk( INPUT ) estimates the fraction of donor fluorescence
%   appearing on the acceptor channel (spectral crosstalk). INPUT may be a
%   TracesFret object or the path to a .traces file.
%
%   ALGORITHM:
%   Spectral crosstalk, where a fraction of donor fluorescence appears on the
%   acceptor channel, is estimated from the elevated acceptor channel baseline
%   after acceptor photobleaching. Bleaching is determined as the last point
%   above 0.2 FRET. This algorithm will not work properly if crosstalk is
%   > 15% or if there are non-zero FRET states with low efficiency also in that
%   range. In such cases, manually correct with an estimated value first.
%
%   See also: crosstalkcorrect, scaleacceptor, gammacorrect, correctTraces.

%   Copyright 2014-2017 Cornell University All Rights Reserved.


%% Process input arguments
narginchk(0,1);
nargoutchk(0,1);

if nargin<1
    input = getFile('*.traces;*.rawtraces');
end
if isempty(input), return; end
if ischar(input)
    data = loadTraces(input);
elseif isa(input,'Traces')
    data = input;
else
    error('Invalid input');
end


%% Calculate crosstalk values
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
crosstalk = crosstalk( crosstalk<1.3 & crosstalk>-1 );

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


