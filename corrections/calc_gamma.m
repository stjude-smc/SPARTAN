function gamma = calc_gamma(input, flag)
% calc_gamma  Estimate donor to acceptor relative apparent brightness.
%
%   GAMMA = calc_gamma( INPUT ) estimates the apparent brightness of donor
%   fluorophore relative to the acceptor (GAMMA). INPUT may be a TracesFret 
%   object or the path to a .traces file.
%
%   ALGORITHM:
%   Ideally, donor and acceptor fluorophores appear equally bright, such that
%   the total fluorescence intensity (D+A) is constant regardless of FRET value.
%   In practice, the fluorophores' quantum yield and detection efficiency
%   ('relative brightness') are not identical, which leads to inaccurate FRET
%   calculations and fluctuating total fluorescence intensity.
%
%   Gamma is the apparent donor brightness relative to the acceptor.
%   Here it is estimated from the ratio of increase in donor intensity to the
%   decrease in acceptor intensity upon acceptor photobleaching.
%   See Ha et al (1999) PNAS 95, p.893. Roy et al (2008), Nat Meth
%   FIXME: 
%
%   CAVEATS:
%   A threshold of 0.2 (before correction) is used to find molecules
%   showing FRET. If all molecules are below this threshold, the algorithm
%   will fail. In this case, you could scale manually first and then
%   fine-tune by running this function.
% 
%   See also: gammacorrect, scaleacceptor, crosstalkcorrect.

%   Copyright 2014-2017 Cornell University All Rights Reserved.


%% Process input arguments
narginchk(0,1);
nargoutchk(0,1);

if nargin<1 || isempty(input),
    % If no files given, obtain a list from the user.
    input = getFile('*.traces');
    if isempty(input), return; end  %user hit cancel.
end

if ischar(input)
    data = loadTraces(input);
elseif isa(input,'TracesFret')
    data = input;
else
    error('Input must be a TracesFret object or path to a .traces file.');
end


%% Estimate crosstalk for each trace
constants = cascadeConstants;

nTraces = size(data.donor,1);
gamma = NaN( nTraces,1 );  %NaN corresponds to no data.

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
        
        % Frames before (pre) and after (post) acceptor bleaching.
        pre_range  = acc_lt-21:acc_lt-1;  %before acceptor bleaching
        post_range = acc_lt+1:donor_lt-1;  %after  acceptor bleaching
        
        % Ignore some corner cases that give strange results, including donor
        % or acceptor blinking during the window.
        if any( fret(pre_range(1):post_range(end))==0 ) || ...
           any( fret(post_range)>0.2 ),
            continue;
        end
        
        % Estimate fold changes in fluorophore intensity upon acceptor
        % bleaching, and gamma as the ratio.
        % NOTE: this is inverted from the Roy 2008 paper, which uses
        % del-A/delD.
        Id = data.donor(i,:);
        Ia = data.acceptor(i,:);
        
        delta_acc   = mean( Ia(pre_range)  )  -  mean( Ia(post_range) );
        delta_donor = mean( Id(post_range) )  -  mean( Id(pre_range)  );
        gamma(i) = delta_donor / delta_acc;
    end

end %for each trace


%% Filter values and get a reasonable population average

% If requested, just return all estimated values
if nargin>2 && ischar(flag) && strcmpi(flag,'all')
    return;
end

% Remove NaN values (from traces where crosstalk could not be estimated).
gamma = gamma( ~isnan(gamma) );

% Remove crosstalk values that are far out of the valid range.
% Without this, the std() may not be useful (outliers increase std).
gamma = gamma( gamma>0 & gamma<10 );

% Show what fraction of traces could be used to calculate gamma.
assert( ~isempty(gamma), 'No useful traces found for calculating gamma' );

percentUsed = 100*numel(gamma)/nTraces;
if percentUsed<5 || numel(gamma)<10,
    warning('Only a few traces (%d, %.0f%%) could be used for correction! May not be accurate.', ...
            numel(gamma),percentUsed );
else
    fprintf('\n%.0f%% of traces were used to calculate gamma.\n',percentUsed);
end

try
    % Gaussian fit to find mean gamma value
    x = 0:0.1:10;
    n = hist( gamma, x );
    f = fit( x', n', 'gauss1');
    gamma = f.b1;
    
    %figure; bar(x,n); hold on; plot(f);
catch
    % Failback if fitting fails: use distribution median.
    gamma = median(gamma);
end



end %function calc_gamma





