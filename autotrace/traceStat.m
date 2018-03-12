function [retval,nTraces] = traceStat( varargin )
% LOADTRACES  Loads fluorescence trace files
%
%   NAMES = traceStat;
%   Returns the names of each trace statistic as a structure.
%
%   STATS = traceStat( DONOR,ACCEPTOR,FRET, const )
%   Calculates metadata from fluorescence/FRET traces that can be
%   used as picking criteria to filter a set of traces.
%   Primarily used in autotrace.m
%
%   STATS = traceStat( FILES, const )
%   Same as above, but loads trace data from FILES, which may be a single
%   filename of a cell array of file names. In the latter case, properties from 
%   all files are combined as a single output result.
%   
%   Names of stats calculated by this function
%     corr          correlation b/t donor and acceptor signals
%     snr           signal-to-noise ratio (over background)
%     snr_s         signal-to-noise ratio (over total intensity)
%     bg            magnitude of background fluctuations/drift
%     t             average total fluorescence intensity
%     d             average donor intensity
%     a             average acceptor intensity
%     maxFRET       highest 1-frame FRET value in trace
%     ncross        number of donor dye blinks (total intensity~0)
%     lifetime      total fluorescence lifetime
%     acclife       number of frames showing FRET (in 5-frame chunks)
%     overlap       0=okay, 1=multiple photobleaching events (overlap)
%     avgfret       average FRET value
%     fretEvents    number of FRET events crossing E=0.14

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% PARAMETERS:
constants = cascadeConstants;  %get from varargin?
chunkSize = 5000;  %maximum number of traces to load at once
useWaitbar = ~( isfield(constants,'quiet') && constants.quiet );
wbh = [];


% If no data given, return names of filtering criteria.
% TODO: also return comments describing each criteria
if nargin<1,
    % Order matters -- will be displayed this way in autotrace!
    ln.t          = 'Mean Total Intensity';
    ln.maxFret    = 'Highest FRET value';
    ln.firstFret  = 'FRET at first frame';
    ln.fretEvents = 'Number of FRET events';
    ln.acclife    = 'FRET Lifetime';
    ln.donorlife  = 'Donor Lifetime';
    ln.lifetime = 'End of trace';
    ln.corr     = 'Correlation of Fluor.';
    ln.corrd    = 'Correlation of Fluor. Derivitive';
    ln.snr      = 'SNR-bg';
    ln.snr_s    = 'SNR-signal';
    ln.nnr      = 'SNR-bg / SNR-signal';
    ln.bg       = 'Background noise';
    ln.ncross   = '# Cy3 Blinks';
    ln.overlap  = 'Multi-step photobleaching';
    ln.safeRegion = 'Single-molecule start';
    ln.a        = 'Mean Acceptor Intensity';
    ln.d        = 'Mean Donor Intensity';
    
    ln.avgFret  = 'Average FRET value';
    
    % These parameters are specific to 3/4-color FRET with a second acceptor.
    % Ideally these should only appear for 3/4-color.
    ln.fret2Lifetime = 'Fret2 Lifetime';
    ln.maxFret2 = 'Highest Fret2 value';
    ln.avgFret2 = 'Average FRET2 value';
    
    retval = ln;
    return;
end


% If a filename, convert it into a list.
if ischar( varargin{1} ),
    varargin{1} = { varargin{1} };
end

% If the user gives a structure, this is a data structure with fluorescence
% data etc from loadTraces. Parse out the fields.
if isstruct( varargin{1} ) || isa(varargin{1},'Traces'),
    data = varargin{1};
    nTraces = size( data.donor, 1 );
    
    % Create a waitbar if this will take awhile.
    if useWaitbar && nTraces>400,
        wbh = parfor_progressbar( sum(nTraces), 'Loading traces and calculating properties...' );
    end
    
    retval = traceStat_data( data,constants );

% If the user gave filenames, load stats from each and combine them.
elseif iscell( varargin{1} ),
    files = varargin{1};
    nTraces = sizeTraces(files);
    retval = struct([]);
    
    % Create a waitbar if this will take awhile.
    if useWaitbar && sum(nTraces)>400,
        wbh = parfor_progressbar( sum(nTraces), 'Loading traces and calculating properties...' );
    end
    
    % Load each file and calculate statistics. If files have >5000 traces, they
    % are loaded in chunks to reduce memory usage.
    for file=1:numel(files),        
        for chunk=1:chunkSize:nTraces(file)
            data = loadTraces( files{file}, (chunk-1)+(1:chunkSize) );
            retval = [ retval traceStat_data( data,constants ) ];
        end
        
    end
   
% Assume the user passed trace data directly.
else
    warning('This syntax is depricated and will not be avilable in future versions');
    
    % Passing in donor/acceptor/fret traces separately.
    if nargin>=3,
        data.donor    = varargin{1};
        data.acceptor = varargin{2};
        data.fret     = varargin{3};
        retval = traceStat_data( data, varargin{4:end} );
    else
        retval = traceStat_data( varargin{:} );
    end
end
    
% Close the waitbar, if it was created.
if ~isempty(wbh),
    close(wbh);
end


% end


%%
function retval = traceStat_data( data, constants )
%

% Handle the single-color case (no acceptors) by creating zero'd data.
if ~data.isChannel('acceptor'),
    data.acceptor = zeros( size(data.donor) );
    data.fret     = zeros( size(data.donor) );
end

[Ntraces,len] = size(data.fret);

% Extract data matrices to variables to make parfor happy.
donorAll    = data.donor;
acceptorAll = data.acceptor;
fretAll     = data.fret;

% Add second acceptor FRET channel if available.
isThreeColor = data.isChannel('fret2');
if isThreeColor,
    fret2All = data.fret2;
    acceptor2All = data.acceptor2;
else
    % Even if unused, parfor requires these to be made.
    % FIXME?: put fret2 elements in a separate loop to avoid parfor.
    fret2All = zeros(size(fretAll));
    acceptor2All = zeros(size(acceptorAll));
end


if nargin<2
    constants = cascadeConstants;
end



% Value are returned into a structure that will be later be saved into 
% Application data as filtering criteria
z = num2cell( zeros(1,Ntraces) );

retval = struct( ...
    'corr',  z, 'corrd', z, ...
    'snr',   z, 'snr_s', z, ...
    'nnr',   z, 'bg',    z, ...
    't', z, 'd', z, 'a', z, ...
    'maxFret',    z, 'ncross',   z, ...
    'lifetime',   z, 'acclife',  z, ...
    'donorlife',  z, 'overlap',  z, ...
    'safeRegion', z, 'avgFret',  z, ...
    'fretEvents', z, 'firstFret', z, ...
    'fret2Lifetime', z, 'maxFret2', z, ...
    'avgFret2', z );
    

%----- CALCULATE STATISTICS FOR EACH TRACE
if numel(fretAll)/2000 > 1000 && constants.enable_parfor,
    % For large computations, parallelize computation across threads.
    pool = gcp;
    M = pool.NumWorkers;
else
    % For small datasets, do everything in the GUI thread (regular for loop).
    M = 0;
end

% parfor (i=1:Ntraces, M)
for i=1:Ntraces
    % Extract trace in a way that makes parpool happy.
    donor    = donorAll(i,:);
    acceptor = acceptorAll(i,:);
    fret     = fretAll(i,:);
    total    = donor+acceptor;
    
    if isThreeColor,
        total = total + acceptor2All(i,:);
    end
    
    %---- Calculate donor lifetime
    % Find donor photobleaching event by finding *last* large drop in total
    % fluorescence. Median-filtering smooths out noise, gradient works
    % better than diff for finding the drops. The value of NSTD is optimal
    % with our data (~300 photons/frame), but another value may be needed
    % with very low intensity data?
    filt_total  = medianfilter(total,constants.TAU);
    dfilt_total = gradient1(filt_total);
    mean_dfilt_total = sum( dfilt_total )/len;
    std_dfilt_total  = std1( dfilt_total );
    
    % Exclude "outliers" from std (including bleaching steps). The std is meant
    % to measure noise, not also signal.
    outliers = abs(dfilt_total) > mean_dfilt_total + 6*std_dfilt_total;    
    thresh = mean_dfilt_total - constants.NSTD*std1( dfilt_total(~outliers) );  %this slows things down a bit...
    
    lt = find( dfilt_total<=thresh, 1,'last' );

    if ~isempty(lt) && lt<len,
        retval(i).lifetime = max(2,lt);
    else
        if mod(i,10)==0 && ~isempty(wbh), wbh.iterate(10); end
        continue; %everything else is hard to calc w/o lifetime!
    end
    
    
    %---- OVERLAP DETECTION: Find multiple photobleaching events
    % Much like the above, we find drops in fluorescence as photobleaching
    % events, but use a more sensitive threshold and only take events where
    % the intensity never returns to its previous level after the "drop",
    % as would be seen with multi-step photobleaching. This is less
    % sensitive than detecting changes in *average* level, but produces way
    % fewer false positives with real data (that have intensity drifting).

    % Find each *run* of points that cross the threshold because the
    % bleaching step may happen over 2 frames because of time averaging.
    thresh = mean_dfilt_total-constants.overlap_nstd*std_dfilt_total;
    events = find( diff(dfilt_total<=thresh)==1 )+1;
    
    % Remove events right at the beginning that could be spurious.
    events = events(events>5 & events<=lt);

    lastPB = 0;

    for j=1:length(events)-1,
        point = events(j);
        
        minb = min( filt_total(2:point-1 ) );
        maxa = max( filt_total(point+1:lt) );

        if maxa < minb,
            lastPB = point;
        end
    end

    retval(i).overlap = lastPB~=0;
    retval(i).safeRegion = max(1,lastPB+2);

    
    %---- Ignore regions where Cy3 is blinking
    s = lt+5;
    bg_range = s:min(s+constants.NBK,len);
    stdbg = std1( total(bg_range) );

    % Calculate background noise over the entire end of the trace.
    if lt+10 < len
        retval(i).bg = std1(donor(s:end))+std1(acceptor(s:end));
    end
    
    % Calculate number of Cy3 PB threshold crossings per frame    
    if lt<3
        donorRange = true(1,lt-1);
        retval(i).donorlife = lt-1;
    else
        % Calculate blinking total fluor cutoff value
        % Find start points where Cy3 blinks
        x = total(1:lt-2) <= constants.blink_nstd*stdbg;  %was lt-3
        retval(i).ncross = sum(  x & ~[0 x(1:end-1)]  );
        
        % Remove from consideration regions where Cy3 is dark
        donorRange = ~x;  %logical index w/o Cy3 blinking region
        
        % Save the length of the "donor-alive" region as the lifetime
        % of the donor fluorophore.
        retval(i).donorlife = sum(donorRange);
        
        % Remove falling edges of the Cy3 blinks
        donorRange = donorRange & [donorRange(2:end) 1] & [1 donorRange(1:end-1)];
    end
    nDonor = sum(donorRange);
    
    
    if nDonor >= 2, %was >2
        % Truncate to regions with donor alive
        donor    = donor(donorRange);
        acceptor = acceptor(donorRange);
    
        % Calculate average amplitudes
        retval(i).d = sum( donor    )/nDonor;
        retval(i).a = sum( acceptor )/nDonor;
        retval(i).t = retval(i).d + retval(i).a;
        
        % Calculate Signal-to-noise ratio and background noise    
        % should be using bg_range for this...
        if lt+10 < len
            retval(i).snr = retval(i).t/stdbg;  % assuming background corrected
        end
        
        % Calculate correlation of derivitive of donor/acceptor signals.
        % This is the method used by {Fei & Gonzalez, 2008}
        ccd = corrcoef( gradient1(donor), gradient1(acceptor) );
        retval(i).corrd = ccd(1,2);
        
        % Calculate correlation coefficient
        cc = corrcoef( donor, acceptor );
        retval(i).corr = cc(1,2);
        
        total_noise = std1( donor+acceptor );
        retval(i).snr_s = retval(i).t ./ total_noise;
        retval(i).nnr = total_noise ./ stdbg;
    end
    
    
    % Find regions (before Cy3 photobleach) where FRET is above a threshold.
    % Properties will be calculated only on these areas.
    if lt>1
        fretRange = fret(1:lt) >= constants.min_fret;

        % Filter the regions so that they must consist of more than 5
        % consecutive points above the threshold
        fretRange = rleFilter( fretRange, constants.rle_min );
        nFret = sum(fretRange);

        retval(i).acclife = nFret;
        
        if nFret>1,
            retval(i).avgFret = sum(fret(fretRange))/nFret;
        end
        
        % Similar calculations when there is a second acceptor (3/4-color).
        if isThreeColor,
            fret2 = fret2All(i,:);
            fretRange = fret2(1:lt) >= constants.min_fret;
            fretRange = rleFilter( fretRange, constants.rle_min );

            retval(i).fret2Lifetime = sum(fretRange);

            if sum(fretRange)>1
                retval(i).avgFret2 = mean( fret2(fretRange) );
            end
            retval(i).maxFret2 = max(fret2);
        end
    
        % Number of events crossing an arbitrary threshold
        % TODO?: additional filtering to detect only anticorrelated events?
        [result] = RLEncode( fret(1:lt) > constants.fretEventTreshold );
        retval(i).fretEvents = sum( result(:,1)==1 );
    end
    
    % Statistics based on FRET distribution
    retval(i).firstFret = fret(1);
    retval(i).maxFret = max(fret);
    
    % Update waitbar, once every 10 traces to reduce overhead.
    if mod(i,10)==0 && ~isempty(wbh),
        wbh.iterate(10);
    end
end


end %function traceStat_data


end
