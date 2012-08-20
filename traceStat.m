function [retval,nTraces] = traceStat( varargin )
% LOADTRACES  Loads fluorescence trace files
%
%   [NAMES] = TRACESTAT;
%   Returns the names of each trace statistic as a structure.
%
%   STATS = TRACESTAT( DONOR,ACCEPTOR,FRET, const )
%   Calculates metadata from fluorescence/FRET traces that can be
%   used as picking criteria to filter a set of traces.
%   Primarily used in autotrace.m
%
%   STATS = TRACESTAT( FILES, const )
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
%
%   snr,snr_s,t are all corrected for gamma (sensitivity/quantum yield ratio)
%

% If no data given, return names of filtering criteria.
% TODO: also return comments describing each criteria
if nargin<1,
    % Order matters -- will be displayed this way in autotrace!
    ln.t        = 'Mean Total Intensity';
    ln.maxFRET  = 'Highest FRET value';
    ln.firstFRET = 'FRET at first frame';
    ln.fretEvents = 'Number of FRET events';
    ln.acclife  = 'FRET Lifetime';
    ln.donorlife  = 'Donor Lifetime';
    ln.lifetime = 'Donor time before bleaching';
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
    
    ln.avgfret  = 'Average FRET value';
    ln.timeSaturated = 'Frames with saturated intensity';
    
    retval = ln;
    return;
end


% If a filename, convert it into a list.
if ischar( varargin{1} ),
    varargin{1} = { varargin{1} };
end

% If the user gives a structure, this is a data structure with fluorescence
% data etc from loadTraces. Parse out the fields.
if isstruct( varargin{1} ),
    data = varargin{1};
    nTraces = size( data.donor, 1 );
    retval = traceStat_data( data.donor, data.acceptor, data.fret );

% If the user gave filenames, load stats from each and combine them.
elseif iscell( varargin{1} ),
    files = varargin{1};
    retval = struct([]);
    nTraces = zeros( numel(files),1 );
    
    h = waitbar(0, 'Loading traces and calculating properties...');
    
    for i=1:numel(files),
        data = loadTraces( files{i} );
        retval = [retval  traceStat_data( data.donor, data.acceptor, data.fret )  ];
        nTraces(i) = size(data.donor,1);

        waitbar(i/numel(files),h);
    end
    close(h);
   
% Assume the user passed trace data directly.
else
    retval = traceStat_data( varargin{:} );
end
    


end



function retval = traceStat_data( donorAll,acceptorAll,fretAll, constants )

[Ntraces,len] = size(fretAll);


if ~exist('constants','var')
    constants = cascadeConstants;
end



% Value are returned into a structure that will be later be saved into 
% Application data as filtering criteria
z = num2cell( zeros(1,Ntraces) );

retval = struct( ...
    'corr',  z, ...
    'corrd', z, ...
    'snr',   z, ...
    'snr_s', z, ...
    'nnr',   z, ...
    'bg',    z, ...
    't',     z, ...
    'd',     z, ...
    'a',     z, ...
    'maxFRET',  z, ...
    'ncross',   z, ...
    'lifetime', z, ...
    'acclife',  z, ...
    'donorlife',  z, ...
    'overlap',  z, ...
    'safeRegion', z, ...
    'avgfret',  z, ...
    'fretEvents', z, ...
    'firstFRET', z, ...
    'timeSaturated', z ...
);



%----- CALCULATE STATISTICS FOR EACH TRACE

for i=1:Ntraces
    
    donor    = donorAll(i,:);
    acceptor = acceptorAll(i,:);
    fret     = fretAll(i,:);
    total    = donor+acceptor;
    
    %---- Calculate donor lifetime
    % FRET is undefined when donor is dark, so it is always set
    % to exactly 0.  The last non-0 point is the photobleaching point.
    % See the end of gettraces. This is a hack. FIXME.
    lt = find( fret~=0, 1,'last' )+1;
    if ~isempty(lt) && lt<len,
        retval(i).lifetime = lt;
    end

    
    %---- OVERLAP DETECTION: Find multiple photobleaching events
    % Ignore first frame of median filtered signal: my be very low b/c of 0
    % padding.
    total2 = constants.gamma*donor + acceptor;
    filt_total  = medianfilter(total2,constants.TAU);
    dfilt_total = gradient(filt_total);
    thresh = mean(dfilt_total) - constants.overlap_nstd*std(dfilt_total);
    
    dips = dfilt_total<=thresh;  %all points beyond threshold
    events = find( ~dips & [dips(2:end) 0] )+1;  %start points of drops in fluor

    lastPB = 0;

    for j=1:length(events)-1, %last event is PB
        point = events(j);
        if point<3, continue; end %ignore first few frames
        
        minb = min( filt_total(2:point-1 ) );
        maxa = max( filt_total(point+1:lt) );

        if maxa < minb,
            lastPB = point+2;
        end
    end

    retval(i).overlap = lastPB~=0;
    retval(i).safeRegion = max(1,lastPB+1);
    
    
    %---- Ignore regions where Cy3 is blinking
    s = lt+5;
    bg_range = s:min(s+constants.NBK,len);
    stdbg = std( total(bg_range) );

    % Calculate number of Cy3 PB threshold crossings per frame    
    if lt<8
        donorRange = false(1,lt);
        retval(i).donorlife = lt;
    else
        % Calculate blinking total fluor cutoff value        
        % Find start points where Cy3 blinks
        x = total(1:lt-3) <= constants.blink_nstd*stdbg;
        retval(i).ncross = sum(  x & ~[0 x(1:end-1)]  );
        
        % Remove from consideration regions where Cy3 is photobleached
        donorRange = ~x;  %logical index w/o Cy3 blinking region
        
        % Save the length of the "donor-alive" region as the lifetime
        % of the donor fluorophore.
        retval(i).donorlife = sum(donorRange);
        
        % Remove falling edges of the Cy3 blinks
        donorRange = donorRange & [donorRange(2:end) 1] & [1 donorRange(1:end-1)];
    end


    %---- 
    if sum(donorRange) > 2
    
        % Calculate average amplitudes
        retval(i).d = mean( donor(donorRange) );
        retval(i).a = mean( acceptor(donorRange) );
        retval(i).t = constants.gamma*retval(i).d + retval(i).a;
        
    end
    
    

    % Calculate Signal-to-noise ratio and background noise    
    % should be using bg_range for this...
    if lt+10 < len
        retval(i).snr = retval(i).t/stdbg;  % assuming background corrected
        
        retval(i).bg = std(donor(s:end))+std(acceptor(s:end));
    end
    
    
    if sum(donorRange) > 2
        % Calculate gradients (derivatives) for correlation.
        del_donor    = gradient(donor);
        del_acceptor = gradient(acceptor);
    
        % Calculate correlation of derivitive of donor/acceptor signals.
        % This is the method used by {Fei & Gonzalez, 2008}
        ccd = corrcoef( del_donor(donorRange), del_acceptor(donorRange) );
        retval(i).corrd = ccd(1,2);
        
        % Calcualte correlation coefficient
        cc = corrcoef( donor(donorRange), acceptor(donorRange) );
        retval(i).corr = cc(1,2);
        
        total_noise = std( constants.gamma*donor(donorRange)+acceptor(donorRange) );
        
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

        retval(i).acclife = sum(fretRange);
        
        if sum(fretRange)>1
            retval(i).avgfret = mean( fret(fretRange) );
        end
    end
    
    % Statistics based on FRET distribution
    retval(i).firstFRET = fret(1);
    retval(i).maxFRET = max(fret);
    
    
    % Number of events crossing an arbitrary threshold
    % TODO?: additional filtering to detect only anticorrelated events?
    [result] = RLEncode( fret > constants.fretEventTreshold );
    retval(i).fretEvents = sum( result(:,1)==1 );
    
    % Number of frames with saturated intensity (max of 16-bit int)
    retval(i).timeSaturated = ...
           sum( donor(donorRange)>30000 | acceptor(donorRange)>30000 );
end


% END FUNCTION infoCalc
end
