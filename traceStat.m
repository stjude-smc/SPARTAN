function retval = traceStat( donor,acceptor,fret, constants )
% LOADTRACES  Loads fluorescence trace files
%
%   [NAMES] = TRACESTAT;
%   Returns the names of each trace statistic as a structure.
%
%   STATS = TRACESTAT( CONST, DONOR,ACCEPTOR,FRET )
%   Calculates metadata from fluorescence/FRET traces that can be
%   used as picking criteria to filter a set of traces.
%   Primarily used in autotrace.m
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
%
%   snr,snr_s,t are all corrected for gamma (sensitivity/quantum yield ratio)
%

% If no data given, return names of filtering criteria.
% TODO: also return comments describing each criteria
if nargin < 3,
    % Order matters -- will be displayed this way in autotrace!
    ln.t        = 'Mean Total Intensity';
    ln.maxFRET  = 'Highest FRET value';
    ln.acclife  = 'FRET Lifetime';
    ln.lifetime = 'Donor Lifetime';
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
    
    retval = ln;
    return;
end


[Ntraces,len] = size(fret);


if ~exist('constants','var')
    constants = cascadeConstants;
end


lifetime = zeros(Ntraces,1);

% Calculate Donor fluor. lifetime;
% FRET is undefined when donor is dark, so it is always set
% to exactly 0.  The last non-0 point is the photobleaching point.
% Bit of a hack, I know, but it works great.
for i=1:Ntraces,
    lt = find( fret(i,:)~=0, 1,'last' )+1;
    if ~isempty(lt) && lt<len,
        lifetime(i) = lt;
    end
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
    'lifetime', num2cell(lifetime'), ...
    'acclife',  z, ...
    'overlap',  z, ...
    'safeRegion', z, ...
    'avgfret',  z ...
);

clear lifetime;



%----- OVERLAP DETECTION: Find multiple photobleaching events
% Ignore first frame of median filtered signal: my be very low b/c of 0
% padding.
total2 = constants.gamma*donor + acceptor;
filt_total  = medianfilter(total2,constants.TAU);
dfilt_total = gradient(filt_total);
thresh = mean(dfilt_total,2) - constants.overlap_nstd*std(dfilt_total,0,2);
    
for i=1:Ntraces,
    lt  = retval(i).lifetime;
    
    dips = dfilt_total(i,1:lt)<=thresh(i);  %all points beyond threshold
    events = find( ~dips & [dips(2:end) 0] )+1;  %start points of drops in fluor

    lastPB = 0;

    for j=1:length(events)-1, %last event is PB
        point = events(j);
        if point<3, continue; end %ignore first few frames
        
        minb = min( filt_total(i,2:point-1 ) );
        maxa = max( filt_total(i,point+1:lt) );

        if maxa < minb,
            lastPB = point+2;
        end
    end

    retval(i).overlap = lastPB~=0;
    retval(i).safeRegion = max(1,lastPB+1);
end



%----- CALCULATE STATISTICS FOR EACH TRACE
total = donor+acceptor;
del_donor    = gradient(donor);
del_acceptor = gradient(acceptor);

for i=1:Ntraces
    lt = retval(i).lifetime;
    retval(i).lifetime = lt;
    
    %---- Ignore regions where Cy3 is blinking
    s = lt+5;
    bg_range = s:min(s+constants.NBK,len);
    stdbg = std( total(i,bg_range) );

    % Calculate number of Cy3 PB threshold crossings per frame    
    if lt<8
        donorRange = logical( zeros(1,lt) );
    else
        % Calculate blinking total fluor cutoff value        
        % Find start points where Cy3 blinks
        x = total(i,1:lt-3) <= constants.blink_nstd*stdbg;
        retval(i).ncross = sum(  x & ~[0 x(1:end-1)]  );
        
        % Remove from consideration regions where Cy3 is photobleached
        donorRange = ~x;  %logical index w/o Cy3 blinking region
        
        % Remove falling edges of the Cy3 blinks
        donorRange = donorRange & [donorRange(2:end) 1] & [1 donorRange(1:end-1)];
        
    end


    %---- 
    if sum(donorRange) > 2
    
        % Calculate average amplitudes
        retval(i).d = mean( donor(i,donorRange) );
        retval(i).a = mean( acceptor(i,donorRange) );
        retval(i).t = constants.gamma*retval(i).d + retval(i).a;
        
    end
    
    

    % Calculate Signal-to-noise ratio and background noise    
    % should be using bg_range for this...
    if lt+10 < len
        retval(i).snr = retval(i).t/stdbg;  % assuming background corrected
        
        retval(i).bg = std(donor(i,s:end))+std(acceptor(i,s:end));
    end
    
    
    if sum(donorRange) > 2
    
        % Calculate correlation of derivitive of donor/acceptor signals.
        % This is the method used by {Fei & Gonzalez, 2008}
        ccd = corrcoef( del_donor(i,donorRange), del_acceptor(i,donorRange) );
        retval(i).corrd = ccd(1,2);
        
        % Calcualte correlation coefficient
        cc = corrcoef( donor(i,donorRange), acceptor(i,donorRange) );
        retval(i).corr = cc(1,2);
        
        total_noise = std( constants.gamma*donor(i,donorRange)+acceptor(i,donorRange) );
        
        retval(i).snr_s = retval(i).t ./ total_noise;
                      
        retval(i).nnr = total_noise ./ stdbg;
    end
    
    
    % Find regions (before Cy3 photobleach) where FRET is above a threshold.
    % Properties will be calculated only on these areas.
    if lt>1
        fretRange = fret(i,1:lt) >= constants.min_fret;

        % Filter the regions so that they must consist of more than 5
        % consecutive points above the threshold
        fretRange = rleFilter( fretRange, constants.rle_min );

        retval(i).acclife = sum(fretRange);
        
        if sum(fretRange)>1
            retval(i).avgfret = mean( fret(i,fretRange) );
        end
    end
    
    % Statistics based on FRET distribution
    %retval(i).maxFRET = fret(i,1);
    retval(i).maxFRET = max(fret(i,:));
end


% END FUNCTION infoCalc
