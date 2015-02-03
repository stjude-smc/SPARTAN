function life = calcLifetime(total,TAU,NSTD)
% CALCLIFETIME  Determine photobleaching lifetimes
% 
%   LIFE = calcLifetime( DATA, TAU, NSTD )
%   calculates the length of each trace before the last photobleaching
%   event, which is assumed to occur when the signal falls below a
%   threshold specified by NSTD. The DATA are first median filtered using a
%   window size of TAU to reduce the effect of noise. DATA is typically the
%   sum of donor and acceptor intensities in FRET experiments.

[Ntraces,nFrames] = size(total);
life = repmat( nFrames, [Ntraces,1] );  %default value = maximum trace length

constants = cascadeConstants();
    
if nargin<3,
    TAU  = constants.TAU;
    NSTD = constants.NSTD;
end

total = double(total'); %makes medianfilter.mex happy.


if numel(total)/2000 > 1000 && constants.enable_parfor,
    % If there is a lot of data and the calculations will take a long time,
    % distribute the work over a number of worker threads.
    pool = gcp('nocreate');
    if isempty(pool),
        pool=parpool('IdleTimeout',360,'SpmdEnabled',false);
    end
    M = pool.NumWorkers;
else
    % For <1000 average length traces, running locally can be faster,
    % particularly if parpool hasn't already been started.
    M = 0;
end


% For each trace, calc Cy3 lifetime
parfor (i=1:Ntraces, M)
    % Median filter traces to remove high frequency noise and find drops in
    % total intensity to detect bleaching steps (or blinking).
    dfilt_total = gradient( medianfilter(total(:,i)',TAU) );
    
    % Exclude "outliers" from std (including bleaching steps).
    % std() should only consider noise, not blinking/bleaching steps.
    outlier_thresh = mean(dfilt_total) + 6*std(dfilt_total);
    outliers = abs(dfilt_total)>outlier_thresh;
    
    thresh = mean( dfilt_total(~outliers) ) - NSTD*std( dfilt_total(~outliers) );
    
    % Find Cy3 photobleaching event by finding *last* large drop in total
    % fluorescence in the median filtered trace. The last frame can have
    % artifacts from medianfiltering, so it is ignored?
    lt = find( dfilt_total(1:end-1)<=thresh, 1,'last' );
    
    if ~isempty(lt)
        life(i) = max(2,lt);
    end
end

end %FUNCTION calcLifetime


