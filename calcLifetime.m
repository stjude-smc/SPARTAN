function life = calcLifetime(total,TAU,NSTD)
% CALCLIFETIME  Determine photobleaching lifetimes
% 
%   LIFE = calcLifetime( DATA, TAU, NSTD )
%   calculates the length of each trace before the last photobleaching
%   event, which is assumed to occur when the signal falls below a
%   threshold specified by NSTD. The DATA are first median filtered using a
%   window size of TAU to reduce the effect of noise. DATA is typically the
%   sum of donor and acceptor intensities in FRET experiments.

[Ntraces,len] = size(total);
life = repmat(len,[Ntraces 1]);  %default value = maximum trace length

if ~exist('TAU','var')
    constants = cascadeConstants();
    TAU  = constants.TAU;
    NSTD = constants.NSTD;
end
    
% Median filter traces to remove high frequency noise.
% SLOWEST STEP IS HERE!
filt_total  = medianfilter(total,TAU);

dfilt_total = gradient(filt_total);
thresh = mean(dfilt_total,2) - NSTD*std(dfilt_total,0,2);


% For each trace, calc Cy3 lifetime
for i=1:Ntraces
    
    % Find Cy3 photobleaching event by finding *last* large drop in total
    % fluorescence in the median filtered trace
    lt = find( dfilt_total(i,:)<=thresh(i), 1,'last' );
    
    if ~isempty(lt)
        life(i) = max(2,lt);
    end
    
end % for each trace

% END FUNCTION CalcLifetime

