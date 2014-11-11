function life = calcLifetime(total,TAU,NSTD)
% CALCLIFETIME  Determine photobleaching lifetimes
% 
%   LIFE = calcLifetime( DATA, TAU, NSTD )
%   calculates the length of each trace before the last photobleaching
%   event, which is assumed to occur when the signal falls below a
%   threshold specified by NSTD. The DATA are first median filtered using a
%   window size of TAU to reduce the effect of noise. DATA is typically the
%   sum of donor and acceptor intensities in FRET experiments.

Ntraces = size(total,1);
life = zeros(Ntraces,1);  %default value = maximum trace length

if nargin<3,
    constants = cascadeConstants();
    TAU  = constants.TAU;
    NSTD = constants.NSTD;
end

total = double(total');

% For each trace, calc Cy3 lifetime
if Ntraces>100,
    parfor i=1:Ntraces
        life(i) = calLifetime2( total(:,i)', TAU, NSTD );
    end
else
    for i=1:Ntraces
        life(i) = calLifetime2( total(:,i)', TAU, NSTD );
    end
end

end %FUNCTION calcLifetime



%%
function lt = calLifetime2( total, TAU, NSTD )


    % Median filter traces to remove high frequency noise.
    dfilt_total = gradient( medianfilter(total,TAU) );
    thresh = mean(dfilt_total) - NSTD*std(dfilt_total);
    
    % Find Cy3 photobleaching event by finding *last* large drop in total
    % fluorescence in the median filtered trace
    lt = find( dfilt_total<=thresh, 1,'last' );
    
    if ~isempty(lt)
        lt = max(2,lt);
    else
        lt = numel(total);
    end
    
end % for each trace

% END FUNCTION CalcLifetime

