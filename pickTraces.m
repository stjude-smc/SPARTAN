function [indexes,values] = pickTraces( stats, criteria  )
% PICKTRACES  Loads fluorescence trace files
%
%   INDEXES = PICK_TRACES( STATS, CRITERIA )
%   Finds INDEXES of traces whose STATS (from traceStat.m) pass the
%   speficied CRITERIA.  If a criteria is not specified, it is not applied.
%   
%   CRITERIA is a structure with the following allowed fields:
%     minTotalIntensity        stats.t
%     maxTotalIntensity
%     minIntensityNoise       stats.snr_s
%     maxIntensityNoise
%     minTotalLifetime        stats.lifetime
%     maxTotalLifetime
%     minFretLifetime         stats.acc_life
%     maxFretLifetime
%     minCorrelation          stats.corr
%     maxCorrelation 
%     minCorrD                stats.corrd
%     maxCorrD 
%     minSNR                  stats.snr
%     maxSNR  (not implemented)
%     maxNNR                  stats.nnr
%     maxBackground
%     maxDonorBlinks          stats.ncross
%     minAverageFret          
%     maxAverageFret          stats.fretEvents
%     minFretEvents
%     minFret                 at least 1 frame above this value
%     overlap                 (0=no filtering, 1=remove overlap, 2=only overlap)
%     random                  remove X percent of traces randomly
%     
%   Using these names is recommended for all filtering code, even if it
%   doesn't use this function.

% Set defaults for criteria not specified
% defaults = constants.defaultCriteria;
% 
% fields = setdiff( fieldnames(defaults), fieldnames(criteria) );
% for i=1:numel(fields),
%     criteria.(fields(i)) = defaults.(fields(i));
% end



% Get trace stat values
t       = [stats.t];
% d     = [stats.d];
% a     = [stats.a];
lifetime = [stats.lifetime];
snr     = [stats.snr];
snr_s = [stats.snr_s];
nnr     = [stats.nnr];
maxFRET = [stats.maxFRET];
bg      = [stats.bg];
corr    = [stats.corr];
corrd   = [stats.corrd];
ncross  = [stats.ncross];
acclife = [stats.acclife];
overlap = [stats.overlap];
avgfret = [stats.avgfret];
Ntraces = numel(t);



% Remove fields with empty criteria values
names = fieldnames(criteria);

for i=1:numel(names)
    name = names{i};
    value = criteria.(name);
    
    if isempty(value),
        criteria = rmfield(criteria,name);
    end
end



% Find traces which fit all the picking criteria (binary array)
picks = logical( ones(1,Ntraces) );

if isfield(criteria,'maxTotalSigma') && ~isempty(criteria.maxTotalSigma)
    % Fit distribution to a Gaussian function
    bins = 0:500:30000;
    [histdata] = hist( t(t>0), bins );
    histdata = histdata / sum(histdata);
    
    f = fit( bins',histdata', 'gauss1' );

%     figure;
%     bar( bins, histdata, 1 ); hold on;
%     plot(f);
    
    mu = f.b1;
    sigma = f.c1;
    
    disp( [mu sigma] );
    
    picks = picks & (t < mu + sigma*criteria.maxTotalSigma);
end

if isfield(criteria,'minFretEvents')
    picks = picks & fretEvents > criteria.minFretEvents;
end

if isfield(criteria,'minTotalIntensity')
    picks = picks & t > criteria.minTotalIntensity;
end    

if isfield(criteria,'maxTotalIntensity')
    picks = picks & t < criteria.maxTotalIntensity;
end

if isfield(criteria,'minIntensityNoise')
    picks = picks & snr_s > criteria.minIntensityNoise;
end    

if isfield(criteria,'maxIntensityNoise')
    picks = picks & snr_s < criteria.maxIntensityNoise;
end

if isfield(criteria,'minTotalLifetime')
    picks = picks & lifetime > criteria.minTotalLifetime;
end

if isfield(criteria,'maxTotalLifetime')
    picks = picks & lifetime < criteria.maxTotalLifetime;
end

if isfield(criteria,'minSNR')
    picks = picks & snr > criteria.minSNR;
end

if isfield(criteria,'maxNNR')
    picks = picks & nnr < criteria.maxNNR;
end

if isfield(criteria,'minFret')
    picks = picks & maxFRET >= criteria.minFret;
end

if isfield(criteria,'maxBackground')
    picks = picks & bg < criteria.maxBackground;
end

if isfield(criteria,'minCorrelation')
    picks = picks & corr > criteria.minCorrelation;
end

if isfield(criteria,'maxCorrelation')
    picks = picks & corr <= criteria.maxCorrelation;
end

if isfield(criteria,'minCorrD')
    picks = picks & corrd > criteria.minCorrD;
end

if isfield(criteria,'maxCorrD')
    picks = picks & corrd <= criteria.maxCorrD;
end

if isfield(criteria,'maxDonorBlinks')
    picks = picks & ncross < criteria.maxDonorBlinks;
end

if isfield(criteria,'minFretLifetime')
    picks = picks & acclife > criteria.minFretLifetime;
end

if isfield(criteria,'overlap')
    % Find which molecules to remove by overlap citeria
    if criteria.overlap == 0,  %disable overlap detection
        overlap = zeros( size(overlap) );
    elseif criteria.overlap == 2,  %select ONLY overlapping molecules
        overlap = ~overlap;
    end
    
    picks = picks & ~overlap;
end

if isfield(criteria,'minAverageFret')
    picks = picks & avgfret > criteria.minAverageFret;
end

if isfield(criteria,'maxAverageFret')
    picks = picks & avgfret < criteria.maxAverageFret;
end

% Return the results: picked molecule indexes, and stat values for picks
indexes = find(picks);
values = stats(picks);


% Random trace removal
if isfield(criteria,'random')
    nTraces = numel(indexes);
    sel = randsample(nTraces, floor(nTraces*criteria.random/100) );
    
    indexes = indexes( sel );
    values  = values( sel );
end











