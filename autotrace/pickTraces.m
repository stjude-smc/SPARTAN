function [indexes,values] = pickTraces( stats, criteria  )
% PICKTRACES  Loads fluorescence trace files
%
%   [INDEXES,VALUES] = PICK_TRACES( STATS, CRITERIA )
%   Finds INDEXES of traces whose STATS (from traceStat.m) pass the
%   speficied CRITERIA.  If a criteria is not specified, it is not applied.
%   The stat VALUES of the selected traces are also returned.
%   
%   CRITERIA is a structure, whose fields are named with a qualifier (min,
%   max, eq) and a stat name (from STATS). Below are some examples:
%     criteria.min_lifetime  = 10; %at least 10 frames of donor lifetime.
%     criteria.max_ncross    = 4;  %no more than 4 donor blinking events.
%     criteria.eq_fretEvents = 1;  %exactly 1 FRET event.
%   
%   If a criteria name is not recognized (because the corrosponding STAT is
%   not recognized), it will is ignored. For a better description of what
%   each criteria will select for, see traceStat.m.
%
%   There is one special criterion that does not follow this pattern. This
%   criterion will select only traces whose total intensity (t) is within 2
%   standard deviations of the mean, calculated using histogram fitting.
%     critiera.maxTotalSigma = 2;

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Set defaults for criteria not specified
% defaults = constants.defaultCriteria;
% 
% fields = setdiff( fieldnames(defaults), fieldnames(criteria) );
% for i=1:numel(fields),
%     criteria.(fields(i)) = defaults.(fields(i));
% end

% Recognized special criteria not using the typical format.
special = {'maxTotalSigma', 'random'};


Ntraces = numel(stats);
picks = ones(1,Ntraces);



% Remove fields with empty criteria values
names = fieldnames(criteria);

for i=1:numel(names)
    name = names{i};
    value = criteria.(name);
    
    if isempty(value),
        criteria = rmfield(criteria,name);
    end
end


% Rename fields with new standard naming
fn = fieldnames(criteria);

for i=1:numel(fn)
    if strfind(fn{i},'min_')==1,
        cn = fn{i}(5:end);
        eq = '>';
    elseif strfind(fn{i},'max_')==1,
        cn = fn{i}(5:end);
        eq = '<';
    elseif strfind(fn{i},'eq_')==1,
        cn = fn{i}(4:end);
        eq = '==';
    elseif any(strcmpi(fn{i},special))
        continue; %these are handled later
    else
        warning( 'pickTraces:badCriteria', 'Filtering criteria %s not recognized. Ignoring.', fn{i} );
        continue;
    end
    picks = picks & eval(['[stats.' cn ']' eq 'criteria.' fn{i}]);
end



% Find traces which fit all the picking criteria (binary array)
if isfield(criteria,'maxTotalSigma')
    t = [stats.t];
    
    % Fit distribution to a Gaussian function
    bins = 0:(max(t)/50):max(t);
    [histdata] = hist( t(t>0), bins );
    histdata = histdata / sum(histdata);
    
    f = fit( double(bins'),double(histdata'), 'gauss1' );
    mu = f.b1;
    sigma = f.c1;
    
    picks = picks & (t < mu + sigma*criteria.maxTotalSigma);
    picks = picks & (t > mu - sigma*criteria.maxTotalSigma);
    
    % Display result of fitting (for debugging code).
    %figure;
    %bar( bins, histdata, 1 ); hold on;
    %plot(f);
    %disp( [mu sigma] );
end


% Return the results: picked molecule indexes, and stat values for picks
indexes = find(picks);
values = stats(indexes);


% Random trace removal
if isfield(criteria,'random')
    nTraces = numel(indexes);
    sel = randsample(nTraces, floor(nTraces*criteria.random/100) );
    
    indexes = indexes( sel );
    values  = values( sel );
end











