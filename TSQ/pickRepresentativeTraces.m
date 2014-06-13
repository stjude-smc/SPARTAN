function pickRepresentativeTraces( filenames )


% 
if nargin<1,
    filenames = getFiles;
end

if nargin>0 && ~iscell(filenames),
    filenames = {filenames};
end

% Load traces
for i=1:numel(filenames),
    data = loadTraces(filenames{i});
    stats = traceStat(data);
    t = [stats.t];
    donorlife = [stats.donorlife];
    bg = [stats.bg];
    
    % Create selection criteria to only keep traces with donor lifetimes,
    % intensity/SNR, etc near the mean, and remove crappy traces.
    criteria.min_t = mean(t)-0.25*std(t);
    criteria.max_t = mean(t)+0.25*std(t);
    criteria.min_donorlife = 0.8*mean(donorlife);
    criteria.max_donorlife = 1.2*mean(donorlife);
    criteria.max_bg = 2*median(bg);
    criteria.eq_overlap = 0;
    
    
    % Select traces and save to a new file.
    ind = pickTraces( stats, criteria  );
    data.donor    = data.donor(ind,:);
    data.acceptor = data.acceptor(ind,:);
    data.fret     = data.fret(ind,:);
    data.traceMetadata = data.traceMetadata(ind);
    
    [p,f] = fileparts(filenames{i});
    saveTraces( fullfile(p, [f '_sel.traces']), 'traces', data );
end


end %function