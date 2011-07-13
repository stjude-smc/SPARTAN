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
    [d,a,f,ids,time] = loadTraces(filenames{i});
    stats = traceStat(d,a,f);
    snr = [stats.snr];
    donorlife = [stats.donorlife];
    bg = [stats.bg];
    
    % Create selection criteria to only keep traces with donor lifetimes,
    % intensity/SNR, etc near the mean, and remove crappy traces.
    criteria.min_snr = mean(snr)-0.5*std(snr);
    criteria.max_snr = mean(snr)+0.5*std(snr);
    criteria.min_donorlife = 0.8*mean(donorlife);
    criteria.max_donorlife = 1.2*mean(donorlife);
    criteria.max_bg = 2*median(bg);
    criteria.overlap = 1;
    
    
    % Select traces and save to a new file.
    ind = pickTraces( stats, criteria  );
    
    fname = strrep(filenames{i},'.txt','_sel.txt');
    saveTraces( fname, 'txt', d(ind,:),a(ind,:),f(ind,:),ids(ind),time );
end


end %function