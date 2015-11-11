function pickRepresentativeTraces( filenames )
% Select traces with average intensity and lifetime.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


%%
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
    criteria.min_t = mean(t)-0.2*std(t);
    criteria.max_t = mean(t)+0.2*std(t);
    criteria.min_donorlife = 0.9*mean(donorlife);
    criteria.max_donorlife = 1.1*mean(donorlife);
    criteria.max_bg = 2*median(bg);
    criteria.eq_overlap = 0;
    
    % Select traces and save to a new file.
    ind = pickTraces(stats, criteria);
    data.subset(ind);
    
    [p,f] = fileparts(filenames{i});
    saveTraces( fullfile(p, [f '_sel.traces']), data );
    
    % Look for an idealization file and select those traces
    dwtfname = fullfile(p, [f '.qub.dwt']);
    
    if exist(dwtfname, 'file')
        [dwt,sampling,offsets,model] = loadDWT(dwtfname);
        idl = dwtToIdl(dwt, offsets, data.nFrames, data.nTraces);
        [dwt,offsets] = idlToDwt( idl(ind,:) );
        
        dwtoutname = fullfile(p, [f '_sel.qub.dwt']);
        saveDWT( dwtoutname, dwt,offsets,model,sampling );
    end
end


end %function
