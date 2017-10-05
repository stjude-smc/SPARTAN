function [trcfiles,dwtfiles] = findTracesDwt( filenames )
% FINDTRACESDWT  Look for .traces and .dwt files associated with input files



%% Process input arguments
narginchk(1,2);

if ischar(filenames),
    filenames={filenames};
end


%%
trcfiles = cell(size(filenames));
dwtfiles = cell(size(filenames));

for i=1:numel(filenames),
    [p,f] = fileparts( filenames{i} );
    if strcmpi(f(end-3:end),'.qub'), f=f(1:end-4); end
    
    % Find traces file
    trcfname = fullfile(p, [f '.traces']);
    if ~exist(trcfname,'file'), error('Unable to find .traces file'); end
    trcfiles{i} = trcfname;
    
    % Find dwell-time file
    dwtfname = fullfile(p, [f '.qub.dwt']);

    if ~exist(dwtfname,'file'),
        dwtfname = fullfile(p, [f '.dwt']);
    end

    if ~exist(dwtfname,'file'), error('Unable to find .dwt file'); end
    dwtfiles{i} = dwtfname;
end



end