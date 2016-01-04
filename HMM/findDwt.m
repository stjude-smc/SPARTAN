function filenames = findDwt( filenames )
% FINDDWT  Look for .dwt files associated with .traces files


for i=1:numel(filenames),
    [p,f,e] = fileparts( filenames{i} );
    
    if strcmpi(e,'.traces'),
        dwtfname = fullfile(p, [f '.qub.dwt']);
        
        if ~exist(dwtfname,'file'),
            dwtfname = fullfile(p, [f '.dwt']);
        end
        
        if ~exist(dwtfname,'file'),
            error('No associated .dwt file found: %s', filenames{i});
        end
        
        filenames{i} = dwtfname;
    end
end



end