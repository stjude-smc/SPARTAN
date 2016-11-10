function output = findDwt( filenames, raiseErrorStr )
% FINDDWT  Look for .dwt files associated with .traces files



%% Process input arguments
narginchk(1,2);

if ischar(filenames),
    filenames={filenames};
end

raiseError = false;
if nargin>1 && strcmpi(raiseErrorStr,'raiseError')
    raiseError = true;
end


%%
output = cell(size(filenames));

for i=1:numel(filenames),
    [p,f,e] = fileparts( filenames{i} );
    
    if strcmpi(e,'.traces'),
        dwtfname = fullfile(p, [f '.qub.dwt']);
        
        if ~exist(dwtfname,'file'),
            dwtfname = fullfile(p, [f '.dwt']);
        end
        
        if exist(dwtfname,'file'),
            output{i} = dwtfname;
        elseif raiseError
            error('No associated .dwt file found: %s', filenames{i});
        end
        
    elseif strcmpi(e,'.dwt'),
        output{i} = filenames{i};
    end
end



end