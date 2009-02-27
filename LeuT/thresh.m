function thresh(fileCatalog)
% KINETIC PARAMETER ESTIMATION  DEC 2, 2008
%  Monitoring conformational change at the intracellular face
%  of LeuT (7/86).
%
%  Daniel Terry

clear;
framerate = 25;

nStates = 3;
mu_start    = [0.01 0.52 0.76];
sigma_start = [0.07 0.07 0.07];
level = mean( mu_start(2:3) );

% Load list of datafiles to process (see distance catalog script)
if ~exist('fileCatalog','var'),
    load('catalog.mat'); %=fileCatalog
end
nFiles = size(fileCatalog,1);


% --- 2. Idealize data using thresholding method

for i=1:nFiles
    
    % Load traces data for next condition
    filename = fileCatalog{i,5};
    if strfind(filename,'.qub.txt'),
        traces = load(filename);
        nd = length(traces);
        traces = reshape(traces,1000,nd/1000)';
    else
        [d,a,traces] = loadTraces( filename );
        clear d;  clear a;
    end
    
    % Idealize data using thresholding method
    [idl,offsets] = tIdealize( traces,level );
    model = [mu_start' sigma_start'];
    
    % Save idealization to file
    dwtFile = strrep(filename,'.qub.txt','.txt');
    dwtFile = strrep(dwtFile,'.txt','.qub.dwt');
    saveDWT( dwtFile, idl, offsets, model, 1000/framerate );
end


