function kineticAnalysis(fileCatalog, options)
%% KINETIC PARAMETER ESTIMATION  DEC 2, 2008
%  Monitoring conformational change at the intracellular face
%  of LeuT (7/86).
%
%  Daniel Terry

% clear;
framerate = 25;


% Define the initial model
nStates = 3;
p0_start = [0.1 0.45 0.45];

A_start       = ones(nStates)/nStates;       % transition probability matrix
x =2*(1/framerate); y=1-((nStates-1)*x);              % initial rate=2/sec
A_start(:) = x;  A_start( find(eye(nStates)) ) = y;

assert( all(sum(A_start,2)-1 <0.0001) );

mu_start    = [0.01 0.52 0.75];
sigma_start = [0.07 0.07 0.07];

% Load list of datafiles to process (see distance catalog script)
if ~exist('fileCatalog','var'),
    load('catalog.mat'); %=fileCatalog
end

% fileCatalog = fileCatalog(1:15,:);

nFiles = size(fileCatalog,1);


%% --- 2. PARAMETER ESTIMATES FOR Baum-Welch TOGETHER
profile clear;

fixMu    = [1 0 1];
fixSigma = [1 0 1];

% If asked to re-idealize using previously optimized parameters,
% skip the Baum-Welch optimization proceedure and just load saved results
if nargin<2,

    % Get parameter estimates many times
    LL = cell(nFiles,1);
    A  = cell(nFiles,1);
    p0 = zeros(nFiles,nStates);
    ps = zeros(nFiles,nStates);
    sigma = zeros(nFiles,nStates);
    mu = zeros(nFiles,nStates);
    
    % Q_parts  = zeros(nRuns-1,nRates);

    % Optimize kinetic parameters
    h = waitbar(0,'Estimating kinetic parameters...');
    for i=1:nFiles

        % Load traces data for next condition
        filename = fileCatalog{i,5};
        
        if ~exist(filename,'file'),
            disp(['Can''t find file ' filename]);
            continue;
        end
        
        if strfind(filename,'.qub.txt'),
            traces = load(filename);
            nd = length(traces);
            traces = reshape(traces,1000,nd/1000)';
        else
            [d,a,traces] = loadTraces( filename );
            clear d;  clear a;
        end

        [LL{i},A{i},mu(i,:),sigma(i,:),p0(i,:),ps(i,:)] = BWoptimize( ...
                    traces, A_start, mu_start, sigma_start, ...
                    p0_start,'FixMu', fixMu, 'FixSigma',fixSigma );

        Qt = framerate.*A{i}';
        disp( [ Qt(2,1) Qt(1,2)+Qt(1,3) Qt(3,2) Qt(2,3) ] );
        
        waitbar(i/nFiles,h);
    end

    close(h);
    
else
    load( 'bw_results.mat' );
end

save( 'bw_results.mat', 'A', 'p0', 'ps', 'mu', 'sigma', 'LL' );


%% --- 3. Save kinetic parameters for later analysis

rates = zeros(nFiles,2);

for i=1:nFiles,
    % column -> row
    rates(i,1) = framerate*A{i}(8);  % k_open->closed (2->3)
    rates(i,2) = framerate*A{i}(6);  % k_closed->open (3->2)
end

% Save data to file
fid = fopen('bwt.txt','w');

for i=1:nFiles,
    title = [fileCatalog{i,1} '/' fileCatalog{i,2} ' ' ...
            fileCatalog{i,3} ' ' fileCatalog{i,4}];
    
    fprintf( fid,'%s\t%f\t%f\n', title,rates(i,1),rates(i,2) );
end

fclose(fid);

return;


%% Idealize

% A{4} = [ ...
%     0.0     0.6    0.05  ;...
%     0.7     0.0    0.5   ;...
%     0.7     5.0    0.0   ] ./framerate;
% 
% A{4}( logical(eye(3)) ) = 1 - sum( A{4}, 2 );
% 
% sum(A{4}')

h = waitbar(0,'Idealizing traces...');
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
    
    % Idealize according to optimized model for each
    model = [mu_start' sigma_start'];
    [idl,offsets,LL2] = idealize( traces, model, p0(i,:), A{i} );
    
    % Save idealization to file
    dwtFile = strrep(filename,'.txt','.qub.dwt');
    saveDWT( dwtFile, idl, offsets, model, 1000/framerate );
    
    waitbar(i/nFiles,h);
end
close(h);


