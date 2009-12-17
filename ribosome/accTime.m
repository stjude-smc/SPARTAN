function output = accTime( tracesFiles, titles )
%ACCPLOT  Plot time-dependant progrssion to tRNA accommodated state
% 
%   OUTPUT = accTime(FILENAMES)
%   Calculates the fraction of molecules that have achieved the
%   accommodated state (tRNA stably bound in the A-site following delivery)
%   as a function of time. FILENAMES is a cell array of filenames to
%   process, where corrosponding DWT files are expected.
%
%   OUTPUT is a set of column vectors with the time axis and the
%   accommodation curves for each of the given files.
%   

% Created by: Daniel Terry (Scott Blanchard Lab)
% Cascade smFRET Analysis Pipeline, Copyright (C) 2009 Daniel Terry
% Date Created: June 18, 2009


%%

sumlen = 120; % in seconds
cutoffTime = 1; % in seconds



%%

output = [];

if nargin<1,
    tracesFiles = getFiles;
end

nFiles = numel(tracesFiles);
if nFiles<1, return; end

% Load model for idealization
constants = cascadeConstants;

model = qub_loadModel( [constants.modelLocation filesep 'tet_selection.qmf'] );
model.fixMu    = ones( model.nStates,1 );
model.fixSigma = ones( model.nStates,1 );
fretModel = [model.mu' model.sigma'];
skmParams.quiet = 1;

for i=1:nFiles,
    
    % Make sure the data has been idealized...
    if ~exist(tracesFiles{i},'file'),
        error('No dwell-times file found!');
    end
    
    % Load FRET data
    [d,a,fret,ids,time] = loadTraces( tracesFiles{i} );
    [nTraces,traceLen] = size(fret);
    clear d; clear a;
    assert( time(1)~=1, 'No time axis found' );
    sampling = time(2)-time(1);

    % Idealize FRET data
    dwtFilename = strrep(tracesFiles{i},'.txt','.qub.dwt');
    [dwt,newModel,LL,offsets] = skm( fret, sampling, model, skmParams );
    saveDWT( dwtFilename, dwt, offsets, fretModel, sampling );
    
    % Load idealization
%     [dwt,sampling,offsets,model] = loadDWT(dwtFilename);
%     assert( sampling==100 ); %ms

    % Find time at which accommodation occurs: estimated as
    % the point  at which a stable ~0.55 FRET state is achieved > 1 sec.
    accTime = repmat(-1,numel(dwt),1);

    for j=1:numel(dwt)

        states = double( dwt{j}(:,1) );
        times  = double( dwt{j}(:,2) ) .* sampling/1000;
        timeline = cumsum( [0; times] );

        % Find the first (if any) dwell in high-FRET longer than cutoffTime sec.
        selection = (states==3) & (times>=cutoffTime);
        idx = find( selection, 1, 'first' );

        if ~isempty(idx),
            accTime(j) = timeline(idx);
        end
    end
    
    % Save accommodation progression curve to output
    [N,X] = hist(accTime(accTime>=0), 0:(sampling/1000):sumlen );    
    NC = cumsum(N);

    N = reshape( N, [1 numel(N)] );
    X = reshape( X, [1 numel(X)] );

    if isempty(output)
        output = [X' (NC/nTraces)'];
    else
        output = [output (NC/nTraces)'];
    end

    
end %for each file

cla;
stairs( output(:,1), output(:,2:end), 'LineWidth',2 );
ylim( [0,1] );
xlim( [0,sumlen] )
set(gca,'xtick',0:30:sumlen);
xlabel('Time (sec)');
ylabel('Fraction Accommodated');

if nargin>1, legend(titles); end

save('accTime.txt','output','-ASCII');


