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

%   Copyright 2007-2015 Cornell University All Rights Reserved.


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

% Construct simple 3-state idealization model
model = QubModel(3);
model.mu    = [0.01 0.25 0.56];
model.sigma = [0.061 0.8 0.15];
model.p0 = [0.8 0.1 0.1];
model.rates = [0    0.1   0     ;
               0.1  0     0.01  ;
               0    1e-5  0    ];
model.fixMu(:)    = true;
model.fixSigma(:) = true;
skmParams.quiet = 1;

for i=1:nFiles,
    
    % Make sure the data has been idealized...
    if ~exist(tracesFiles{i},'file'),
        error('No dwell-times file found!');
    end
    
    % Load FRET data
    data = loadTraces( tracesFiles{i} );
    assert( data.time(1)~=1, 'No time axis found' );

    % Idealize FRET data
    [p,n] = fileparts(tracesFiles{i});
    dwtFilename = fullfile( p, [n '.qub.dwt'] );
    idl = skm( data.fret, data.sampling, model, skmParams );
    [dwt,offsets] = idlToDwt(idl);
    saveDWT( dwtFilename, dwt, offsets, model, data.sampling );
    
    % Load idealization
%     [dwt,sampling,offsets,model] = loadDWT(dwtFilename);
%     assert( sampling==100 ); %ms

    % Find time at which accommodation occurs: estimated as
    % the point  at which a stable ~0.55 FRET state is achieved > 1 sec.
    accTime = repmat(-1,numel(dwt),1);

    for j=1:numel(dwt)

        states = double( dwt{j}(:,1) );
        times  = double( dwt{j}(:,2) ) .* data.sampling/1000;
        timeline = cumsum( [0; times] );

        % Find the first (if any) dwell in high-FRET longer than cutoffTime sec.
        selection = (states==3) & (times>=cutoffTime);
        idx = find( selection, 1, 'first' );

        if ~isempty(idx),
            accTime(j) = timeline(idx);
        end
    end
    
    % Save accommodation progression curve to output
    [N,X] = hist(accTime(accTime>=0), 0:(data.sampling/1000):sumlen );
    NC = cumsum(N);

    N = reshape( N, [1 numel(N)] );
    X = reshape( X, [1 numel(X)] );

    if isempty(output)
        output = [X' (NC/data.nTraces)'];
    else
        output = [output (NC/data.nTraces)'];
    end

    
end %for each file

figure;  cax = axes;
stairs( cax, output(:,1), output(:,2:end), 'LineWidth',2 );
ylim( cax, [0,1] );
xlim( cax, [0,sumlen] )
set(cax,'xtick',0:30:sumlen);
xlabel(cax, 'Time (sec)');
ylabel(cax, 'Fraction Accommodated');

if nargin>1, legend(cax, titles); end

save('accTime.txt','output','-ASCII');


