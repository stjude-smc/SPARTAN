function output = accTime( tracesFiles )
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

output = [];

if nargin<1,
    tracesFiles = getFiles;
end

nFiles = numel(tracesFiles);

for i=1:nFiles,
    
    % Make sure the data has been idealized...
    if ~exist(tracesFiles{i},'file'),
        error('No dwell-times file found!');
    end
    
    % Load FRET data
    [d,a,fret] = loadTraces( tracesFiles{i} );
    [nTraces,traceLen] = size(fret);
    clear d; clear a;

    % Load idealization
    dwtFilename = strrep(tracesFiles{i},'.txt','.qub.dwt');
    [dwt,sampling,offsets,model] = loadDWT(dwtFilename);
    % nTraces = numel(dwt);
    % idl = dwtToIdl( dwt, traceLen, offsets );
    assert( sampling==100 ); %ms

    % Find time at which accommodation occurs: estimated as
    % the point  at which a stable ~0.55 FRET state is achieved > 1 sec.
    accTime = zeros(nTraces,1);

    for j=1:nTraces

        states = double( dwt{j}(:,1) );
        times  = double( dwt{j}(:,2) );
        timeline = cumsum( [1; times] );

        % Find the first (if any) dwell in high-FRET longer than 1 sec
        selection = (states==3) & (times>=30);
        idx = find( selection, 1, 'first' );

        if ~isempty(idx),
            accTime(j) = timeline(idx)*(sampling/1000);
        end
    end
    
    % Save accommodation progression curve to output
    [N,X] = hist(accTime(accTime>0),0:1:600);    
    NC = cumsum(N);
%     stairs(X,NC/nTraces)
%     ylim([0 1]);

    if isempty(output)
        output = [X' (NC/nTraces)'];
    else
        output = [output (NC/nTraces)'];
    end

    
end %for each file

cla;
stairs( output(:,1), output(:,2:end) );
ylim( [0,1] );
xlabel('Time (sec)');
ylabel('Fraction Accommodated');


save('accTime.txt','output','-ASCII');


