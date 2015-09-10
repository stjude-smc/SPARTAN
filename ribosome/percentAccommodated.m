function [fractionAccommodated,stdFA] = percentAccommodated()
% percentAccommodated    Calculates the extent of acceptor labeling
%
%   ACC = percentAccommodated( FILENAMES )
%   Returns the fraction of all single, donor-labeled particles that also
%   have an acceptor dye. This is calculated as the fraction of traces
%   with an average FRET value > 0.2.
%   If no files are specified, the user will be prompted for them.
%

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% For each file...
fractionAccommodated = [];
stdFA = [];

while 1,
    % Get filenames (could also do this by selecting directories...
    [f,p] = uigetfile('*.traces', 'Select traces files', 'MultiSelect','on');
    if p==0, return; end %user hit cancel
    filenames = strcat(p,f);
    if ~iscell(filenames),  filenames = {filenames};  end
    nFiles = numel(filenames);
    
    % Calculate fraction accommodated for each file
    acc = zeros(nFiles,1);
    for i=1:nFiles
        acc(i) = calcFractionAccommodated( filenames{i} );
    end
    
    % Average over multiple equilibrium movies
    fractionAccommodated(end+1) = mean(acc);
    stdFA(end+1) = std(acc);
    
    % Output results
    [p,name] = fileparts( filenames{1} );
    disp( sprintf('%-30s: %.1f%% (%.1f%%)', ...
            name,100*fractionAccommodated(end),100*stdFA(end)) );
end





function fractionAccommodated = calcFractionAccommodated( filename )

% Define filtering criteria
selectionCriteria.min_snr       = 6;
selectionCriteria.eq_overlap    = 0;
selectionCriteria.maxTotalSigma = 2;
selectionCriteria.max_ncross    = 4;

% Calculate trace statistics
stats = traceStat( filename );

% Filter traces based on standard criteria
idxSelected = pickTraces( stats, selectionCriteria );
nTraces = numel(idxSelected);

% Calculate percent accommodated
avgFret = [stats.avgFret];
avgFret = avgFret(idxSelected);

fractionAccommodated = sum(avgFret>0.2)/nTraces;




