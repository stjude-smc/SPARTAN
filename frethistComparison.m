function f = frethistComparison(files, titles)
% FRETHISTCOMPARISON  Overlays multiple 1D FRET histograms for comparison
%
%   Prompts user for locations of traces files to load.  Creates 1D
%   FRET histograms.  Plots each seperately as stair plots, with
%   different colors and a legend.
%   


% OPTIONS:
nBootstrap = 100;
constants = cascadeConstants();
sampling = 40; %doesn't really matter, as long as its close
contour_bin_size = 0.03;
pophist_sumlen = 50;
fretaxis = (-0.1:contour_bin_size:1.0)';
nbins = length(fretaxis);

removeDarkState = 1;
darkModel = [constants.modelLocation 'darkstate.qmf'];


colors = [ 0  0      0   ; ...                  % black
           0.75   0        0.75 ; ... % purple
           0  0.75   0.75 ; ... % cyan
           0 0.5 0 ;                     % green
           0.75   0.75   0 ; ...      % yellow
           1 0 0   ; ...                  % red
           0.6 0 0   ];                   % dark red



% Prompt user for filenames if not supplied
if nargin < 1,
    disp('Select traces files, hit cancel when finished');
    files = getFiles('*.txt','Choose a traces file:');
end

nFiles = numel(files);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end


%
figure;

%% 
frethist = zeros(nbins,2*nFiles);

model = qub_loadModel( darkModel );
model.mu       = [0.01 0.3];
model.sigma    = [0.061 0.074];
model.fixMu    = ones(1,2);
model.fixSigma = ones(1,2);
skmParams.quiet = 1;

% Load FRET histograms
for i=1:nFiles
    
    % Load FRET data
    [d,a,fret] = loadTraces( files{i} );
    fret = fret(:,1:pophist_sumlen);
    [nTraces,nFrames] = size(fret);
    
    
    % Idealize data to 2-state model to eliminate dark-state dwells
    if removeDarkState
        
        % Idealize to two-state model to select dark state
        dwt = skm( fret, sampling, model, skmParams );

        % Expand DWT into an idealization: 2=non-zero FRET state
        idl = zeros(nTraces,nFrames);

        for dwtID=1:nTraces,
            states = dwt{dwtID}(:,1);
            times  = double(dwt{dwtID}(:,2));

            ends = cumsum( times );
            starts = cumsum( [1; times(1:end-1)] );
            for j=1:numel(states),
                idl(dwtID,starts(j):ends(j)) = states(j);
            end
        end        
    else
        % Use all datapoints for histogram otherwise
        idl = repmat(2,size(fret));
    end %if removeDarkState
    
        
    % Calculate FRET histograms from many bootstrap datasets
    pophist = zeros(nbins,nBootstrap);
    fret = fret(:,1:pophist_sumlen);
    
    for s=1:nBootstrap,
        % Construct bootstrap datasets
        if s==1,
            idxBootstrap = 1:nTraces;
        else
            idxBootstrap = floor(rand(nTraces,1)*nTraces)+1;
        end
        
        data = fret( idxBootstrap, : );     %bootstrap traces
        data = data( idl(idxBootstrap,:)==2 ); %non-zero FRET only
        
        % Create FRET histogram from the bootstrapped dataset
        histdata  = hist( data, fretaxis );
        pophist(:,s) = histdata/sum(histdata);   %normalization

        % Plot FRET histograms
        lw = 4;
        if s~=1,
            lw = 1;
        end
        
        plot( fretaxis, pophist(:,s), 'Color',colors(i,:), 'LineWidth',lw );
        hold on;
        drawnow;
    end
    
    pophistErrors = std(pophist,[],2);
        
    % Add histogram from current dataset to output
    frethist(:,2*i-1) = pophist(:,1);
    
    % Add bootstrapped errors from current dataset
    frethist(:,2*i) = pophistErrors;

    %
%     X = [fret_axis ; fret_axis];
%     Y = [pophist(:,1)-pophistErrors ; pophist(:,1)+pophistErrors];
%     
%     patch( X,Y,colors(i,:)); hold on;
%     alpha(0.3);
%     drawnow;
end

hold off;
ylabel( 'Percent of total time' );

xlabel( 'FRET Efficiency' );
xlim( [0.1 1.0] );


if nargin>=2,
    legend( titles );
end

f = [fretaxis frethist];
save( 'pophist.txt', 'f', '-ASCII' );


end



