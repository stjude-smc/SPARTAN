function f = frethistComparison(files, titles, settings)
% FRETHISTCOMPARISON  Overlays multiple 1D FRET histograms for comparison
%
%   Prompts user for locations of traces files to load.  Creates 1D
%   FRET histograms.  Plots each seperately as stair plots, with
%   different colors and a legend.
%   


if nargin<3,
    constants = cascadeConstants();
    settings = constants.defaultMakeplotsOptions;
end

% Display settings:
% Most of these are the same as makeplots so that they look the same, but
% you can override those defaults below.
pophist_sumlen = settings.contour_length; % how many frames to use.
pophist_offset = settings.pophist_offset; % first N frames to throw out.
fretaxis = settings.fret_axis';
nbins = length(fretaxis);

colors = [ 0      0      0    ; ...  % black
           0.75   0      0.75 ; ...  % purple
           0      0.75   0.75 ; ...  % cyan
           0      0.5    0    ; ...  % green
           0.75   0.75   0    ; ...  % yellow
           1      0      0    ; ...  % red
           0.6    0      0    ];     % dark red

% Settings for removing the zero-FRET state.
% TURN THIS OFF if you want the full histogram.
removeDarkState = false;
darkModel = [constants.modelLocation 'darkstate.qmf'];
sampling = 20; %ms. unused unless time axis is undefined.

% Settings for error bar calculation.
calcErrorBars = false;  % if true, use bootstrapping to display error bars

if calcErrorBars
    nBootstrap = 100;  %number of bootstrap samples to make.
else 
    nBootstrap = 1;
end




% Prompt user for filenames if not supplied
if nargin < 1 || isempty(files),
    disp('Select traces files, hit cancel when finished');
    files = getFiles([],'Choose a traces file:');
end

nFiles = numel(files);

if nFiles == 0,
    disp('No files specified, exiting.');
    return;
end


% If there are a ton of datapoints, we can't use the simple colors above.
% Instead, just use a simple blue-to-red gradient
if nFiles>size(colors,1),
    colors = zeros(nFiles,3);
    interval = (1/(nFiles-1));
    colors(:,1) = 0:interval:1;
    colors(:,3) = 1:-interval:0;
end

%
figure;


% If the dark state fitting model is not found, ask the user for it.
if removeDarkState && ~exist(darkModel,'file')
    warning('The model used for removing the dark state was not found');

    [f,p] = uigetfile(darkModel,'Select model used for removing the dark state');
    if f==0,
        removeDarkState = 0;
        disp('No model. Giving up on trying to remove the dark state.');
    else
        darkModel = [p f];
    end
end

% Load the dark state model. The parameters can be adjusted, especially if
% the lowest (non-zero) state is lower than 0.4.
if removeDarkState
    model = qub_loadModel( darkModel );
    model.mu       = [0.01 0.3];
    model.sigma    = [0.061 0.074];
    model.fixMu    = ones(1,2);
    model.fixSigma = ones(1,2);
    skmParams.quiet = 1;
end


%% 
frethist = zeros(nbins,2*nFiles);

% Load FRET histograms
for i=1:nFiles
    
    % Load FRET data
    data = loadTraces( files{i} );
    fret = data.fret( :, pophist_offset+(1:pophist_sumlen) );
    [nTraces,nFrames] = size(fret);
    
    
    % Idealize data to 2-state model to eliminate dark-state dwells
    if removeDarkState    
        if data.time(1)==0,
            sampling = diff( data.time(1:2) );
        end
        
        % Idealize to two-state model to select dark state
        disp('Removing dark state using hidden Markov modeling...');
        [dwt,~,~,offsets] = skm( fret, sampling, model, skmParams );
        idl = dwtToIdl( dwt, nFrames, offsets );
    else
        % Use all datapoints for histogram otherwise
        idl = repmat(2,size(fret));
    end
    
        
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
    end
    
    % Spline interpolate the data so it's easier to see (but not saved that way)
    sx = fretaxis(1):0.001:fretaxis(end);
    sy = spline( fretaxis, pophist(:,1), sx );
    plot( sx, sy, 'Color',colors(i,:), 'LineWidth',3 );
    hold on;
    
    % Calculate and plot error bars
    if calcErrorBars
        pophistErrors = std(pophist,[],2);
        frethist(:,2*i) = pophistErrors;
        errorbar( fretaxis, pophist(:,1), pophistErrors, 'x', 'Color',colors(i,:), 'LineWidth',1 );
    end
    
    % Add histogram from current dataset to output
    % Add bootstrapped errors from current dataset
    frethist(:,2*i-1) = pophist(:,1); %use the first set: all traces.

    %
%     X = [fret_axis ; fret_axis];
%     Y = [pophist(:,1)-pophistErrors ; pophist(:,1)+pophistErrors];
%     
%     patch( X,Y,colors(i,:)); hold on;
%     alpha(0.3);
    drawnow;
end


% Decorate the plot with axes etc.
hold off;
ylabel( 'Percent of total time' );

xlabel( 'FRET Efficiency' );
xlim( [0.1 1.0] );
yl = ylim;
ylim( [0 yl(2)] );

if nargin>=2,
    legend( titles );
end


% Save the results. If calcErrorBars is true, every other column has the
% error bars of the associated histogram.
if ~calcErrorBars,
    frethist = frethist(:,1:2:end); %remove error bar columns.
end
    
f = [fretaxis frethist];
save( 'pophist.txt', 'f', '-ASCII' );


end



