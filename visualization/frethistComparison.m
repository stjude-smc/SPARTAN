function f = frethistComparison(files, titles, settings)
% FRETHISTCOMPARISON  Overlays multiple 1D FRET histograms for comparison
%
%   Prompts user for locations of traces files to load.  Creates 1D
%   FRET histograms.  Plots each seperately as stair plots, with
%   different colors and a legend.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


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
removeDarkState = true;

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

if ~iscell(files),
    files = {files};
end

nFiles = numel(files);
if nFiles == 0,  return;  end


% If there are a ton of datapoints, we can't use the simple colors above.
% Instead, just use a simple blue-to-red gradient
if nFiles>size(colors,1),
    colors = zeros(nFiles,3);
    interval = (1/(nFiles-1));
    colors(:,1) = 0:interval:1;
    colors(:,3) = 1:-interval:0;
end

%
hf = figure;
cax = axes('Parent',hf);


% Model for removing dark state noise. Adjust if any state is < 0.4.
if removeDarkState,
    model.p0    = [0.01 0.99]';
    model.rates = [0 5; 1 0];
    model.mu       = [0.01 0.3];
    model.sigma    = [0.061 0.074];
    model.fixMu    = true(1,2);
    model.fixSigma = true(1,2);
    model = QubModel(model);
    
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
    if removeDarkState,
        [dwt,~,~,offsets] = skm( fret, data.sampling, model, skmParams );
        idl = dwtToIdl( dwt, offsets, nFrames, nTraces );
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
    plot( cax, sx, sy, 'Color',colors(i,:), 'LineWidth',3 );
    hold(cax,'on');
    
    % Calculate and plot error bars
    if calcErrorBars
        pophistErrors = std(pophist,[],2);
        frethist(:,2*i) = pophistErrors;
        errorbar( cax, fretaxis, pophist(:,1), pophistErrors/2, 'x', 'Color',colors(i,:), 'LineWidth',1 );
    end
    
    % Add histogram from current dataset to output
    % Add bootstrapped errors from current dataset
    frethist(:,2*i-1) = pophist(:,1); %use the first set: all traces.

    drawnow;
end


% Decorate the plot with axes etc.
hold(cax,'off');
ylabel( cax, 'Percent of total time' );
xlabel( cax, 'FRET Efficiency' );
xlim( cax, [0.1 1.0] );
yl = ylim(cax);
ylim( cax, [0 yl(2)] );

if nargin>=2,
    legend( cax, titles );
end


% Save the results. If calcErrorBars is true, every other column has the
% error bars of the associated histogram.
if ~calcErrorBars,
    frethist = frethist(:,1:2:end); %remove error bar columns.
end
    
f = [fretaxis frethist];
save( 'pophist.txt', 'f', '-ASCII' );


end



