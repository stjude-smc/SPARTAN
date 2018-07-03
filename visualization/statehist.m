function [shist,histmax] = statehist(varargin)
% STATEHIST  state occupancy FRET histogram
%
%   S = STATEHIST( DWT, DATA, options )   where
%   S is a collection (in cols) of FRET histograms for each state, as
%   calculated from the DWT file.  The first col specifies the bins.
%   DWT is the filename of the idealization file from QuB.
%   DATA is the auto.txt filename containing raw Fluorescence/FRET data.
%
%   [...] = STATEHIST(AX, ...) plots the state FRET histograms in the axes AX.
%
%   OPTIONS (optional), can have any of the following fields:
%    - pophist_sumlen:    number of frames to consider when summing
%         histograms. Used to insure population histograms (e.g., 50 frames)
%         match the state histograms by not considering later data.
% 
%    - fret_axis:  specifies FRET histogram bin centers.
%    - fretField:  name of the field in the Traces object to use.

%   Copyright 2007-2016 Cornell University All Rights Reserved.



%% Process input parameters, setting defaults if not provided.

[ax,args] = axescheck(varargin{:});

% Get filenames from user if not passed
if numel(args)<2
    args{1} = getFile('*.dwt','Choose QuB dwt file:');
    if isempty(args{1}), return; end
    
    args{2} = getFile('*.traces','Choose traces file:');
    if isempty(args{2}), return; end
end

options = cascadeConstants('defaultMakeplotsOptions');
if numel(args)>=3,
    options = mergestruct( options, args{3} );
end

[dwt_input,traces_input] = args{1:2};
colors = QubModelViewer.colors;



%% Load data
if ischar(traces_input) || isa(traces_input,'TracesFret'),
    % Load Traces object from input or file and select appropriate FRET field.
    if ischar(traces_input),
        traces_input = loadTraces(traces_input);
    end

    if isfield(options,'fretField') && ~isempty(options.fretField),
        fret = traces_input.(options.fretField);
    else
        fret = traces_input.fret;
    end
    
elseif isnumeric(traces_input),
    % FRET data provided directly.
    fret = traces_input;

else
    error('statehist: Invalid traces input');
end

[nTraces,traceLen] = size(fret);


% --- Load the dwell-time data and create an idealization
if ischar(dwt_input),
    [dwt,~,offsets] = loadDWT(dwt_input);
    idl = dwtToIdl(dwt,offsets,traceLen,nTraces);
    
elseif isnumeric(dwt_input),
    idl = dwt_input;

else
    error('statehist: Invalid idealization input');
end
nStates = max(idl(:));


% --- Truncate traces so they match the display length in makeplots
if options.truncate_statehist,
    traceLen = options.contour_length;
    idl  = idl(:,1:traceLen);
    fret = fret(:,1:traceLen);
end


%% Create normalized histograms for each state.
fret_axis = options.fret_axis;
shist = zeros( numel(fret_axis), nStates+1 );
shist(:,1) = fret_axis;

for j=1:nStates,
    newdata = hist( fret(idl==j), fret_axis ) /sum(idl(:)~=0);
    shist(:,j+1) = newdata;
end

% Normalize histogram height to plot with most molecules (see makeplots)
if options.cplot_normalize_to_max,
    shist(:,2:end) = shist(:,2:end) * nTraces/options.cplot_normalize_to_max;
end


%% Plot, if applicable
if isempty(ax)
    if nargin>0
        % Do not make plots if raw data is requested
        return;
    else
        ax = gca;
    end
end
cla(ax,'reset');


bins = shist(:,1);
histdata = shist(:,2:end)*100;
[~,nStates] = size(histdata);

% Pad with empty bins for display
df = mean(diff(bins));
bins = [bins(1)-df; bins; bins(end)+df];
histdata = [zeros(1,nStates); histdata; zeros(1,nStates)];

% If requested, remove 0-FRET peak and renormalize
if options.ignoreState0
    histdata = histdata(:,2:end);
    histdata = 100*histdata ./ sum(histdata(:));
    colors = colors(2:end);
end

% If the option is set, rescale so that plots with only a few
% molecules show low occupancy in the statehist.
if options.cplot_normalize_to_max,
    histdata = histdata * nTraces/options.cplot_normalize_to_max;
end

hold(ax,'on');
set(ax,'ColorOrder',vertcat(colors{:}));

% Draw translucent, filled area underneath curves
for j=1:size(histdata,2)
    patch( histdata(:,j), bins, colors{j}, ...
            'EdgeColor','none','FaceAlpha',0.25, 'Parent',ax );
end

% Draw state FRET histograms as solid lines
plot( ax, histdata, bins, '--', 'LineWidth',1.5 );

% Add a line with total occupancy.
totalHist = sum( histdata, 2 );
plot( ax, totalHist, bins, 'k-', 'LineWidth',1.5 );

% Formatting -- duplicated in makeplots...
ylim(ax, options.fretRange);
xlabel(ax,'Counts (%)' );
ylabel(ax, 'FRET');  %fixme
set(ax, 'YGrid','on', 'Box','on');

histmax = max( totalHist(bins>0.05) );


end %FUNCTION statehist


