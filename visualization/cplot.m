function output = cplot( varargin )
% CPLOT  Draw 2D FRET contour plot
%
%   CPLOT( HIST, BOUNDS, OPTIONS )
%   Draw a contour plot from the FRET histogram (HIST)
%   BOUNDS(1,2) define X-axis limits; BOUNDS(3,4) define Y-axis limits.
%   OPTIONS is a struct array for additional options -- see cascadeConstants.m.
%   
%   CPLOT( AX, ... ) plots into AX instead of the current axes (GCA).
%
%   HIST = CPLOT(...) returns the histogram as plotted.
%
%   See also: makeplots, makecplot.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


if nargout>0, output = []; end



%% Process input arguments
[cax,args] = axescheck(varargin{:});
isNewAx = isempty(cax);

% Load FRET data from input
if numel(args)>=1,
    if isnumeric(args{1}),
        hist2d = args{1};
    elseif isa(args{1},'Traces')
        hist2d = makecplot(args{1});
    elseif ischar(args{1})
        hist2d = makecplot(loadTraces(args{1}));
    else
        error('Invalid input');
    end
    
% No input data given; ask the user for a .traces file
else
    fname = getFile;
    if isempty(fname), return; end
    hist2d = makecplot(loadTraces(fname));
end


% Parse optional arguments
if numel(args)>=3,
    options = args{3};
else
    c = cascadeConstants();
    options = c.defaultMakeplotsOptions;
end

if numel(args)>=2
    bounds = args{2};
else
    bounds = [1 options.contour_length options.fretRange];
end



%% Bin FRET data into 50 x-axis (time) bins.
% Larger plots (100s to 1000s of points) can cause MATLAB to slow or crash.
binFactor = ceil( (bounds(2)-bounds(1)+1)/50 );

if binFactor>1,
    % This method splits the traces up into binFactor strides for summing
    % over time bins. For this to work, the trace length must be divisble
    % by binFactor. Truncate here to make this so.
    len = size(hist2d,2)-1;
    len = len-mod(len,binFactor);  %truncated size
    time_axis = hist2d(1,2:len+1);
    fret_axis = hist2d(2:end,1);
    histdata  = hist2d(2:end,2:len+1);
    
    time_axis = time_axis(1:binFactor:end); %new (binned) time axis
    nTime = numel(time_axis);
    
    % Average over each time window to bin consecutively.
    histdata_b = zeros( numel(fret_axis), nTime );
    
    for i=1:binFactor,
        histdata_b = histdata_b + histdata(:,i:binFactor:end);
    end
    
    % Rebuild FRET histogram matrix with axes.
    hist2d = zeros( numel(fret_axis)+1, nTime+1 );
    hist2d(2:end,1) = fret_axis;
    hist2d(1,2:end) = time_axis;
    hist2d(2:end,2:end) = histdata_b./binFactor;
end



%% Construct histogram for plotting
time_axis = hist2d(1,2:end);
fret_axis = hist2d(2:end,1);

% Define time region to plot as frame numbers in the time-binned histogram.
lims = ceil( bounds(1:2)./binFactor );
lims = lims(1):lims(2);

% Setup plot axes and contour levels
% bounds(2) = min(bounds(2),time_axis(end));

% max_mol = sum( hist2d(2:end,2) )/options.cplot_scale_factor;  %red=?% of total
max_mol = 1/options.cplot_scale_factor;
nl = size(options.cmap,1)-1;       %number of contour levels
con = 0:(max_mol/nl):max_mol;      %contour levels

% Truncate the plot to the display window.
hist2d = hist2d(:,[1 lims+1]);
if nargout>0, output = hist2d; end

% If the top contour levels are not filled, the levels get distorted.
% Add a permanent, very high peak in the corner to prevent this.
hist2d(end,end) = max_mol*2;



%% Draw the filled contour plot in current axis
cax = newplot(cax);
[~,hand] = contourf(cax, time_axis(lims), fret_axis, hist2d(2:end,2:end), con);
colormap(cax, options.cmap);
set(hand, 'LineColor', 'none');
set(cax,'ytick', 0:0.2:bounds(4));
axis(cax, bounds);

% Add all appearance details only if this is an independent plot.
if isNewAx,
    set(gca, 'PlotBoxAspectRatio', [1.5 2 1]);
    xlabel(cax,'Time (frames)');
    ylabel(cax,'FRET');
    zoom(cax,'on');
end



