function hist2d = cplot( varargin )
% CPLOT  Draws FRET 2D contour plot
%
%   CPLOT( HIST, BOUNDS, OPTIONS )
%   Draws a contour plot from the FRET histogram (HIST)
%   BOUNDS(1,2) define X-axis limits; BOUNDS(3,4) define Y-axis limits
%   
%   If no histogram specified, the user will be propted for a cplot file.
%   Used by: autotrace, makeplots.
%
%   DT 090619: added time-binning for very large histograms
%

[cax,args] = axescheck(varargin{:});
cax = newplot(cax);

% Parse input arguments
if numel(args)>=1,
    hist2d = args{1};
else
    [histfile,histpath] = uigetfile('*.txt','Choose a contour plot:');
    infile=strcat(histpath,histfile);
    hist2d=load(infile);
    
    figure();
end

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

scale = options.cplot_scale_factor;


% Time-bin the data if histogram axes are very long.
% 50 datapoints (frames) gives good resolution; any more is excessive.
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
    
    % Average over each time window to bin consecutively. Since the number
    % of frames per bin will always be much less than the number of
    % timepoints, iterate over each frame in the bin instead of every
    % window (as before). This is much faster.
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

time_axis = hist2d(1,2:end);
fret_axis = hist2d(2:end,1);


% Define time region to plot as frame numbers in the time-binned histogram.
lims = ceil( bounds(1:2)./binFactor );
lims = lims(1):lims(2);

% Setup plot axes and contour levels
% bounds(2) = min(bounds(2),time_axis(end));

% max_mol = sum( hist2d(2:end,2) )/scale;     %red=?% of total
max_mol = 1/scale;
nl = size(options.cmap,1)-1;                        %number of contour levels
con = 0:(max_mol/nl):max_mol;               %contour levels


% Truncate the plot to the display window.
hist2d = hist2d(:,[1 lims+1]);

% If the top contour levels are not filled, the levels get distorted.
% Adds a permanent, very high peak in the corner to prevent this.
hist2d_n = hist2d;
hist2d_n(end,end) = max_mol*2;


% Draw the filled contour plot in current axis
[~,hand] = contourf( cax, ...
        time_axis(lims), fret_axis, ...
        hist2d_n( 2:end, 2:end ), con );

colormap(cax,options.cmap);
set(hand, 'LineColor', 'none');

% set(gca, 'PlotBoxAspectRatio', [1.5 2 1]);
% set(gca,'xtick',bounds(1):10:bounds(2))
set(cax,'ytick', 0:0.2:bounds(4))
axis( cax, bounds );

xlabel(cax,'Time (frames)');
ylabel(cax,'FRET');

zoom(cax,'on');



