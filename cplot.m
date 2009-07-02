% function cplot(hist2d, bounds, constants)
function cplot( varargin )
% CPLOT  Draws FRET 2D contour plot
%
%   CPLOT( HIST, BOUNDS, CONST )
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
    [histfile histpath]=uigetfile('*.txt','Choose a contour plot:');
    infile=strcat(histpath,histfile);
    hist2d=load(infile);
    
    figure();
end

if numel(args)>=2
    bounds = args{2};
else
    bounds = [1 constants.contour_length -0.1 1.0];
end

if numel(args)>=3,
    constants = args{3};
else
    constants = cascadeConstants();
end
scale = constants.cplot_scale_factor;


% Time-bin the data if histogram axes are very long
binFactor = ceil( (bounds(2)-bounds(1)+1)/50 );

time_axis = hist2d(1,2:end);
fret_axis = hist2d(2:end,1);
histdata  = hist2d(2:end,2:end);

if binFactor>1,
    [nRows,nCols] = size(histdata);
    
    time_axis = time_axis(1:binFactor:end);
    nTime = numel(time_axis);
    
    histdata_b  = zeros( numel(fret_axis), nTime );
    
    for i=1:nRows
        for j=1:nTime,
            idx = (1:binFactor) + binFactor*(j-1);
            idx = idx(idx<=nCols);
            histdata_b(i,j) = mean( histdata(i,idx) );
        end
    end
    
    hist2d = zeros(nRows+1,nTime+1);
    hist2d(2:end,1) = fret_axis;
    hist2d(1,2:end) = time_axis;
    hist2d(2:end,2:end) = histdata_b;
end


% Load Stanford colormap
cmap = dlmread( 'frethist_colormap.txt' )/255;


% Setup plot axes and contour levels
% bounds(2) = min(bounds(2),time_axis(end));

max_mol = sum( hist2d(2:end,2) )/scale;     %red=?% of total
nl = size(cmap,1)-1;                        %number of contour levels
con = 0:(max_mol/nl):max_mol;               %contour levels

% If the top contour levels are not filled, the levels get distorted.
% Adds a permanent, very high peak in the corner to prevent this.
hist2d(end,bounds(2)) = max_mol*2;


% Draw the filled contour plot in current axis
lims = ceil( (bounds(1):bounds(2))/binFactor );

[C,hand] = contourf( cax, ...
        time_axis(lims), fret_axis, ...
        hist2d( 2:end, lims+1 ), con );

colormap(cax,cmap);
set(hand, 'LineColor', 'none');

% set(gca, 'PlotBoxAspectRatio', [1.5 2 1]);
% set(gca,'xtick',bounds(1):10:bounds(2))
set(cax,'ytick', 0:0.2:bounds(4))
axis( cax, bounds );

xlabel(cax,'Time (frames)');
ylabel(cax,'FRET');

zoom(cax,'on');



