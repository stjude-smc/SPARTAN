function cplot(hist2d, bounds, constants)
% CPLOT  Draws FRET 2D contour plot
%
%   CPLOT( HIST, BOUNDS, CONST )
%   Draws a contour plot from the FRET histogram (HIST)
%   BOUNDS(1,2) define X-axis limits; BOUNDS(3,4) define Y-axis limits
%   
%   If no histogram specified, the user will be propted for a cplot file.
%   Used by: autotrace, makeplots.

% Parse input arguments
if nargin<1
    [histfile histpath]=uigetfile('*.txt','Choose a contour plot:');
    infile=strcat(histpath,histfile);
    hist2d=load(infile);
    
    figure();
end

if nargin<3
    constants = cascadeConstants();
end
scale = constants.cplot_scale_factor;

if nargin<2
    bounds = [1 constants.contour_length -0.1 1.0];
end



% Load Stanford colormap
cmap = dlmread( 'frethist_colormap.txt' )/255;


% Setup plot axes and contour levels
% time_axis = hist2d(1,2:end);
fret_axis = hist2d(2:end,1);
% bounds(2) = min(bounds(2),time_axis(end));

max_mol = sum( hist2d(2:end,2) )/scale;     %red=?% of total
nl = size(cmap,1)-1;                        %number of contour levels
con = 0:(max_mol/nl):max_mol;               %contour levels

% If the top contour levels are not filled, the levels get distorted.
% Adds a permanent, very high peak in the corner to prevent this.
hist2d(end,bounds(2)) = max_mol*2;


% Draw the filled contour plot in current axis
[C,hand] = contourf( ...
        bounds(1):bounds(2), fret_axis, ...
        hist2d( 2:end, (bounds(1):bounds(2))+1 ), con );

colormap(cmap);
set(hand, 'LineColor', 'none');

% set(gca, 'PlotBoxAspectRatio', [1.5 2 1]);
% set(gca,'xtick',bounds(1):10:bounds(2))
set(gca,'ytick', 0:0.2:bounds(4))
axis( bounds );

xlabel('Time (frames)');
ylabel('FRET');

zoom on;



