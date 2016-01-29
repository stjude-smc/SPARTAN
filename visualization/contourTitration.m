function hist_out = contourTitration( files )
% contourTitration   
%
%   Draws a contour plot with the titrated condition on the x axis and the
%   FRET histogram on the y axis.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


%=============  Process input arguments

if nargin<1,
    files = getFiles;
end

% Make sure file list is always a structure array.
if nargin==1 && ischar(files),
    files = {files};
end


%=============  Setup parameters.
constants = cascadeConstants();
options = constants.defaultMakeplotsOptions;


%=============  1. Make histograms for each condition.
% Easiest way to do this is to run frethistComparison. This script will ask
% the user for the filenames itself, if not specified.
frethist = frethistComparison( files );

fret_axis = frethist(:,1);
frethist = frethist(:,2:end); %remove fret axis


%=============  2. Build the contour plot.

% Specify how many pixels wide each condition will be in the x-axis.
% This only changes how the plot is displayed. Change the value to 0 to
% disable.
% options.plotWidth = 10;
width = 10;  %number of pixels on either side


%
X = linspace(0,2,1+2*width);
taper = normpdf(X,1,0.4);
taper = taper./max(taper);  %normalize so maximum is unity

nX = numel(X);
nY = size(frethist,1);

hist2d = zeros( nY, nX*numel(files) );

for i=1:numel(files),
    % Fill the area for this file with the tapering curve.
    % Distribute the 1D histograms across the full area.
    hist2d( :, (i-1)*nX+(1:nX) ) = ...
        repmat(taper,nY,1)  .*  repmat( frethist(:,i), 1,nX );
end

% Add the x-axis. Without more information, we don't know what it really
% is, so just number them for now.
time_axis = linspace(1-0.5,numel(files)+0.5,numel(files)*nX);

% Save the histogram for importing into Origin.
hist_out = zeros( 1+size(hist2d) );
hist_out(1,2:end) = time_axis;
hist_out(2:end,1) = fret_axis;
hist_out(2:end,2:end) = hist2d;

save( 'contourTitration.txt', 'hist_out', '-ASCII' );


%=============  3. Display the plot

max_mol = 1/options.cplot_scale_factor;
nl = size(options.cmap,1)-1;                        %number of contour levels
con = 0:(max_mol/nl):max_mol;               %contour levels

% If the top contour levels are not filled, the levels get distorted.
% Adds a permanent, very high peak in the corner to prevent this.
hist2d(end,end) = max_mol*2;

% Draw the filled contour plot in current axis
% [cax,args] = axescheck(varargin{:});
% cax = newplot(cax);
figure;
cax = axes;

[~,hand] = contourf( cax, time_axis, fret_axis, hist2d, con );

colormap(cax,options.cmap);
set(hand, 'LineColor', 'none');

set(cax,'xtick', 1:numel(files));
set(cax,'ytick', 0:0.2:1)
ylim(cax, [0.2,1] );

xlabel(cax,'Condition');
ylabel(cax,'FRET');

zoom(cax,'on');





end %function contourTitration




